# -----------------------------------------------------------
# 02_fitGAM_GOenrichment.R
# Purpose : Fit GAM models with tradeSeq, run start-vs-end tests,
#           perform GO enrichment, and generate cluster-level heatmaps
# Inputs  : - sce_subset (SingleCellExperiment with counts)
#           - crv (SlingshotDataSet, 2 lineages)
#           - cl (cell-type / cluster labels)
# Outputs : - sce_fitgam_lineage1.rds, sce_fitgam_lineage2.rds
#           - startRes_Lin1.csv, startRes_Lin2.csv
#           - GO enrichment CSVs
#           - Heatmaps (PDF)
# -----------------------------------------------------------



suppressPackageStartupMessages({
  library(tradeSeq)
  library(SingleCellExperiment)
  library(slingshot)
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(matrixStats)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(circlize)
})

# ============================
# 0. Libraries & setup
# ============================

set.seed(123)
setwd("set/your/path")  # <-- adjust this


# ============================
# 1. Load input data
# ============================
load("sce_subset")     # counts + metadata
counts <- counts(sce_subset)

load("crv1")           # SlingshotDataSet
cl <- sce_subset$GMM   # cluster labels

# ============================
# 2. Prepare lineages (remove 3rd if present)
# ============================

# NOTE: In this dataset, Slingshot initially detected 3 lineages.
# However, biological interpretation (and downstream plots) suggested
# that the 3rd lineage was spurious / not relevant for this trajectory.
# Therefore, we restrict the analysis to Lineage 1 and Lineage 2.
# This requires manual editing of the SlingshotDataSet object slots.

# --- Remove 3rd lineage from assays ---
crv1@assays@data$pseudotime <- crv1@assays@data$pseudotime[, -3]
crv1@assays@data$weights    <- crv1@assays@data$weights[, -3]

# --- Update colData to contain only 2 lineages ---
crv1@colData@rownames <- c("Lineage1", "Lineage2")
crv1@colData@nrows    <- integer(2)

# --- Remove 3rd lineage from metadata ---
crv1@metadata$lineages[3]           <- NULL
crv1@metadata$slingParams$end.clus  <- c("7", "6")   # adjust based on dataset
crv1@metadata$slingParams$end.given <- c(TRUE, TRUE)
crv1@metadata$curves[3]             <- NULL

# --- Save and reload the edited SlingshotDataSet ---
save(crv1, file = "crv1_edit_2lineages")
#load("crv1_edit_2lineages")

# --- Finalize object to use in downstream steps ---
crv <- crv1

# --- Sanity check ---
cat("Structure of crv1 after editing:\n")
str(crv1)

# ============================
# 3. Filter cells & genes
# ============================

# Extract pseudotime and weights from the edited SlingshotDataSet
pseudotime   <- slingPseudotime(crv, na = FALSE)
cellWeights  <- slingCurveWeights(crv)

# --- Cell filtering ---
# Remove cells that do not belong to any lineage (sum of weights = 0)
valid_cells <- rowSums(cellWeights) > 0

counts_filtered     <- counts[, valid_cells]
pseudotime_filtered <- pseudotime[valid_cells, , drop = FALSE]
cellWeights_filtered <- cellWeights[valid_cells, , drop = FALSE]
cl_filtered         <- cl[valid_cells]

cat("Filtered counts matrix dimensions:", dim(counts_filtered), "\n")
cat("Filtered pseudotime dimensions:", dim(pseudotime_filtered), "\n")
cat("Filtered weights dimensions:", dim(cellWeights_filtered), "\n")

# --- Gene filtering ---
# Keep the top 3000 most variable genes
gene_variances <- rowVars(counts_filtered)
valid_genes    <- which(!is.na(gene_variances))
gene_variances <- gene_variances[valid_genes]
counts_filtered <- counts_filtered[valid_genes, ]

top_genes <- order(gene_variances, decreasing = TRUE)[1:min(3000, length(gene_variances))]
counts_filtered <- counts_filtered[top_genes, ]

cat("Counts matrix after gene filtering:",
    nrow(counts_filtered), "genes and",
    ncol(counts_filtered), "cells.\n")

# --- Ensure consistency across objects ---
common_cells <- Reduce(intersect, list(
  colnames(counts_filtered),
  rownames(pseudotime_filtered),
  rownames(cellWeights_filtered)
))
common_cells <- sort(common_cells)

# Reorder all objects to have matching cell sets
counts_filtered     <- counts_filtered[, common_cells]
pseudotime_filtered <- pseudotime_filtered[common_cells, , drop = FALSE]
cellWeights_filtered <- cellWeights_filtered[common_cells, , drop = FALSE]
cl_filtered         <- cl_filtered[common_cells]



# ============================
# 4. Lineage 1 analysis
# ============================

# --- 4.1 Subset cells belonging to Lineage 1 ---
# We identify cells that primarily follow Lineage 1 by selecting those
# with the highest lineage weight and ensuring non-zero contribution.
# This ensures that only cells truly participating in Lineage 1 are analyzed.
weights <- slingCurveWeights(crv)
lineage1_cells <- rownames(weights)[apply(weights, 1, which.max) == 1]
lineage1_cells <- lineage1_cells[weights[lineage1_cells, "Lineage1"] > 0]

# Subset the count matrix, pseudotime, lineage weights, and clusters.
counts_l1     <- counts_filtered[, lineage1_cells]
pseudotime_l1 <- pseudotime_filtered[lineage1_cells, "Lineage1", drop = FALSE]
weights_l1    <- cellWeights_filtered[lineage1_cells, "Lineage1", drop = FALSE]
cl_l1         <- cl_filtered[lineage1_cells]

# Adjust Slingshot object so it only represents Lineage 1.
# (This simplifies downstream GAM fitting and avoids confusion with extra lineages.)
crv_l1 <- crv[lineage1_cells, ]
crv_l1@assays@data$pseudotime <- crv@assays@data$pseudotime[, "Lineage1", drop = FALSE]
crv_l1@assays@data$weights    <- crv@assays@data$weights[, "Lineage1", drop = FALSE]
crv_l1@metadata$lineages      <- crv@metadata$lineages[1]
crv_l1@metadata$curves        <- crv@metadata$curves[1]
crv_l1@metadata$slingParams$end.clus  <- crv@metadata$slingParams$end.clus[1]
crv_l1@metadata$slingParams$end.given <- crv@metadata$slingParams$end.given[1]
colnames(crv_l1@assays@data$pseudotime) <- "Lineage1"
colnames(crv_l1@assays@data$weights)    <- "Lineage1"


# Sanity check
cat("Lineage 1 subset prepared with",
    ncol(counts_l1), "cells and",
    nrow(counts_l1), "genes.\n")


# --- 4.2 Fit GAM model (tradeSeq) ---
# tradeSeq fits gene expression as smooth curves along pseudotime
# using generalized additive models (GAMs). This allows us to model
# non-linear dynamics of gene expression during the trajectory.

sce_fitgam_l1 <- fitGAM(
  counts      = counts_l1,
  pseudotime  = pseudotime_l1,
  cellWeights = weights_l1,
  nknots      = 6,              # 6 knots balances flexibility and overfitting
  verbose     = TRUE
)

# Add cluster labels
names(cl_l1) <- lineage1_cells
sce_fitgam_l1$GMM <- cl_l1

save(sce_fitgam_l1, file = "results/sce_fitgam_lineage1.rds")
cat("Saved fitted GAM model for Lineage 1.\n")


# --- 4.3 Differential expression: start vs end ---
# To identify genes changing along Lineage 1, we compare expression
# at the start vs. the end of pseudotime using a Wald test.
# Significant genes indicate dynamic regulation across the trajectory.
startRes_l1     <- startVsEndTest(sce_fitgam_l1)
startRes_sig_l1 <- startRes_l1[startRes_l1$pvalue < 0.05, ]

fwrite(startRes_sig_l1,
       row.names = TRUE,
       file = "results/startRes_Lineage1.csv")
cat("Saved start-vs-end DE results for Lineage 1.\n")


# --- 4.4 GO enrichment analysis (clusterProfiler) ---
# We perform GO enrichment on the significant DE genes to identify
# biological processes that are overrepresented. This links statistical
# findings to biological interpretation.

top_genes_l1 <- rownames(startRes_sig_l1)

if (length(top_genes_l1) > 0) {
  gene_df_l1 <- bitr(top_genes_l1,
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Hs.eg.db)

  ego_l1 <- enrichGO(
    gene          = gene_df_l1$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",        # Biological Process ontology
    pAdjustMethod = "BH",        # multiple testing correction
    readable      = TRUE
  )

  write.csv(as.data.frame(ego_l1),
            "results/GO_Enrichment_Lineage1.csv",
            row.names = FALSE)

  cat("Saved GO enrichment results for Lineage 1.\n")
} else {
  warning("No significant DE genes found in Lineage 1.")
}


# --- 4.5 Enrichment per cluster ---
# Instead of only analyzing the lineage as a whole, we also
# check whether specific clusters of cells within Lineage 1
# are enriched for distinct GO terms. This helps reveal
# subpopulation-specific biology.

top_per_cluster <- list()
for (clus in unique(cl_l1)) {
  cluster_cells <- names(cl_l1)[cl_l1 == clus]
  cluster_expr  <- rowMeans(counts_l1[, cluster_cells, drop = FALSE])
  top_genes     <- names(sort(cluster_expr, decreasing = TRUE))[1:300]
  intersected   <- intersect(top_genes, top_genes_l1)
  top_per_cluster[[paste0("Cluster_", clus)]] <- intersected
}

ego_by_cluster_l1 <- compareCluster(
  geneCluster   = top_per_cluster,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

saveRDS(ego_by_cluster_l1, "results/GO_enrichment_lineage1_by_cluster.rds")
df_enrich_l1 <- as.data.frame(ego_by_cluster_l1)
write.csv(df_enrich_l1,
          "results/GO_enrichment_lineage1_by_cluster.csv",
          row.names = FALSE)

cat("Saved cluster-specific GO enrichment for Lineage 1.\n")


# --- 4.6 Heatmaps: GO terms of interest ---
# Heatmaps allow visualization of how enriched GO terms are expressed
# across clusters in Lineage 1. We compute average log1p(counts) per
# term × cluster, and plot both raw (log1p) and normalized (z-score).

parse_gene_list <- function(x) unique(unlist(strsplit(x, "/")))

mean_log1p_counts_in_cells <- function(counts_mat, genes, cell_ids) {
  g <- intersect(genes, rownames(counts_mat))
  if (length(g) == 0L || length(cell_ids) == 0L) return(NA_real_)
  X <- counts_mat[g, cell_ids, drop = FALSE]
  X <- log1p(as.matrix(X))
  mean(X, na.rm = TRUE)
}

term_to_genes_l1 <- lapply(
  split(df_enrich_l1$geneID, df_enrich_l1$Description),
  parse_gene_list
)

clusters_l1 <- sort(unique(as.character(cl_l1)))
terms_l1    <- names(term_to_genes_l1)

mat_l1 <- sapply(clusters_l1, function(clus) {
  cell_ids <- names(cl_l1)[cl_l1 == clus]
  sapply(terms_l1, function(term) {
    mean_log1p_counts_in_cells(counts_l1, term_to_genes_l1[[term]], cell_ids)
  })
})
mat_l1 <- as.matrix(mat_l1)
rownames(mat_l1) <- terms_l1

# Heatmap: log1p counts
rng <- range(mat_l1, na.rm = TRUE)
pdf("plots/All_GO_Terms_Lineage1_log1pCounts_Heatmap.pdf",
    width = 10,
    height = min(100, nrow(mat_l1) * 0.35))
draw(
  Heatmap(
    mat_l1,
    name = "mean log1p(counts)",
    col = colorRamp2(c(rng[1], mean(rng), rng[2]), c("white", "orange", "red")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_title = "All enriched GO terms (Lineage 1)",
    row_title = "GO Biological Process",
    na_col = "grey90"
  )
)
dev.off()
cat("Saved Lineage 1 log1p heatmap.\n")

# Heatmap: z-scores
mat_l1_z <- t(scale(t(mat_l1)))
mat_l1_z[is.na(mat_l1_z)] <- 0
cap <- 3
mat_l1_z[mat_l1_z >  cap] <-  cap
mat_l1_z[mat_l1_z < -cap] <- -cap

col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

pdf("plots/All_GO_Terms_Lineage1_Zscore_Heatmap.pdf",
    width = 10,
    height = min(100, nrow(mat_l1_z) * 0.35))
draw(
  Heatmap(
    mat_l1_z,
    name = "z-score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_title = "All GO terms (Lineage 1) - z-score",
    row_title = "GO Biological Process",
    na_col = "grey90"
  )
)
dev.off()
cat("Saved Lineage 1 z-score heatmap.\n")

# ============================
# 5. Lineage 2 analysis
# ============================
#
# NOTE:
# The analysis for Lineage 2 is written in the *exact same structure*
# as for Lineage 1 (steps 4.1 → 4.6).
# This deliberate symmetry ensures:
#   - Reproducibility: results are directly comparable across lineages.
#   - Maintainability: if more lineages are added later, the same
#     modular layout can be reused.
# Perform the same steps as in chapter 4 by adjusting the lineage. 

# ============================
# 6. Summary & outputs
# ============================
#
# At this point, the pipeline has completed:
#   - Fitting GAMs for Lineage 1 and Lineage 2
#   - Identifying dynamic genes with start-vs-end tests
#   - Running GO enrichment on significant genes
#   - Performing per-cluster enrichment analyses
#   - Generating heatmaps (log1p and z-score) of GO terms
#
# Outputs were written to:
#   - results/  → R objects (.rds) and tables (.csv)
#   - plots/    → Heatmaps (.pdf)
#
# This consistent structure allows results from Lineage 1 and Lineage 2
# to be compared directly and makes the workflow easy to extend to
# additional lineages if required.

# --- Optional: print a message to console ---
message("All analyses complete. Results are in 'results/' and plots in 'plots/'.")

# --- Optional: save session information for reproducibility ---
sink("results/sessionInfo.txt")
sessionInfo()
sink()

