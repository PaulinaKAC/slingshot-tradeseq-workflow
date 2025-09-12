# -----------------------------------------------------------
# 01_slingshot_pseudotime.R
# Purpose:  Convert Seurat object to SCE, run PCA/DiffusionMap,
#           infer lineages & pseudotime with Slingshot
# Input :   processed Seurat .rds
# Output:   SingleCellExperiment .rds with pseudotime, lineage info
# -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(destiny)
  library(mclust)
  library(RColorBrewer)
  library(ggplot2)
})

set.seed(123)
setwd("set/your/path")  # <-- change to your folder

# --------- 1. Load data ----------
input_file <- "yourfile.rds"     # <-- replace with actual file name
seurat_obj <- readRDS(input_file)
message("Loaded Seurat object with ", ncol(seurat_obj), " cells")
sce <- as.SingleCellExperiment(seurat_obj) # Convert to SingleCellExperiment

# --------- 2. Gene filtering ----------
keep_genes <- rowSums(counts(sce) >= 3) >= 10 # Keep genes expressed in ≥10 cells with ≥3 counts
sce <- sce[keep_genes, ]
message("Genes kept after filtering: ", nrow(sce))

# --------- 3. Quantile normalization ----------
FQnorm <- function(counts) {
  rk <- apply(counts, 2, rank, ties.method = "min")   # Compute within-cell ranks for each gene
  counts.sort <- apply(counts, 2, sort)               # Sort expression values within each cell
  refdist <- apply(counts.sort, 1, median)            # Compute the median expression at each rank across all cells
  norm <- apply(rk, 2, function(r) refdist[r])        # Replace each gene’s value by the reference value of its rank 
  rownames(norm) <- rownames(counts)                  # Preserve original gene names
  norm
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)        # Apply quantile normalization to the raw count

# --------- 4. Dimensionality reduction ----------
#       We compute:
#         (a) Principal Components (PCA)
#         (b) Diffusion Map (preferred for smooth trajectories)

## PCA on log-transformed normalized counts
# Transpose so cells are rows and genes are columns
pca <- prcomp(
  t(log1p(assays(sce)$norm)),  # log1p: log(1 + x) for variance stabilization
  scale. = FALSE               # no additional scaling; data already normalized
)
# Keep first two principal components for quick plotting
rd_pca <- pca$x[, 1:2]

## Diffusion Map
library(destiny)
dm <- DiffusionMap(t(log1p(assays(sce)$norm)))  # Captures continuous cell-state transitions
rd_dm <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)      # Store first two diffusion components

## Add both low-dimensional representations to the SCE object
# reducedDims() holds matrices of dimension-reduced embeddings
reducedDims(sce) <- SimpleList(
  PCA     = rd_pca,
  DiffMap = rd_dm
)

## Quick diagnostic plots (saved to PDF)
pdf("plots/dimred_pca.pdf", width = 6, height = 6)
plot(rd_pca, col = "grey", pch = 16, asp = 1,
     main = "PCA: first two components")
dev.off()

pdf("plots/dimred_diffmap.pdf", width = 6, height = 6)
plot(rd_dm, col = "grey", pch = 16, asp = 1,
     main = "Diffusion Map: DC1 vs DC2")
dev.off()

# --------- 5. Add clustering information ----------
                
sce$Seurat_clusters <- as.character(Idents(seurat_obj)) # Copy Seurat cluster assignments
num_clusters <- length(unique(sce$Seurat_clusters))     # Determine the number of unique clusters

# Generate a color palette with one color per cluster
cluster_cols <- setNames(
  rainbow(num_clusters),
  unique(sce$Seurat_clusters)
)

# Plot Diffusion Map colored by Seurat clusters
pdf("plots/clustering_diffmap_seurat.pdf", width = 7, height = 7)
plot(
  reducedDims(sce)$DiffMap,
  col = cluster_cols[sce$Seurat_clusters],
  pch = 20,
  asp = 1,
  main = "Diffusion Map colored by Seurat clusters"
)
legend(
  "bottomright",
  legend = names(cluster_cols),
  col    = cluster_cols,
  pch    = 20,
  cex    = 0.8,
  title  = "Seurat clusters"
)
dev.off()

# --------- 6. Additional clustering ---------- 
## 6a) Gaussian Mixture Model (GMM) clustering
# Fits a mixture of multivariate Gaussians to the Diffusion Map
# and automatically selects the optimal number of clusters
# based on Bayesian Information Criterion (BIC).
## Gaussian Mixture Model (Mclust)

# Store GMM cluster assignments in the SCE object
gmm_result <- Mclust(rd_dm)                     # run the model
sce$GMM <- as.factor(gmm_result$classification) # save clusters
                
# Plot Diffusion Map colored by GMM clusters
pdf("plots/clustering_diffmap_gmm.pdf", width = 7, height = 7)
plot(
  reducedDims(sce)$DiffMap,
  col = brewer.pal(max(as.numeric(sce$GMM)), "Set1")[as.numeric(sce$GMM)],
  pch = 20,
  asp = 1,
  main = "Diffusion Map colored by GMM clusters"
)
legend(
  "bottomright",
  legend = levels(sce$GMM),
  col    = brewer.pal(max(as.numeric(sce$GMM)), "Set1"),
  pch    = 20,
  cex    = 0.8,
  title  = "GMM clusters"
)
dev.off()


# --------- 7. Save final SCE ----------
saveRDS(sce, file = output_file)
message("Saved SCE with reducedDims and clustering to: ", output_file)

# End of script

