# slingshot-tradeseq-workflow
This repository contains an R-based workflow for pseudotime and trajectory inference on single-cell RNA-seq data. It starts from a processed Seurat object and uses Slingshot for lineage detection and pseudotime estimation, followed by tradeSeq for differential expression along trajectories.
## Overview
This workflow takes a processed Seurat object (after QC, normalization, clustering)
and performs:
1. Conversion to `SingleCellExperiment`
2. Dimensionality reduction (PCA and DiffusionMap)
3. Trajectory inference with Slingshot
4. Differential expression along pseudotime with tradeSeq

The repo is dataset-agnostic: provide any processed Seurat object.

## Example Visualizations
<img width="300" height="500" alt="image" src="https://github.com/user-attachments/assets/8f4ff7c9-e9e9-4993-9fb1-deffb67be1cb" />
<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/e7aeee75-7b50-421f-947b-cb984967a56e" />




## Requirements
R ≥ 4.2 and the following R packages:
`slingshot`, `SingleCellExperiment`, `destiny`, `tradeSeq`, `Seurat`, `mclust`, `ggplot2`, `ggthemes`, `RColorBrewer`

## Usage
```bash
# Run trajectory inference
Rscript scripts/01_slingshot_pseudotime.R \
       --input path/to/seurat_object.rds

# Run tradeSeq differential expression
Rscript scripts/02_tradeseq_differential.R \
       --input results/sce_after_slingshot.rds
```
## Project Structure
```
scRNAseq-trajectory-analysis/
├─ README.md                 # Overview, usage instructions, citations
├─ LICENSE                   # e.g. MIT
├─ .gitignore                # Ignore large data and outputs
├─ scripts/
│   ├─ 01_slingshot_pseudotime.R   # Convert Seurat → SCE, run PCA/DiffusionMap, Slingshot
│   └─ 02_tradeseq_differential.R  # Differential expression along pseudotime with tradeSeq
└─ images/                   # Example figures for the README
    ├─ diffusion_map_example.png
    └─ tradeseq_example.png
```

