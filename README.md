# Slingshot + tradeSeq Workflow

Trajectory inference and pseudotime-based differential expression analysis of single-cell RNA-seq data using **Slingshot** and **tradeSeq**.  
This workflow is intended to run **after Seurat pre-processing** (QC, normalization, clustering) and identifies lineage relationships and genes that change along pseudotime.

---

## 🚀 Overview

This repository provides:

* Conversion of a processed Seurat object to `SingleCellExperiment`
* Dimensionality reduction (PCA and Diffusion Map)
* Trajectory inference and pseudotime ordering with **Slingshot**
* Differential expression along lineages with **tradeSeq**
* Example plots of trajectories and gene-expression trends

The code is dataset-agnostic: supply any single-cell dataset that has been cleaned and clustered.

---

## 📁 Project Structure
```text
slingshot-tradeseq-workflow/
├─ README.md                     
├─ LICENSE                        
├─ .gitignore                      
├─ scripts/
    ├─ 01_slingshot_pseudotime.R     # Convert Seurat → SCE, PCA/DiffusionMap, Slingshot
    └─ 02_tradeseq_differential.R    # tradeSeq differential expression along pseudotime

```

## Example Visualizations
These figures are illustrative only; results will vary depending on your data.

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/8f4ff7c9-e9e9-4993-9fb1-deffb67be1cb" />
<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/e7aeee75-7b50-421f-947b-cb984967a56e" />



## Requirements
R ≥ 4.2 and the following R packages:
`slingshot`, `SingleCellExperiment`, `destiny`, `tradeSeq`, `Seurat`, `mclust`, `ggplot2`, `ggthemes`, `RColorBrewer`
Enough computing RAM/cores especially for tradeSeq and diffusion mapping

## 📚 Citations

If you use this workflow or any of its code, please cite the relevant software and methods:

* **Slingshot**  
  Street K, et al. *Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics.*  
  BMC Genomics. 2018;19:477.  
  [https://doi.org/10.1186/s12864-018-4772-0](https://doi.org/10.1186/s12864-018-4772-0)

* **tradeSeq**  
  Van den Berge K, et al. *Trajectory-based differential expression analysis for single-cell sequencing data.*  
  Nature Communications. 2020;11:1201.  
  [https://doi.org/10.1038/s41467-020-14766-3](https://doi.org/10.1038/s41467-020-14766-3)

* **Seurat**  
  Hao Y, et al. *Integrated analysis of multimodal single-cell data.*  
  Cell. 2021;184(13):3573–3587.  
  [https://doi.org/10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048)

* **destiny (Diffusion Maps)**  
  Angerer P, et al. *destiny: diffusion maps for large-scale single-cell data analysis.*  
  Bioinformatics. 2016;32(8):1241–1243.  
  [https://doi.org/10.1093/bioinformatics/btv715](https://doi.org/10.1093/bioinformatics/btv715)

* **mclust**  
  Scrucca L, et al. *mclust 5: Clustering, Classification and Density Estimation Using Gaussian Finite Mixture Models.*  
  The R Journal. 2016;8(1):289–317.  
  [https://doi.org/10.32614/RJ-2016-021](https://doi.org/10.32614/RJ-2016-021)


---

## 💡 Notes & Tips

* Clustering resolution in Seurat strongly influences inferred trajectories—experiment with different resolutions.
* Diffusion Map often captures smooth biological processes better than PCA for Slingshot.
* Ensure enough cells per cluster and along each lineage for robust tradeSeq testing.
* Keep a log of parameters and package versions for reproducibility.
* Large datasets may require substantial RAM/cores, especially for tradeSeq fits.

