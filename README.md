# Inhibition of Aggression by Social Touch — Rat PIL snRNA-seq

This repository contains R code for subclustering and annotation of neurons from rat posterior intralaminar thalamic nucleus (PIL) single-nucleus RNA sequencing (snRNA-seq) data. The analysis covers preprocessing/QC, integration and clustering, neuron subclustering, marker visualization, and CellChat-based cell–cell communication summaries. Figures are saved in TIFF and PDF formats to comply with journal requirements.

## Study

**Inhibition of Aggression by Social Touch**

Tamás Láng, Botond B. Drahos, Dávid Keller, Fanni Dóra, Ingrid Csordás, Vivien Szendi, Gina Puska, Valery Grinevich, and Árpád Dobolyi

*Abstract:* TBD (manuscript in preparation).

**One-sentence summary:** Social touch–activated PIL neurons project to the MPOA and acutely inhibit aggression; chemogenetic/optogenetic manipulation confirms causality.

## Raw Data

- **Raw FASTQ** files and **10x feature–barcode matrices** for this study will be deposited to **NCBI GEO**: **GSEXXXXXX** (to be provided).


Until accessions are assigned, place local 10x outputs under `data/raw/<sample>/outs/` (see below). These paths are ignored by Git and will not be uploaded.

## System requirements

- **R**: 4.4.3 (tested), RStudio recommended  
- **OS**: Windows 11 and Ubuntu 24.04 (WSL) (tested)  
- **Packages**:  
  - CRAN: Seurat (v5.3.1), SeuratDisk, SoupX (v1.6.2), tidyverse, dplyr, tidyr, Matrix, data.table, ggplot2, patchwork, cowplot, future, readr, yaml  
  - Bioconductor: scDblFinder (v1.18.0), SingleCellExperiment, scater  
  - Other: CellChat (v1.6.1), NMF, ggalluvial  
- **Hardware**: ≥16 GB RAM (32 GB recommended for full dataset)

## Installation

```r
# CRAN
install.packages(c("Seurat","SeuratDisk","SoupX","tidyverse","dplyr","tidyr",
                   "Matrix","data.table","ggplot2","patchwork","cowplot","future","readr","yaml"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("scDblFinder","SingleCellExperiment","scater"))

# CellChat and dependencies
install.packages(c("CellChat","NMF","ggalluvial"))
```

## Instructions for use

### Clone this repository

```bash
git clone https://github.com/DobolyiLab/Rat_PIL_snRNAseq.git
cd Rat_PIL_snRNAseq
```

### Place input data

- Download from GEO (GSE305279; reviewer token required until publication),  
  or use local 10x outputs under `data/raw/<sample>/outs/`.  
- Sample metadata and paths must be defined in `config/params.yml`.

### Run scripts in order

1. `scripts/01_qc.R` – QC, Seurat object creation, violin plots  
2. `scripts/02_integration.R` – SCTransform normalization, integration, clustering, marker dotplots  
3. `scripts/03_soupx.R` – Ambient RNA correction with SoupX, reintegration  
4. `scripts/04_doublet.R` – DoubletFinder classification and singlet filtering  
5. `scripts/05_annotation.R` – Cluster annotation, dotplots, renaming via CSV mapping  
6. `scripts/06_subneurons_analysis.R` – Neuron subset analysis, subclustering, Fos/Calb1, neurotransmitter and neuropeptide markers  
7. `scripts/07_cellchat.R` – CellChat analysis (communication networks, scatter, heatmap, river, dot plots)

### Outputs

- **Processed Seurat objects** → `data/processed/`  
- **Figures** → `figures/` (TIFF, PDF)  
- **Tables** → `data/processed/` (CSV)  
- **Session info logs** → `data/processed/`

## Demo

You can test the workflow on a small synthetic dataset:

```r
library(Seurat)
set.seed(1)
mat <- matrix(rpois(2000*300, lambda=1), nrow=2000, ncol=300,
              dimnames=list(paste0("Gene",1:2000), paste0("Cell",1:300)))
obj <- CreateSeuratObject(mat)
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA() |> RunUMAP(dims=1:10)
saveRDS(obj, file="data/processed/demo_seurat.rds")
```

## Data & Code availability

- **Data**: Raw sequencing data deposited in GEO: **GSEXXXXXX** (reviewer token provided; public upon publication).  
- **Code**: This GitHub repository: https://github.com/DobolyiLab/Rat_PIL_snRNAseq  
- **License**: MIT (see LICENSE).

## Reproducibility

Record R environment with:

```r
sink("data/processed/sessionInfo.txt"); sessionInfo(); sink()
```

- Key parameters (filters, PCs, clustering resolutions) are documented in script headers.

