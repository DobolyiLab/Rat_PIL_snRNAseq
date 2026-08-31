# Inhibition of Aggression by Social Touch — Rat PIL snRNA-seq

Single-nucleus RNA-seq analysis of the rat posterior intralaminar thalamic
nucleus (PIL): preprocessing and QC, ambient RNA removal, doublet detection,
integration and clustering, cell-type labeling, and neuronal subclustering.

## Study

**Inhibition of Aggression by Social Touch**

Tamás Láng, Botond B. Drahos, Dávid Keller, Fanni Dóra, Ingrid Csordás,
Vivien Szendi, Gina Puska, Valery Grinevich, and Árpád Dobolyi

**One-sentence summary:** Social touch–activated PIL neurons project to the MPOA
and acutely inhibit aggression; chemogenetic/optogenetic manipulation confirms
causality.

Twelve adult male Wistar rats were assigned to three experimental groups
(n = 4 per group). Nuclei from animals within each group were pooled into one
sequencing library per group, giving three datasets: `R1_Cont`, `R2_Affi`, `R3_Aggr`.

## Data

Raw FASTQ files and 10x feature–barcode matrices are deposited at NCBI GEO:
**GSE305279**.

## Contents

```
scripts/
  PIL_snRNAseq_pipeline.R      -- raw data through neuronal subclusters and
                                  neuropeptide markers
  celltype_labeling.R          -- names the clusters (run after the pipeline)
  Fos_Calb1_quantification.R   -- Fos/Calb1 detection per subcluster and sample
config/
  paths.example.R              -- copy to paths.R and edit
  neuropeptide_genes.txt       -- one gene symbol per line
  markers_used_in_figures.csv  -- one "gene" column
utils/
  plotting.R             -- plot_theme() and write_plot()
```

## You need to add

- `config/paths.R` -- copy from `paths.example.R`, fill in your data paths
- `config/cluster_rename.csv` -- `celltype_labeling.R` writes a template on
  its first run; fill it in and re-run

## Running

```r
source("scripts/PIL_snRNAseq_pipeline.R")
source("scripts/celltype_labeling.R")
source("scripts/Fos_Calb1_quantification.R")
```

Every stochastic step is seeded (`set.seed(42)`), so re-running the pipeline
reproduces the same clusters.

## Requirements

CRAN: Seurat (v5), SeuratObject, SoupX, Matrix, data.table, dplyr, tidyr,
ggplot2, future.

Bioconductor: scDblFinder, SingleCellExperiment, BiocParallel, MAST.
