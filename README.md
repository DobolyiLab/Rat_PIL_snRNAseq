# Inhibition of Aggression by Social Touch — Rat PIL snRNA-seq

R code for the single-nucleus RNA sequencing analysis of the rat posterior intralaminar thalamic nucleus (PIL): preprocessing and quality control, ambient RNA removal, doublet detection, integration and clustering, cell-type annotation, neuronal subclustering, marker visualisation, and CellChat-based cell–cell communication inference.

Publication-quality figures are exported in both TIFF and PDF format.

## Study

**Inhibition of Aggression by Social Touch**

Tamás Láng, Botond B. Drahos, Dávid Keller, Fanni Dóra, Ingrid Csordás, Vivien Szendi, Gina Puska, Valery Grinevich, and Árpád Dobolyi

*Abstract:* TBD (manuscript in preparation).

**One-sentence summary:** Social touch–activated PIL neurons project to the MPOA and acutely inhibit aggression; chemogenetic and optogenetic manipulation confirms causality.

Twelve adult male Wistar rats were assigned to three experimental groups (n = 4 animals per group). For snRNA-seq, nuclei from animals within each experimental group were pooled into one sequencing library, resulting in three datasets:

- `R1_Cont` — control
- `R2_Affi` — affiliative/social-touch condition
- `R3_Aggr` — aggressive condition

Because there is one pooled sequencing library per experimental condition, sequencing-library identity is completely confounded with experimental condition. Consequently, comparisons between conditions are interpreted descriptively rather than as replicate-level statistical tests of condition effects.

## Data

Raw FASTQ files and 10x feature–barcode matrices will be deposited at NCBI GEO:

**GSE305279**

Until then, place the 10x outputs locally and point `config/params.yml` to them.

The pipeline expects the filtered matrices in one flat folder with a per-sample prefix, for example:

```text
R1_Cont_matrix.mtx.gz
R1_Cont_barcodes.tsv.gz
R1_Cont_features.tsv.gz

R2_Affi_matrix.mtx.gz
R2_Affi_barcodes.tsv.gz
R2_Affi_features.tsv.gz

R3_Aggr_matrix.mtx.gz
R3_Aggr_barcodes.tsv.gz
R3_Aggr_features.tsv.gz
```

Each sample's Cell Ranger `raw_feature_bc_matrix` folder should also be available in the standard 10x layout.

The raw droplet matrices are required by SoupX. If raw matrices are unavailable, step 03 falls back to decontX.

## Repository layout

```text
config/
  params.yml                       analysis parameters and input paths
  markers_used_in_figures.csv      marker panel used in figures
  neuropeptide_genes.txt           neuropeptide gene list
  clusters_to_drop.csv             optional manual cluster exclusions
  cluster_rename.csv               optional manual annotation overrides

scripts/
  00_setup.R                       packages, options, shared helpers
  01_qc_and_merge.R                per-sample QC and merge
  02_integration_SCT.R             preliminary SCT integration and clustering
  03_ambient_removal.R             SoupX/decontX and re-integration
  04_doublets_and_recluster.R      scDblFinder, doublet removal, reclustering
  05_annotation_and_renaming.R     cluster exclusion, annotation and naming
  06_neuron_subset.R               neuronal subset, composition filter, figures
  07_cellchat.R                    secreted-signalling analysis

scripts/analyses/
  clustered_dotplot.R              scCustomize clustered dotplot
  oligo_maturation.R               oligodendrocyte maturation staging

utils/
  plotting.R                       plot_theme() and write_plot()
```

Outputs are written to:

```text
data/processed/
data/processed/qc/
figures/
```

These directories are created automatically when required.

## Running the pipeline

Before the first run, edit `matrix_dir` and the per-sample `raw_dir` entries at the top of `config/params.yml` so that they point to the local data.

All remaining entries in `params.yml` are analysis parameters rather than paths.

Set the working directory to the project root — the directory containing `scripts/` and `config/` — and run the steps in order:

```r
source("scripts/01_qc_and_merge.R")
source("scripts/02_integration_SCT.R")
source("scripts/03_ambient_removal.R")
source("scripts/04_doublets_and_recluster.R")
source("scripts/05_annotation_and_renaming.R")
source("scripts/06_neuron_subset.R")
source("scripts/07_cellchat.R")
```

Each script sources `scripts/00_setup.R` automatically, so `00_setup.R` does not need to be run separately.

The scripts under `scripts/analyses/` are optional and can be run after the required upstream objects have been generated.

The main pipeline steps are intended to be run one at a time, with diagnostics inspected between steps. In particular, inspect:

```text
data/processed/qc/04_cluster_diagnostics.csv
```

after step 04 before defining any manual exclusions in:

```text
config/clusters_to_drop.csv
```

Each step reads the `.rds` object written by the previous step and produces a new output object. Re-running a step therefore recomputes the analysis from its saved input rather than relying on an in-memory object from an earlier R session.

The computationally intensive steps are primarily ambient-RNA correction/re-integration, doublet detection, and neuronal re-integration/subclustering.

## Cluster identifiers change between steps

Clustering is performed multiple times during the pipeline, and each clustering pass can assign new numerical cluster identifiers.

Clustering occurs in steps 02, 04, 05, and 06.

Step 02 clusters are provisional. They are generated primarily to provide cluster structure for ambient-RNA estimation and are stored as:

```text
clusters_prelim
```

The editable cluster configuration files therefore refer to different clustering stages:

| File | Refers to | Evidence table with matching IDs |
|---|---|---|
| `config/clusters_to_drop.csv` | step-04 clusters | `qc/04_cluster_diagnostics.csv` |
| `config/cluster_rename.csv` | step-05 clusters after exclusion and reclustering | `qc/05a_cluster_annotation_evidence.csv` |

Cell-type names are initially proposed in:

```text
config/cluster_rename_suggested.csv
```

This file is regenerated by step 05 using lineage-marker evidence from the clusters present in that particular run.

If manual annotation overrides are required, create:

```text
config/cluster_rename.csv
```

Manual entries take precedence over suggested names.

Step 05 checks that the supplied cluster IDs match the clustering produced by the current run and stops rather than silently applying a renaming table to incompatible cluster IDs.

In step 06, neuronal subclusters are renumbered after filtering to produce consecutive cluster identifiers.

The corresponding pre-renumbering IDs are retained in the metadata column:

```text
cluster_original
```

and the mapping between original and final neuronal cluster IDs is written to:

```text
data/processed/06_subneuron_cluster_id_map.csv
```

## Analysis decisions relevant to the Methods

### Ambient RNA removal

`params$ambient$method` is set to `auto`.

SoupX is used when raw droplet matrices are available for all samples. If the required raw matrices are unavailable, the pipeline falls back to decontX.

The method actually used and the estimated contamination for each sample are recorded in:

```text
data/processed/qc/03_ambient_method.csv
```

### Cluster exclusion

No cluster is excluded unless:

```text
config/clusters_to_drop.csv
```

exists.

Each exclusion must include a written reason.

Step 05 cross-checks requested exclusions against the diagnostics generated in step 04 and warns if a cluster is removed despite not being flagged by the quantitative diagnostics.

The final exclusion decisions together with the supporting evidence are written to:

```text
data/processed/qc/05_dropped_clusters_with_evidence.csv
```

This preserves the rationale for manual exclusions and makes the decision process auditable.

### Cross-sample neuronal composition filter

Step 06 applies an optional composition filter to neuronal subclusters.

A subcluster is removed if any sample contributes less than:

```text
subneurons$composition_filter$min_sample_pct
```

of the size-adjusted sample composition.

Sample contributions are first normalised by the total number of neuronal nuclei recovered from that sample. Under proportional parity, each of the three samples would therefore contribute approximately 33.3% of a subcluster.

This size adjustment prevents differences in the total number of recovered nuclei per sample from directly determining the composition metric.

However, there is one pooled sequencing library per experimental condition. Library identity is therefore completely confounded with condition, and condition-restricted populations cannot be distinguished statistically from library-specific technical effects.

For this reason, exclusion based on sample composition must be interpreted cautiously.

Step 06 records markers for the excluded neuronal subclusters in:

```text
data/processed/06_dropped_clusters_markers.csv
```

and saves the complete neuronal object before application of the composition filter as:

```text
data/processed/06_neurons_unfiltered.rds
```

This allows sensitivity analyses to be performed with the unfiltered neuronal population if required.

### Marker p-values

`FindAllMarkers` performs tests at the nucleus level rather than at the level of independent biological replicates.

Nuclei originating from the same pooled sequencing library are not independent experimental units. Marker-test p-values should therefore be interpreted as descriptive evidence supporting cluster identity rather than as biological-replicate-level inferential statistics.

These p-values do not support statistical claims about differential expression between experimental conditions.

### Cell–cell communication

CellChat is used to infer putative ligand–receptor signalling networks from gene-expression patterns.

Because there is one pooled sequencing library per experimental condition, comparisons of inferred signalling between `R1_Cont`, `R2_Affi`, and `R3_Aggr` are descriptive and hypothesis-generating.

They do not constitute replicate-level statistical evidence for condition-dependent changes in cell–cell communication.

### Mitochondrial genes

Rat transcriptome references may represent mitochondrial genes using several naming conventions, including:

```text
Mt-nd1
AY172581.x
ND1
COX1
CYTB
```

`scripts/00_setup.R` tests the supported mitochondrial-gene naming conventions and raises an error if none are detected.

This prevents the mitochondrial QC filter from silently assigning zero mitochondrial content because of an incompatible gene-name convention.

## Reproducibility

Random-number generation is controlled for stochastic analysis steps.

| Step | Source of randomness | Controlled by |
|---|---|---|
| 02, 03 | CCA/IRLBA during `FindIntegrationAnchors` | `set.seed(42)` in `00_setup.R` |
| 03 | decontX fallback | `seed = 42` |
| 04 | scDblFinder simulated doublets | `doublets$seed` and BiocParallel RNG seed |
| 04, 05, 06 | PCA, UMAP and clustering | explicit seeds in `params.yml` where applicable |
| 07 | CellChat permutation procedures / NMF | `cellchat$seed` |

Two reproducibility considerations are important.

First, `set.seed()` guarantees identical stochastic behaviour only when the sequence of random-number-consuming operations is unchanged. For this reason, the main pipeline steps should preferably be run in fresh R sessions rather than interleaved with unrelated commands.

Second, package versions can alter analysis results even when random seeds are fixed.

`sessionInfo()` is therefore written alongside outputs, and:

```text
renv.lock
```

records the package environment used for the analysis.

The complete analysis can be reproduced from the raw count matrices deposited in GEO together with this repository and the package versions recorded in `renv.lock`.

Processed `.rds` objects are additionally retained to facilitate exact verification of the published analysis and figures if future package versions change numerical results or implementation details.

## Requirements

### CRAN

- Seurat v5
- SeuratObject
- SoupX
- Matrix
- data.table
- dplyr
- tidyr
- ggplot2
- patchwork
- future
- readr
- yaml
- NMF

### Bioconductor

- scDblFinder
- SingleCellExperiment
- BiocParallel
- ComplexHeatmap

### GitHub

- CellChat (`jinworks/CellChat`)

### Optional

- scCustomize — used by `scripts/analyses/clustered_dotplot.R`
- celda — provides decontX for the ambient-RNA fallback
- ggalluvial — used for optional CellChat river/alluvial plots

## Citation

Citation information will be added after publication of the associated manuscript.

## Data availability

Raw sequencing data and count matrices will be available through NCBI GEO under accession:

**GSE305279**

Analysis code, configuration files, and documentation required to reproduce the computational workflow are provided in this repository.
