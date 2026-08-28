# Inhibition of Aggression by Social Touch — Rat PIL snRNA-seq

R code for the single-nucleus RNA sequencing analysis of the rat posterior
intralaminar thalamic nucleus (PIL): preprocessing and QC, ambient RNA removal,
doublet detection, integration and clustering, cell-type annotation, neuronal
subclustering, marker visualisation, and CellChat-based cell–cell communication
inference. Publication-quality figures are exported in both TIFF and PDF format.

## Study

**Inhibition of Aggression by Social Touch**

Tamás Láng, Botond B. Drahos, Dávid Keller, Fanni Dóra, Ingrid Csordás,
Vivien Szendi, Gina Puska, Valery Grinevich, and Árpád Dobolyi

*Abstract:* TBD (manuscript in preparation).

**One-sentence summary:** Social touch–activated PIL neurons project to the MPOA
and acutely inhibit aggression; chemogenetic/optogenetic manipulation confirms
causality.

Twelve adult male Wistar rats were assigned to three experimental groups (n = 4
animals per group). For snRNA-seq, nuclei from animals within each experimental group
were pooled into one sequencing library, resulting in three datasets: R1_Cont,
R2_Affi, and R3_Aggr.

## Data

Raw FASTQ files and 10x feature–barcode matrices will be deposited at NCBI GEO:
**GSEXXXXXX** (accession to be added).

Until then, place the 10x outputs locally and point `config/params.yml` at them.
Data directories are listed in `.gitignore` and are never uploaded.

The pipeline expects the filtered matrices in one flat folder with a per-sample
prefix (`R1_Cont_matrix.mtx.gz`, `R1_Cont_barcodes.tsv.gz`,
`R1_Cont_features.tsv.gz`, …), and each sample's cellranger
`raw_feature_bc_matrix` folder in the standard layout. The raw matrices are
required by SoupX; without them step 03 falls back to decontX.

## Layout

```
config/
  params.yml                    every threshold and path; no numbers in the code
  markers_used_in_figures.csv   marker panel for the figures
  neuropeptide_genes.txt        neuropeptide gene list
  cluster_rename.csv            cluster_id -> cell-type name (optional override)
  clusters_to_drop.csv          cluster_id -> written reason for exclusion
scripts/
  00_setup.R                    packages, options, shared helpers (sourced by all)
  01_qc_and_merge.R             per-sample QC, merge
  02_integration_SCT.R          preliminary SCT integration and clustering
  03_ambient_removal.R          SoupX (or decontX) + re-integration
  04_doublets_and_recluster.R   scDblFinder, removal, reclustering, diagnostics
  05_annotation_and_renaming.R  exclusion, reclustering, lineage evidence, naming
  06_neuron_subset.R            neuronal subset, composition filter, figures
  07_cellchat.R                 secreted-signalling analysis
  analyses/
    clustered_dotplot.R         scCustomize clustered dotplot (optional figure)
    oligo_maturation.R          oligodendrocyte maturation staging (optional)
utils/
  plotting.R                    plot_theme() and write_plot()
```

Outputs go to `data/processed/` (objects and tables), `data/processed/qc/`
(diagnostics) and `figures/`. All three are created on first run.

## Running it

Before the first run, edit `matrix_dir` and the per-sample `raw_dir` entries at
the top of `config/params.yml` so they point at your data. Everything else in
that file is an analysis parameter, not a path.

Set the working directory to the project root — the folder holding `scripts/`
and `config/` — and run the steps in order:

```r
source("scripts/01_qc_and_merge.R")
source("scripts/02_integration_SCT.R")
source("scripts/03_ambient_removal.R")
source("scripts/04_doublets_and_recluster.R")
source("scripts/05_annotation_and_renaming.R")
source("scripts/06_neuron_subset.R")
source("scripts/07_cellchat.R")
```

Each script sources `scripts/00_setup.R` itself, so it does not need to be run
separately. The two scripts under `scripts/analyses/` are optional and can be
run at any point after step 05.

The steps are meant to be run one at a time, with the diagnostics inspected in
between — in particular, `qc/04_cluster_diagnostics.csv` after step 04, which is
what `config/clusters_to_drop.csv` refers to.

Every step reads the `.rds` written by the previous one; nothing is cached, so a
re-run always recomputes from its input. The slow steps are 03 and 06
(integration) and 04 (scDblFinder).

## Cluster identifiers change between steps

Clustering is performed four times — in steps 02, 04, 05 and 06 — and each pass
assigns new cluster IDs. Step 02's are provisional, exist only to give the
ambient-removal step its cluster structure, and are stored as `clusters_prelim`.

The two editable config files therefore refer to different numbering:

| file | keyed to | evidence table with matching IDs |
|---|---|---|
| `config/clusters_to_drop.csv` | step-04 clusters | `qc/04_cluster_diagnostics.csv` |
| `config/cluster_rename.csv` | step-05 clusters (after exclusion and reclustering) | `qc/05a_cluster_annotation_evidence.csv` |

Cell-type names are assigned from `config/cluster_rename_suggested.csv`, which
step 05 regenerates from the lineage marker panels against the clusters that
exist in that run. Write `config/cluster_rename.csv` only to override those
names by hand; it takes precedence, and step 05 stops rather than apply a table
that no longer matches the current clustering.

In step 06 the neuronal subclusters are renumbered once more, to close the gaps
left by the composition filter. The original IDs are kept in the
`cluster_original` metadata column and in
`data/processed/06_subneuron_cluster_id_map.csv`.

## Analysis decisions that belong in the Methods

**Ambient RNA.** `params$ambient$method` is `auto`: SoupX when raw droplet
matrices are available for every sample, decontX otherwise. The method actually
used and the estimated contamination per sample are recorded in
`qc/03_ambient_method.csv`.

**Cluster exclusion.** Nothing is excluded unless `config/clusters_to_drop.csv`
exists, and every entry requires a written reason. Step 05 cross-checks each
entry against the step-04 diagnostics, warns when a cluster is dropped that the
diagnostics do not flag, and writes `qc/05_dropped_clusters_with_evidence.csv` —
the decision together with its supporting numbers.

**Neuronal composition filter.** Step 06 removes subclusters in which any sample
contributes less than `subneurons$composition_filter$min_sample_pct` of the
size-adjusted composition. Size-adjusted means each sample's contribution is
first divided by its own neuron total, so 33.3% per sample represents
proportional parity; without this the raw split is confounded by differing
sample sizes.

With one pooled sequencing library per condition, library identity is completely
confounded with experimental condition. Although each library represents nuclei pooled
from multiple animals, animal-level biological replication is not retained in the
sequencing data. Therefore, an excluded cluster cannot be attributed confidently to
technical variation rather than to a condition-restricted population. Step 06 records
the excluded clusters' markers in `data/processed/06_dropped_clusters_markers.csv` and
saves the unfiltered object as `data/processed/06_neurons_unfiltered.rds`, so a
sensitivity analysis can be run by pointing step 07 at that file instead.

**Marker p-values.** `FindAllMarkers` tests nuclei rather than independent biological
replicates. Nuclei originating from the same pooled sequencing library are not
independent experimental units, so marker-test p-values should be interpreted only as
descriptive evidence for cluster identity. They must not be used to infer differential
expression between experimental conditions.
This is the accepted convention for describing cluster markers, but it does not
support claims about differences between conditions.

**Mitochondrial genes.** Rat references name these in three different ways
(`Mt-nd1`, `AY172581.x`, or bare `ND1`/`COX1`/`CYTB`). `00_setup.R` matches all
three and raises an error if none match, so the QC filter cannot silently pass
every cell.

**Cell–cell communication.** CellChat is used to infer putative ligand–receptor
signalling networks from gene-expression patterns. Because there is one pooled
sequencing library per experimental condition, comparisons between conditions are
descriptive and hypothesis-generating; they do not constitute replicate-level
statistical evidence for condition-dependent changes in signalling.

## Reproducibility

Stochastic steps and how they are pinned:

| step | source of randomness | pinned by |
|---|---|---|
| 02, 03 | CCA/irlba inside `FindIntegrationAnchors` | `set.seed(42)` in `00_setup.R` |
| 03 | decontX (fallback path only) | `seed = 42` in the call |
| 04 | scDblFinder's simulated doublets | `doublets$seed` and the `BPPARAM` RNG seed |
| 04, 05, 06 | PCA, UMAP, Louvain | Seurat defaults plus the explicit seeds in `params.yml` |
| 07 | CellChat permutation test, NMF | `cellchat$seed` |

Two caveats:

* `set.seed()` in `00_setup.R` determines the integration steps only as long as
  the amount of RNG consumed before them is identical. Run each step in a fresh
  R session rather than interleaving other commands.
* Package versions matter. `sessionInfo()` is written next to the outputs at
  every step; `renv.lock` records the exact versions used.

The analysis can be reproduced from the raw count matrices deposited in GEO together
with this repository and the package versions recorded in renv.lock. Processed .rds
objects are additionally retained to facilitate exact verification of the published
results and figures across future software-version changes.

## Requirements

CRAN:
Seurat (v5), SeuratObject, SoupX, Matrix, data.table, dplyr, tidyr,
ggplot2, patchwork, future, readr, yaml, NMF.

Bioconductor:
scDblFinder, SingleCellExperiment, BiocParallel, ComplexHeatmap.

GitHub:
CellChat (`jinworks/CellChat`).

Optional:
scCustomize (`analyses/clustered_dotplot.R`),
celda (Bioconductor; decontX fallback),
ggalluvial (CellChat river plots).
