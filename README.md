# Rat PIL snRNA-seq pipeline

Single-nucleus RNA-seq of the rat posterior intralaminar thalamic nucleus (PIL)
across three samples: `R1_Cont`, `R2_Affi`, `R3_Aggr`.

## Layout

```
config/
  params.yml                    every threshold and path; no numbers in the code
  markers_used_in_figures.csv   marker panel for the figures
  neuropeptide_genes.txt        neuropeptide gene list
  cluster_rename.csv            cluster_id -> cell-type name
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
run_all.R                       runs steps 01-07
```

Outputs go to `data/processed/` (objects and tables), `data/processed/qc/`
(diagnostics) and `figures/`. All three are created on first run.

## Running it

```r
setwd("E:/Rat_PIL/github")
source("run_all.R")            # or: Rscript run_all.R
```

Single steps:

```r
source("scripts/00_setup.R")
source("scripts/05_annotation_and_renaming.R")
```

## No checkpoints

Every step runs from the `.rds` written by the previous one. Earlier versions
cached the slow parts of steps 04 and 06; that cache is gone, because a cached
file that outlives a change upstream silently overrides the fresh computation.
Both failure modes actually occurred during development: a stale doublet call
and a stale neuron subset each propagated a wrong clustering into every figure
downstream, visible only as biologically impossible plots.

The slow steps are 03 (integration), 04 (scDblFinder) and 06 (subset
integration). All are seeded, so re-running them reproduces the same result.

## Cluster IDs are renumbered twice — this matters

Clustering happens four times -- in steps 02, 04, 05 and 06 -- and each pass
creates new cluster IDs. Step 02's are provisional (they exist only to give the
ambient-removal step its cluster structure) and are stored as `clusters_prelim`.
The two config files answer to different numbering:

| file | keyed to | evidence table with matching IDs |
|---|---|---|
| `config/clusters_to_drop.csv` | step-04 clusters | `qc/04_cluster_diagnostics.csv` |
| `config/cluster_rename.csv` | step-05 clusters (after exclusion + reclustering) | `qc/05a_cluster_annotation_evidence.csv` |

A rename table from an earlier run cannot be reused: cluster 3 today is not
cluster 3 from before. Whenever step 04 or the reclustering in step 05 is
re-run, the names come from `config/cluster_rename_suggested.csv`, which step 05
regenerates against the current clusters on every run. Write
`config/cluster_rename.csv` only when you want to override those names by hand;
it takes precedence, and step 05 refuses to apply it if it no longer matches.

In step 06 the neuronal subclusters are renumbered once more, to close the gaps
left by the composition filter. The original IDs are kept in the
`cluster_original` metadata column and in
`data/processed/06_subneuron_cluster_id_map.csv`.

## Decisions that belong in the Methods

**Ambient RNA.** `params$ambient$method` is `auto`: SoupX when raw droplet
matrices are available for every sample, otherwise decontX. The method actually
used and the estimated contamination per sample are written to
`qc/03_ambient_method.csv`. Do not write "SoupX" in the Methods without checking
that file.

**Cluster exclusion.** Nothing is excluded unless `config/clusters_to_drop.csv`
exists, and every entry needs a written reason. Step 05 cross-checks each entry
against the step-04 diagnostics, warns when a cluster is dropped that nothing
flags, and writes `qc/05_dropped_clusters_with_evidence.csv` — the decision plus
its supporting numbers, ready for the supplement.

**Neuronal composition filter.** Step 06 removes subclusters where any sample
contributes less than `subneurons$composition_filter$min_sample_pct` of the
size-adjusted composition. Size-adjusted means each sample's contribution is
first divided by its own neuron total, so 33.3% each is proportional parity;
without this the raw split is confounded by differing sample sizes.

With one animal per condition, sample identity is confounded with experimental
condition, so an excluded cluster cannot be attributed to technical variation
rather than to a condition-restricted population. Step 06 warns about this,
writes the excluded clusters' markers to
`data/processed/06_dropped_clusters_markers.csv`, and saves the unfiltered
object as `data/processed/06_neurons_unfiltered.rds` so a sensitivity analysis
can be run by pointing step 07 at that file instead.

**Marker p-values.** `FindAllMarkers` tests nuclei, not animals. Thousands of
nuclei from one animal are not independent observations, so the p-values are
anticonservative. This is the accepted convention for describing cluster
markers, but it does not support statements about differences between
conditions.

**Mitochondrial genes.** Rat references name them three different ways
(`Mt-nd1`, `AY172581.x`, or bare `ND1`/`COX1`/`CYTB`). `00_setup.R` matches all
three and raises an error if none match, so the QC filter can never silently
pass everything.

## Requirements

CRAN: Seurat (v5), SeuratObject, SoupX, Matrix, data.table, dplyr, tidyr,
ggplot2, patchwork, future, readr, yaml.
Bioconductor: scDblFinder, SingleCellExperiment, BiocParallel, ComplexHeatmap,
CellChat, NMF.
Optional: scCustomize (`analyses/clustered_dotplot.R`), celda (decontX fallback),
ggalluvial (CellChat river plots).

`sessionInfo()` is written next to the outputs at every step.

## Reproducibility

Stochastic steps and how they are pinned:

| step | source of randomness | pinned by |
|---|---|---|
| 02, 03 | CCA/irlba inside `FindIntegrationAnchors` | `set.seed(42)` in `00_setup.R` |
| 03 | decontX (only on the fallback path) | `seed = 42` in the call |
| 04 | scDblFinder's simulated doublets | `doublets$seed` + `BPPARAM` RNGseed |
| 04, 05, 06 | PCA, UMAP, Louvain | Seurat defaults (`seed.use = 42`, `random.seed = 0`) plus explicit seeds in `params.yml` |
| 07 | CellChat permutation test, NMF | `cellchat$seed` |

Two caveats worth stating in a data-availability section rather than glossing
over:

* `set.seed()` at the top of `00_setup.R` only determines the integration steps
  as long as the amount of RNG consumed before them is identical. Running the
  scripts in a different order, or interactively with extra commands in between,
  can shift it. Run each step in a fresh session, or via `run_all.R`.
* Package versions matter. `sessionInfo()` is written next to the outputs at
  every step; consider `renv::init()` and depositing `renv.lock` with the code.

For a reader to reproduce the figures exactly, deposit: the raw count matrices,
this repository including `config/`, and the processed objects from
`data/processed/`. The processed objects are what make the result verifiable
even if a package version changes.
