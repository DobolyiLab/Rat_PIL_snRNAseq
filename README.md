# Inhibition of Aggression by Social Touch — Rat PIL snRNA-seq

R code for the single-nucleus RNA sequencing analysis of the rat posterior
intralaminar thalamic nucleus (PIL): preprocessing and QC, ambient RNA removal,
doublet detection, integration and clustering, cell-type annotation, neuronal
subclustering, marker visualisation, and CellChat-based cell–cell communication
inference. Figures are written as TIFF and PDF to meet journal requirements.

## Study

**Inhibition of Aggression by Social Touch**

Tamás Láng, Botond B. Drahos, Dávid Keller, Fanni Dóra, Ingrid Csordás,
Vivien Szendi, Gina Puska, Valery Grinevich, and Árpád Dobolyi

*Abstract:* TBD (manuscript in preparation).

**One-sentence summary:** Social touch–activated PIL neurons project to the MPOA
and acutely inhibit aggression; chemogenetic/optogenetic manipulation confirms
causality.

Three samples are analysed: `R1\_Cont`, `R2\_Affi`, `R3\_Aggr`.

## Data

Raw FASTQ files and 10x feature–barcode matrices will be deposited at NCBI GEO:
**GSEXXXXXX** (accession to be added).

Until then, place the 10x outputs locally and point `config/params.yml` at them.
Data directories are listed in `.gitignore` and are never uploaded.

The pipeline expects the filtered matrices in one flat folder with a per-sample
prefix (`R1\_Cont\_matrix.mtx.gz`, `R1\_Cont\_barcodes.tsv.gz`,
`R1\_Cont\_features.tsv.gz`, …), and each sample's cellranger
`raw\_feature\_bc\_matrix` folder in the standard layout. The raw matrices are
required by SoupX; without them step 03 falls back to decontX.

## Layout

```
config/
  params.yml                    every threshold and path; no numbers in the code
  markers\_used\_in\_figures.csv   marker panel for the figures
  neuropeptide\_genes.txt        neuropeptide gene list
scripts/
  00\_setup.R                    packages, options, shared helpers (sourced by all)
  01\_qc\_and\_merge.R             per-sample QC, merge
  02\_integration\_SCT.R          preliminary SCT integration and clustering
  03\_ambient\_removal.R          SoupX (or decontX) + re-integration
  04\_doublets\_and\_recluster.R   scDblFinder, removal, reclustering, diagnostics
  05\_annotation\_and\_renaming.R  exclusion, reclustering, lineage evidence, naming
  06\_neuron\_subset.R            neuronal subset, composition filter, figures
  07\_cellchat.R                 secreted-signalling analysis
  analyses/
    clustered\_dotplot.R         scCustomize clustered dotplot (optional figure)
    oligo\_maturation.R          oligodendrocyte maturation staging (optional)
utils/
  plotting.R                    plot\_theme() and write\_plot()
```

Outputs go to `data/processed/` (objects and tables), `data/processed/qc/`
(diagnostics) and `figures/`. All three are created on first run.

## Running it

Before the first run, edit `matrix\_dir` and the per-sample `raw\_dir` entries at
the top of `config/params.yml` so they point at your data. Everything else in
that file is an analysis parameter, not a path.

Set the working directory to the project root — the folder holding `scripts/`
and `config/` — and run the steps in order:

```r
source("scripts/01\_qc\_and\_merge.R")
source("scripts/02\_integration\_SCT.R")
source("scripts/03\_ambient\_removal.R")
source("scripts/04\_doublets\_and\_recluster.R")
source("scripts/05\_annotation\_and\_renaming.R")
source("scripts/06\_neuron\_subset.R")
source("scripts/07\_cellchat.R")
```

Each script sources `scripts/00\_setup.R` itself, so it does not need to be run
separately. The two scripts under `scripts/analyses/` are optional and can be
run at any point after step 05.

The steps are meant to be run one at a time, with the diagnostics inspected in
between — in particular, `qc/04\_cluster\_diagnostics.csv` after step 04, which is
what `config/clusters\_to\_drop.csv` refers to.

Every step reads the `.rds` written by the previous one; nothing is cached, so a
re-run always recomputes from its input. The slow steps are 03 and 06
(integration) and 04 (scDblFinder).

## Cluster identifiers change between steps

Clustering is performed four times — in steps 02, 04, 05 and 06 — and each pass
assigns new cluster IDs. Step 02's are provisional, exist only to give the
ambient-removal step its cluster structure, and are stored as `clusters\_prelim`.

The two editable config files therefore refer to different numbering:

|file|keyed to|evidence table with matching IDs|
|-|-|-|
|`config/clusters\_to\_drop.csv`|step-04 clusters|`qc/04\_cluster\_diagnostics.csv`|
|`config/cluster\_rename.csv`|step-05 clusters (after exclusion and reclustering)|`qc/05a\_cluster\_annotation\_evidence.csv`|

Cell-type names are assigned from `config/cluster\_rename\_suggested.csv`, which
step 05 regenerates from the lineage marker panels against the clusters that
exist in that run. Write `config/cluster\_rename.csv` only to override those
names by hand; it takes precedence, and step 05 stops rather than apply a table
that no longer matches the current clustering.

In step 06 the neuronal subclusters are renumbered once more, to close the gaps
left by the composition filter. The original IDs are kept in the
`cluster\_original` metadata column and in
`data/processed/06\_subneuron\_cluster\_id\_map.csv`.

## Analysis decisions that belong in the Methods

**Ambient RNA.** `params$ambient$method` is `auto`: SoupX when raw droplet
matrices are available for every sample, decontX otherwise. The method actually
used and the estimated contamination per sample are recorded in
`qc/03\_ambient\_method.csv`.

**Cluster exclusion.** Nothing is excluded unless `config/clusters\_to\_drop.csv`
exists, and every entry requires a written reason. Step 05 cross-checks each
entry against the step-04 diagnostics, warns when a cluster is dropped that the
diagnostics do not flag, and writes `qc/05\_dropped\_clusters\_with\_evidence.csv` —
the decision together with its supporting numbers.

**Neuronal composition filter.** Step 06 removes subclusters in which any sample
contributes less than `subneurons$composition\_filter$min\_sample\_pct` of the
size-adjusted composition. Size-adjusted means each sample's contribution is
first divided by its own neuron total, so 33.3% per sample represents
proportional parity; without this the raw split is confounded by differing
sample sizes.

With one animal per condition, sample identity is confounded with experimental
condition, so an excluded cluster cannot be attributed to technical variation
rather than to a condition-restricted population. Step 06 records the excluded
clusters' markers in `data/processed/06\_dropped\_clusters\_markers.csv` and saves
the unfiltered object as `data/processed/06\_neurons\_unfiltered.rds`, so a
sensitivity analysis can be run by pointing step 07 at that file instead.

**Marker p-values.** `FindAllMarkers` tests nuclei, not animals. Nuclei from one
animal are not independent observations, so the p-values are anticonservative.
This is the accepted convention for describing cluster markers, but it does not
support claims about differences between conditions.

**Mitochondrial genes.** Rat references name these in three different ways
(`Mt-nd1`, `AY172581.x`, or bare `ND1`/`COX1`/`CYTB`). `00\_setup.R` matches all
three and raises an error if none match, so the QC filter cannot silently pass
every cell.

## Reproducibility

Stochastic steps and how they are pinned:

|step|source of randomness|pinned by|
|-|-|-|
|02, 03|CCA/irlba inside `FindIntegrationAnchors`|`set.seed(42)` in `00\_setup.R`|
|03|decontX (fallback path only)|`seed = 42` in the call|
|04|scDblFinder's simulated doublets|`doublets$seed` and the `BPPARAM` RNG seed|
|04, 05, 06|PCA, UMAP, Louvain|Seurat defaults plus the explicit seeds in `params.yml`|
|07|CellChat permutation test, NMF|`cellchat$seed`|

Two caveats:

* `set.seed()` in `00\_setup.R` determines the integration steps only as long as
the amount of RNG consumed before them is identical. Run each step in a fresh
R session rather than interleaving other commands.
* Package versions matter. `sessionInfo()` is written next to the outputs at
every step; `renv.lock` records the exact versions used.

To reproduce the figures exactly, the raw count matrices (GEO), this repository,
and the processed objects from `data/processed/` are all needed. The processed
objects are what keep the result verifiable if a package version changes.

## Requirements

CRAN: Seurat (v5), SeuratObject, SoupX, Matrix, data.table, dplyr, tidyr,
ggplot2, patchwork, future, readr, yaml.

Bioconductor: scDblFinder, SingleCellExperiment, BiocParallel, ComplexHeatmap,
CellChat, NMF.

Optional: scCustomize (`analyses/clustered\_dotplot.R`), celda (decontX
fallback), ggalluvial (CellChat river plots).

