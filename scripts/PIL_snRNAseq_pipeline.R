# =============================================================================
# PIL_snRNAseq_pipeline.R
#
# Rat posterior intralaminar thalamic nucleus (PIL) single-nucleus RNA-seq
# analysis: raw data through the final neuronal subclusters and their
# neuropeptide marker genes. This is the code underlying the cell numbers,
# cluster counts, and marker tables reported in the manuscript's Results
# and Methods (QC, ambient RNA correction, doublet detection, integration,
# cell-type marker identification, neuronal subclustering, neuropeptide
# quantification).
#
# NOT YET IN THIS FILE (described in the Methods but verified separately --
# see README.md for status): CellChat cell-cell communication analysis,
# Fos/Calb1 quantification across experimental samples, assignment of
# cell-type / neurotransmitter-identity labels to clusters, generation of
# the QC metrics table (Table S1).
#
# Every stochastic step is seeded (STAGE 0). Read README.md before running.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(future)
})

# =============================================================================
# STAGE 0 -- setup
# =============================================================================
set.seed(42)
future::plan("sequential")
options(future.globals.maxSize = 20 * 1024^3, future.seed = TRUE)

source("config/paths.R")   # defines r1_path, r2_path, r3_path, sample.dirs, out_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

N_FEATURES        <- 3000   # integration features, used at every integration step
DOUBLET_THRESHOLD <- 60     # % scDblFinder doublet rate above which a cluster is excluded
MITO_TOP_N        <- 20     # how many top marker genes per cluster are screened for
                            # mitochondrial genes, when flagging low-quality clusters


# =============================================================================
# STAGE 1 -- raw data load and per-sample QC
# =============================================================================
seu1 <- CreateSeuratObject(counts = Read10X(data.dir = r1_path), project = "R1_Cont", min.cells = 1, min.features = 1)
seu2 <- CreateSeuratObject(counts = Read10X(data.dir = r2_path), project = "R2_Affi", min.cells = 1, min.features = 1)
seu3 <- CreateSeuratObject(counts = Read10X(data.dir = r3_path), project = "R3_Aggr", min.cells = 1, min.features = 1)
seurat.list <- list(seu1, seu2, seu3)
sample.size <- length(seurat.list)

seurat.list[[1]]$condition <- "control"
seurat.list[[2]]$condition <- "affiliative"
seurat.list[[3]]$condition <- "aggressive"

# Mitochondrial content: two naming conventions are summed because this rat
# reference mixes them (lower-case "mt-" gene symbols and "AY"-accession
# mitochondrial contigs).
seurat.list <- lapply(seurat.list, function(x) {
  x <- PercentageFeatureSet(x, pattern = "^mt-",    col.name = "percent.mt", assay = "RNA")
  x <- PercentageFeatureSet(x, pattern = "^AY[0-9]", col.name = "percent.ay", assay = "RNA")
  x@meta.data$percent.mt <- as.numeric(x@meta.data$percent.mt) + as.numeric(x@meta.data$percent.ay)
  x
})

seurat.list <- lapply(seurat.list, function(x) subset(x, subset = nFeature_RNA > 200 & percent.mt < 5))

message("--- CHECK: nuclei per sample after stage-1 QC ---")
for (x in seurat.list) message(sprintf("  %s: %d nuclei", unique(x$orig.ident), ncol(x)))


# =============================================================================
# STAGE 2 -- preliminary integration (provides cluster labels for SoupX)
# =============================================================================
seurat.list <- lapply(seurat.list, function(x) {
  SCTransform(x, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE) %>%
    RunPCA(npcs = 20, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = N_FEATURES)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT",
                                  anchor.features = features)
sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

sct <- RunPCA(sct, npcs = 20, verbose = FALSE)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:20, verbose = FALSE)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:20)
sct <- FindClusters(sct, resolution = 0.8)

message("--- CHECK: preliminary clusters (stage 2) ---")
print(table(sct$seurat_clusters))


# =============================================================================
# STAGE 3 -- barcode / sample bookkeeping, needed to map SoupX corrections
#            back onto the stage-2 cluster labels
# =============================================================================
sct$orig_cellid <- colnames(sct)
sct$cell_short  <- sub("_\\d+$", "", sct$orig_cellid)
sct$sample_num  <- sapply(colnames(sct), function(x) as.numeric(sub(".*_(\\d+)$", "\\1", x)))


# =============================================================================
# STAGE 4 -- SoupX ambient RNA correction, using stage-2 clusters
# =============================================================================
adjusted_matrices <- vector("list", sample.size)
seurat_objects    <- vector("list", sample.size)

for (i in seq_len(sample.size)) {
  filtered_mat <- Read10X(file.path(sample.dirs[[i]], "outs/filtered_feature_bc_matrix/"))
  raw_mat      <- Read10X(file.path(sample.dirs[[i]], "outs/raw_feature_bc_matrix/"))
  sc <- SoupChannel(tod = raw_mat, toc = filtered_mat)
  sc$metaData$cell <- rownames(sc$metaData)
  sc$metaData$celltype <- apply(sc$metaData, 1, function(x)
    as.character(subset(sct@meta.data, cell_short == x["cell"] & sample_num == i)$seurat_clusters))
  sc  <- setClusters(sc, clusters = as.character(sc$metaData$celltype))
  sc  <- autoEstCont(sc, forceAccept = TRUE)
  out <- adjustCounts(sc, roundToInt = TRUE)
  adjusted_matrices[[i]] <- out

  adj.cells <- colnames(out)
  meta <- sct@meta.data[sct$cell_short %in% adj.cells & sct$sample_num == i, ]
  rownames(meta) <- sub("_\\d+$", "", meta$orig_cellid)
  seurat_objects[[i]] <- CreateSeuratObject(out, meta.data = meta)
  seurat_objects[[i]] <- seurat_objects[[i]][, colnames(seurat_objects[[i]]) %in% rownames(meta)]
}


# =============================================================================
# STAGE 5 -- QC and integration of the SoupX-corrected counts -> soupx.sct
# =============================================================================
soupx <- seurat_objects
soupx <- lapply(soupx, function(x) PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt", assay = "RNA"))
soupx <- lapply(soupx, function(x) subset(x, features = rownames(x)[!grepl("^MT-", rownames(x))]))
soupx <- lapply(soupx, function(x) subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5))
soupx <- lapply(soupx, function(x) SCTransform(x, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE) %>%
                  RunPCA(npcs = 20, verbose = FALSE))

features <- SelectIntegrationFeatures(object.list = soupx, nfeatures = N_FEATURES)
soupx <- PrepSCTIntegration(object.list = soupx, anchor.features = features)
soupx.anchors <- FindIntegrationAnchors(object.list = soupx, normalization.method = "SCT",
                                        anchor.features = features)
soupx.sct <- IntegrateData(anchorset = soupx.anchors, normalization.method = "SCT")

soupx.sct <- RunPCA(soupx.sct, npcs = 20, verbose = FALSE)
soupx.sct <- RunUMAP(soupx.sct, reduction = "pca", dims = 1:20, verbose = FALSE)
soupx.sct <- FindNeighbors(soupx.sct, reduction = "pca", dims = 1:20)
soupx.sct <- FindClusters(soupx.sct, resolution = 0.8)

message("--- CHECK: clusters after SoupX + re-integration (stage 5) ---")
print(table(soupx.sct$seurat_clusters))

saveRDS(soupx.sct, file.path(out_dir, "soupx_sct.rds"))


# =============================================================================
# STAGE 6 -- doublet detection and cluster exclusion
#
# Two independent exclusion criteria, both computed rather than eyeballed:
#   (a) clusters where >60% of nuclei are scDblFinder doublets
#   (b) clusters whose top marker genes are dominated by mitochondrial genes
#       (mitochondrial/ambient-RNA signature, not a real cell type)
#
# (b) is deliberately marker-based, not an average-expression z-score: with
# several genuinely elevated clusters, a mean/sd z-score on raw expression
# gets pulled up by those same clusters and under-flags all but the most
# extreme one. Screening the TOP MARKER GENES (i.e. what actually defines
# each cluster relative to the others) is the criterion that was actually
# used, and is robust to that masking effect.
# =============================================================================
DefaultAssay(soupx.sct) <- "RNA"
sce <- as.SingleCellExperiment(soupx.sct, assay = "SCT")
sce <- scDblFinder(sce, samples = "sample_num", BPPARAM = SnowParam(16), includePCs = 20)
soupx.sct$scDblFinder.class <- sce$scDblFinder.class[match(colnames(soupx.sct), colnames(sce))]

Idents(soupx.sct) <- "seurat_clusters"

# -- (a) doublet-dominated clusters --------------------------------------------
tbl <- table(soupx.sct$seurat_clusters, soupx.sct$scDblFinder.class)
doublet_percent <- prop.table(tbl, margin = 1)[, "doublet"] * 100
doublet_clusters <- names(doublet_percent)[doublet_percent > DOUBLET_THRESHOLD]
message("--- CHECK: doublet % per cluster ---")
print(round(doublet_percent, 1))
message(sprintf("Doublet-dominated clusters (>%d%%): %s", DOUBLET_THRESHOLD,
                paste(doublet_clusters, collapse = ", ")))

# -- (b) mitochondrial-marker-dominated clusters -------------------------------
DefaultAssay(soupx.sct) <- "SCT"
soupx.sct <- PrepSCTFindMarkers(soupx.sct, verbose = FALSE)

cluster_markers <- FindAllMarkers(
  soupx.sct,
  only.pos            = TRUE,
  min.pct             = 0.25,
  logfc.threshold     = 0.5,
  max.cells.per.ident  = 500,   # diagnostic only -- speeds this up, does not
  verbose             = FALSE   # affect the final reported cell-type markers
)

mito_pattern <- "^[Mm][Tt]-|^AY172581|^AY[0-9]{6}|^(ND[1-6]|ND4L|COX[1-3]|ATP6|ATP8|CYTB)$"

top_markers_per_cluster <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = MITO_TOP_N) %>%
  summarise(
    n_mito_in_top = sum(grepl(mito_pattern, gene)),
    top_genes     = paste(gene, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_mito_in_top))

message("--- CHECK: top marker genes per cluster, mitochondrial-gene count ---")
print(top_markers_per_cluster, n = Inf)

mito_cluster <- as.character(top_markers_per_cluster$cluster[
  top_markers_per_cluster$n_mito_in_top >= MITO_TOP_N / 2
])
message("Mitochondrial-marker-dominated cluster(s): ", paste(mito_cluster, collapse = ", "))

# -- apply both exclusions ------------------------------------------------------
exclude_clusters <- union(doublet_clusters, mito_cluster)
message("Excluding: ", paste(exclude_clusters, collapse = ", "))

Idents(soupx.sct) <- "scDblFinder.class"
soupx.sct <- subset(soupx.sct, idents = "singlet")
Idents(soupx.sct) <- "seurat_clusters"
sct <- if (length(exclude_clusters)) {
  subset(soupx.sct, idents = exclude_clusters, invert = TRUE)
} else {
  soupx.sct
}

message(sprintf("--- CHECK: %d nuclei remain after doublet + cluster exclusion ---", ncol(sct)))


# =============================================================================
# STAGE 7 -- final integration and clustering
#
# This is the integration and clustering that all downstream cell-type and
# neuronal-subcluster analyses are built on: rpca-based anchor finding
# (canonical correlation analysis was used at the earlier, provisional
# integration steps in stages 2 and 5; rpca is used here for the final,
# reported clustering), 50 PCs, UMAP/neighbours on the top 15 PCs with
# min.dist = 0.5, and Louvain clustering at resolution 0.5.
# =============================================================================
seurat.list <- SplitObject(sct, split.by = "orig.ident")
seurat.list <- lapply(seurat.list, function(x) {
  SCTransform(x, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = N_FEATURES)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = "rpca")

sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

sct <- RunPCA(sct, npcs = 50, verbose = FALSE)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:15, min.dist = 0.5, verbose = FALSE)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:15)
sct <- FindClusters(sct, resolution = 0.5)

message(sprintf("--- CHECK: %d final clusters after integration (stage 7) ---",
                nlevels(sct$seurat_clusters)))
print(table(sct$seurat_clusters))

saveRDS(sct, file.path(out_dir, "sct_final_integrated.rds"))


# =============================================================================
# STAGE 8 -- cell-type marker genes
#
# Markers are computed on the RNA assay (log-normalized, not the batch-
# corrected "integrated" assay -- integrated values are for clustering/
# visualization, not differential expression). MAST is used for cell-type
# marker identification; neuropeptide quantification (stage 11) uses the
# Wilcoxon rank-sum test instead -- this is intentional, not an
# inconsistency, and both are described separately in the Methods.
# =============================================================================
DefaultAssay(sct) <- "RNA"
sct <- JoinLayers(sct)
sct <- NormalizeData(sct, verbose = FALSE)
Idents(sct) <- "seurat_clusters"

sct_joined <- sct

sct_markers <- FindAllMarkers(
  sct_joined,
  assay           = "RNA",
  only.pos        = TRUE,
  min.pct         = 0.1,
  logfc.threshold = 0.5,
  test.use        = "MAST"
)
sct_markers <- sct_markers %>% filter(p_val_adj < 0.05)

write.csv(sct_markers, file.path(out_dir, "celltype_markers_MAST.csv"), row.names = FALSE)
message(sprintf("--- CHECK: %d significant cell-type markers across %d clusters ---",
                nrow(sct_markers), length(unique(sct_markers$cluster))))

saveRDS(sct_joined, file.path(out_dir, "sct_joined_annotated.rds"))


# =============================================================================
# STAGE 9 -- neuron subset and re-clustering
#
# Cluster IDs identified as neuronal by marker expression (Snap25, Rbfox3)
# on sct_joined -- see figures/celltype_markers_MAST.csv and the UMAP/
# DotPlot for this assignment.
# =============================================================================
neuron_cluster_ids <- c("2", "4", "8", "11", "14", "16", "17")

Idents(sct_joined) <- "seurat_clusters"
neurons <- subset(sct_joined, idents = neuron_cluster_ids)
neurons <- SCTransform(neurons, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE)
neurons <- RunPCA(neurons, npcs = 50, verbose = FALSE)
neurons <- RunUMAP(neurons, dims = 1:20, min.dist = 0.5)
neurons <- FindNeighbors(neurons, dims = 1:20)
neurons <- FindClusters(neurons, resolution = 0.6)

message(sprintf("--- CHECK: %d neuronal nuclei in %d subclusters before composition filter ---",
                ncol(neurons), nlevels(neurons$seurat_clusters)))


# =============================================================================
# STAGE 10 -- composition filter
#
# A subcluster is retained only if EACH of the three experimental samples
# (Cont, Affi, Aggr) contributes at least 10% of that subcluster's nuclei.
# =============================================================================
tab <- as.data.frame(table(neurons$seurat_clusters, neurons$orig.ident))
colnames(tab) <- c("cluster", "sample", "cell_count")

cluster_pct <- tab %>%
  group_by(cluster) %>%
  mutate(pct = 100 * cell_count / sum(cell_count)) %>%
  summarise(min_pct = min(pct), max_pct = max(pct), .groups = "drop")

message("--- CHECK: per-sample composition of each pre-filter subcluster ---")
print(cluster_pct, n = Inf)

keep_clusters <- as.character(cluster_pct$cluster[cluster_pct$min_pct >= 10])

message(sprintf("--- CHECK: %d of %d clusters retained (all three samples >=10%%) ---",
                length(keep_clusters), nrow(cluster_pct)))

subneurons <- subset(neurons, idents = keep_clusters)
subneurons <- SCTransform(subneurons, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE)
subneurons <- RunPCA(subneurons, npcs = 50, verbose = FALSE)
subneurons <- RunUMAP(subneurons, dims = 1:20)
subneurons <- FindNeighbors(subneurons, dims = 1:20)
subneurons <- FindClusters(subneurons, resolution = 0.6)

message(sprintf("--- CHECK: %d nuclei in %d final neuronal subclusters ---",
                ncol(subneurons), nlevels(subneurons$seurat_clusters)))

saveRDS(subneurons, file.path(out_dir, "subneurons.rds"))


# =============================================================================
# STAGE 11 -- neuropeptide gene quantification
#
# Separate FindAllMarkers() call, deliberately with no test.use argument
# (-> Seurat's Wilcoxon rank-sum default), distinct from stage 8's
# MAST-based cell-type markers.
# =============================================================================
neuropep_path    <- file.path("config", "neuropeptide_genes.txt")  # one gene symbol per line
neuropep_list    <- readLines(neuropep_path)
neuropep_list    <- unique(neuropep_list[nzchar(neuropep_list)])
neuropep_present <- intersect(neuropep_list, rownames(subneurons))

Idents(subneurons) <- "seurat_clusters"
neuropep_markers <- FindAllMarkers(
  subneurons,
  features = neuropep_present,
  only.pos = TRUE
)
neuropep_markers <- neuropep_markers %>% filter(p_val_adj < 0.05)

write.csv(neuropep_markers, file.path(out_dir, "neuropeptide_marker_quantification.csv"),
         row.names = FALSE)
message(sprintf("--- CHECK: %d of %d listed neuropeptides are significant markers of >=1 subcluster ---",
                length(unique(neuropep_markers$gene)), length(neuropep_present)))

message("\nDone through subneurons + neuropeptide markers. CellChat, Fos/Calb1 ",
        "quantification, and cell-type/neurotransmitter labelling are separate ",
        "scripts -- see README.md.")
