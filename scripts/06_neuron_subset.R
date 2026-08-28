# =============================================================================
# 06_neuron_subset.R  -- neuronal subset, subclustering, figures
#
# Single file, and no checkpoint: the script always runs from
# 05_annotated_all.rds. Caching the integrated subset would save a few minutes,
# but a cached file that outlives a change in step 04 or 05 silently overrides
# the fresh computation -- which is how a neuron subset built from mislabelled
# clusters once propagated into every figure below.
#
# FIXES vs the original version of this script:
#  * reads 05_annotated_all.rds, so the curated cell-type labels are used
#  * neurons are selected at CLUSTER level. The old per-cell rule
#    "Snap25 > 0 OR Rbfox3 > 0" keeps almost every nucleus in snRNA-seq,
#    because ambient Snap25 is present nearly everywhere -- oligodendrocytes
#    included
#  * batch integration is REDONE on the subset. The old code ran SCTransform +
#    PCA on the merged subset, silently discarding the integration from 02/03
#  * the composition filter is mathematically fixed. Percentages used to be
#    normalized WITHIN cluster, so with 3 samples the maximum was always >=33%
#    and the filter could never remove anything
#  * FindAllMarkers runs on `subneurons`, not on the all-cell object;
#    latent.vars dropped (invalid with the default Wilcoxon test)
#  * Fos / Calb1 expression values are no longer overwritten with zeros;
#    thresholding is recorded as an annotation, and co-expressing cells get
#    their own "Both" category instead of being silently called Fos
#  * UMAP column names and gene lists resolved defensively
# =============================================================================
source("scripts/00_setup.R")

npcs_sub <- params$subneurons$npcs
dims_sub <- params$subneurons$dims

obj <- readRDS("data/processed/05_annotated_all.rds")
DefaultAssay(obj) <- "RNA"
obj <- join_layers_if_split(obj, "RNA")
obj <- ensure_rna_normalized(obj)

# =============================================================================
# 1) Neuron selection at cluster level
# =============================================================================
neuro <- present_features(obj, c("Snap25", "Rbfox3", "Syt1"), "RNA", "neuronal markers")
glial <- present_features(obj, c("Plp1", "Mbp", "Aqp4", "Csf1r", "Pdgfra"), "RNA", "glial markers")
cts   <- get_expr(obj, "RNA", "counts")

cluster_var <- if ("celltype_detailed" %in% colnames(obj@meta.data)) "celltype_detailed" else "seurat_clusters"
cl <- as.character(obj@meta.data[[cluster_var]])

neuron_score <- data.frame(
  cluster    = cl,
  neuro_pos  = Matrix::colSums(cts[neuro, , drop = FALSE] > 0) >= 2,
  glial_pos  = Matrix::colSums(cts[glial, , drop = FALSE] > 0) >= 2
) %>%
  group_by(cluster) %>%
  summarise(n_cells = n(),
            neuronal_pct = round(100 * mean(neuro_pos), 1),
            glial_pct    = round(100 * mean(glial_pos), 1),
            .groups = "drop") %>%
  mutate(label_is_neuronal = grepl(params$subneurons$neuron_label_regex, cluster, ignore.case = TRUE),
         data_is_neuronal  = neuronal_pct >= params$subneurons$min_neuronal_pct &
           glial_pct    <= params$subneurons$max_glial_pct,
         keep = label_is_neuronal | (!any(label_is_neuronal) & data_is_neuronal))

readr::write_csv(neuron_score, "data/processed/qc/06_neuron_cluster_selection.csv")
print(as.data.frame(neuron_score))

neuron_clusters <- neuron_score$cluster[neuron_score$keep]
stopifnot("no neuronal cluster identified -- check neuron_label_regex in params.yml" =
            length(neuron_clusters) > 0)
disagree <- neuron_score$cluster[neuron_score$label_is_neuronal != neuron_score$data_is_neuronal]
if (length(disagree)) message("[06] label and marker evidence disagree for: ",
                              paste(disagree, collapse = ", "), " -- inspect before proceeding.")

# -------------------------------------------------------------------------
# Refuse to subset on names the markers contradict.
#
# Selection is driven by the cluster NAME (neuron_label_regex), so a rename
# table applied to the wrong clusters in step 05 silently produces a "neuron"
# subset made of oligodendrocytes. Every downstream figure then looks broken
# in a puzzling way rather than failing: the neurotransmitter dotplot comes
# out nearly empty, because Gad1 and Slc17a6 are genuinely absent.
#
# A cluster named GABA/GLU/DOPA must look neuronal by its markers too.
# -------------------------------------------------------------------------
bad <- neuron_score %>%
  filter(keep, neuronal_pct < params$subneurons$min_neuronal_pct)
if (nrow(bad)) {
  stop("Clusters selected as neuronal are not neuronal by their markers:\n",
       paste(sprintf("  %s: %.1f%% neuronal, %.1f%% glial (n=%d)",
                     bad$cluster, bad$neuronal_pct, bad$glial_pct, bad$n_cells),
             collapse = "\n"),
       "\n\nThe cell-type labels in 05_annotated_all.rds do not match the clusters ",
       "they are attached to. Re-run 05 with a cluster_rename.csv regenerated ",
       "for the current clustering, then re-run this script.",
       call. = FALSE)
}

Idents(obj) <- cluster_var
neurons <- subset(obj, idents = neuron_clusters)
log_step("06/neurons", neurons, sprintf("from clusters: %s", paste(neuron_clusters, collapse = ", ")))

# =============================================================================
# 2) Re-integrate the subset (do NOT lose the batch correction)
# =============================================================================
integrate_subset <- function(x, npcs, nfeat) {
  lst <- SplitObject(x, split.by = "orig.ident")
  # ncol() on a Seurat v5 object returns a double, so the vapply template must
  # be numeric(1); integer(1) fails with "values must be type 'integer'".
  sizes <- vapply(lst, function(y) as.integer(ncol(y)), integer(1))
  message("[06] cells per sample: ", paste(sprintf("%s=%d", names(sizes), sizes), collapse = ", "))
  
  if (length(lst) < 2 || min(sizes) < params$subneurons$min_cells_for_integration) {
    warning("Too few cells per sample for anchor-based integration; ",
            "falling back to a single SCTransform. Batch effects will NOT be corrected -- ",
            "state this in the Methods.")
    x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
    return(RunPCA(x, npcs = npcs, verbose = FALSE))
  }
  
  lst <- lapply(lst, function(y) {
    y <- SCTransform(y, vst.flavor = "v2", verbose = FALSE)
    RunPCA(y, npcs = npcs, verbose = FALSE)
  })
  feats <- SelectIntegrationFeatures(lst, nfeatures = nfeat)
  lst   <- PrepSCTIntegration(lst, anchor.features = feats)
  k <- max(5, min(100, min(sizes) - 1))
  anch <- FindIntegrationAnchors(object.list = lst, normalization.method = "SCT",
                                 anchor.features = feats, dims = 1:npcs)
  out <- IntegrateData(anchorset = anch, normalization.method = "SCT",
                       dims = 1:npcs, k.weight = k)
  DefaultAssay(out) <- "integrated"
  RunPCA(out, npcs = npcs, verbose = FALSE)
}

neurons <- integrate_subset(neurons, npcs_sub, params$integration$nfeatures)
neurons <- RunUMAP(neurons, dims = 1:dims_sub, min.dist = 0.5, verbose = FALSE)
neurons <- FindNeighbors(neurons, dims = 1:dims_sub, verbose = FALSE)
neurons <- FindClusters(neurons, resolution = params$subneurons$resolution, verbose = FALSE)

write_plot(DimPlot(neurons, label = TRUE, label.size = params$plots$label_size,
                   repel = TRUE) + ggtitle("Neuron subclustering") + plot_theme(),
           "06_umap_neuron_subclustering", pdf = TRUE)

# =============================================================================
# 3) Cluster composition -- BOTH normalizations, and an honest filter
# -----------------------------------------------------------------------------
# pct_of_cluster : how a cluster splits across samples (columns sum to 100 per
#                  cluster). With 3 samples the max is ALWAYS >= 33.3%, so a
#                  ">=10%" cut on this quantity removes nothing. This is the
#                  bug that made the original "data-driven" subsetting a no-op.
# pct_of_sample  : how much of a sample's neurons a cluster represents. This is
#                  the quantity a size-based inclusion rule should use.
# =============================================================================
tab <- as.data.frame(table(cluster = neurons$seurat_clusters, sample = neurons$orig.ident),
                     responseName = "cell_count")

comp <- tab %>%
  group_by(cluster) %>% mutate(pct_of_cluster = 100 * cell_count / sum(cell_count)) %>%
  group_by(sample)  %>% mutate(pct_of_sample  = 100 * cell_count / sum(cell_count)) %>%
  ungroup() %>%
  # Size-adjusted composition. The raw pct_of_cluster is confounded by how many
  # neurons each sample contributed overall: if one sample yields twice as many
  # nuclei, an entirely unbiased cluster still looks skewed toward it. Dividing
  # by each sample's own total first, then renormalizing within the cluster,
  # removes that confound -- 33.3% each means "proportionally equal".
  group_by(cluster) %>%
  mutate(pct_of_cluster_adj = 100 * pct_of_sample / sum(pct_of_sample)) %>%
  ungroup()
readr::write_csv(comp, "data/processed/06_cluster_sample_percent.csv")

cf <- params$subneurons$composition_filter
comp_col <- if (isTRUE(cf$size_adjusted)) "pct_of_cluster_adj" else "pct_of_cluster"

# The upper bound is optional and off by default. With three samples the
# percentages sum to 100, so a cluster where one sample exceeds 80% almost
# always has another below 10% and is already caught by the lower bound. Only
# an exact 80/10/10 split would escape, which does not occur in practice.
max_thr <- if (is.null(cf$max_sample_pct)) Inf else cf$max_sample_pct

cluster_summary <- comp %>%
  group_by(cluster) %>%
  summarise(n_cells           = sum(cell_count),
            max_pct_of_sample = max(pct_of_sample),
            min_contribution  = min(.data[[comp_col]]),
            max_contribution  = max(.data[[comp_col]]),
            dominant_sample   = sample[which.max(.data[[comp_col]])],
            .groups = "drop") %>%
  mutate(
    fails_size        = max_pct_of_sample < params$subneurons$min_pct_of_sample,
    fails_underrep    = isTRUE(cf$enabled) & min_contribution < cf$min_sample_pct,
    fails_dominated   = isTRUE(cf$enabled) & max_contribution > max_thr,
    keep              = !(fails_size | fails_underrep | fails_dominated))

readr::write_csv(cluster_summary, "data/processed/06_cluster_composition_filter.csv")
print(as.data.frame(cluster_summary))

keep_clusters <- as.character(cluster_summary$cluster[cluster_summary$keep])
dropped <- setdiff(as.character(levels(neurons$seurat_clusters)), keep_clusters)

if (isTRUE(cf$enabled)) {
  message(sprintf("[06] composition filter on %s: every sample must contribute >= %.0f%%%s",
                  comp_col, cf$min_sample_pct,
                  if (is.finite(max_thr)) sprintf(" and <= %.0f%%", max_thr) else ""))
  for (i in which(cluster_summary$fails_underrep | cluster_summary$fails_dominated)) {
    message(sprintf("      cluster %s (n=%d): range %.1f-%.1f%%, dominant %s -> dropped",
                    cluster_summary$cluster[i], cluster_summary$n_cells[i],
                    cluster_summary$min_contribution[i], cluster_summary$max_contribution[i],
                    cluster_summary$dominant_sample[i]))
  }
  if (any(cluster_summary$fails_underrep | cluster_summary$fails_dominated)) {
    warning("Clusters were removed for being sample-skewed. With one animal per ",
            "condition, a sample-specific cluster cannot be distinguished from a ",
            "condition-specific one, so this filter also removes any genuinely ",
            "condition-restricted population. State the criterion and the number ",
            "of clusters it removed in the Methods, and report what those clusters ",
            "were in a supplementary table.")
  }
}
message(sprintf("[06] keeping %d cluster(s), dropping %d (%s)",
                length(keep_clusters), length(dropped),
                if (length(dropped)) paste(dropped, collapse = ", ") else "none"))

# Keep the unfiltered object so a sensitivity analysis never requires redoing
# the integration: step 07 can simply be pointed at this file.
saveRDS(neurons, "data/processed/06_neurons_unfiltered.rds")

# --- characterise the clusters that are about to be removed ------------------
# A cluster excluded on composition grounds still has to be described, so a
# reader can judge whether a real population was discarded. Markers are computed
# on the UNFILTERED object, before the cluster disappears.
if (length(dropped)) {
  neurons <- ensure_rna_normalized(neurons)
  dropped_markers <- lapply(dropped, function(cl) {
    m <- FindMarkers(neurons, assay = "RNA", ident.1 = cl,
                     only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
    if (!nrow(m)) return(NULL)
    m$gene <- rownames(m); m$cluster <- cl; m
  })
  dropped_markers <- dplyr::bind_rows(dropped_markers)
  readr::write_csv(dropped_markers, "data/processed/06_dropped_clusters_markers.csv")
  
  top_dropped <- dropped_markers %>%
    filter(!grepl("^LOC|^AABR|^ENSRNOG", gene)) %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
  print(as.data.frame(top_dropped))
  
  readr::write_csv(
    cluster_summary %>% mutate(cluster = as.character(cluster)) %>%
      filter(cluster %in% dropped) %>%
      dplyr::left_join(top_dropped, by = "cluster"),
    "data/processed/06_dropped_clusters_summary.csv")
  message("[06] dropped-cluster markers written for the supplementary table.")
}

# NO re-clustering after the filter, deliberately.
# Whole clusters are removed here, not scattered cells, so the remaining
# clusters are structurally intact -- unlike doublet removal in 04, which
# hollowed clusters out and did require reclustering. Re-clustering now would
# also invalidate the stated criterion: nothing guarantees that clusters
# derived AFTER the filter still satisfy the composition rule, yet the Methods
# would claim they do.
subneurons <- if (length(dropped)) subset(neurons, idents = keep_clusters) else neurons
subneurons$seurat_clusters <- droplevels(factor(subneurons$seurat_clusters))
Idents(subneurons) <- "seurat_clusters"

# The criterion must still hold for exactly the clusters that were kept.
final_check <- comp %>%
  filter(as.character(cluster) %in% keep_clusters) %>%
  group_by(cluster) %>%
  summarise(min_contribution = min(.data[[comp_col]]), .groups = "drop")
stopifnot("composition criterion violated by a retained cluster" =
            !isTRUE(cf$enabled) || all(final_check$min_contribution >= cf$min_sample_pct))

# --- renumber the retained clusters consecutively ----------------------------
# Dropping clusters leaves gaps in the IDs (0..20 with three missing), which
# looks like an error on a figure. The original ID is kept in `cluster_original`
# and the mapping is written out, so any cluster on a figure can still be traced
# back to the composition table and the marker lists.
if (isTRUE(params$subneurons$renumber)) {
  old_levels <- levels(droplevels(factor(subneurons$seurat_clusters)))
  
  if (identical(params$subneurons$renumber_order, "size")) {
    sizes_old  <- table(as.character(subneurons$seurat_clusters))[old_levels]
    old_levels <- old_levels[order(sizes_old, decreasing = TRUE)]
  } else {
    old_levels <- old_levels[order(as.integer(old_levels))]
  }
  
  new_levels <- as.character(seq_along(old_levels) - 1)
  id_map <- setNames(new_levels, old_levels)
  
  # Assign via explicitly named vectors. id_map[x] returns a vector whose names
  # are the OLD CLUSTER IDS, and Seurat matches metadata names against cell
  # names -- hence "No cell overlap between new meta data and Seurat object".
  orig_vec <- factor(as.character(subneurons$seurat_clusters), levels = old_levels)
  names(orig_vec) <- colnames(subneurons)
  subneurons$cluster_original <- orig_vec
  
  new_vec <- factor(unname(id_map[as.character(orig_vec)]), levels = new_levels)
  names(new_vec) <- colnames(subneurons)
  subneurons$seurat_clusters <- new_vec
  Idents(subneurons) <- "seurat_clusters"
  stopifnot(!anyNA(subneurons$seurat_clusters))
  
  map_tbl <- data.frame(cluster_original = old_levels,
                        cluster_new      = new_levels,
                        n_cells          = as.integer(table(subneurons$seurat_clusters)[new_levels]))
  map_tbl <- dplyr::left_join(
    map_tbl,
    cluster_summary %>% mutate(cluster = as.character(cluster)) %>%
      select(cluster, min_contribution, max_contribution, dominant_sample),
    by = c("cluster_original" = "cluster"))
  readr::write_csv(map_tbl, "data/processed/06_subneuron_cluster_id_map.csv")
  
  message(sprintf("[06] renumbered %d clusters to 0-%d (%s order); original IDs kept in ",
                  length(new_levels), length(new_levels) - 1,
                  if (identical(params$subneurons$renumber_order, "size")) "size" else "original"),
          "`cluster_original`, mapping in 06_subneuron_cluster_id_map.csv")
}

log_step("06/subneurons", subneurons,
         sprintf("(%d clusters)", nlevels(subneurons$seurat_clusters)))

write_plot(DimPlot(subneurons, label = TRUE, label.size = params$plots$label_size,
                   repel = TRUE) + ggtitle("Subneurons") + plot_theme(),
           "06_umap_subneurons", pdf = TRUE)

# =============================================================================
# 4) Markers -- from THIS object, on normalized RNA
# =============================================================================
DefaultAssay(subneurons) <- "RNA"
subneurons <- join_layers_if_split(subneurons, "RNA")   # re-split by integrate_subset()
subneurons <- ensure_rna_normalized(subneurons)

markers <- FindAllMarkers(subneurons, assay = "RNA", only.pos = TRUE,
                          min.pct         = params$subneurons$marker_min_pct,
                          logfc.threshold = params$subneurons$marker_logfc)
readr::write_csv(markers, "data/processed/06_subneurons_markers_all.csv")

top3_table <- markers %>%
  filter(!grepl("^LOC|^AABR|^ENSRNOG", gene), pct.1 >= 0.1) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 3) %>% ungroup()
readr::write_csv(top3_table, "data/processed/06_subneurons_markers_top3.csv")

gene_levels <- present_features(subneurons, rev(unique(top3_table$gene)), "RNA", "top3 markers")
p_top3 <- DotPlot(subneurons, features = gene_levels, dot.scale = 8, dot.min = params$plots$dot_min) +
  coord_flip() + plot_theme() + rotated_x()
write_plot(p_top3, "06_dotplot_top3_markers", pdf = TRUE)

p_top3 <- DotPlot(subneurons, features = gene_levels) +
  coord_flip() + plot_theme() + RotatedAxis() +
  scale_radius(range = c(0, 10))
ggsave("figures/06_dotplot_top3_markers2.tiff", p_top3, width = 10, height = 12, units = "in", dpi = 300)


p_top3 <- DotPlot(
  subneurons,
  features = gene_levels
) +
  coord_flip() +
  plot_theme() +
  RotatedAxis() +
  scale_radius(range = c(0, 10)) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(
  "figures/06_dotplot_top3_markers2.tiff",
  p_top3,
  width = 12,
  height = 14,
  units = "in",
  dpi = 300
)

# =============================================================================
# 5) Fos / Calb1
# -----------------------------------------------------------------------------
# The original code set expression to 0 for every cell of a cluster that had
# fewer than 8 cells above threshold. That silently edits the measurements and
# is hard to defend in review. Here the real values are plotted, the threshold
# is recorded as an annotation, and clusters failing it are listed in a table.
# Fos is an IEG, so it is also split by sample -- that is where the
# Cont/Affi/Aggr contrast lives.
# =============================================================================
genes_fc <- present_features(subneurons, c("Fos", "Calb1"), "RNA", "activity markers")
thr      <- params$subneurons$fos_calb_threshold
min_pos  <- params$subneurons$fos_calb_min_cells

expr_data <- FetchData(subneurons, vars = c(genes_fc, "seurat_clusters", "orig.ident"))
expr_data$cell <- rownames(expr_data)

cluster_pass <- expr_data %>%
  pivot_longer(all_of(genes_fc), names_to = "gene", values_to = "expr") %>%
  group_by(seurat_clusters, gene) %>%
  summarise(n_pos = sum(expr > thr), passes = n_pos >= min_pos, .groups = "drop")
readr::write_csv(cluster_pass, "data/processed/qc/06_Fos_Calb1_cluster_threshold.csv")

write_plot(VlnPlot(subneurons, features = genes_fc, group.by = "seurat_clusters",
                   pt.size = 0, flip = TRUE, stack = length(genes_fc) > 1) + plot_theme(),
           "06_violin_Fos_Calb1", pdf = TRUE)
write_plot(VlnPlot(subneurons, features = genes_fc, group.by = "seurat_clusters",
                   split.by = "orig.ident", pt.size = 0) + plot_theme(),
           "06_violin_Fos_Calb1_by-sample", pdf = TRUE)

# UMAP highlight, with an explicit "Both" class
pass_map <- cluster_pass %>% filter(passes) %>%
  transmute(key = paste(seurat_clusters, gene)) %>% pull(key)

status <- expr_data %>%
  pivot_longer(all_of(genes_fc), names_to = "gene", values_to = "expr") %>%
  mutate(pos = expr > thr & paste(seurat_clusters, gene) %in% pass_map) %>%
  select(cell, gene, pos) %>%
  pivot_wider(names_from = gene, values_from = pos)

status$marker_group <- with(status, dplyr::case_when(
  Fos & Calb1 ~ "Both",
  Fos         ~ "Fos",
  Calb1       ~ "Calb1",
  TRUE        ~ "None"))

uc <- umap_cols(subneurons)
df_umap <- as.data.frame(Embeddings(subneurons, "umap"))
df_umap$cell <- rownames(df_umap)
df_umap <- dplyr::left_join(df_umap, status[, c("cell", "marker_group")], by = "cell")

pal <- c(None = "lightgrey", Calb1 = "deepskyblue3", Fos = "deeppink2", Both = "darkorchid3")
p_hl <- ggplot(df_umap, aes(x = .data[[uc[1]]], y = .data[[uc[2]]])) +
  geom_point(data = subset(df_umap, marker_group == "None"),
             aes(color = marker_group), size = 0.6, alpha = 0.3) +
  geom_point(data = subset(df_umap, marker_group != "None"),
             aes(color = marker_group), size = 1, alpha = 0.9) +
  scale_color_manual(values = pal) +
  labs(title = sprintf("Fos / Calb1 (expr > %.2f, clusters with >= %d positive cells)", thr, min_pos)) +
  plot_theme()
write_plot(p_hl, "06_umap_Fos_Calb1", pdf = TRUE)

# =============================================================================
# 6) Neurotransmitter and neuropeptide dotplots
# =============================================================================
nt_list <- list(
  Glutamatergic = c("Slc17a6", "Slc17a7"),
  GABAergic     = c("Gad1", "Gad2", "Slc32a1"),
  Dopaminergic  = c("Th", "Slc6a3")
)
nt_list <- lapply(nt_list, function(g) present_features(subneurons, g, "RNA", "neurotransmitter"))
nt_list <- nt_list[lengths(nt_list) > 0]
write_plot(
  DotPlot(
    subneurons,
    features = nt_list,
    group.by = "seurat_clusters",
    dot.scale = 8,
    dot.min = 0.05
  ) +
    scale_radius(range = c(1, 10)) +
    scale_y_discrete(limits = rev) +
    RotatedAxis() +
    plot_theme() +
    theme(
      strip.text = element_text(
        size = 15,
        face = "bold"
      )
    ),
  "06_dotplot_neurotransmitters",
  pdf = TRUE
)

neuropep <- readLines("config/neuropeptide_genes.txt")
neuropep <- unique(trimws(neuropep[nzchar(trimws(neuropep))]))
np_cfg <- params$neuropeptides
n_listed <- length(neuropep)
neuropep <- present_features(subneurons, neuropep, "RNA", "neuropeptides")

if (identical(np_cfg$selection, "expressed") ||
    (identical(np_cfg$selection, "all") && isTRUE(params$plots$drop_unexpressed))) {
  neuropep <- expressed_features(subneurons, neuropep, "seurat_clusters")
  
} else if (identical(np_cfg$selection, "significant")) {
  # Keep only neuropeptides that are significant POSITIVE markers of at least
  # one subcluster, using the FindAllMarkers result computed in section 4 --
  # one-vs-rest across clusters, not a comparison between conditions.
  np_sig <- markers %>%
    filter(gene %in% neuropep,
           p_val_adj < np_cfg$padj,
           avg_log2FC >= np_cfg$min_log2fc) %>%
    group_by(gene) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  readr::write_csv(
    np_sig %>% select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val, p_val_adj),
    "data/processed/06_neuropeptides_significant.csv")
  
  if (isTRUE(np_cfg$order_by_cluster)) {
    np_sig <- np_sig %>% arrange(as.integer(as.character(cluster)), desc(avg_log2FC))
  }
  neuropep <- as.character(np_sig$gene)
  
  message(sprintf(paste0("[06] neuropeptides: %d listed -> %d present -> %d significant ",
                         "markers (padj < %.3g, log2FC >= %.2f)"),
                  n_listed, length(present_features(subneurons, readLines("config/neuropeptide_genes.txt"),
                                                    "RNA", "neuropeptides")),
                  length(neuropep), np_cfg$padj, np_cfg$min_log2fc))
  if (!length(neuropep)) stop("No neuropeptide passed the significance filter. ",
                              "Relax neuropeptides$padj / min_log2fc, or set selection: expressed.",
                              call. = FALSE)
}
write_plot(DotPlot(subneurons, features = neuropep, group.by = "seurat_clusters", dot.scale = 8, dot.min = params$plots$dot_min) +
             plot_theme() + rotated_x() + scale_y_discrete(limits = rev),
           "06_dotplot_neuropeptides", pdf = TRUE)

subneurons <- trim_scale_data(subneurons)
saveRDS(subneurons, "data/processed/06_subneurons.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_06.txt")
