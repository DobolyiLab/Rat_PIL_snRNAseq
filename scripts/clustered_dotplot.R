# =============================================================================
# fig_clustered_dotplot.R  -- scCustomize::Clustered_DotPlot of the all-cell object
#
# Replaces the older standalone snippet, which referred to `sct_joined` and
# re-applied a cluster_rename map keyed to the PRE-reclustering IDs. Step 05
# already writes `celltype_detailed`, so no renaming happens here -- doing it
# again with a stale map would mislabel the plot.
# =============================================================================
source("scripts/00_setup.R")

if (!requireNamespace("scCustomize", quietly = TRUE)) {
  stop("scCustomize is not installed. Install it with:\n",
       "  install.packages('scCustomize')\n",
       "It also needs ComplexHeatmap:\n",
       "  BiocManager::install('ComplexHeatmap')", call. = FALSE)
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  stop("ComplexHeatmap is not installed: BiocManager::install('ComplexHeatmap')",
       call. = FALSE)
}

obj <- readRDS("data/processed/05_annotated_all.rds")
DefaultAssay(obj) <- "RNA"
obj <- join_layers_if_split(obj, "RNA")
obj <- ensure_rna_normalized(obj)

# Use the curated labels written by 05; fall back to numeric clusters.
ident_col <- if ("celltype_detailed" %in% colnames(obj@meta.data)) {
  "celltype_detailed"
} else {
  message("[clustered_dotplot] celltype_detailed not found -- falling back to seurat_clusters. ",
          "Run 05 with config/cluster_rename.csv first for labelled axes.")
  "seurat_clusters"
}
Idents(obj) <- ident_col
message("[clustered_dotplot] grouping by '", ident_col, "': ",
        paste(levels(Idents(obj)), collapse = ", "))

genes <- readr::read_csv("config/markers_used_in_figures.csv", show_col_types = FALSE)$gene
genes <- present_features(obj, genes, "RNA", "clustered dotplot markers")
stopifnot(length(genes) > 1)

# plot_km_elbow = FALSE makes this return the heatmap directly. With TRUE it
# returns a list(elbow_plot, heatmap) instead, which will not print to a device.
p <- scCustomize::Clustered_DotPlot(
  obj,
  features       = genes,
  plot_km_elbow  = FALSE,
  x_lab_rotate   = TRUE,
  row_label_size = params$plots$axis_text_size
)

# A ComplexHeatmap object does NOT auto-print inside a script or a device block.
# The original snippet worked only because the console auto-printed it; wrapped
# in tiff()/dev.off() it would have written an empty file. draw() is explicit.
save_heatmap <- function(hm, path, width, height, res = 400) {
  ext <- tools::file_ext(path)
  if (ext == "tiff") {
    tiff(path, width = width, height = height, res = res, compression = "lzw")
  } else {
    pdf(path, width = width / res, height = height / res)
  }
  on.exit(dev.off())
  ComplexHeatmap::draw(hm)
  invisible(path)
}

save_heatmap(p, "figures/fig_clustered_dotplot.tiff", width = 3200, height = 2400)
save_heatmap(p, "figures/fig_clustered_dotplot.pdf",  width = 2300, height = 2000)
message("[clustered_dotplot] written: figures/fig_clustered_dotplot.{tiff,pdf}")

writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_clustered_dotplot.txt")
