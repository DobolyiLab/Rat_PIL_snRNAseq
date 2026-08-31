# =============================================================================
# celltype_labeling.R
#
# Applies cell-type / neurotransmitter-identity labels to the clusters in
# sct_joined_annotated.rds (produced by stage 8 of PIL_snRNAseq_pipeline.R).
#
# Workflow:
#   1. Run this script once with no config/cluster_rename.csv present. It
#      writes config/cluster_rename_template.csv listing every current
#      cluster ID with an empty name column.
#   2. Inspect figures/celltype_dotplot_all.pdf and
#      figures/celltype_umap_by_cluster.pdf (written below) together with
#      celltype_markers_MAST.csv (from stage 8) to decide each cluster's
#      identity.
#   3. Fill in the "new_name" column of the template, save it as
#      config/cluster_rename.csv, and re-run this script.
#
# The cluster IDs here are specific to the seeded run that produced
# sct_joined_annotated.rds. If the upstream pipeline is re-run with different
# parameters, the cluster numbering changes and cluster_rename.csv must be
# redone against the new clustering -- it is not portable across runs.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("utils/plotting.R")   # defines plot_theme() and write_plot(plot, name, pdf = TRUE)
source("config/paths.R")     # defines out_dir

sct_joined <- readRDS(file.path(out_dir, "sct_joined_annotated.rds"))
Idents(sct_joined) <- "seurat_clusters"

markers <- read.csv("config/markers_used_in_figures.csv", stringsAsFactors = FALSE)$gene
markers <- intersect(markers, rownames(sct_joined))

if (!dir.exists(file.path(out_dir, "figures"))) dir.create(file.path(out_dir, "figures"), recursive = TRUE)

write_plot(
  DotPlot(sct_joined, features = markers, cols = c("lightgrey", "midnightblue"),
          col.min = -1, scale.min = 0) + RotatedAxis() + plot_theme(),
  "celltype_dotplot_all", pdf = TRUE
)
write_plot(
  DimPlot(sct_joined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + plot_theme(),
  "celltype_umap_by_cluster", pdf = TRUE
)
write_plot(
  DimPlot(sct_joined, reduction = "umap", group.by = "orig.ident") + plot_theme(),
  "celltype_umap_by_sample", pdf = TRUE
)

cur <- sort(as.character(unique(sct_joined$seurat_clusters)))
rename_path <- "config/cluster_rename.csv"

if (!file.exists(rename_path)) {

  write.csv(data.frame(cluster_id = cur, new_name = ""),
           "config/cluster_rename_template.csv", row.names = FALSE)
  message("config/cluster_rename.csv not found.\n",
          "  1. Inspect ", file.path(out_dir, "figures/celltype_dotplot_all.pdf"),
          " and celltype_markers_MAST.csv\n",
          "  2. Fill in config/cluster_rename_template.csv\n",
          "  3. Save it as config/cluster_rename.csv and re-run this script.")

} else {

  ren <- read.csv(rename_path, stringsAsFactors = FALSE)
  stopifnot("cluster_rename.csv needs cluster_id and new_name columns" =
              all(c("cluster_id", "new_name") %in% names(ren)))
  ren$cluster_id <- as.character(ren$cluster_id)
  ren <- ren[ren$cluster_id %in% cur & nzchar(ren$new_name), , drop = FALSE]

  unmapped <- setdiff(cur, ren$cluster_id)
  if (length(unmapped)) {
    message("Clusters with no name in cluster_rename.csv (left as numeric labels): ",
            paste(unmapped, collapse = ", "))
  }

  map <- setNames(ren$new_name, ren$cluster_id)
  sct_joined$celltype_detailed <- factor(ifelse(
    as.character(Idents(sct_joined)) %in% names(map),
    map[as.character(Idents(sct_joined))],
    as.character(Idents(sct_joined))
  ))
  Idents(sct_joined) <- "celltype_detailed"

  write_plot(
    DimPlot(sct_joined, reduction = "umap", label = TRUE) + plot_theme(),
    "celltype_umap_labeled", pdf = TRUE
  )
  write_plot(
    DotPlot(sct_joined, features = markers, group.by = "celltype_detailed",
            cols = c("lightgrey", "midnightblue"), col.min = -1, scale.min = 0) +
      RotatedAxis() + plot_theme(),
    "celltype_dotplot_labeled", pdf = TRUE
  )

  saveRDS(sct_joined, file.path(out_dir, "sct_joined_labeled.rds"))
  message(sprintf("--- CHECK: %d clusters labeled, %d left numeric ---",
                  length(unique(ren$new_name)), length(unmapped)))
  print(table(sct_joined$celltype_detailed))
}
