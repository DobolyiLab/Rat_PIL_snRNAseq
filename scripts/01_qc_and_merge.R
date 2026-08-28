# =============================================================================
# 01_qc_and_merge.R  -- per-sample QC, merge
# FIXES vs original:
#  * mito pattern is asserted (was "^mt-", which matches nothing in most rat refs)
#  * percent.mt is no longer built by adding two overlapping percentages
#  * upper nFeature bound applied here as well (was only applied in step 03)
#  * merge() uses add.cell.ids -> deterministic "SAMPLE_BARCODE" cell names,
#    instead of Reduce()'s cumulative _1/_2 suffixes
#  * per-sample QC table written out for the Methods / reviewer response
# =============================================================================
source("scripts/00_setup.R")
check_sample_paths(need_raw = TRUE)   # warns if SoupX will not be possible

qc_rows <- list()

objs <- lapply(params$samples, function(s) {
  so <- CreateSeuratObject(counts = read_sample_matrix(s, "filtered"),
                           project = s$id, min.cells = 1, min.features = 1)
  so <- add_percent_mt(so, strict = TRUE)
  n_before <- ncol(so)
  
  so <- subset(so, subset = nFeature_RNA > params$qc$min_features &
                 nFeature_RNA < params$qc$max_features &
                 percent.mt   < params$qc$max_percent_mt)
  
  qc_rows[[s$id]] <<- data.frame(
    sample          = s$id,
    cells_cellranger = n_before,
    cells_after_qc   = ncol(so),
    median_nFeature  = median(so$nFeature_RNA),
    median_nCount    = median(so$nCount_RNA),
    median_percent_mt = round(median(so$percent.mt), 3)
  )
  log_step(paste0("01/", s$id), so, sprintf("(%d -> %d after QC)", n_before, ncol(so)))
  so
})
names(objs) <- SAMPLE_IDS

# Deterministic cell names: "R1_Cont_AAACCCAAGAAACACT-1"
combined <- merge(objs[[1]], y = objs[-1], add.cell.ids = SAMPLE_IDS)
stopifnot(!anyDuplicated(colnames(combined)))
stopifnot(all(sub("_[ACGT]+-\\d+$", "", colnames(combined)) %in% SAMPLE_IDS))

readr::write_csv(dplyr::bind_rows(qc_rows), "data/processed/qc/01_qc_per_sample.csv")

p_vln <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "orig.ident", pt.size = 0, ncol = 3) + plot_theme()
write_plot(p_vln, "01_violin_feature-count-mt", pdf = TRUE)

saveRDS(combined, "data/processed/01_combined_qc.rds")
writeLines(capture.output(sessionInfo()), "data/processed/sessionInfo_01.txt")

