# =============================================================================
# config/paths.R
#
# Edit this file to point at your own data. Nothing in scripts/ should
# contain a hardcoded path -- if you find one, it is a bug, please report it.
# =============================================================================

r1_path <- "<path-to-data>/R1_Cont_with_WPRE/outs/filtered_feature_bc_matrix"
r2_path <- "<path-to-data>/R2_Affi_with_WPRE/outs/filtered_feature_bc_matrix"
r3_path <- "<path-to-data>/R3_Aggr_with_WPRE/outs/filtered_feature_bc_matrix"

sample.dirs <- list(
  "<path-to-data>/R1_Cont_with_WPRE",
  "<path-to-data>/R2_Affi_with_WPRE",
  "<path-to-data>/R3_Aggr_with_WPRE"
)

out_dir <- "<path-to-output>/outputs"
