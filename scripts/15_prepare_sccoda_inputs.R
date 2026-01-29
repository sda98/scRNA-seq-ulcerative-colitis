#!/usr/bin/env Rscript

# =============================================================================
# 15_prepare_sccoda_inputs.R
#
# Purpose:
#   1) Load final cleaned/annotated Seurat object (celltype_final present).
#   2) Build per-sample cell count matrix (samples x cell types).
#   3) Export scCODA inputs:
#        - sccoda_counts.csv   (samples x cell types; no sample/condition columns)
#        - sccoda_metadata.csv (sample_id, condition, condition_code)
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_clean.rds
#     - Must contain meta.data columns:
#         - celltype_final
#         - sample_id (or gsm_id)
#         - condition (or group)
#
# Outputs:
#   results/03_compositional_analysis/sessionInfo.txt
#   results/03_compositional_analysis/tables/sccoda_counts.csv
#   results/03_compositional_analysis/tables/sccoda_metadata.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(
  root_dir, "results", "02_clustering_analysis", "objects",
  "seurat_integrated_clean.rds"
)
if (!file.exists(in_obj)) {
  stop(
    "Input object not found: ", in_obj,
    "\nExpected a cleaned object from your final annotation step (e.g., Script 13)."
  )
}

out_dir   <- file.path(root_dir, "results", "03_compositional_analysis")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# Load object
# =========================
cat("\nLoading cleaned Seurat object...\n")
seurat_clean <- readRDS(in_obj)

if (!"celltype_final" %in% colnames(seurat_clean@meta.data)) {
  stop("Missing 'celltype_final' in metadata.")
}

sample_col <- if ("sample_id" %in% colnames(seurat_clean@meta.data)) {
  "sample_id"
} else if ("gsm_id" %in% colnames(seurat_clean@meta.data)) {
  "gsm_id"
} else {
  stop("No sample identifier found (expected 'sample_id' or 'gsm_id').")
}

condition_col <- if ("condition" %in% colnames(seurat_clean@meta.data)) {
  "condition"
} else if ("group" %in% colnames(seurat_clean@meta.data)) {
  "group"
} else {
  stop("No condition column found (expected 'condition' or 'group').")
}

cat("✓ Object loaded\n")
cat("Sample column:", sample_col, "\n")
cat("Condition column:", condition_col, "\n")

# ============================================
# Prepare data for scCODA
# ============================================
cat("\n=== Preparing data for scCODA ===\n")

# Per-sample counts (samples x cell types)
count_matrix <- seurat_clean@meta.data %>%
  dplyr::group_by(
    sample_id = .data[[sample_col]],
    condition = .data[[condition_col]],
    celltype_final
  ) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = celltype_final,
    values_from = n_cells,
    values_fill = 0
  )

cat("Count matrix created\n")
cat("Samples:", nrow(count_matrix), "\n")
cat("Cell types:", ncol(count_matrix) - 2, "\n")

# Split into metadata + counts-only (scCODA expects separate files)
count_df <- as.data.frame(count_matrix)

sample_info <- count_df %>% dplyr::select(sample_id, condition)
counts_only <- count_df %>% dplyr::select(-sample_id, -condition)

# Simplify condition names for Python
sample_info$condition_code <- dplyr::case_when(
  sample_info$condition == "Healthy Control" ~ "HC",
  sample_info$condition == "UC Self-Control" ~ "UCSC",
  sample_info$condition == "Ulcerative Colitis" ~ "UC",
  TRUE ~ as.character(sample_info$condition)
)

write.csv(counts_only, file.path(table_dir, "sccoda_counts.csv"), row.names = FALSE)
write.csv(sample_info, file.path(table_dir, "sccoda_metadata.csv"), row.names = FALSE)

cat("✓ scCODA input files saved:\n")
cat("  -", file.path(table_dir, "sccoda_counts.csv"), "\n")
cat("  -", file.path(table_dir, "sccoda_metadata.csv"), "\n")
