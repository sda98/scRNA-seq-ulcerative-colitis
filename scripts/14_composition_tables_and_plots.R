#!/usr/bin/env Rscript

# =============================================================================
# 14_composition_tables_and_plots.R
#
# Purpose:
#   1) Load final cleaned/annotated Seurat object.
#   2) Compute per-sample cell counts and proportions by celltype_final.
#   3) Produce descriptive summaries by condition and fold-change tables.
#   4) Export compositional analysis tables.
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_clean.rds
#     - Must contain meta.data columns:
#         - celltype_final
#         - sample_id (or gsm_id)
#         - condition (or group)
#
# Outputs:
#   results/03_compositional_analysis/sessionInfo.txt
#   results/03_compositional_analysis/tables/cell_counts_per_sample.csv
#   results/03_compositional_analysis/tables/composition_summary_by_condition.csv
#   results/03_compositional_analysis/tables/composition_fold_changes.csv
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
  "seurat_clean.rds"
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
# COMPOSITIONAL ANALYSIS
# ============================================
cat("\n=== COMPOSITIONAL ANALYSIS ===\n")

# ============================================
# Step 1: Calculate cell counts per sample
# ============================================
count_matrix <- seurat_clean@meta.data %>%
  dplyr::group_by(
    sample_id   = .data[[sample_col]],
    condition   = .data[[condition_col]],
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
print(count_matrix)

write.csv(count_matrix, file.path(table_dir, "cell_counts_per_sample.csv"), row.names = FALSE)

# ============================================
# Step 2: Calculate proportions
# ============================================
sample_condition_map <- seurat_clean@meta.data %>%
  dplyr::distinct(
    sample_id = .data[[sample_col]],
    condition = .data[[condition_col]]
  )

composition_data <- seurat_clean@meta.data %>%
  dplyr::count(sample_id = .data[[sample_col]], celltype_final, name = "n_cells") %>%
  tidyr::complete(sample_id, celltype_final, fill = list(n_cells = 0)) %>%
  dplyr::left_join(sample_condition_map, by = "sample_id") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    total_cells = sum(n_cells),
    proportion  = n_cells / total_cells,
    percentage  = proportion * 100
  ) %>%
  dplyr::ungroup()

# Keep your preferred condition order (only applied if those levels exist)
composition_data$condition <- factor(
  composition_data$condition,
  levels = c("Healthy Control", "UC Self-Control", "Ulcerative Colitis")
)

# ============================================
# Step 3: Descriptive summary
# ============================================
composition_summary <- composition_data %>%
  dplyr::group_by(condition, celltype_final) %>%
  dplyr::summarise(
    n_samples         = dplyr::n_distinct(sample_id),
    mean_percentage   = mean(percentage, na.rm = TRUE),
    sd_percentage     = sd(percentage, na.rm = TRUE),
    median_percentage = median(percentage, na.rm = TRUE),
    .groups           = "drop"
  ) %>%
  dplyr::arrange(celltype_final, condition)

cat("\n--- Composition Summary (Descriptive) ---\n")
print(composition_summary, n = 60)

write.csv(
  composition_summary,
  file.path(table_dir, "composition_summary_by_condition.csv"),
  row.names = FALSE
)

# ============================================
# Step 4: Calculate fold changes
# ============================================
composition_wide <- composition_summary %>%
  dplyr::select(condition, celltype_final, mean_percentage) %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_percentage)

fold_changes <- composition_wide %>%
  dplyr::mutate(
    FC_UCSC_vs_HC    = `UC Self-Control`    / `Healthy Control`,
    FC_UC_vs_HC      = `Ulcerative Colitis` / `Healthy Control`,
    FC_UC_vs_UCSC    = `Ulcerative Colitis` / `UC Self-Control`,
    log2FC_UCSC_vs_HC = log2(FC_UCSC_vs_HC),
    log2FC_UC_vs_HC   = log2(FC_UC_vs_HC),
    log2FC_UC_vs_UCSC = log2(FC_UC_vs_UCSC)
  ) %>%
  dplyr::arrange(desc(abs(log2FC_UC_vs_HC)))

cat("\n--- Fold Changes (sorted by |log2FC| UC vs HC) ---\n")
print(
  fold_changes %>%
    dplyr::select(celltype_final, `Healthy Control`, `Ulcerative Colitis`, FC_UC_vs_HC, log2FC_UC_vs_HC)
)

write.csv(
  fold_changes,
  file.path(table_dir, "composition_fold_changes.csv"),
  row.names = FALSE
)

# ============================================
# Step 5: Colonocyte-specific summary
# ============================================
cat("\n--- Colonocyte Composition ---\n")

colonocyte_summary <- composition_summary %>%
  dplyr::filter(celltype_final %in% c("Colonocytes*", "Colonocytes"))

print(colonocyte_summary)

colonocyte_fc <- fold_changes %>%
  dplyr::filter(celltype_final %in% c("Colonocytes*", "Colonocytes"))

cat("\n--- Colonocyte Fold Changes ---\n")
print(
  colonocyte_fc %>%
    dplyr::select(celltype_final, `Healthy Control`, `Ulcerative Colitis`, FC_UC_vs_HC, log2FC_UC_vs_HC)
)

cat("\n✓✓✓ COMPOSITIONAL ANALYSIS COMPLETE ✓✓✓\n\n")
