# ============================================
# COMPOSITIONAL ANALYSIS
# ============================================

cat("\n=== COMPOSITIONAL ANALYSIS ===\n")

# Note: seurat_clean already has standardized celltype_final names
# (Colonocytes* and Colonocytes)

# ============================================
# Step 1: Calculate cell counts per sample
# ============================================

count_matrix <- seurat_clean@meta.data %>%
  dplyr::group_by(sample_id, condition, celltype_final) %>%
  dplyr::summarise(n_cells = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = celltype_final,
    values_from = n_cells,
    values_fill = 0
  )

cat("Count matrix created\n")
cat("Samples:", nrow(count_matrix), "\n")
cat("Cell types:", ncol(count_matrix) - 2, "\n")

print(count_matrix)

# Save count matrix
write.csv(count_matrix, file.path(table_dir, "cell_counts_per_sample.csv"), row.names = FALSE)

# ============================================
# Step 2: Calculate proportions
# ============================================

# Get condition mapping for each sample (needed after complete())
sample_condition_map <- seurat_clean@meta.data %>%
  dplyr::distinct(sample_id, condition)

composition_data <- seurat_clean@meta.data %>%
  dplyr::count(sample_id, celltype_final, name = "n_cells") %>%
  tidyr::complete(sample_id, celltype_final, fill = list(n_cells = 0)) %>%
  dplyr::left_join(sample_condition_map, by = "sample_id") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells,
    percentage = proportion * 100
  ) %>%
  dplyr::ungroup()

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
    n_samples = dplyr::n_distinct(sample_id),
    mean_percentage = mean(percentage, na.rm = TRUE),
    sd_percentage = sd(percentage, na.rm = TRUE),
    median_percentage = median(percentage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(celltype_final, condition)

cat("\n--- Composition Summary (Descriptive) ---\n")
print(composition_summary, n = 60)

write.csv(composition_summary, file.path(table_dir, "composition_summary_by_condition.csv"), row.names = FALSE)

# ============================================
# Step 4: Calculate fold changes
# ============================================

composition_wide <- composition_summary %>%
  dplyr::select(condition, celltype_final, mean_percentage) %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_percentage)

fold_changes <- composition_wide %>%
  dplyr::mutate(
    FC_UCSC_vs_HC = `UC Self-Control` / `Healthy Control`,
    FC_UC_vs_HC = `Ulcerative Colitis` / `Healthy Control`,
    FC_UC_vs_UCSC = `Ulcerative Colitis` / `UC Self-Control`,
    log2FC_UCSC_vs_HC = log2(FC_UCSC_vs_HC),
    log2FC_UC_vs_HC = log2(FC_UC_vs_HC),
    log2FC_UC_vs_UCSC = log2(FC_UC_vs_UCSC)
  ) %>%
  dplyr::arrange(desc(abs(log2FC_UC_vs_HC)))

cat("\n--- Fold Changes (sorted by |log2FC| UC vs HC) ---\n")
print(fold_changes %>% dplyr::select(celltype_final, `Healthy Control`, `Ulcerative Colitis`, FC_UC_vs_HC, log2FC_UC_vs_HC))

write.csv(fold_changes, file.path(table_dir, "composition_fold_changes.csv"), row.names = FALSE)

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
print(colonocyte_fc %>% dplyr::select(celltype_final, `Healthy Control`, `Ulcerative Colitis`, FC_UC_vs_HC, log2FC_UC_vs_HC))
