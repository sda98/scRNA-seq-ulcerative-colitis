#!/usr/bin/env Rscript

# =============================================================================
# 09_clustering_qc_cluster4_investigation.R
#
# Purpose (defensible QC investigation of a problematic cluster):
#   1) Summarize QC metrics per cluster (res=0.4) to flag suspicious clusters.
#   2) Evaluate Cluster 4 distribution across samples (composition).
#   3) Compare Cluster 4 vs all other cells using per-sample paired medians
#      (defensible; avoids cell-level p-values / pseudoreplication).
#   4) Visualize:
#       - Pie chart: Cluster 4 composition across samples (white background).
#       - Spaghetti plots: per-sample medians (Cluster 4 vs Others) with p-values.
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_clustered.rds
#
# Outputs (tables):
#   results/02_clustering_analysis/tables/qc_by_cluster_res0.4.csv
#   results/02_clustering_analysis/tables/cluster4_sample_composition.csv
#   results/02_clustering_analysis/tables/cluster4_vs_others_by_sample.csv
#   results/02_clustering_analysis/tables/cluster4_vs_others_paired_tests.csv
#
# Outputs (plots):
#   results/02_clustering_analysis/plots/cluster4_composition_pie.png
#   results/02_clustering_analysis/plots/qc_cluster4_vs_others_by_sample.png
#
# Notes:
#   - Requires meta.data columns: nFeature_RNA, nCount_RNA, percent.mt
#   - Requires clustering column: RNA_snn_res.0.4
#   - Uses sample_id (preferred) or gsm_id for per-sample pairing
#   - Reports raw p-values (no multiple testing correction; small targeted QC set)
# =============================================================================

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(
  root_dir, "results", "02_clustering_analysis", "objects",
  "seurat_integrated_clustered.rds"
)
if (!file.exists(in_obj)) {
  stop(
    "Input object not found: ", in_obj,
    "\nRun 07_clustering.R first (save seurat_integrated_clustered.rds)."
  )
}

out_dir   <- file.path(root_dir, "results", "02_clustering_analysis")
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_09.txt"))

# =========================
# Load object
# =========================
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  CLUSTER 4 QC INVESTIGATION (DEFENSIBLE: SAMPLE-LEVEL)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

seurat_obj <- readRDS(in_obj)

required_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "RNA_snn_res.0.4")
missing <- setdiff(required_cols, colnames(seurat_obj@meta.data))
if (length(missing) > 0) {
  stop("Missing metadata columns: ", paste(missing, collapse = ", "))
}

# Set identities
Idents(seurat_obj) <- "RNA_snn_res.0.4"
DefaultAssay(seurat_obj) <- "RNA"

# Confirm cluster 4 exists
if (!"4" %in% levels(Idents(seurat_obj))) {
  stop("Cluster '4' not found in RNA_snn_res.0.4 identities.")
}

# Sample column for pairing
sample_col <- if ("sample_id" %in% colnames(seurat_obj@meta.data)) {
  "sample_id"
} else if ("gsm_id" %in% colnames(seurat_obj@meta.data)) {
  "gsm_id"
} else {
  stop("No sample identifier found (expected 'sample_id' or 'gsm_id').")
}

cat("Sample column used for pairing:", sample_col, "\n")
cat("Total cells:", ncol(seurat_obj), "\n")
cat("Total clusters (res=0.4):", length(unique(Idents(seurat_obj))), "\n\n")

# Cluster 4 flag (for plotting)
seurat_obj$cluster4_flag <- ifelse(as.character(Idents(seurat_obj)) == "4", "Cluster 4", "Other clusters")
seurat_obj$cluster4_flag <- factor(seurat_obj$cluster4_flag, levels = c("Other clusters", "Cluster 4"))

# =========================
# Plot theme (your style)
# =========================
theme_pub <- function(base_size = 30) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", size = 34),
      axis.text  = element_text(size = 28),
      legend.title = element_text(face = "bold", size = 26),
      legend.text  = element_text(size = 20),
      plot.title   = element_text(face = "bold", size = 34)
    )
}

# =============================================================================
# 1) QC SUMMARY BY CLUSTER
# =============================================================================
cat("=== 1) QC Summary by Cluster (res=0.4) ===\n")

global_nFeature <- median(seurat_obj$nFeature_RNA, na.rm = TRUE)
global_nCount   <- median(seurat_obj$nCount_RNA,   na.rm = TRUE)
global_mt       <- median(seurat_obj$percent.mt,   na.rm = TRUE)

qc_summary <- seurat_obj@meta.data %>%
  dplyr::mutate(cluster = as.character(RNA_snn_res.0.4)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    pct_total = 100 * n_cells / ncol(seurat_obj),
    median_nFeature = median(nFeature_RNA, na.rm = TRUE),
    median_nCount   = median(nCount_RNA,   na.rm = TRUE),
    median_pct_mt   = median(percent.mt,   na.rm = TRUE),
    pct_mt_gt10     = 100 * mean(percent.mt > 10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    nFeature_pct_diff = 100 * (median_nFeature - global_nFeature) / global_nFeature,
    nCount_pct_diff   = 100 * (median_nCount   - global_nCount)   / global_nCount,
    mt_diff           = median_pct_mt - global_mt
  ) %>%
  dplyr::arrange(as.numeric(cluster))

print(qc_summary, n = 50)

write.csv(qc_summary, file.path(table_dir, "qc_by_cluster_res0.4.csv"), row.names = FALSE)
cat("✓ Saved qc_by_cluster_res0.4.csv\n\n")

c4 <- qc_summary %>% dplyr::filter(cluster == "4")
cat("CLUSTER 4 SUMMARY:\n")
cat("  Cells:", c4$n_cells, "(", round(c4$pct_total, 2), "% of dataset)\n")
cat("  Median nFeature:", c4$median_nFeature,
    sprintf("(%+.1f%% vs global median %d)", c4$nFeature_pct_diff, global_nFeature), "\n")
cat("  Median nCount  :", c4$median_nCount,
    sprintf("(%+.1f%% vs global median %d)", c4$nCount_pct_diff, global_nCount), "\n")
cat("  Median %MT     :", round(c4$median_pct_mt, 3),
    sprintf("(Δ %+.3f vs global median %.3f)", c4$mt_diff, global_mt), "\n\n")

# =============================================================================
# 2) SAMPLE COMPOSITION (Cluster 4)
# =============================================================================
cat("=== 2) Cluster 4 Sample Composition ===\n")

composition <- seurat_obj@meta.data %>%
  dplyr::mutate(
    cluster = as.character(RNA_snn_res.0.4),
    sample  = .data[[sample_col]]
  ) %>%
  dplyr::group_by(sample, cluster) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(pct_of_sample = 100 * n / sum(n)) %>%
  dplyr::ungroup()

c4_comp <- composition %>%
  dplyr::filter(cluster == "4") %>%
  dplyr::arrange(desc(n))

write.csv(c4_comp, file.path(table_dir, "cluster4_sample_composition.csv"), row.names = FALSE)
cat("✓ Saved cluster4_sample_composition.csv\n")

cat("Cluster 4 cells by sample:\n")
print(c4_comp, n = 50)

c4_total <- sum(c4_comp$n)
top_sample_contribution <- 100 * max(c4_comp$n) / c4_total
cat("\nCluster 4 is present in", nrow(c4_comp), "samples.\n")
cat("Top sample contributes", round(top_sample_contribution, 1), "% of Cluster 4 cells.\n")

if (top_sample_contribution > 80) {
  cat("⚠ Note: Cluster 4 is highly concentrated in a single sample (possible batch/sample-specific effect).\n\n")
} else {
  cat("✓ Cluster 4 is distributed across samples (less likely to be single-sample artifact).\n\n")
}

# =============================================================================
# 2A) PIE CHART: Cluster 4 composition by sample (white background, readable labels)
# =============================================================================
cat("Creating pie chart...\n")

sample_plot_col <- if ("sample_label" %in% colnames(seurat_obj@meta.data)) "sample_label" else sample_col

c4_comp_plot <- seurat_obj@meta.data %>%
  dplyr::mutate(
    cluster = as.character(RNA_snn_res.0.4),
    sample_plot = .data[[sample_plot_col]]
  ) %>%
  dplyr::filter(cluster == "4") %>%
  dplyr::count(sample_plot, name = "n") %>%
  dplyr::mutate(pct_of_c4 = 100 * n / sum(n)) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(
    label = ifelse(pct_of_c4 >= 5, paste0(sample_plot, "\n", round(pct_of_c4, 1), "%"), "")
  )

base_sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"
)
samples_unique <- c4_comp_plot$sample_plot
sample_colors <- rep(base_sample_colors, length.out = length(samples_unique))
names(sample_colors) <- samples_unique

.hex_luma <- function(hex) {
  rgb <- grDevices::col2rgb(hex) / 255
  0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
}
c4_comp_plot$fill <- sample_colors[c4_comp_plot$sample_plot]
c4_comp_plot$label_color <- ifelse(.hex_luma(c4_comp_plot$fill) < 0.55, "white", "black")

p_pie <- ggplot(c4_comp_plot, aes(x = "", y = n, fill = sample_plot)) +
  geom_col(width = 1, color = "black", linewidth = 0.6) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label, color = label_color),
    position = position_stack(vjust = 0.5),
    size = 5,
    fontface = "bold",
    lineheight = 0.9,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = sample_colors) +
  scale_color_identity() +
  labs(fill = "Sample") +
  theme_void(base_size = 30) +
  theme(
    legend.title = element_text(face = "bold", size = 26),
    legend.text  = element_text(size = 20),
    legend.position = "right",
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(plot_dir, "cluster4_composition_pie.png"),
  p_pie,
  width = 14, height = 10, dpi = 600,
  bg = "white"
)
cat("✓ Saved cluster4_composition_pie.png\n\n")

# =============================================================================
# 3) PER-SAMPLE PAIRED COMPARISON (DEFENSIBLE)
# =============================================================================
cat("=== 3) Per-sample paired comparison (Cluster 4 vs Others) ===\n")

min_c4_cells <- 20

by_sample <- seurat_obj@meta.data %>%
  dplyr::mutate(
    sample = .data[[sample_col]],
    group  = ifelse(as.character(RNA_snn_res.0.4) == "4", "C4", "Other")
  ) %>%
  dplyr::group_by(sample, group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    med_nFeature = median(nFeature_RNA, na.rm = TRUE),
    med_nCount   = median(nCount_RNA,   na.rm = TRUE),
    med_mt       = median(percent.mt,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = group, values_from = c(n, med_nFeature, med_nCount, med_mt)) %>%
  dplyr::filter(!is.na(n_C4), n_C4 >= min_c4_cells)

cat("Samples with >=", min_c4_cells, "Cluster 4 cells:", nrow(by_sample), "\n")

by_sample <- by_sample %>%
  dplyr::mutate(
    diff_nFeature = med_nFeature_C4 - med_nFeature_Other,
    diff_nCount   = med_nCount_C4   - med_nCount_Other,
    diff_mt       = med_mt_C4       - med_mt_Other,
    pct_diff_nFeature = 100 * diff_nFeature / med_nFeature_Other,
    pct_diff_nCount   = 100 * diff_nCount   / med_nCount_Other
  )

write.csv(by_sample, file.path(table_dir, "cluster4_vs_others_by_sample.csv"), row.names = FALSE)
cat("✓ Saved cluster4_vs_others_by_sample.csv\n")

paired_test_results <- NULL
if (nrow(by_sample) >= 4) {

  test_nf <- wilcox.test(by_sample$med_nFeature_C4, by_sample$med_nFeature_Other, paired = TRUE)
  test_nc <- wilcox.test(by_sample$med_nCount_C4,   by_sample$med_nCount_Other,   paired = TRUE)
  test_mt <- wilcox.test(by_sample$med_mt_C4,       by_sample$med_mt_Other,       paired = TRUE)

  pvals <- c(nFeature_RNA = test_nf$p.value, nCount_RNA = test_nc$p.value, percent.mt = test_mt$p.value)

  effects <- c(
    nFeature_RNA = median(by_sample$diff_nFeature, na.rm = TRUE),
    nCount_RNA   = median(by_sample$diff_nCount,   na.rm = TRUE),
    percent.mt   = median(by_sample$diff_mt,       na.rm = TRUE)
  )

  pct_effects <- c(
    nFeature_RNA = median(by_sample$pct_diff_nFeature, na.rm = TRUE),
    nCount_RNA   = median(by_sample$pct_diff_nCount,   na.rm = TRUE),
    percent.mt   = NA_real_
  )

  paired_test_results <- data.frame(
    metric = names(pvals),
    n_samples = nrow(by_sample),
    effect_median_delta = round(effects, 2),
    effect_median_pct   = round(pct_effects, 1),
    p_value = signif(pvals, 3),
    stringsAsFactors = FALSE
  )

  write.csv(
    paired_test_results,
    file.path(table_dir, "cluster4_vs_others_paired_tests.csv"),
    row.names = FALSE
  )
  cat("✓ Saved cluster4_vs_others_paired_tests.csv\n\n")

  cat("Paired Wilcoxon tests (per-sample medians; raw p-values):\n")
  print(paired_test_results)

} else {
  cat("Insufficient samples (need >=4) for paired Wilcoxon testing.\n")
}

# =============================================================================
# 4) PLOT: Paired spaghetti plot with p-values on facets
# =============================================================================
cat("\n=== Plotting paired spaghetti (per-sample medians) ===\n")

if (nrow(by_sample) >= 4) {

  spaghetti_df <- by_sample %>%
    dplyr::select(
      sample,
      med_nFeature_C4, med_nFeature_Other,
      med_nCount_C4,   med_nCount_Other,
      med_mt_C4,       med_mt_Other
    ) %>%
    tidyr::pivot_longer(-sample, names_to = "var", values_to = "value") %>%
    tidyr::separate(var, into = c("stat", "metric", "group"), sep = "_", extra = "merge") %>%
    dplyr::mutate(
      group = factor(group, levels = c("Other", "C4"), labels = c("Other clusters", "Cluster 4")),
      metric = factor(metric, levels = c("nFeature", "nCount", "mt"),
                      labels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
    )

  ann <- paired_test_results %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt")),
      label = dplyr::case_when(
        metric == "percent.mt" ~ paste0("p = ", p_value, "\nΔmedian = ", effect_median_delta),
        TRUE ~ paste0("p = ", p_value, "\nΔmedian = ", effect_median_delta,
                      " (", effect_median_pct, "%)")
      )
    ) %>%
    dplyr::select(metric, label)

  y_max <- spaghetti_df %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(y = max(value, na.rm = TRUE) * 1.05, .groups = "drop")

  ann <- dplyr::left_join(ann, y_max, by = "metric") %>%
    dplyr::mutate(x = 1.5)

  p_paired <- ggplot(spaghetti_df, aes(x = group, y = value, group = sample)) +
    geom_line(alpha = 0.55, linewidth = 1, color = "grey40") +
    geom_point(aes(color = group), size = 3.5) +
    scale_color_manual(values = c("Other clusters" = "grey50", "Cluster 4" = "firebrick")) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    geom_text(
      data = ann,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold",
      hjust = 0.5,
      vjust = 1
    ) +
    labs(x = NULL, y = "Per-sample median", color = NULL) +
    theme_pub() +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 30)
    )

  ggsave(
    file.path(plot_dir, "qc_cluster4_vs_others_by_sample.png"),
    p_paired,
    width = 22, height = 7, dpi = 600,
    bg = "white"
  )

  cat("✓ Saved qc_cluster4_vs_others_by_sample.png\n")

} else {
  cat("Not enough samples for paired spaghetti plot.\n")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  CLUSTER 4 QC INVESTIGATION COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
