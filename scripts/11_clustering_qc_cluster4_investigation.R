#!/usr/bin/env Rscript

# =============================================================================
# 11_clustering_qc_cluster4_investigation.R
#
# Purpose:
#   Investigate whether Cluster 4 (RNA_snn_res.0.4) represents low-quality cells
#   that should be removed before downstream analysis.
#
# Overview (what this script does):
#   1) Load clustered Seurat object (Harmony-integrated; res = 0.4).
#   2) Summarize QC metrics per cluster to flag suspicious clusters.
#   3) Read precomputed top-20 marker table and subset Cluster 4 markers.
#      - Categorize Cluster 4 markers by functional type and export table.
#      - Plot Cluster 4 top markers by type.
#   4) Assess Cluster 4 distribution across samples (composition) and plot a pie chart.
#   5) Compare Cluster 4 vs all other cells using per-sample paired medians
#      (sample-based paired Wilcoxin test; avoids cell-level p-values / pseudoreplication).
#      - Export per-sample table and paired Wilcoxon test results.
#      - Plot bar comparison (mean of per-sample medians ± SD) with p-values.
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_clustered.rds
#   results/02_clustering_analysis/tables/top20_markers_per_cluster_res0.4.csv
#
# Outputs (tables):
#   results/02_clustering_analysis/sessionInfo_11.txt
#   results/02_clustering_analysis/tables/qc_by_cluster_res0.4.csv
#   results/02_clustering_analysis/tables/cluster4_top20_markers_by_type.csv
#   results/02_clustering_analysis/tables/cluster4_vs_others_by_sample.csv
#   results/02_clustering_analysis/tables/cluster4_vs_others_paired_tests.csv
#
# Outputs (plots):
#   results/02_clustering_analysis/plots/cluster4_top20_markers_by_type.png
#   results/02_clustering_analysis/plots/cluster4_composition_pie.png
#   results/02_clustering_analysis/plots/qc_cluster4_vs_others_by_sample.png
# =============================================================================

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggh4x)
  library(ggtext)
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

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_11.txt"))

# =========================
# 1) LOAD SEURAT OBJECT
# =========================
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  CLUSTER 4 QC INVESTIGATION\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

seurat_integrated <- readRDS(in_obj)

required_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "RNA_snn_res.0.4")
missing <- setdiff(required_cols, colnames(seurat_integrated@meta.data))
if (length(missing) > 0) {
  stop("Missing metadata columns: ", paste(missing, collapse = ", "))
}

Idents(seurat_integrated) <- "RNA_snn_res.0.4"
DefaultAssay(seurat_integrated) <- "RNA"

if (!"4" %in% levels(Idents(seurat_integrated))) {
  stop("Cluster '4' not found in RNA_snn_res.0.4 identities.")
}

sample_col <- if ("sample_id" %in% colnames(seurat_integrated@meta.data)) {
  "sample_id"
} else if ("gsm_id" %in% colnames(seurat_integrated@meta.data)) {
  "gsm_id"
} else {
  stop("No sample identifier found (expected 'sample_id' or 'gsm_id').")
}

cat("Sample column used for pairing:", sample_col, "\n")
cat("Total cells:", ncol(seurat_integrated), "\n")
cat("Total clusters (res=0.4):", length(unique(Idents(seurat_integrated))), "\n\n")

seurat_integrated$cluster4_flag <- ifelse(as.character(Idents(seurat_integrated)) == "4", "Cluster 4", "Other clusters")
seurat_integrated$cluster4_flag <- factor(seurat_integrated$cluster4_flag, levels = c("Other clusters", "Cluster 4"))

# =============================================================================
# 2) QC SUMMARY BY CLUSTER
# =============================================================================
cat("=== 1) QC Summary by Cluster (res=0.4) ===\n")

global_nFeature <- median(seurat_integrated$nFeature_RNA, na.rm = TRUE)
global_nCount   <- median(seurat_integrated$nCount_RNA,   na.rm = TRUE)
global_mt       <- median(seurat_integrated$percent.mt,   na.rm = TRUE)

qc_summary <- seurat_integrated@meta.data %>%
  dplyr::mutate(cluster = as.character(RNA_snn_res.0.4)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    pct_total = 100 * n_cells / ncol(seurat_integrated),
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
    sprintf("(%+.3f vs global median %.3f)", c4$mt_diff, global_mt), "\n\n")

# =============================================================================
# 3) CLUSTER 4 TOP MARKERS BAR PLOT - ANNOTATION
# =============================================================================
cat("=== 3) Cluster 4 top markers grouped by type ===\n")

markers_path <- file.path(table_dir, "top20_markers_per_cluster_res0.4.csv")
if (!file.exists(markers_path)) {
  stop("Missing markers file: ", markers_path)
}

markers_tbl <- read.csv(markers_path, stringsAsFactors = FALSE)

top_c4 <- markers_tbl %>%
  filter(cluster == 4) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  mutate(
    type = case_when(
      # ---- Immune ----
      gene %in% c(
        "LTB", "CORO1A", "RAC2", "KLRB1", "CD52", "ARHGDIB",
        "CD3D", "CD3E", "CD3G", "CD2", "CD7",
        "IL7R", "IL32", "PTPRC",
        "TRAC", "TRBC1", "TRBC2"
      ) ~ "Immune",
      
      # ---- Ubiquitous: Cytoskeleton ----
      gene %in% c(
        "ACTB", "ACTG1", "PFN1", "ARPC3", "COTL1",
        "TMSB4X", "TMSB10", "CFL1", "MYL6"
      ) ~ "Cytoskeleton (ubiquitous)",
      
      # ---- Ubiquitous: RNA Translation ----
      grepl("^EIF", gene) |
        grepl("^RPL", gene) |
        grepl("^RPS", gene) |
        gene %in% c(
          "SNRPD2", "SNRPG", "SNRPE",
          "EEF1A1", "EEF2", "PABPC1"
        ) ~ "RNA Translation (ubiquitous)",
      
      # ---- Ubiquitous: Metabolism ----
      gene %in% c(
        "LDHB", "LDHA", "GAPDH", "PKM",
        "ENO1", "PGAM1", "ATP5F1E", "ATP5MC3"
      ) ~ "Metabolism (ubiquitous)",
      
      # ---- Catch-all ----
      TRUE ~ "Other / unclear"
    ),
    type = factor(type, levels = c(
      "Immune",
      "Cytoskeleton (ubiquitous)",
      "RNA Translation (ubiquitous)",
      "Metabolism (ubiquitous)",
      "Other / unclear"
    ))
  )

write.csv(top_c4, file.path(table_dir, "cluster4_top20_markers_by_type.csv"), row.names = FALSE)

# ---- Summary stats ----
cat("\nCluster 4 marker composition:\n")
type_summary <- top_c4 %>% 
  dplyr::count(type) %>% 
  dplyr::mutate(pct = round(100 * n / sum(n), 1))
print(type_summary)

# Colors
type_cols <- c(
  "Immune"                       = "#F8766D",
  "Cytoskeleton (ubiquitous)"    = "#A3A500",
  "RNA Translation (ubiquitous)"     = "#00BF7D",
  "Metabolism (ubiquitous)"      = "#00B0F6",
  "Other / unclear"              = "#E76BF3"
)

# Order genes within each type by lfc
top_c4 <- top_c4 %>%
  arrange(type, avg_log2FC) %>%
  mutate(gene = factor(gene, levels = gene))

# Panel heights proportional to gene count
panel_heights <- top_c4 %>% count(type) %>% pull(n)

# Plot
p_markers <- ggplot(top_c4, aes(x = gene, y = avg_log2FC, fill = type)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 2.5) +
  geom_col(color = "black", linewidth = 0.4, width = 0.8) +
  coord_flip() +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = type_cols) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5), expand = c(0, 0)) +
  labs(title = "Cluster 4 top-20 markers", x = NULL, y = expression(bold("Average Log"[2]*"FC"))) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 60, face = "bold", margin = margin(b = 40, unit = "pt")),
    axis.text.x = element_text(size = 52, color = "black"),
    axis.text.y = element_text(size = 28, face = "bold.italic", color = "black"),
    axis.title.x = element_text(size = 40, face = "bold", margin = margin(t = 30, unit = "pt")),
    strip.text = element_text(size = 33, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.6),
    legend.position = "none"
  ) +
  ggh4x::force_panelsizes(rows = panel_heights) +
  ggh4x::facetted_pos_scales(y = list(scale_y_continuous(limits = c(0, 1.55), breaks = seq(0, 1.5, 0.5))))

ggsave(
  file.path(plot_dir, "cluster4_top20_markers_by_type.png"),
  p_markers,
  width = 13, height = 20, dpi = 600,
  bg = "white"
)

cat("✓ Saved cluster4_top_markers_by_type.png\n")

# =============================================================================
# 4) PIE CHART: Cluster 4 composition by sample
# =============================================================================
cat("Creating pie chart...\n")

sample_plot_col <- if ("sample_label" %in% colnames(seurat_integrated@meta.data)) "sample_label" else sample_col

c4_comp_plot <- seurat_integrated@meta.data %>%
  mutate(
    cluster = as.character(RNA_snn_res.0.4),
    sample_plot = .data[[sample_plot_col]]
  ) %>%
  filter(cluster == "4") %>%
  count(sample_plot, name = "n") %>%
  mutate(pct_of_c4 = 100 * n / sum(n)) %>%
  arrange(desc(n)) %>%
  mutate(
    sample_plot = factor(sample_plot, levels = sample_plot),
    legend_label = paste0(sample_plot, " (", round(pct_of_c4, 1), "%)"),
    legend_label = factor(legend_label, levels = legend_label)
  )

base_sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"
)
samples_unique <- levels(c4_comp_plot$sample_plot)
sample_colors <- setNames(
  rep(base_sample_colors, length.out = length(samples_unique)),
  samples_unique
)

p_pie <- ggplot(c4_comp_plot, aes(x = "", y = n, fill = sample_plot)) +
  geom_col(width = 1, color = "black", linewidth = 0.6) +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = sample_colors,
    labels = paste0(c4_comp_plot$sample_plot, " (**", round(c4_comp_plot$pct_of_c4, 1), "%**)"),
    breaks = c4_comp_plot$sample_plot
  ) +
  labs(
    title = "Cluster 4 sample composition",
    fill = "Sample"
  ) +
  theme_void(base_size = 30) +
  theme(
    plot.title = element_text(face = "bold", size = 50, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 35),
    legend.text = ggtext::element_markdown(size = 25),
    legend.key.size = unit(1, "cm"),
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(byrow = TRUE))

ggsave(
  file.path(plot_dir, "cluster4_composition_pie.png"),
  p_pie,
  width = 16, height = 12, dpi = 600,
  bg = "white"
)
cat("✓ Saved cluster4_composition_pie.png\n\n")

# =============================================================================
# 5) PER-SAMPLE PAIRED COMPARISON (DEFENSIBLE)
# =============================================================================
cat("=== 3) Per-sample paired comparison (Cluster 4 vs Others) ===\n")

min_c4_cells <- 20

by_sample <- seurat_integrated@meta.data %>%
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
# 6) PLOT: Bar plot comparing Cluster 4 vs Others (per-sample medians)
# =============================================================================
cat("\n=== Plotting bar comparison (per-sample medians) ===\n")
if (nrow(by_sample) >= 4) {
  
  # Calculate mean of per-sample medians for each group
  bar_df <- data.frame(
    metric = factor(c("nFeature_RNA", "nCount_RNA", "percent.mt",
                      "nFeature_RNA", "nCount_RNA", "percent.mt"),
                    levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    labels = c("nFeature_RNA", "nCount_RNA", "Mitochondrial reads (%)")),
    group = factor(rep(c("Other clusters", "Cluster 4"), each = 3),
                   levels = c("Other clusters", "Cluster 4")),
    value = c(
      mean(by_sample$med_nFeature_Other),
      mean(by_sample$med_nCount_Other),
      mean(by_sample$med_mt_Other),
      mean(by_sample$med_nFeature_C4),
      mean(by_sample$med_nCount_C4),
      mean(by_sample$med_mt_C4)
    ),
    sd = c(
      sd(by_sample$med_nFeature_Other),
      sd(by_sample$med_nCount_Other),
      sd(by_sample$med_mt_Other),
      sd(by_sample$med_nFeature_C4),
      sd(by_sample$med_nCount_C4),
      sd(by_sample$med_mt_C4)
    )
  )
  
  # Annotation with paired sample-based Wilcoxin p-values
  ann <- paired_test_results %>%
    mutate(
      metric = factor(metric, 
                      levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      labels = c("nFeature_RNA", "nCount_RNA", "Mitochondrial reads (%)")),
      label = case_when(
        metric == "Mitochondrial reads (%)" ~ paste0("p = ", p_value, "\nΔ = ", effect_median_delta),
        TRUE ~ paste0("p = ", p_value, "\nΔ = ", effect_median_delta, " (", effect_median_pct, "%)")
      )
    )
  
  # Y position for annotations
  y_max <- bar_df %>%
    group_by(metric) %>%
    summarise(y = max(value + sd) * 1.15, .groups = "drop")
  
  ann <- left_join(ann, y_max, by = "metric")
  
  p_bar <- ggplot(bar_df, aes(x = group, y = value, fill = group)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    geom_col(color = "black", linewidth = 0.5, width = 0.7) +
    geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.2, linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("Other clusters" = "grey70", "Cluster 4" = "firebrick")) +
    geom_text(
      data = ann,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold",
      hjust = 0.5,
      vjust = 0.8
    ) +
    labs(x = NULL, y = "Mean of per-sample medians (± SD)") +
    theme_classic(base_size = 24) +
    theme(
      axis.title.y = element_text(face = "bold", size = 22),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 24),
      axis.text.y = element_text(size = 25),
      strip.text = element_text(face = "bold", size = 26),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.6),
      legend.position = "none"
    )
  
  ggsave(
    file.path(plot_dir, "qc_cluster4_vs_others_by_sample.png"),
    p_bar,
    width = 16, height = 8, dpi = 600,
    bg = "white"
  )
  
cat("✓ Saved qc_cluster4_vs_others_by_sample.png\n")
  
} else {
  cat("Not enough samples for bar plot.\n")
}
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  CLUSTER 4 QC INVESTIGATION COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
