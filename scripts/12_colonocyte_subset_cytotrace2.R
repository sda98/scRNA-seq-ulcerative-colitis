#!/usr/bin/env Rscript

# =============================================================================
# 12_colonocyte_subset_cytotrace2.R
#
# Purpose:
#   1) Load annotated Seurat object (manual labels already added).
#   2) Create a colonocyte subset using res=0.4 clusters 6 and 7.
#   3) Run CytoTRACE2 on colonocyte RNA counts (genes x cells) to obtain an
#      independent differentiation signal.
#   4) Summarize CytoTRACE2 results:
#        - Cell-level summary (descriptive): mean/median/SD by colonocyte cluster.
#        - Sample-level summary: per-sample medians for clusters 6 vs 7
#            (paired within the same sample when both clusters are present).
#   5) Visualize CytoTRACE2 results:
#        - Boxplot of per-cell CytoTRACE2 scores (Cluster 6 vs 7) with a sample-level
#            paired Wilcoxon p-value label.
#        - Spaghetti plot of per-sample medians (paired Cluster 6 vs 7) to show
#            per-sample direction consistency.
#   6) Save outputs:
#         - Tables (cell-level summary, per-sample medians, Wilcoxon test results).
#         - Plots (boxplot, spaghetti plot).
#         - Colonocyte subset Seurat object with CytoTRACE2 metadata.
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_annotated.rds
#
# Outputs:
#   results/02_clustering_analysis/sessionInfo.txt
#   results/02_clustering_analysis/tables/cytotrace2_by_colonocyte_cluster.csv
#   results/02_clustering_analysis/tables/cytotrace2_persample_medians.csv
#   results/02_clustering_analysis/tables/cytotrace2_wilcox_samplelevel.csv
#   results/02_clustering_analysis/plots/cytotrace2_boxplot.png
#   results/02_clustering_analysis/plots/cytotrace2_spaghetti_persample.png
#   results/02_clustering_analysis/objects/seurat_colonocytes_cytotrace2.rds
#
# Notes:
#   - Uses RNA assay counts as input; filters genes expressed in >= 1% of colonocyte cells.
# =============================================================================

# =========================
# Libraries 
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(CytoTRACE2)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "02_clustering_analysis", "objects",
                    "seurat_integrated_annotated.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 09_celltype_annotation_manual.R first (or ensure annotated object exists).")
}

out_dir    <- file.path(root_dir, "results", "02_clustering_analysis")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# 1) LOAD SEURAT OBJECT
# =========================
cat("\nLoading annotated Seurat object...\n")
seurat_integrated <- readRDS(in_obj)

if (!"RNA_snn_res.0.4" %in% colnames(seurat_integrated@meta.data)) {
  stop("Missing 'RNA_snn_res.0.4' in metadata.")
}
if (!"celltype_final" %in% colnames(seurat_integrated@meta.data)) {
  stop("Missing 'celltype_final' in metadata (manual annotations expected).")
}

cat("✓ Object loaded\n")

# =========================
# 2) SUBSET COLONOCYTES
# =========================
cat("\n=== Creating colonocyte subset ===\n")

# Check available cell types
cat("Available cell types:\n")
print(table(seurat_integrated$celltype_final))

# Detect sample and condition columns
sample_col <- if ("sample_id" %in% colnames(seurat_integrated@meta.data)) {
  "sample_id"
} else if ("gsm_id" %in% colnames(seurat_integrated@meta.data)) {
  "gsm_id"
} else {
  stop("No sample identifier found (expected 'sample_id' or 'gsm_id').")
}

condition_col <- if ("condition" %in% colnames(seurat_integrated@meta.data)) {
  "condition"
} else if ("group" %in% colnames(seurat_integrated@meta.data)) {
  "group"
} else {
  NULL
}

cat("Sample column:", sample_col, "\n")
cat("Condition column:", ifelse(is.null(condition_col), "not found", condition_col), "\n")

# Report sample counts
cat("\nSamples per condition:\n")
if (!is.null(condition_col)) {
  print(table(seurat_integrated@meta.data[[condition_col]]))
}

# Subset colonocytes (clusters 6 and 7)
Idents(seurat_integrated) <- "RNA_snn_res.0.4"
seurat_colonocytes <- subset(seurat_integrated, idents = c("6", "7"))

# Add colonocyte type labels
seurat_colonocytes$colonocyte_type <- ifelse(
  seurat_colonocytes$RNA_snn_res.0.4 == "6",
  "Colonocytes (less mature)",
  "Colonocytes (more mature)"
)

cat("\nColonocyte subset created\n")
cat("Total colonocytes:", ncol(seurat_colonocytes), "\n")
cat("Cluster 6 (less mature):", sum(seurat_colonocytes$colonocyte_type == "Colonocytes (less mature)"), "\n")
cat("Cluster 7 (more mature):", sum(seurat_colonocytes$colonocyte_type == "Colonocytes (more mature)"), "\n")

# Check samples with colonocytes
cat("\nSamples contributing colonocytes:\n")
print(table(seurat_colonocytes@meta.data[[sample_col]]))

# ============================================
# 3) CytoTRACE2 ANALYSIS
# ============================================
cat("\n=== CytoTRACE2 Analysis ===\n")

# Get expression matrix for colonocytes (genes x cells)
expr <- GetAssayData(seurat_colonocytes, assay = "RNA", layer = "counts")

# Keep genes expressed in at least 1% of cells
keep_genes <- Matrix::rowSums(expr > 0) >= ceiling(0.01 * ncol(expr))
expr <- expr[keep_genes, ]

cat("Genes before filtering:", nrow(GetAssayData(seurat_colonocytes, assay = "RNA", layer = "counts")), "\n")
cat("Genes after filtering:", nrow(expr), "\n")

# Densify after filtering
expr_matrix <- as.matrix(expr)

cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "cells\n")

# Run CytoTRACE2
cat("\nRunning CytoTRACE2 (this may take a few minutes)...\n")

set.seed(42)
cytotrace2_result <- cytotrace2(expr_matrix, species = "human", ncores = 4)

cat("✓ CytoTRACE2 complete\n")

# Check result structure
cat("\nResult structure:\n")
print(names(cytotrace2_result))

# ============================================
# 4) EXTRACT AND ANALYZE CytoTRACE2 SCORES
# ============================================
cat("\n=== Analyzing CytoTRACE2 Results ===\n")

# Extract scores
cytotrace_scores  <- cytotrace2_result$CytoTRACE2_Score
cytotrace_potency <- cytotrace2_result$CytoTRACE2_Potency

# Ensure scores match input matrix column order
if (!is.null(names(cytotrace_scores))) {
  cytotrace_scores  <- cytotrace_scores[colnames(expr_matrix)]
  cytotrace_potency <- cytotrace_potency[colnames(expr_matrix)]
}
stopifnot(length(cytotrace_scores) == ncol(seurat_colonocytes))

cat("Score structure:\n")
print(head(cytotrace_scores))
cat("\nLength of scores:", length(cytotrace_scores), "\n")
cat("Number of colonocytes:", ncol(seurat_colonocytes), "\n")

# Add scores to Seurat object
seurat_colonocytes$CytoTRACE2_Score   <- cytotrace_scores
seurat_colonocytes$CytoTRACE2_Potency <- cytotrace_potency

cat("✓ Scores added to Seurat object\n")

# ============================================
# 5) SUMMARY STATISTIC BY CLUSTER
# ============================================
cat("\n--- CytoTRACE2 Scores by Cluster ---\n")

cytotrace_by_cluster <- seurat_colonocytes@meta.data %>%
  group_by(colonocyte_type) %>%
  dplyr::summarise(
    n_cells      = dplyr::n(),
    mean_score   = mean(CytoTRACE2_Score, na.rm = TRUE),
    median_score = median(CytoTRACE2_Score, na.rm = TRUE),
    sd_score     = sd(CytoTRACE2_Score, na.rm = TRUE),
    .groups      = "drop"
  )

print(cytotrace_by_cluster)

write.csv(cytotrace_by_cluster, 
          file.path(table_dir, "cytotrace2_by_colonocyte_cluster.csv"), 
          row.names = FALSE)

# ============================================
# 6) Per-sample medians for statistical test
# ============================================
cat("\n--- Computing per-sample medians ---\n")

meta_ct <- seurat_colonocytes@meta.data

# Per-sample medians within each cluster
idx6 <- meta_ct$RNA_snn_res.0.4 == "6"
idx7 <- meta_ct$RNA_snn_res.0.4 == "7"

med6 <- tapply(meta_ct$CytoTRACE2_Score[idx6], meta_ct[[sample_col]][idx6], median, na.rm = TRUE)
med7 <- tapply(meta_ct$CytoTRACE2_Score[idx7], meta_ct[[sample_col]][idx7], median, na.rm = TRUE)

# Samples with cells in both clusters
common_samples <- intersect(names(med6), names(med7))

cat("Samples in cluster 6:", length(med6), "\n")
cat("Samples in cluster 7:", length(med7), "\n")
cat("Samples in BOTH clusters:", length(common_samples), "\n")

# Build per-sample data frame
sample_pairs <- data.frame(
  sample = common_samples,
  cluster6_median = as.numeric(med6[common_samples]),
  cluster7_median = as.numeric(med7[common_samples])
)

# Add condition if available
if (!is.null(condition_col)) {
  sample_condition <- meta_ct %>%
    select(all_of(c(sample_col, condition_col))) %>%
    distinct()
  colnames(sample_condition) <- c("sample", "condition")
  sample_pairs <- left_join(sample_pairs, sample_condition, by = "sample")
}

# Calculate difference (Cluster 6 - Cluster 7)
sample_pairs$difference <- sample_pairs$cluster6_median - sample_pairs$cluster7_median

# Direction check
sample_pairs$direction <- ifelse(sample_pairs$difference > 0, 
                                 "Cluster 6 > Cluster 7", 
                                 "Cluster 7 > Cluster 6")

cat("\nPer-sample medians:\n")
print(sample_pairs)

write.csv(sample_pairs, 
          file.path(table_dir, "cytotrace2_persample_medians.csv"), 
          row.names = FALSE)

# Report consistency
n_consistent <- sum(sample_pairs$difference > 0)
cat("\nConsistency check:\n")
cat("  Samples where Cluster 6 > Cluster 7:", n_consistent, "of", nrow(sample_pairs), "\n")
cat("  Samples where Cluster 7 > Cluster 6:", nrow(sample_pairs) - n_consistent, "of", nrow(sample_pairs), "\n")

# ============================================
# 7) Statistical test (paired Wilcoxon)
# ============================================
cat("\n--- Statistical Comparison (paired Wilcoxon, sample-level) ---\n")

if (length(common_samples) >= 4) {
  wilcox_result <- wilcox.test(
    x      = sample_pairs$cluster6_median,
    y      = sample_pairs$cluster7_median,
    paired = TRUE,
    exact  = FALSE
  )
  test_note <- paste0("Paired Wilcoxon signed-rank test on per-sample medians (n = ", 
                      length(common_samples), " samples)")
} else {
  wilcox_result <- wilcox.test(
    x      = as.numeric(med6),
    y      = as.numeric(med7),
    paired = FALSE,
    exact  = FALSE
  )
  test_note <- paste0("Unpaired Wilcoxon test on per-sample medians (n6 = ", 
                      length(med6), ", n7 = ", length(med7), ")")
}

cat(test_note, "\n")
cat("p-value:", wilcox_result$p.value, "\n")

# Effect size: median of differences
median_diff <- median(sample_pairs$difference)
cat("Median difference (Cluster 6 - Cluster 7):", round(median_diff, 4), "\n")

# Save test results
test_results <- data.frame(
  test                    = test_note,
  n_samples               = length(common_samples),
  p_value                 = wilcox_result$p.value,
  median_cluster6         = median(sample_pairs$cluster6_median),
  median_cluster7         = median(sample_pairs$cluster7_median),
  median_difference       = median_diff,
  n_cluster6_higher       = n_consistent,
  n_cluster7_higher       = nrow(sample_pairs) - n_consistent
)

write.csv(test_results, 
          file.path(table_dir, "cytotrace2_wilcox_samplelevel.csv"), 
          row.names = FALSE)

# Build p-value label for plot
p_label <- if (is.na(wilcox_result$p.value)) {
  "p = NA"
} else if (wilcox_result$p.value < 0.001) {
  paste0("p = ", formatC(wilcox_result$p.value, format = "e", digits = 2))
} else {
  paste0("p = ", round(wilcox_result$p.value, 3))
}

# ============================================
# 8) Boxplot visualization
# ============================================
cat("\n=== Creating visualizations ===\n")

Idents(seurat_colonocytes) <- "RNA_snn_res.0.4"

# Cell-level medians for display
median_6 <- median(seurat_colonocytes$CytoTRACE2_Score[Idents(seurat_colonocytes) == "6"], na.rm = TRUE)
median_7 <- median(seurat_colonocytes$CytoTRACE2_Score[Idents(seurat_colonocytes) == "7"], na.rm = TRUE)
n_6 <- sum(Idents(seurat_colonocytes) == "6")
n_7 <- sum(Idents(seurat_colonocytes) == "7")

# Create plot data frame
plot_df <- data.frame(
  cluster   = as.character(Idents(seurat_colonocytes)),
  cytotrace = seurat_colonocytes$CytoTRACE2_Score
)
plot_df$x_pos <- ifelse(plot_df$cluster == "6", 1, 3)

# Max for annotation positioning
y_max <- max(plot_df$cytotrace, na.rm = TRUE)

# CytoTRACE2 Boxplot
p_cytotrace_box <- ggplot(plot_df, aes(x = x_pos, y = cytotrace, fill = cluster, group = cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey40") +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, linewidth = 1, width = 0.6) +
  scale_fill_manual(values = c("6" = "#FFFF33", "7" = "#A65628")) +
  scale_x_continuous(
    breaks = c(1, 3),
    labels = c("Cluster 6", "Cluster 7"),
    limits = c(0, 4.5)
  ) +
  labs(
    title = NULL,
    x     = NULL,
    y     = "CytoTRACE2 Score\n(higher = less differentiated)"
  ) +
  # Median labels
  
  annotate("text", x = 1.38, y = median_6,
           label = paste0("median = ", round(median_6, 3)),
           size = 8, fontface = "bold", hjust = 0) +
  annotate("text", x = 3.38, y = median_7,
           label = paste0("median = ", round(median_7, 3)),
           size = 8, fontface = "bold", hjust = 0) +
  # N labels
  annotate("text", x = 1, y = -0.03,
           label = paste0("n = ", format(n_6, big.mark = ",")),
           size = 10, fontface = "bold") +
  annotate("text", x = 3, y = -0.03,
           label = paste0("n = ", format(n_7, big.mark = ",")),
           size = 10, fontface = "bold") +
  # Significance bar
  annotate("segment", x = 1, xend = 3,
           y = y_max * 1.0 - 0.03, yend = y_max * 1.0 - 0.03,
           linewidth = 1) +
  annotate("text", x = 2, y = y_max * 1.05,
           label = paste0(p_label, "\n(paired, n=", length(common_samples), " samples)"),
           size = 10, fontface = "bold") +
  coord_cartesian(ylim = c(-0.07, y_max * 1.15), clip = "off") +
  theme_classic(base_size = 32) +
  theme(
    axis.title.y     = element_text(size = 32, face = "bold", margin = margin(r = 15)),
    axis.text.x      = element_text(face = "bold", size = 28),
    axis.text.y      = element_text(size = 38),
    axis.ticks       = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm"),
    legend.position  = "none",
    plot.margin      = margin(t = 20, r = 100, b = 20, l = 20, unit = "pt")
  )

ggsave(file.path(plot_dir, "cytotrace2_boxplot.png"), 
       p_cytotrace_box, 
       width = 14, height = 10, dpi = 600, bg = "white")

cat("✓ CytoTRACE2 boxplot saved\n")

# ============================================
# 9) Spaghetti plot (per-sample paired comparison)
# ============================================
cat("\n--- Creating spaghetti plot ---\n")

# Reshape for plotting
spaghetti_df <- sample_pairs %>%
  pivot_longer(
    cols      = c(cluster6_median, cluster7_median),
    names_to  = "cluster",
    values_to = "cytotrace_median"
  ) %>%
  mutate(
    cluster_label = ifelse(cluster == "cluster6_median", "Cluster 6", "Cluster 7"),
    x_pos         = ifelse(cluster == "cluster6_median", 1, 2)
  )

# Define condition colors if available
if ("condition" %in% colnames(spaghetti_df)) {
  condition_colors <- c(
    "Ulcerative Colitis" = "#E41A1C",
    "UC Self-Control"    = "#377EB8",
    "Healthy Control"    = "#4DAF4A"
  )
  # Fallback for any missing conditions
  unique_conditions <- unique(spaghetti_df$condition)
  missing_conditions <- setdiff(unique_conditions, names(condition_colors))
  if (length(missing_conditions) > 0) {
    extra_colors <- setNames(scales::hue_pal()(length(missing_conditions)), missing_conditions)
    condition_colors <- c(condition_colors, extra_colors)
  }
}

p_spaghetti <- ggplot(spaghetti_df, aes(x = x_pos, y = cytotrace_median, group = sample)) +
  geom_line(aes(color = if ("condition" %in% colnames(spaghetti_df)) condition else NULL),
            linewidth = 1.2, alpha = 0.8) +
  geom_point(aes(color = if ("condition" %in% colnames(spaghetti_df)) condition else NULL),
             size = 5) +
  {if ("condition" %in% colnames(spaghetti_df)) 
    scale_color_manual(values = condition_colors, name = "Condition")} +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Cluster 6", "Cluster 7"),
    limits = c(0.5, 2.5)
  ) +
  labs(
    title    = "CytoTRACE2 scores: per-sample comparison",
    subtitle = paste0(n_consistent, " of ", nrow(sample_pairs), 
                      " samples show Cluster 6 > Cluster 7 (Wilcoxin ", p_label, ")"),
    x        = NULL,
    y        = "CytoTRACE2 Score\n(per sample median)"
  ) +
  theme_classic(base_size = 24) +
  theme(
    plot.title       = element_text(size = 28, face = "bold"),
    plot.subtitle    = element_text(size = 20, face = "italic", color = "grey30"),
    axis.title.y     = element_text(size = 24, face = "bold", margin = margin(r = 15)),
    axis.text.x      = element_text(face = "bold", size = 28),
    axis.text.y      = element_text(size = 30),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 22),
    legend.text      = element_text(size = 19),
    legend.key.spacing.y = unit(0.3, units = "cm")
  )

ggsave(file.path(plot_dir, "cytotrace2_spaghetti_persample.png"),
       p_spaghetti,
       width = 12, height = 9, dpi = 600, bg = "white")

cat("✓ Spaghetti plot saved\n")

# ============================================
# 10) Save Seurat object
# ============================================
saveRDS(seurat_colonocytes, file.path(object_dir, "seurat_colonocytes_cytotrace2.rds"))
cat("✓ Colonocyte Seurat object with CytoTRACE2 scores saved\n")

# ============================================
# Summary
# ============================================
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("Total colonocytes analyzed:", ncol(seurat_colonocytes), "\n")
cat("  - Cluster 6 (less mature):", n_6, "\n")
cat("  - Cluster 7 (more mature):", n_7, "\n")
cat("\n")
cat("Samples with cells in both clusters:", length(common_samples), "\n")
if (!is.null(condition_col)) {
  cat("  - By condition:\n")
  print(table(sample_pairs$condition))
}
cat("\n")
cat("CytoTRACE2 validation:\n")
cat("  - Cluster 6 median (cell-level):", round(median_6, 4), "\n")
cat("  - Cluster 7 median (cell-level):", round(median_7, 4), "\n")
cat("  - Expected direction (Cluster 6 > Cluster 7):", 
    ifelse(median_6 > median_7, "CONFIRMED", "NOT CONFIRMED"), "\n")
cat("\n")
cat("Statistical test:\n")
cat("  -", test_note, "\n")
cat("  - p-value:", wilcox_result$p.value, "\n")
cat("  - Consistency:", n_consistent, "of", nrow(sample_pairs), "samples\n")
cat("\n")
cat("=== Script complete ===\n")
