#!/usr/bin/env Rscript

# =============================================================================
# 12_colonocyte_subset_cytotrace2.R
#
# Purpose:
#   1) Load annotated Seurat object (manual labels already added).
#   2) Subset colonocytes (clusters 6 and 7 at res=0.4).
#   3) Run CytoTRACE2 on colonocyte counts (genes x cells) as an independent
#      differentiation validation signal.
#   4) Save:
#        - Colonocyte subset object with CytoTRACE2 scores
#        - Summary table + stats
#        - Boxplot comparing Cluster 6 vs 7 CytoTRACE2 scores
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_annotated.rds
#
# Outputs:
#   results/04_colonocyte_cytotrace2/sessionInfo.txt
#   results/04_colonocyte_cytotrace2/tables/cytotrace2_by_colonocyte_cluster.csv
#   results/04_colonocyte_cytotrace2/tables/cytotrace2_wilcox_celllevel.csv
#   results/04_colonocyte_cytotrace2/plots/cytotrace2_boxplot_cluster6_vs_7.png
#   results/04_colonocyte_cytotrace2/objects/seurat_colonocytes_cytotrace2.rds
#
# Notes:
#   - Uses RNA assay counts as input; filters genes expressed in >=1% of cells.
#   - CytoTRACE2 returns per-cell scores; these are added to the subset metadata.
#   - Cell-level Wilcoxon is shown as a descriptive test (cells are not independent).
#     If you want a defensible inferential test, do sample-level pairing instead.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  # CytoTRACE2 package must be installed/loaded elsewhere; keep explicit if needed:
  # library(CytoTRACE2)
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

out_dir    <- file.path(root_dir, "results", "04_colonocyte_cytotrace2")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# 1) Load object
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
# 2) Subset colonocytes (clusters 6 and 7)
# =========================
cat("\n=== Creating colonocyte subset ===\n")
cat("Available cell types:\n")
print(table(seurat_integrated$celltype_final, useNA = "ifany"))

Idents(seurat_integrated) <- "RNA_snn_res.0.4"
seurat_colonocytes <- subset(seurat_integrated, idents = c("6", "7"))

seurat_colonocytes$colonocyte_type <- ifelse(
  seurat_colonocytes$RNA_snn_res.0.4 == "6",
  "Colonocytes*",
  "Colonocytes"
)

cat("\nColonocyte subset created\n")
cat("Total colonocytes:", ncol(seurat_colonocytes), "\n")
cat("Colonocytes* (cluster 6):", sum(seurat_colonocytes$colonocyte_type == "Colonocytes*"), "\n")
cat("Colonocytes  (cluster 7):", sum(seurat_colonocytes$colonocyte_type == "Colonocytes"), "\n")

# =========================
# 3) CytoTRACE2
# =========================
cat("\n=== CytoTRACE2 Analysis ===\n")

DefaultAssay(seurat_colonocytes) <- "RNA"

expr <- GetAssayData(seurat_colonocytes, assay = "RNA", layer = "counts")

# Filter genes expressed in >=1% of cells (reduces memory + noise)
keep_genes <- Matrix::rowSums(expr > 0) >= ceiling(0.01 * ncol(expr))
cat("Genes before filtering:", nrow(expr), "\n")
expr <- expr[keep_genes, , drop = FALSE]
cat("Genes after filtering :", nrow(expr), "\n")

# Densify AFTER filtering
expr_matrix <- as.matrix(expr)
cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "cells\n")

cat("\nRunning CytoTRACE2 (this may take a few minutes)...\n")
set.seed(42)
cytotrace2_result <- cytotrace2(expr_matrix, species = "human", ncores = 4)
cat("✓ CytoTRACE2 complete\n")

cat("\nResult structure:\n")
print(names(cytotrace2_result))

# Extract scores
cat("\n=== Extracting CytoTRACE2 Results ===\n")
cytotrace_scores  <- cytotrace2_result$CytoTRACE2_Score
cytotrace_potency <- cytotrace2_result$CytoTRACE2_Potency

# Align to cell order (safety)
if (!is.null(names(cytotrace_scores))) {
  cytotrace_scores  <- cytotrace_scores[colnames(expr_matrix)]
  cytotrace_potency <- cytotrace_potency[colnames(expr_matrix)]
}
stopifnot(length(cytotrace_scores) == ncol(seurat_colonocytes))

seurat_colonocytes$CytoTRACE2_Score   <- cytotrace_scores
seurat_colonocytes$CytoTRACE2_Potency <- cytotrace_potency

cat("✓ Scores added to colonocyte subset\n")

# =========================
# 4) Summaries + stats
# =========================
cat("\n=== CytoTRACE2 Scores by Colonocyte Type ===\n")

cytotrace_by_cluster <- seurat_colonocytes@meta.data %>%
  group_by(colonocyte_type) %>%
  summarise(
    n_cells = n(),
    mean_score   = mean(CytoTRACE2_Score, na.rm = TRUE),
    median_score = median(CytoTRACE2_Score, na.rm = TRUE),
    sd_score     = sd(CytoTRACE2_Score, na.rm = TRUE),
    .groups = "drop"
  )

print(cytotrace_by_cluster)

write.csv(
  cytotrace_by_cluster,
  file.path(table_dir, "cytotrace2_by_colonocyte_cluster.csv"),
  row.names = FALSE
)

# Cell-level Wilcoxon (descriptive; not strictly independent)
immature_scores <- seurat_colonocytes$CytoTRACE2_Score[seurat_colonocytes$colonocyte_type == "Colonocytes*"]
mature_scores   <- seurat_colonocytes$CytoTRACE2_Score[seurat_colonocytes$colonocyte_type == "Colonocytes"]

wilcox_result <- wilcox.test(immature_scores, mature_scores)

wilcox_tbl <- data.frame(
  test = "Wilcoxon rank-sum (cell-level; descriptive)",
  p_value = wilcox_result$p.value,
  mean_Colonocytes_star = mean(immature_scores, na.rm = TRUE),
  mean_Colonocytes      = mean(mature_scores, na.rm = TRUE),
  median_Colonocytes_star = median(immature_scores, na.rm = TRUE),
  median_Colonocytes      = median(mature_scores, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cat("\n--- Statistical Comparison (descriptive) ---\n")
print(wilcox_tbl)

write.csv(
  wilcox_tbl,
  file.path(table_dir, "cytotrace2_wilcox_celllevel.csv"),
  row.names = FALSE
)

# =========================
# 5) Plot (boxplot; your style)
# =========================
cat("\n=== Plotting CytoTRACE2 boxplot ===\n")

Idents(seurat_colonocytes) <- "RNA_snn_res.0.4"

median_6 <- median(seurat_colonocytes$CytoTRACE2_Score[Idents(seurat_colonocytes) == "6"], na.rm = TRUE)
median_7 <- median(seurat_colonocytes$CytoTRACE2_Score[Idents(seurat_colonocytes) == "7"], na.rm = TRUE)
n_6 <- sum(Idents(seurat_colonocytes) == "6")
n_7 <- sum(Idents(seurat_colonocytes) == "7")

plot_df <- data.frame(
  cluster   = as.character(Idents(seurat_colonocytes)),
  cytotrace = seurat_colonocytes$CytoTRACE2_Score,
  stringsAsFactors = FALSE
)
plot_df$x_pos <- ifelse(plot_df$cluster == "6", 1, 3)

y_max <- max(plot_df$cytotrace, na.rm = TRUE)

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
    x = NULL,
    y = "CytoTRACE2 Score\n(higher = less differentiated)"
  ) +
  annotate("text", x = 1.38, y = median_6, label = paste0("median = ", round(median_6, 3)),
           size = 8, fontface = "bold", hjust = 0) +
  annotate("text", x = 3.38, y = median_7, label = paste0("median = ", round(median_7, 3)),
           size = 8, fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = -0.03, label = paste0("n = ", format(n_6, big.mark = ",")),
           size = 10, fontface = "bold") +
  annotate("text", x = 3, y = -0.03, label = paste0("n = ", format(n_7, big.mark = ",")),
           size = 10, fontface = "bold") +
  annotate("segment", x = 1, xend = 3, y = y_max * 1.0, yend = y_max * 1.0, linewidth = 1) +
  annotate("text", x = 2, y = y_max * 1.05,
           label = ifelse(wilcox_result$p.value < 2.2e-16, "p < 2.2e-16",
                          paste0("p = ", signif(wilcox_result$p.value, 3))),
           size = 12, fontface = "bold") +
  coord_cartesian(ylim = c(-0.07, y_max * 1.12), clip = "off") +
  theme_classic(base_size = 32) +
  theme(
    axis.title.y = element_text(size = 32, face = "bold", margin = margin(r = 15)),
    axis.text.x  = element_text(face = "bold", size = 32),
    axis.text.y  = element_text(size = 38),
    axis.ticks   = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm"),
    legend.position = "none",
    plot.margin = margin(t = 20, r = 100, b = 20, l = 20, unit = "pt")
  )

ggsave(
  file.path(plot_dir, "cytotrace2_boxplot_cluster6_vs_7.png"),
  p_cytotrace_box, width = 14, height = 10, dpi = 600
)

cat("✓ Saved cytotrace2_boxplot_cluster6_vs_7.png\n")

# =========================
# 6) Save subset object
# =========================
out_obj <- file.path(object_dir, "seurat_colonocytes_cytotrace2.rds")
saveRDS(seurat_colonocytes, out_obj)

cat("\n✓ Saved colonocyte subset with CytoTRACE2 metadata: ", out_obj, "\n")
cat("\n✓✓✓ CytoTRACE2 colonocyte validation complete ✓✓✓\n\n")

