#!/usr/bin/env Rscript

# =============================================================================
# 09_clustering_qc.R
#
# Purpose:
#   1) Summarize QC metrics by cluster (res = 0.4) as a diagnostic
#   2) Save QC summary table
#   3) Create violin plots of key QC metrics by cluster:
#       - percent.mt
#       - nFeature_RNA
#       - nCount_RNA
#      and save as a combined figure
#   4) Compare Cluster 4 vs all other clusters (medians + Wilcoxon tests)
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_clustered.rds
#
# Outputs:
#   results/02_clustering_analysis/tables/qc_metrics_by_cluster.csv
#   results/02_clustering_analysis/plots/qc_metrics_by_cluster.png
#
# Notes:
#   - Assumes active clustering exists at RNA_snn_res.0.4
#   - Assumes meta.data includes: nFeature_RNA, nCount_RNA, percent.mt
# =============================================================================

# =========================
# Libraries (required)
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "02_clustering_analysis", "objects", "seurat_integrated_clustered.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 07_clustering.R first (or ensure the clustered object is saved).")
}

out_dir    <- file.path(root_dir, "results", "02_clustering_analysis")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# Load object
# =========================
cat("\nLoading clustered Seurat object...\n")
seurat_integrated <- readRDS(in_obj)

# Basic sanity checks
needed_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "RNA_snn_res.0.4")
missing_cols <- setdiff(needed_cols, colnames(seurat_integrated@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Set active identities to res 0.4 (required for WhichCells(idents=...))
Idents(seurat_integrated) <- "RNA_snn_res.0.4"

# ============================================
# QC METRICS BY CLUSTER - DIAGNOSTIC
# ============================================
cat("\n=== Checking QC metrics by cluster ===\n")

# Summary statistics per cluster
qc_by_cluster <- seurat_integrated@meta.data %>%
  dplyr::group_by(RNA_snn_res.0.4) %>%
  dplyr::summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_pct_mt = median(percent.mt),
    mean_pct_mt = mean(percent.mt),
    pct_high_mt = sum(percent.mt > 10) / n() * 100,  # % cells with >10% mito
    .groups = "drop"
  ) %>%
  dplyr::arrange(as.numeric(as.character(RNA_snn_res.0.4)))

cat("\nQC metrics by cluster:\n")
print(qc_by_cluster, n = 25)

write.csv(
  qc_by_cluster,
  file.path(table_dir, "qc_metrics_by_cluster.csv"),
  row.names = FALSE
)

cat("✓ QC metrics table saved\n")

# ============================================
# Violin plots - QC by cluster
# ============================================
cat("\nCreating QC violin plots by cluster...\n")

p_mt_by_cluster <- VlnPlot(
  seurat_integrated,
  features = "percent.mt",
  group.by = "RNA_snn_res.0.4",
  pt.size = 0
) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
  ggtitle("Mitochondrial content by cluster") +
  xlab("Cluster") +
  ylab("Mitochondrial %") +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

p_nfeature_by_cluster <- VlnPlot(
  seurat_integrated,
  features = "nFeature_RNA",
  group.by = "RNA_snn_res.0.4",
  pt.size = 0
) +
  ggtitle("Genes detected by cluster") +
  xlab("Cluster") +
  ylab("nFeature_RNA") +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

p_ncount_by_cluster <- VlnPlot(
  seurat_integrated,
  features = "nCount_RNA",
  group.by = "RNA_snn_res.0.4",
  pt.size = 0
) +
  ggtitle("UMI counts by cluster") +
  xlab("Cluster") +
  ylab("nCount_RNA") +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

# Combined plot
p_qc_combined <- p_mt_by_cluster / p_nfeature_by_cluster / p_ncount_by_cluster

ggsave(
  file.path(plot_dir, "qc_metrics_by_cluster.png"),
  p_qc_combined, width = 14, height = 16, dpi = 600
)

cat("✓ QC by cluster plots saved\n")

# ============================================
# Statistical comparison: Cluster 4 vs others
# ============================================
cat("\n=== Cluster 4 vs all other clusters ===\n")

cluster4_cells <- WhichCells(seurat_integrated, idents = "4")
other_cells <- setdiff(colnames(seurat_integrated), cluster4_cells)

if (length(cluster4_cells) == 0) {
  cat("No cells found for cluster 4. Skipping cluster-4 comparison.\n")
} else {

  comparison_stats <- data.frame(
    metric = c("percent.mt", "nFeature_RNA", "nCount_RNA"),
    cluster4_median = c(
      median(seurat_integrated$percent.mt[cluster4_cells]),
      median(seurat_integrated$nFeature_RNA[cluster4_cells]),
      median(seurat_integrated$nCount_RNA[cluster4_cells])
    ),
    others_median = c(
      median(seurat_integrated$percent.mt[other_cells]),
      median(seurat_integrated$nFeature_RNA[other_cells]),
      median(seurat_integrated$nCount_RNA[other_cells])
    )
  )

  comparison_stats$fold_change <- comparison_stats$cluster4_median / comparison_stats$others_median

  cat("\nCluster 4 vs all others:\n")
  print(comparison_stats)

  # Wilcoxon tests
  cat("\nStatistical tests (Wilcoxon):\n")
  cat(
    "percent.mt p-value:",
    wilcox.test(seurat_integrated$percent.mt[cluster4_cells],
                seurat_integrated$percent.mt[other_cells])$p.value, "\n"
  )
  cat(
    "nFeature_RNA p-value:",
    wilcox.test(seurat_integrated$nFeature_RNA[cluster4_cells],
                seurat_integrated$nFeature_RNA[other_cells])$p.value, "\n"
  )
  cat(
    "nCount_RNA p-value:",
    wilcox.test(seurat_integrated$nCount_RNA[cluster4_cells],
                seurat_integrated$nCount_RNA[other_cells])$p.value, "\n"
  )
}

cat("\n✓✓✓ CLUSTER QC COMPLETE ✓✓✓\n\n")
