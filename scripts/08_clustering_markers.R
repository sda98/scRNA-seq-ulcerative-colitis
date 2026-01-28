#!/usr/bin/env Rscript

# =============================================================================
# 08_clustering_markers.R
#
# Purpose:
#   1) Identify cluster marker genes using FindAllMarkers (Wilcoxon)
#   2) Filter out nuisance genes (mitochondrial, ribosomal, hemoglobin)
#   3) Export:
#       - all markers (filtered)
#       - top 5 / top 10 / top 20 markers per cluster
#       - per-cluster marker summary stats
#
# Inputs:
#   results/02_clustering_analysis/objects/seurat_integrated_clustered.rds
#
# Outputs:
#   results/02_clustering_analysis/tables/all_markers_res0.4.csv
#   results/02_clustering_analysis/tables/top5_markers_per_cluster_res0.4.csv
#   results/02_clustering_analysis/tables/top10_markers_per_cluster_res0.4.csv
#   results/02_clustering_analysis/tables/top20_markers_per_cluster_res0.4.csv
#   results/02_clustering_analysis/tables/marker_summary_per_cluster_res0.4.csv
#
# Notes:
#   - Uses clustering resolution 0.4 (RNA_snn_res.0.4)
#   - Assay set to RNA for marker detection
#   - For Seurat v5: joins RNA layers if multiple layers exist
# =============================================================================

# =========================
# Libraries 
# =========================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
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
if (!"RNA" %in% names(seurat_integrated@assays)) {
  stop("RNA assay not found in object.")
}
if (!"RNA_snn_res.0.4" %in% colnames(seurat_integrated@meta.data)) {
  stop("Clustering column 'RNA_snn_res.0.4' not found in metadata.\n",
       "Make sure clustering at res=0.4 exists.")
}

# ============================================
# 1. FIND CLUSTER MARKERS
# ============================================
cat("\n=== FINDING CLUSTER MARKERS ===\n")
cat("Using resolution 0.4 for marker identification\n\n")

# Set resolution 0.4 as active identity
Idents(seurat_integrated) <- "RNA_snn_res.0.4"

# Switch to RNA assay for marker detection
DefaultAssay(seurat_integrated) <- "RNA"

# Join RNA layers for Seurat v5
if (exists("Layers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
  lyr <- SeuratObject::Layers(seurat_integrated[["RNA"]])
  if (length(lyr) > 1) {
    cat("Joining RNA layers...\n")
    seurat_integrated <- JoinLayers(seurat_integrated)
    cat("✓ RNA layers joined\n")
  } else {
    cat("✓ No layer joining needed\n")
  }
} else {
  cat("✓ SeuratObject::Layers() not available; skipping JoinLayers\n")
}

cat("Finding markers for", length(unique(Idents(seurat_integrated))), "clusters...\n\n")

# Find markers for all clusters
all_markers <- FindAllMarkers(
  seurat_integrated,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.75,
  return.thresh = 0.01,
  assay = "RNA",
  verbose = TRUE
)

cat("\n✓ Marker detection complete\n")
cat("Total markers found:", nrow(all_markers), "\n")

# ============================================
# 2. REMOVE NUISANCE GENES
# ============================================

cat("Filtering out ribosomal, mitochondrial, and hemoglobin genes...\n")
all_markers_filtered <- all_markers %>%
  dplyr::filter(
    !grepl("^MT-", gene),
    !grepl("^RPL", gene),
    !grepl("^RPS", gene),
    !grepl("^HBA|^HBB", gene)
  )

cat("Markers after filtering:", nrow(all_markers_filtered), "\n\n")

# ============================================
# 3. EXPORTING RESULTS
# ============================================

# Save all markers
write.csv(
  all_markers_filtered,
  file.path(table_dir, "all_markers_res0.4.csv"),
  row.names = FALSE
)
cat("✓ All markers saved\n")

# Top 5 markers per cluster
top5 <- all_markers_filtered %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

write.csv(
  top5,
  file.path(table_dir, "top5_markers_per_cluster_res0.4.csv"),
  row.names = FALSE
)
cat("✓ Top 5 markers saved\n\n")

# Print top 5 markers
cat("Top 5 markers per cluster (by avg_log2FC):\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
print(top5 %>% dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj), n = 100)

# Top 10 markers per cluster
top10 <- all_markers_filtered %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

write.csv(
  top10,
  file.path(table_dir, "top10_markers_per_cluster_res0.4.csv"),
  row.names = FALSE
)
cat("\n✓ Top 10 markers also saved for reference\n")

# Top 20 markers per cluster
top20 <- all_markers_filtered %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  dplyr::ungroup()

write.csv(
  top20,
  file.path(table_dir, "top20_markers_per_cluster_res0.4.csv"),
  row.names = FALSE
)
cat("\n✓ Top 20 markers also saved for reference\n")

# Summary statistics per cluster
marker_summary <- all_markers_filtered %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_markers = n(),
    mean_log2FC = mean(avg_log2FC),
    max_log2FC = max(avg_log2FC),
    .groups = "drop"
  )

cat("\nMarker summary per cluster:\n")
print(marker_summary, n = 50)

write.csv(
  marker_summary,
  file.path(table_dir, "marker_summary_per_cluster_res0.4.csv"),
  row.names = FALSE
)

cat("\n✓✓✓ MARKER IDENTIFICATION COMPLETE ✓✓✓\n\n")
