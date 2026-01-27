#!/usr/bin/env Rscript

# =============================================================================
# 02_qc_metrics_and_plots.R
#
# Purpose:
#   - Load Seurat object from Script 01
#   - Plot QC metrics BEFORE filtering (per sample)
#   - Apply robust per-sample lower-bound filters (MAD-based) + mito cutoff
#   - Plot QC metrics AFTER filtering (per sample)
#   - Scatter plot (nCount vs nFeature) faceted by sample with cutoff lines
#   - Save filtered Seurat object for downstream scripts
#
# Inputs:
#   results/01_quality_control/objects/seurat_raw.rds
#
# Outputs:
#   results/01_quality_control/sessionInfo.txt
#   results/01_quality_control/plots/qc_violin_before_filtering.png
#   results/01_quality_control/plots/qc_violin_after_filtering.png
#   results/01_quality_control/plots/qc_scatter_by_sample.png
#   results/01_quality_control/tables/qc_bounds_by_sample.csv
#   results/01_quality_control/objects/seurat_qc_filtered.rds
# =============================================================================

# =========================
# Libraries (required)
# =========================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggh4x)
  library(scales)
  library(grid)   
})

# =========================
# Setting up directions and files
# =========================

root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "01_quality_control", "objects", "seurat_raw.rds")
if (!file.exists(in_obj)) {
  stop("Input Seurat object not found: ", in_obj,
       "\nRun 01_load_data_create_seurat.R first.")
}

out_dir    <- file.path(root_dir, "results", "01_quality_control")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
object_dir <- file.path(out_dir, "objects")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# =========================
# Load object
# =========================
cat("\nLoading Seurat object...\n")
seurat_obj <- readRDS(in_obj)

# Ensure QC metrics exist
if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}
if (!"percent.ribo" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
}

# Basic sanity
if (!"sample_id" %in% colnames(seurat_obj@meta.data)) {
  stop("metadata column 'sample_id' not found. Script 01 should create it.")
}

# =============================================================================
# 1) QC FILTERING - ROBUST BOUNDS METHOD
# =============================================================================
cat("\n=== QC FILTERING ===\n")
cat("Using per-sample robust bounds method\n")
cat("This accounts for technical variation between samples\n\n")

# QC parameters
mt_cutoff <- 15  # mitochondrial cutoff (%)
k <- 3           # MAD multiplier
q_low <- 0.001   # lower quantile threshold

robust_bounds <- function(x, k = 3, q = 0.001, min_scale = 1) {
  x <- x[is.finite(x)]
  med <- median(x)
  sc  <- mad(x)
  if (is.na(sc) || sc == 0) sc <- IQR(x) / 1.349
  if (is.na(sc) || sc == 0) sc <- sd(x)
  if (is.na(sc) || sc == 0) sc <- min_scale
  lower <- max(med - k * sc, as.numeric(quantile(x, q)))
  c(lower = lower, median = med, scale = sc)
}

# =============================================================================
# 2) QC PLOTS - BEFORE FILTERING
# =============================================================================
cat("\nGenerating QC plots before filtering...\n")

new_titles <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "Mitochondrial\ncontent (%)"
)

new_xlabs <- c("nCount_RNA", "nFeature_RNA", "percent.mt")

p_list <- VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0,
  combine = FALSE
)

for (i in seq_along(p_list)) {
  p_list[[i]] <- p_list[[i]] +
    ggtitle(new_titles[i]) +
    xlab(new_xlabs[i]) +
    labs(fill = "Sample") +
    theme_classic(base_size = 20) +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
      axis.ticks.x = element_blank(),
      plot.title   = element_text(size = 18, face = "bold"),
      legend.text  = element_text(face = "bold", size = 14)
    )
}

p_vln_before <- (wrap_plots(p_list, ncol = 3) + plot_layout(guides = "collect")) &
  theme(legend.position = "right")

p_vln_before <- p_vln_before +
  plot_annotation(
    title = "QC metrics before filtering (per sample)",
    theme = theme(plot.title = element_text(face = "bold", size = 24))
  )

ggsave(
  file.path(plot_dir, "qc_violin_before_filtering.png"),
  p_vln_before, width = 24, height = 5, dpi = 600
)

# =============================================================================
# 3) APPLY FILTERING BY SAMPLE
# =============================================================================
cat("\nApplying robust QC filters per sample...\n")
cat("Filters: lower bounds (MAD-based) + mito% cutoff\n\n")

cells_by_sample <- split(colnames(seurat_obj), seurat_obj$sample_id)

keep <- rep(TRUE, ncol(seurat_obj))
names(keep) <- colnames(seurat_obj)

bounds_df <- data.frame()

for (s in names(cells_by_sample)) {
  cells <- cells_by_sample[[s]]
  md <- seurat_obj@meta.data[cells, , drop = FALSE]
  
  b_feat <- robust_bounds(md$nFeature_RNA, k = k, q = q_low)
  b_cnt  <- robust_bounds(md$nCount_RNA,   k = k, q = q_low)
  
  keep[cells] <- (
    md$nFeature_RNA > b_feat["lower"] &
      md$nCount_RNA  > b_cnt["lower"] &
      md$percent.mt  < mt_cutoff
  )
  
  bounds_df <- rbind(bounds_df, data.frame(
    sample    = s,
    feat_lower = b_feat["lower"],
    cnt_lower  = b_cnt["lower"],
    mt_cutoff  = mt_cutoff,
    n_before   = length(cells),
    n_after    = sum(keep[cells]),
    pct_kept   = round(100 * sum(keep[cells]) / length(cells), 1)
  ))
}

write.csv(bounds_df, file.path(table_dir, "qc_bounds_by_sample.csv"), row.names = FALSE)

cat("\nQC filtering summary:\n")
cat("Total cells before:", ncol(seurat_obj), "\n")
cat("Total cells after :", sum(keep), "\n")
cat("Cells removed     :", ncol(seurat_obj) - sum(keep),
    "(", round(100 * (ncol(seurat_obj) - sum(keep)) / ncol(seurat_obj), 2), "%)\n")

cat("\nPer-sample filtering:\n")
print(bounds_df[, c("sample", "n_before", "n_after", "pct_kept")])

seurat_obj_filtered <- subset(seurat_obj, cells = names(keep)[keep])

cat("\nCells per condition after filtering:\n")
if ("condition_short" %in% colnames(seurat_obj_filtered@meta.data)) {
  print(table(seurat_obj_filtered$condition_short))
} else {
  cat("condition_short not found (ok, but you probably want it from Script 01).\n")
}

cat("\nCells per sample after filtering:\n")
print(table(seurat_obj_filtered$sample_id))

# =============================================================================
# 4) QC PLOTS - AFTER FILTERING
# =============================================================================
cat("\nGenerating QC plots after filtering...\n")

p_list <- VlnPlot(
  seurat_obj_filtered,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0,
  combine = FALSE
)

for (i in seq_along(p_list)) {
  p_list[[i]] <- p_list[[i]] +
    ggtitle(new_titles[i]) +
    xlab(new_xlabs[i]) +
    labs(fill = "Sample") +
    theme_classic(base_size = 20) +
    theme(
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
      axis.ticks.x = element_blank(),
      plot.title   = element_text(size = 18, face = "bold"),
      legend.text  = element_text(face = "bold", size = 14)
    )
}

p_vln_after <- (wrap_plots(p_list, ncol = 3) + plot_layout(guides = "collect")) &
  theme(legend.position = "right")

p_vln_after <- p_vln_after +
  plot_annotation(
    title = "QC metrics after filtering (per sample)",
    theme = theme(plot.title = element_text(face = "bold", size = 24))
  )

ggsave(
  file.path(plot_dir, "qc_violin_after_filtering.png"),
  p_vln_after, width = 24, height = 5, dpi = 600
)

# =============================================================================
# 5) SCATTER PLOT - PER SAMPLE (WITH THRESHOLD LINES)
# =============================================================================

# Add readable sample labels (optional but helps faceting titles)
# (kept on BOTH objects so downstream scripts can use it)
make_sample_label <- function(gsm_id, sample_id) {
  dplyr::case_when(
    gsm_id == "GSM7307094" ~ "UC Self-Control 1",
    gsm_id == "GSM7307095" ~ "UC Self-Control 2",
    gsm_id == "GSM7307096" ~ "UC Self-Control 3",
    gsm_id == "GSM7307097" ~ "UC Self-Control 4",
    gsm_id == "GSM7307098" ~ "Ulcerative Colitis 1",
    gsm_id == "GSM7307099" ~ "Ulcerative Colitis 2",
    gsm_id == "GSM7307100" ~ "Ulcerative Colitis 3",
    gsm_id == "GSM7307101" ~ "Ulcerative Colitis 4",
    gsm_id == "GSM7307102" ~ "Healthy Control 5",
    gsm_id == "GSM7307103" ~ "Healthy Control 6",
    gsm_id == "GSM7307104" ~ "Healthy Control 7",
    gsm_id == "GSM7307105" ~ "Healthy Control 8",
    TRUE ~ sample_id
  )
}

if ("gsm_id" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$sample_label <- make_sample_label(seurat_obj$gsm_id, seurat_obj$sample_id)
  seurat_obj_filtered$sample_label <- make_sample_label(seurat_obj_filtered$gsm_id, seurat_obj_filtered$sample_id)
} else {
  # Fallback: just reuse sample_id
  seurat_obj$sample_label <- seurat_obj$sample_id
  seurat_obj_filtered$sample_label <- seurat_obj_filtered$sample_id
}

cat("\nGenerating scatter plot...\n")

# Use ORIGINAL object (before filtering) to show all cells
df_scatter <- seurat_obj@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, percent.mt, sample_label, sample_id) %>%
  dplyr::mutate(sample_label = factor(sample_label))

# Map bounds_df sample_id -> sample_label (to join cutoffs correctly)
sample_map <- data.frame(
  sample = unique(seurat_obj$sample_id),
  sample_label = unique(seurat_obj$sample_label),
  stringsAsFactors = FALSE
)

bounds_plot <- bounds_df %>%
  dplyr::left_join(sample_map, by = c("sample" = "sample")) %>%
  dplyr::transmute(
    sample_label = factor(sample_label, levels = levels(df_scatter$sample_label)),
    cnt_lower  = as.numeric(cnt_lower),
    feat_lower = as.numeric(feat_lower)
  )

df_scatter2 <- df_scatter %>%
  dplyr::left_join(bounds_plot, by = "sample_label") %>%
  dplyr::mutate(
    in_bounds = nCount_RNA > cnt_lower &
      nFeature_RNA > feat_lower &
      percent.mt < mt_cutoff
  )

p_scatter_counts_features <- ggplot(df_scatter2, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(
    aes(fill = percent.mt),
    shape = 21, color = "black",
    size = 4, stroke = 0.25
  ) +
  geom_vline(
    data = bounds_plot,
    aes(xintercept = cnt_lower),
    linetype = "dashed",
    linewidth = 3,
    color = "red"
  ) +
  geom_hline(
    data = bounds_plot,
    aes(yintercept = feat_lower),
    linetype = "dashed",
    linewidth = 3,
    color = "red"
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 100),
    breaks = c(0, 20, 40, 60, 80, 100),
    oob = scales::squish,
    name = "Mitochondrial reads (%)"
  ) +
  guides(
    fill = guide_colorbar(
      ticks = TRUE,
      ticks.linewidth = 2.8,
      ticks.colour = "black",
      frame.colour = "black",
      frame.linewidth = 0.8,
      barheight = unit(150, "mm"),
      barwidth = unit(18, "mm")
    )
  ) +
  labs(
    x = "UMI counts (nCount_RNA)",
    y = "Genes detected (nFeature_RNA)"
  ) +
  ggh4x::facet_wrap2(~ sample_label, scales = "fixed", axes = "all", ncol = 4) +
  coord_fixed(ratio = 1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic(base_size = 38) +
  theme(
    plot.title = element_text(face = "bold", size = 42),
    legend.position = "right",
    axis.title = element_text(face = "bold", size = 70),
    axis.title.x = element_text(margin = margin(t = 50, r = 0, b = 0, l = 0, unit = "pt")),
    axis.title.y = element_text(margin = margin(t = 0, r = 50, b = 0, l = 0, unit = "pt")),
    axis.text.y = element_text(size = 38),
    axis.text.x = element_text(size = 36),
    strip.text = element_text(face = "bold", size = 33),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.title = element_text(face = "bold", size = 41, lineheight = 2.1, margin = margin(b = 50)),
    legend.text = element_text(face = "bold", size = 37),
    legend.ticks = element_line(linewidth = 3.2, color = "black"),
    legend.ticks.length = unit(10, "pt"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 50, unit = "pt")
  )

ggsave(
  file.path(plot_dir, "qc_scatter_by_sample.png"),
  p_scatter_counts_features,
  width = 35,
  height = 18,
  dpi = 600
)

cat("✓ Scatter plot saved\n")

# =============================================================================
# Save filtered object
# =============================================================================
saveRDS(seurat_obj_filtered, file.path(object_dir, "seurat_qc_filtered.rds"))
cat("\n✓ Saved filtered object: ", file.path(object_dir, "seurat_qc_filtered.rds"), "\n")
cat("\n✓✓✓ QC COMPLETE ✓✓✓\n\n")
