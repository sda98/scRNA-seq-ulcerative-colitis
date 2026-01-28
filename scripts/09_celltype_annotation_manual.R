#!/usr/bin/env Rscript

# ==============================================================================
# 09_celltype_annotation_manual.R
# ------------------------------------------------------------------------------
# Purpose
#   1) Load clustered Seurat object (Harmony-integrated; RNA_snn_res.0.4; UMAP present)
#   2) Add manual, marker-based final cell type labels (cluster -> celltype_final)
#   3) Export final cell type distribution table
#   4) Export DotPlot of epithelial + colonocyte differentiation markers
#   5) Save annotated Seurat object for downstream scripts
#
# Input
#   - results/02_clustering_analysis/objects/seurat_integrated_clustered.rds
#
# Outputs
#   - results/02_clustering_analysis/sessionInfo.txt
#   - results/02_clustering_analysis/tables/celltype_final_counts.csv
#   - results/02_clustering_analysis/plots/dotplot_colonocyte_differentiation_markers.png
#   - results/02_clustering_analysis/objects/seurat_integrated_annotated.rds
#
# Notes
#   - Cluster 4 is retained but labeled "Low-quality (excluded)".
#   - DotPlot uses the RNA assay and Seurat's avg.exp.scaled (scaled expression).
#   - ggnewscale is used to display 3 independent color scales (epi / early / late) with
#     separate colorbars in the legend.
#
# Usage
#   Rscript scripts/09_celltype_annotation_manual.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(grid)
  library(ggnewscale)
})

# =========================
# Paths and I/O
# =========================
root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

in_obj <- file.path(root_dir, "results", "02_clustering_analysis", "objects", "seurat_integrated_clustered.rds")
if (!file.exists(in_obj)) {
  stop("Input object not found: ", in_obj,
       "\nRun 07_clustering.R first (or ensure clustered object exists).")
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
# 1) LOAD SEURAT OBJECT
# =========================
cat("\nLoading clustered Seurat object...\n")
seurat_integrated <- readRDS(in_obj)

if (!"RNA_snn_res.0.4" %in% colnames(seurat_integrated@meta.data)) {
  stop("Clustering column 'RNA_snn_res.0.4' not found. Make sure res=0.4 exists.")
}
if (!"umap" %in% Reductions(seurat_integrated)) {
  stop("UMAP reduction not found in object. Make sure UMAP was computed in Script 06/07.")
}

cat("✓ Object loaded\n")

# =============================================================================
# 2) MANUAL CELL TYPE ANNOTATION (RES=0.4)
# =============================================================================
cat("\n=== CELL TYPE ANNOTATION (MANUAL) ===\n")

Idents(seurat_integrated) <- "RNA_snn_res.0.4"

final_annotations <- c(
  "0"  = "CD4+ T cells",
  "1"  = "Plasma cells (IgA+)",
  "2"  = "Cytotoxic lymphocytes (CD8+ T / NK)",
  "3"  = "B cells (naive/memory)",
  "4"  = "Low-quality (excluded)",
  "5"  = "Plasma cells (IgG+)",
  "6"  = "Colonocytes (less mature)",
  "7"  = "Colonocytes (more mature)",
  "8"  = "Inflammatory fibroblasts",
  "9"  = "Macrophages",
  "10" = "Cycling B cells",
  "11" = "Goblet cells",
  "12" = "Endothelial cells",
  "13" = "Mast cells",
  "14" = "Enteric glia",
  "15" = "Tuft cells",
  "16" = "Pericytes / SMCs",
  "17" = "Germinal center B cells",
  "18" = "Enteroendocrine cells (L-cell enriched)"
)

cluster_ids <- as.character(Idents(seurat_integrated))
missing_clusters <- setdiff(sort(unique(cluster_ids)), names(final_annotations))
if (length(missing_clusters) > 0) {
  stop("Annotation map is missing these cluster IDs: ", paste(missing_clusters, collapse = ", "))
}

seurat_integrated$celltype_final <- unname(final_annotations[cluster_ids])
seurat_integrated$celltype_final <- factor(seurat_integrated$celltype_final)

cat("✓ Final annotations added (celltype_final)\n")
cat("\nCell type distribution:\n")
print(table(seurat_integrated$celltype_final, useNA = "ifany"))

# Export counts table
celltype_counts <- as.data.frame(table(seurat_integrated$celltype_final, useNA = "ifany"))
colnames(celltype_counts) <- c("celltype_final", "n_cells")
write.csv(celltype_counts, file.path(table_dir, "celltype_final_counts.csv"), row.names = FALSE)


# =============================================================================
# 3) DOTPLOT: COLONOCYTE DIFFERENTIATION MARKERS (CLUSTERS 0-18)
# =============================================================================
cat("\n=== DotPlot: colonocyte differentiation markers ===\n")

Idents(seurat_integrated) <- "RNA_snn_res.0.4"
DefaultAssay(seurat_integrated) <- "RNA"

cluster_numbers <- as.character(0:18)
cluster_labels  <- paste0("Cluster ", cluster_numbers)

seurat_integrated$cluster_ordered <- factor(
  paste0("Cluster ", as.character(Idents(seurat_integrated))),
  levels = cluster_labels
)
Idents(seurat_integrated) <- "cluster_ordered"

# ---- marker sets ----
# Pan-epithelial markers (keep as-is)
genes_epithelial <- c("EPCAM", "LGALS4", "CLDN4", "AGR2", "CA1", "CA2")

# Early/immature colonocytes - REVISED
genes_early <- c(
  "SLC12A2",
  "EPHB2",
  "EPHB3",
  "OLFM4",    
  "MKI67",  
  "TOP2A",   
  "PCNA"     
)

# Late/mature colonocytes (mostly keep, minor tweaks)
genes_late <- c(
  "GUCA2A", "GUCA2B",    # Surface/crypt-top specific 
  "SLC26A3",             # Surface colonocyte 
  "CA4",                 # Mature BEST4+ colonocyte 
  "CEACAM7",             # Mature colonocyte 
  "KRT20",               # Differentiation marker 
  "BEST4"               # Specialized mature colonocytes
)

features_use <- c(genes_epithelial, genes_early, genes_late)
features_use <- features_use[features_use %in% rownames(seurat_integrated)]
cat("Using", length(features_use), "markers\n")

# ---- DotPlot data ----
p_dot_raw <- DotPlot(
  seurat_integrated,
  features  = features_use,
  dot.scale = 10
) + RotatedAxis()

df_dot <- p_dot_raw$data %>%
  dplyr::mutate(
    gene = as.character(features.plot),
    group = dplyr::case_when(
      gene %in% genes_epithelial ~ "Epithelial\ncells",
      gene %in% genes_early      ~ "Colonocytes\n(less differentiated)",
      gene %in% genes_late       ~ "Colonocytes\n(more differentiated)",
      TRUE ~ "Other"
    ),
    gene = factor(gene, levels = features_use),
    id   = factor(id, levels = rev(cluster_labels))
  ) %>%
  dplyr::filter(group != "Other")

# ---- x positions for separators + label placement ----
gene_pos <- tibble::tibble(
  gene  = factor(features_use, levels = features_use),
  x_pos = seq_along(features_use)
) %>%
  dplyr::mutate(
    gene_chr = as.character(gene),
    group = dplyr::case_when(
      gene_chr %in% genes_epithelial ~ "Epithelial\ncells",
      gene_chr %in% genes_early      ~ "Colonocytes\n(less differentiated)",
      gene_chr %in% genes_late       ~ "Colonocytes\n(more differentiated)",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(group != "Other")

sep_1 <- max(gene_pos$x_pos[gene_pos$group == "Epithelial\ncells"]) + 0.5
sep_2 <- max(gene_pos$x_pos[gene_pos$group == "Colonocytes\n(less differentiated)"]) + 0.5

group_centers <- gene_pos %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(x = mean(x_pos), .groups = "drop")

y_top <- length(cluster_labels) + 2.0

# ---- section header label fills (light tints matched to palettes) ----
label_df <- group_centers %>%
  dplyr::mutate(
    y = y_top,
    fill_col = dplyr::case_when(
      group == "Epithelial\ncells" ~ "#D6ECFF",
      group == "Colonocytes\n(less differentiated)" ~ "#FDE0D2",
      group == "Colonocytes\n(more differentiated)"  ~ "#D9F2E3",
      TRUE ~ "grey95"
    )
  )

# ---- distinct palettes ----
pal_epi   <- c("white", "#2C7FB8")
pal_early <- c("white", "#F46D43")
pal_late  <- c("white", "#31A354")

# --- legend title ---
expr_legend_title <- "Average \nrelative expression \nwith respect to other clusters"

# --- colorbar sizing ---
cb_h <- unit(62, "mm")
cb_w <- unit(14, "mm")

# --- add tick marks on the bars ---
cb_ticks     <- TRUE
cb_tick_col  <- "black"
cb_tick_lw   <- 0.9
cb_frame_lw  <- 0.6

p_dot_final <- ggplot() +
  
  # ---- Epithelial markers layer + colorbar 1 ----
geom_point(
  data = df_dot[df_dot$group == "Epithelial\ncells", ],
  aes(x = gene, y = id, size = pct.exp, color = avg.exp.scaled)
) +
  scale_color_gradient(
    low = pal_epi[1], high = pal_epi[2],
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    oob = scales::squish,
    name = expr_legend_title,
    guide = guide_colorbar(
      order = 1,
      barheight = cb_h, barwidth = cb_w,
      ticks = cb_ticks,
      ticks.colour = cb_tick_col,
      ticks.linewidth = cb_tick_lw,
      frame.colour = "black",
      frame.linewidth = cb_frame_lw
    )
  ) +
  
  ggnewscale::new_scale_color() +
  
  # ---- LESS differentiated (early) markers layer + colorbar 2 ----
geom_point(
  data = df_dot[df_dot$group == "Colonocytes\n(less differentiated)", ],
  aes(x = gene, y = id, size = pct.exp, color = avg.exp.scaled)
) +
  scale_color_gradient(
    low = pal_early[1], high = pal_early[2],
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    oob = scales::squish,
    name = expr_legend_title,
    guide = guide_colorbar(
      order = 2,
      barheight = cb_h, barwidth = cb_w,
      ticks = cb_ticks,
      ticks.colour = cb_tick_col,
      ticks.linewidth = cb_tick_lw,
      frame.colour = "black",
      frame.linewidth = cb_frame_lw
    )
  ) +
  
  ggnewscale::new_scale_color() +
  
  # ---- MORE differentiated (late) markers layer + colorbar 3 ----
geom_point(
  data = df_dot[df_dot$group == "Colonocytes\n(more differentiated)", ],
  aes(x = gene, y = id, size = pct.exp, color = avg.exp.scaled)
) +
  scale_color_gradient(
    low = pal_late[1], high = pal_late[2],
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    oob = scales::squish,
    name = expr_legend_title,
    guide = guide_colorbar(
      order = 3,
      barheight = cb_h, barwidth = cb_w,
      ticks = cb_ticks,
      ticks.colour = cb_tick_col,
      ticks.linewidth = cb_tick_lw,
      frame.colour = "black",
      frame.linewidth = cb_frame_lw
    )
  ) +
  # ---- shared size legend ----
scale_size_continuous(
  range  = c(0.5, 14),
  limits = c(0, 100),
  breaks = c(10, 30, 50, 70, 90),
  name   = "Percent expressed",
  guide  = guide_legend(
    order = 4,
    override.aes = list(size = 8),
    keyheight = unit(10, "mm"),
    keywidth  = unit(10, "mm")
  )
) +
  # ---- separators ----
geom_vline(
  xintercept = c(sep_1, sep_2),
  linetype = "solid", linewidth = 2.5, color = "black"
) +
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = group, fill = I(fill_col)),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 10,
    color = "black",
    label.size = 0.8,
    label.r = unit(0.15, "lines"),
    label.padding = unit(0.35, "lines"),
    show.legend = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 28) +
  theme(
    axis.text.x = element_text(face = "bold.italic", size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 23),
    axis.ticks  = element_line(linewidth = 0.6),
    legend.title = element_text(face = "bold", size = 24, margin = margin(b = 20, unit = "pt")),
    legend.text  = element_text(face = "bold", size = 25),
    legend.box   = "vertical",
    legend.spacing.y = unit(14, "mm"),
    legend.box.spacing = unit(6, "mm"),
    legend.margin      = margin(0, 0, 0, 0, unit = "mm"),
    plot.margin  = margin(t = 100, r = 20, b = 10, l = 20, unit = "pt")
  )


ggsave(
  file.path(plot_dir, "dotplot_colonocyte_differentiation_markers.png"),
  p_dot_final, width = 22, height = 15, dpi = 600
)


# =============================================================================
# 4) SAVE ANNOTATED OBJECT
# =============================================================================
Idents(seurat_integrated) <- "RNA_snn_res.0.4" 
out_obj <- file.path(object_dir, "seurat_integrated_annotated.rds")
saveRDS(seurat_integrated, out_obj)

cat("\n✓ Saved annotated object: ", out_obj, "\n")
cat("\n✓✓✓ MANUAL ANNOTATION COMPLETE ✓✓✓\n\n")
