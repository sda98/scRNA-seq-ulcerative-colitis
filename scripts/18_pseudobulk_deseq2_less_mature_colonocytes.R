#!/usr/bin/env Rscript
# =============================================================================
# 18_pseudobulk_less_mature_colonocytes.R
#
# Purpose:
#   Downstream analysis for pseudobulk DESeq2 (continuation of Script 17):
#     1) Load dds object from Script 17
#     2) PCA on VST (PC1/2 and PC3/4)
#     3) DE results: UC vs Healthy Control
#     4) Annotate results: UC vs HC
#     5) Volcano plot: UC vs HC
#     6) DE results: UC vs Self-Control
#     7) Annotate results: UC vs UCSC
#     8) Volcano plot: UC vs UCSC
#
# Inputs:
#   results/04_pseudobulk_deseq2/objects/dds_<SAFE_CELLTYPE>.rds
#   results/04_pseudobulk_deseq2/objects/pseudobulk_<SAFE_CELLTYPE>.rds
#
# Outputs:
#   results/04_pseudobulk_deseq2/tables/pseudobulk_DESeq2/<SAFE_CELLTYPE>/
#     - DESeq2_UC_vs_HC_annotated.csv
#     - DESeq2_UC_vs_UCSC_annotated.csv
#   results/04_pseudobulk_deseq2/plots/pseudobulk_DESeq2/<SAFE_CELLTYPE>/
#     - PCA_PC1_vs_PC2.png
#     - PCA_PC3_vs_PC4.png
#     - MA_plot_UC_vs_HC_unshrunk.png
#     - MA_plot_UC_vs_HC_shrunk.png
#     - Volcano_UC_vs_HC.png
#     - MA_plot_UC_vs_UCSC_unshrunk.png
#     - MA_plot_UC_vs_UCSC_shrunk.png
#     - Volcano_UC_vs_UCSC.png
#
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ashr)
  library(matrixStats)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# =============================================================================
# Paths and I/O
# =============================================================================

root_dir <- getwd()
if (basename(root_dir) == "scripts") root_dir <- dirname(root_dir)

ct <- "Colonocytes (less mature)"
safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)

out_dir    <- file.path(root_dir, "results", "04_pseudobulk_deseq2")
object_dir <- file.path(out_dir, "objects")
table_dir  <- file.path(out_dir, "tables")
plot_dir   <- file.path(out_dir, "plots")

dds_path <- file.path(object_dir, paste0("dds_", safe_ct, ".rds"))
pb_path  <- file.path(object_dir, paste0("pseudobulk_", safe_ct, ".rds"))

if (!file.exists(dds_path)) {
  stop("Missing dds object: ", dds_path, "\nRun Script 17 first.")
}
if (!file.exists(pb_path)) {
  stop("Missing pseudobulk bundle: ", pb_path, "\nRun Script 17 first.")
}

table_dir_ct <- file.path(table_dir, "pseudobulk_DESeq2", safe_ct)
plot_dir_ct  <- file.path(plot_dir,  "pseudobulk_DESeq2", safe_ct)

dir.create(table_dir_ct, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_ct,  recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  PSEUDOBULK DESEQ2 - LESS MATURE COLONOCYTES\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Cell type  :", ct, "\n")
cat("dds path   :", dds_path, "\n")
cat("Tables dir :", table_dir_ct, "\n")
cat("Plots dir  :", plot_dir_ct, "\n\n")

# =============================================================================
# 1) Load objects from Script 17
# =============================================================================

cat("Loading dds and pseudobulk bundle...\n")

dds <- readRDS(dds_path)
pb_bundle <- readRDS(pb_path)
sample_md <- pb_bundle$sample_md

cat("dds:", nrow(dds), "genes x", ncol(dds), "samples\n")
cat("✓ Objects loaded\n")

# =============================================================================
# 2) PCA
# =============================================================================

cat("\n=== PCA ===\n")

vsd <- vst(dds, blind = TRUE)

ntop <- 500
rv <- rowVars(assay(vsd))
top_idx <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

pca <- prcomp(t(assay(vsd)[top_idx, ]), center = TRUE, scale. = FALSE)

pca_data <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4]
) %>%
  left_join(sample_md %>% rownames_to_column("Sample"), by = "Sample")

pc_variance <- (pca$sdev^2) / sum(pca$sdev^2) * 100

cond_cols <- c(
  "Healthy Control"    = "blue",
  "UC Self-Control"    = "green",
  "Ulcerative Colitis" = "red"
)

# PC1 vs PC2
p_pca1 <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = condition), shape = 21,
             stroke = 0.5, size = 8.5, color = "black", alpha = 0.7) +
  scale_fill_manual(values = cond_cols, name = NULL) +
  labs(
    x = paste0("PC1: ", formatC(pc_variance[1], format = "f", digits = 1), "%"),
    y = paste0("PC2: ", formatC(pc_variance[2], format = "f", digits = 1), "%")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 32, color = "black"),
    axis.title.y = element_text(size = 32, face = "bold"),
    axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 15, unit = "pt")),
    legend.position = "right",
    legend.text = element_text(size = 18, color = "black", face = "bold"),
    legend.key.height = unit(1, "cm"),
    legend.margin = margin(l = 25, unit = "pt")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 8.5))) +
  coord_fixed(ratio = 1)

ggsave(file.path(plot_dir_ct, "PCA_PC1_vs_PC2.png"), 
       plot = p_pca1, width = 9, height = 8, units = "in", dpi = 600)

cat("✓ PCA (PC1 vs PC2) saved\n")

# PC3 vs PC4
p_pca2 <- ggplot(pca_data, aes(x = PC3, y = PC4)) +
  geom_point(aes(fill = condition), shape = 21,
             stroke = 0.5, size = 8.5, color = "black", alpha = 0.7) +
  scale_fill_manual(values = cond_cols, name = NULL) +
  labs(
    x = paste0("PC3: ", formatC(pc_variance[3], format = "f", digits = 1), "%"),
    y = paste0("PC4: ", formatC(pc_variance[4], format = "f", digits = 1), "%")
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 32, color = "black"),
    axis.title.y = element_text(size = 32, face = "bold"),
    axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 15, unit = "pt")),
    legend.position = "right",
    legend.text = element_text(size = 18, color = "black", face = "bold"),
    legend.key.height = unit(1, "cm"),
    legend.margin = margin(l = 25, unit = "pt")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 8.5))) +
  coord_fixed(ratio = 1)

ggsave(file.path(plot_dir_ct, "PCA_PC3_vs_PC4.png"), 
       plot = p_pca2, width = 9, height = 8, units = "in", dpi = 600)

cat("✓ PCA (PC3 vs PC4) saved\n")

# =============================================================================
# 3) DE results: UC vs Healthy Control
# =============================================================================

cat("\n=== UC vs Healthy Control ===\n")

res_uc_hc <- results(dds, 
                     contrast = c("condition", "Ulcerative Colitis", "Healthy Control"), 
                     alpha = 0.05)

# MA plot unshrunk
png(file.path(plot_dir_ct, "MA_plot_UC_vs_HC_unshrunk.png"), 
    width = 8, height = 7, units = "in", res = 600)
plotMA(res_uc_hc, ylim = c(-10, 10))
dev.off()

# LFC shrinkage
res_uc_hc_shrunk <- lfcShrink(dds, 
                              contrast = c("condition", "Ulcerative Colitis", "Healthy Control"), 
                              res = res_uc_hc, 
                              type = "ashr")

# MA plot shrunk
png(file.path(plot_dir_ct, "MA_plot_UC_vs_HC_shrunk.png"), 
    width = 8, height = 7, units = "in", res = 600)
plotMA(res_uc_hc_shrunk, ylim = c(-10, 10))
dev.off()

cat("✓ MA plots saved\n")

# =============================================================================
# 4) Annotate results: UC vs HC
# =============================================================================

res_tbl <- as.data.frame(res_uc_hc_shrunk) %>%
  rownames_to_column("gene_id")

gene_ids_raw <- res_tbl$gene_id

# Detect ID type
is_ensembl <- all(grepl("^ENSG\\d+", gene_ids_raw))
is_numeric <- all(grepl("^\\d+$", gene_ids_raw))
keytype <- if (is_ensembl) "ENSEMBL" else if (is_numeric) "ENTREZID" else "SYMBOL"

cat("Detected keytype:", keytype, "\n")

gene_key <- if (keytype == "ENSEMBL") sub("\\..*$", "", gene_ids_raw) else gene_ids_raw

# Annotation lookup
anno_raw <- tryCatch(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(gene_key),
    keytype = keytype,
    columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
  ) %>% distinct(),
  error = function(e) {
    message("Annotation lookup failed: ", e$message)
    NULL
  }
)

if (!is.null(anno_raw) && nrow(anno_raw) > 0) {
  if (keytype == "SYMBOL") {
    anno <- anno_raw %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarise(
        ENTREZID = paste(unique(na.omit(ENTREZID)), collapse = ";"),
        ENSEMBL  = paste(unique(na.omit(ENSEMBL)), collapse = ";"),
        GENENAME = paste(unique(na.omit(GENENAME)), collapse = ";"),
        .groups = "drop"
      ) %>%
      dplyr::rename(gene_key = SYMBOL)
  } else {
    anno <- anno_raw %>%
      dplyr::group_by(.data[[keytype]]) %>%
      dplyr::summarise(
        ENTREZID = paste(unique(na.omit(ENTREZID)), collapse = ";"),
        ENSEMBL  = paste(unique(na.omit(ENSEMBL)), collapse = ";"),
        SYMBOL   = paste(unique(na.omit(SYMBOL)), collapse = ";"),
        GENENAME = paste(unique(na.omit(GENENAME)), collapse = ";"),
        .groups = "drop"
      ) %>%
      dplyr::rename(gene_key = !!keytype)
  }
  
  deseq_tbl <- res_tbl %>%
    mutate(gene_key = gene_key) %>%
    left_join(anno, by = "gene_key") %>%
    mutate(SYMBOL_clean = gene_id)
} else {
  message("No annotation found - using gene_id as SYMBOL_clean")
  deseq_tbl <- res_tbl %>%
    mutate(
      ENSEMBL = NA_character_,
      ENTREZID = NA_character_,
      GENENAME = NA_character_,
      SYMBOL_clean = gene_id
    )
}

deseq_tbl <- deseq_tbl %>%
  dplyr::select(gene_id, any_of(c("ENSEMBL", "ENTREZID", "GENENAME")), SYMBOL_clean,
                baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  arrange(padj)

write_csv(deseq_tbl, file.path(table_dir_ct, "DESeq2_UC_vs_HC_annotated.csv"))

# Summary
cat("\nDE Summary (UC vs HC):\n")
cat("  Total genes tested:", nrow(deseq_tbl), "\n")
cat("  Significant (padj < 0.05):", sum(deseq_tbl$padj < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated:", sum(deseq_tbl$padj < 0.05 & deseq_tbl$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Downregulated:", sum(deseq_tbl$padj < 0.05 & deseq_tbl$log2FoldChange < 0, na.rm = TRUE), "\n")

# =============================================================================
# 5) Volcano plot: UC vs HC
# =============================================================================

deseq_volcano <- deseq_tbl %>%
  filter(!is.na(padj)) %>%
  mutate(
    diffexpressed = case_when(
      log2FoldChange >= 1 & padj <= 0.05  ~ "Upregulated in UC",
      log2FoldChange <= -1 & padj <= 0.05 ~ "Downregulated in UC",
      TRUE ~ "No significant change"
    ),
    diffexpressed = factor(diffexpressed, 
                           levels = c("Upregulated in UC", "No significant change", "Downregulated in UC"))
  )

# Functional categories for labeling
gene_categories <- data.frame(
  gene = c(
    "REG1A", "LCN2",
    "XBP1", "MANF",
     "PI3",
    "CXCL1", "CCL20",
    "HMGCS2", "ZG16", "OLFM4"
  ),
  category = c(
    rep("Anti-Apoptosis Defense", 2),
    rep("ER Stress Response", 2),
    rep("Intestinal Epithelium Regeneration", 1),
    rep("Immune Cells Recruitment", 2),
    rep("Intestinal Mucus Barrier", 3)
  )
)

category_colors <- c(
  "Anti-Apoptosis Defense"              = "cornsilk1",
  "ER Stress Response"                  = "lightblue1",
  "Intestinal Epithelium Regeneration"  = "pink",
  "Immune Cells Recruitment"            = "lightgreen",
  "Intestinal Mucus Barrier"            = "gold"
)

deseq_volcano <- deseq_volcano %>%
  left_join(gene_categories, by = c("SYMBOL_clean" = "gene")) %>%
  mutate(
    delabel = ifelse(
      !is.na(category) & padj < 0.05 & abs(log2FoldChange) > 1,
      SYMBOL_clean, NA
    )
  )

label_df <- subset(deseq_volcano, !is.na(delabel))
cat("Genes being labeled:", nrow(label_df), "\n")

volcano_plot <- ggplot(deseq_volcano, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(deseq_volcano, diffexpressed == "No significant change"),
             size = 7, shape = 21, color = "black", fill = "grey50", alpha = 1) +
  geom_point(data = subset(deseq_volcano, diffexpressed == "Upregulated in UC"),
             size = 7, shape = 21, color = "black", fill = "#b2182b", alpha = 1) +
  geom_point(data = subset(deseq_volcano, diffexpressed == "Downregulated in UC"),
             size = 7, shape = 21, color = "black", fill = "#2166ac", alpha = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 1) +
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 2)) +
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2), expand = c(0, 0)) +
  labs(
    x = expression("Log"[2]*"FC"),
    y = expression("-log"[10]*"(P"[adj]*")"),
    title = "Differential Expression: UC vs Healthy Control",
    subtitle = ct
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 45, face = "bold", color = "black"),
    axis.title.x = element_text(size = 45, face = "bold", color = "black", 
                                margin = margin(t = 20, unit = "pt")),
    axis.text.x = element_text(size = 52, color = "black", margin = margin(t = 5, unit = "pt")),
    axis.text.y = element_text(size = 52, color = "black", margin = margin(r = 5, unit = "pt")),
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black", size = 52,
                              margin = margin(b = 5, unit = "pt")),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 45,
                                 margin = margin(b = 15, unit = "pt")),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(t = 10, r = 10, b = 10, l = 10),
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 27),
    legend.key.size = unit(1.5, "lines")
  )

if (nrow(label_df) > 0) {
  volcano_plot <- volcano_plot +
    ggrepel::geom_label_repel(
      data = label_df,
      aes(label = delabel, fill = category),
      color = "black",
      box.padding = 2.5,
      point.padding = 0.8,
      segment.color = "black",
      segment.size = 0.9,
      max.overlaps = Inf,
      max.iter = 10000,
      size = 9,
      fontface = "bold.italic",
      show.legend = TRUE
    ) +
    scale_fill_manual(values = category_colors, name = "Category") +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 10))) +
    theme(
      legend.key.height = unit(1.9, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 37, margin = margin(b = 30, unit = "pt"))
    )
}

ggsave(file.path(plot_dir_ct, "Volcano_UC_vs_HC.png"), 
       plot = volcano_plot, width = 30, height = 18, units = "in", dpi = 600)

cat("✓ Volcano plot (UC vs HC) saved\n")

# =============================================================================
# 6) DE results: UC vs Self-Control
# =============================================================================

cat("\n=== UC vs Self-Control ===\n")

res_uc_ucsc <- results(dds, 
                       contrast = c("condition", "Ulcerative Colitis", "UC Self-Control"), 
                       alpha = 0.05)

# MA plot unshrunk
png(file.path(plot_dir_ct, "MA_plot_UC_vs_UCSC_unshrunk.png"), 
    width = 8, height = 7, units = "in", res = 600)
plotMA(res_uc_ucsc, ylim = c(-10, 10))
dev.off()

# LFC shrinkage
res_uc_ucsc_shrunk <- lfcShrink(dds, 
                                contrast = c("condition", "Ulcerative Colitis", "UC Self-Control"), 
                                res = res_uc_ucsc, 
                                type = "ashr")

# MA plot shrunk
png(file.path(plot_dir_ct, "MA_plot_UC_vs_UCSC_shrunk.png"), 
    width = 8, height = 7, units = "in", res = 600)
plotMA(res_uc_ucsc_shrunk, ylim = c(-10, 10))
dev.off()

cat("✓ MA plots saved\n")

# =============================================================================
# 7) Annotate results: UC vs UCSC
# =============================================================================

res_tbl2 <- as.data.frame(res_uc_ucsc_shrunk) %>%
  rownames_to_column("gene_id")

gene_key2 <- if (keytype == "ENSEMBL") sub("\\..*$", "", res_tbl2$gene_id) else res_tbl2$gene_id

if (!is.null(anno_raw) && nrow(anno_raw) > 0) {
  deseq_tbl2 <- res_tbl2 %>%
    mutate(gene_key = gene_key2) %>%
    left_join(anno, by = "gene_key") %>%
    mutate(SYMBOL_clean = gene_id)
} else {
  deseq_tbl2 <- res_tbl2 %>%
    mutate(
      ENSEMBL = NA_character_,
      ENTREZID = NA_character_,
      GENENAME = NA_character_,
      SYMBOL_clean = gene_id
    )
}

deseq_tbl2 <- deseq_tbl2 %>%
  dplyr::select(gene_id, any_of(c("ENSEMBL", "ENTREZID", "GENENAME")), SYMBOL_clean,
                baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  arrange(padj)

write_csv(deseq_tbl2, file.path(table_dir_ct, "DESeq2_UC_vs_UCSC_annotated.csv"))

# Summary
cat("\nDE Summary (UC vs UCSC):\n")
cat("  Total genes tested:", nrow(deseq_tbl2), "\n")
cat("  Significant (padj < 0.05):", sum(deseq_tbl2$padj < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated:", sum(deseq_tbl2$padj < 0.05 & deseq_tbl2$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Downregulated:", sum(deseq_tbl2$padj < 0.05 & deseq_tbl2$log2FoldChange < 0, na.rm = TRUE), "\n")

# =============================================================================
# 8) Volcano plot: UC vs UCSC
# =============================================================================

deseq_volcano2 <- deseq_tbl2 %>%
  filter(!is.na(padj)) %>%
  mutate(
    diffexpressed = case_when(
      log2FoldChange >= 1 & padj <= 0.05  ~ "Upregulated in UC",
      log2FoldChange <= -1 & padj <= 0.05 ~ "Downregulated in UC",
      TRUE ~ "No significant change"
    ),
    diffexpressed = factor(diffexpressed, 
                           levels = c("Upregulated in UC", "No significant change", "Downregulated in UC"))
  )

# Functional categories for UC vs UCSC
gene_categories2 <- data.frame(
  gene = c("MCL1", "REG1A", "LCN2", 
           "DERL3",
           "TIMP1", "PI3", 
           "CXCL1", "CCL20",
           "HMGCS2"),
  category = c(
    rep("Anti-Apoptosis Defense", 3),
    rep("ER Stress Response", 1),
    rep("Intestinal Epithelium Regeneration", 2),
    rep("Immune Cells Recruitment", 2),
    rep("Intestinal Mucus Barrier", 1)
  )
)

category_colors2 <- c(
  "Anti-Apoptosis Defense"              = "cornsilk1",
  "ER Stress Response"                  = "lightblue1",
  "Intestinal Epithelium Regeneration"  = "pink",
  "Immune Cells Recruitment"            = "lightgreen",
  "Intestinal Mucus Barrier"            = "gold"
)

deseq_volcano2 <- deseq_volcano2 %>%
  left_join(gene_categories2, by = c("SYMBOL_clean" = "gene")) %>%
  mutate(
    delabel = ifelse(
      !is.na(category) & padj < 0.05 & abs(log2FoldChange) > 1,
      SYMBOL_clean, NA
    )
  )

label_df2 <- subset(deseq_volcano2, !is.na(delabel))
cat("Genes being labeled:", nrow(label_df2), "\n")

volcano_plot2 <- ggplot(deseq_volcano2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(deseq_volcano2, diffexpressed == "No significant change"),
             size = 7, shape = 21, color = "black", fill = "grey50", alpha = 1) +
  geom_point(data = subset(deseq_volcano2, diffexpressed == "Upregulated in UC"),
             size = 7, shape = 21, color = "black", fill = "#b2182b", alpha = 1) +
  geom_point(data = subset(deseq_volcano2, diffexpressed == "Downregulated in UC"),
             size = 7, shape = 21, color = "black", fill = "#2166ac", alpha = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 1) +
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 2)) +
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2), expand = c(0, 0)) +
  labs(
    x = expression("Log"[2]*"FC"),
    y = expression("-log"[10]*"(P"[adj]*")"),
    title = "Differential Expression: UC vs Self-Control",
    subtitle = ct
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 45, face = "bold", color = "black"),
    axis.title.x = element_text(size = 45, face = "bold", color = "black", 
                                margin = margin(t = 20, unit = "pt")),
    axis.text.x = element_text(size = 52, color = "black", margin = margin(t = 5, unit = "pt")),
    axis.text.y = element_text(size = 52, color = "black", margin = margin(r = 5, unit = "pt")),
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black", size = 52,
                              margin = margin(b = 5, unit = "pt")),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 45,
                                 margin = margin(b = 15, unit = "pt")),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(t = 10, r = 10, b = 10, l = 10),
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 27),
    legend.key.size = unit(1.5, "lines")
  )

if (nrow(label_df2) > 0) {
  volcano_plot2 <- volcano_plot2 +
    ggrepel::geom_label_repel(
      data = label_df2,
      aes(label = delabel, fill = category),
      color = "black",
      box.padding = 2.2,
      nudge_y = 0,
      point.padding = 0.8,
      segment.color = "black",
      segment.size = 0.9,
      max.overlaps = Inf,
      max.iter = 10000,
      size = 9,
      fontface = "bold.italic",
      show.legend = TRUE
    ) +
    scale_fill_manual(values = category_colors2, name = "Category") +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 10))) +
    theme(
      legend.key.height = unit(1.9, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 37, margin = margin(b = 30, unit = "pt"))
    )
}

ggsave(file.path(plot_dir_ct, "Volcano_UC_vs_UCSC.png"), 
       plot = volcano_plot2, width = 30, height = 18, units = "in", dpi = 600)

cat("✓ Volcano plot (UC vs UCSC) saved\n")

# =============================================================================
# Done
# =============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Tables:", table_dir_ct, "\n")
cat("Plots:", plot_dir_ct, "\n")
