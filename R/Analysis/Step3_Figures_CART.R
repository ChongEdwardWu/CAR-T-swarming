#!/usr/bin/env Rscript
# CAR T scRNA-seq — Step 3: Figures & summary exports

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
  library(clustree)
  library(openxlsx2)
  library(future)
  library(viridis)
  library(Nebulosa)
  library(UCell)
  library(ggpubr)
  library(EnhancedVolcano)
})

set.seed(123)
# ------------------------- Config -----------------------------------------
NWORKERS              <- as.integer(Sys.getenv("NWORKERS", unset = "8"))
PLAN_STRATEGY         <- Sys.getenv("PLAN", unset = "sequential")  # or "multicore"
PROJECT_BASE          <- Sys.getenv("PROJECT_BASE", unset = ".")
WORKDIR               <- Sys.getenv("WORKDIR", unset = file.path(PROJECT_BASE, "03_R"))
RES_DIR               <- file.path(WORKDIR, "results")
FIG_DIR               <- file.path(WORKDIR, "figures")
FIG_SUBDIR            <- Sys.getenv("FIG_SUBDIR", unset = "CART")
ANNOTATED_RDS         <- Sys.getenv("ANNOTATED_RDS", unset = file.path(RES_DIR, "02_Annotation_CART.rds"))
INCLUDE_DATE_IN_NAMES <- as.logical(Sys.getenv("INCLUDE_DATE_IN_NAMES", unset = "FALSE"))
SPECIES               <- Sys.getenv("SPECIES", unset = "hs")
GROUP_LEVELS_ENV      <- Sys.getenv("GROUP_LEVELS", unset = "PostCAR,PostCC5")

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
if (!dir.exists(file.path(FIG_DIR, FIG_SUBDIR))) dir.create(file.path(FIG_DIR, FIG_SUBDIR), recursive = TRUE)
setwd(WORKDIR)
plan(PLAN_STRATEGY, workers = NWORKERS)

# ------------------------- Load object ------------------------------------
seu <- readRDS(ANNOTATED_RDS)

# Consistent grouping & l2 order (if available)
if ("group" %in% colnames(seu@meta.data)) {
  group_levels <- trimws(strsplit(GROUP_LEVELS_ENV, ",")[[1]])
  seu$group <- factor(seu$group, levels = group_levels)
}

l2_order <- c("CD8T_Cycling_MKI67", "CD8T_PolyCTL_IFNG", "CD8T_TCM_CCR7",
              "CD8T_HyperAct_MHCII", "CD8T_Tex_TOX", "CD4T_Helper_TNFRSF4")
if ("CellType_l2" %in% colnames(seu@meta.data)) {
  seu$CellType_l2 <- factor(seu$CellType_l2, levels = l2_order)
  Idents(seu) <- seu$CellType_l2
}

# Color palette for l2 classes
colorpanel <- c(
  "CD8T_Cycling_MKI67"  = "#fdbc22",
  "CD8T_PolyCTL_IFNG"   = "#d62e2d",
  "CD8T_TCM_CCR7"       = "#45c5f8",
  "CD8T_HyperAct_MHCII" = "#4FBF8C",
  "CD8T_Tex_TOX"        = "#A56BFF",
  "CD4T_Helper_TNFRSF4" = "#8C7A61"
)

# Workbook to collect summary tables
wb <- wb_workbook()

# ------------------------- Fig 1: UMAP (l2) --------------------------------
p1 <- DimPlot(seu, reduction = "umap", group.by = "CellType_l2", label = FALSE) +
  coord_fixed(ratio = 1) + scale_color_manual(values = colorpanel)

if (!inherits(p1, "try-error")) {
  ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig1_CART_UMAP_L2.png"), p1, width = 250, height = 150, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig1_CART_UMAP_L2.pdf"), p1, width = 250, height = 150, units = "mm", device = "pdf", bg = "transparent")
}

# ------------------------- Fig 2: l2 marker DotPlot -----------------------
DefaultAssay(seu) <- "SCT"
feature_genes <- list(
  "1"    = c("MKI67","STMN1","MCM7","TOP2A","PCNA"),
  "4"    = c("IFNG","GZMB","PRF1","CCL3","CCL4"),
  "6"    = c("CCR7","TCF7","CCL22","CCL1","SPIB","FSCN1"),
  "2"    = c("HLA-DRA","HLA-DRB1","HLA-DQA1","CD74"),
  "3579" = c("TOX","CD96","TIGIT","HAVCR2","LAG3","CTLA4","KLRB1","KLRD1","NR4A3"),
  "8"    = c("CD4","TNFRSF4","TNFRSF18")
)

gene_panel <- unique(unlist(feature_genes))
gene_panel <- intersect(gene_panel, rownames(seu))

p2 <- DotPlot(seu, group.by = "CellType_l2", assay = "SCT", features = rev(gene_panel)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90"), legend.direction = "horizontal", legend.position = "bottom") +
  scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))

ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig2_CART_L2_CellTypeFeat_dotplot.png"), p2, width = 100, height = 300, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig2_CART_L2_CellTypeFeat_dotplot.pdf"), p2, width = 100, height = 300, units = "mm", device = "pdf", bg = "transparent")

# ------------------------- Fig 3: cluster abundance -----------------------
source_cluster <- tibble(Cluster = seu$CellType_l2, Group = seu$group) %>%
  group_by(Cluster, Group) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no)

wb <- wb_add_worksheet(wb, sheet = "group_cluster")
wb <- wb_add_data(wb, sheet = "group_cluster", x = as.data.frame(source_cluster))

p3_bar <- ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") + scale_fill_manual(values = colorpanel) +
  coord_fixed(ratio = 1/10) + theme_bw() + xlab("group") + ylab("%")

# Donut plot (faceted by group)
donut_df <- source_cluster %>% mutate(x_pos = as.numeric(factor(Group)))
p3_donut <- ggplot(donut_df, aes(x = x_pos, y = perc, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1, colour = "black", size = 0.2) +
  coord_polar(theta = "y", direction = -1) + xlim(0.5, 2.5) + ylim(0, 100) +
  scale_fill_manual(values = colorpanel) + theme_void() + facet_wrap(~ Group, nrow = 1)

# Save abundance plots
for (nm in c("fig3_CART_L2_cluster_abundance.png", "fig3_CART_L2_cluster_abundance.pdf")) {
  ggsave(file.path(FIG_DIR, FIG_SUBDIR, nm), p3_bar, width = 100, height = 100, units = "mm", dpi = 300, device = tools::file_ext(nm))
}
# Donut
for (nm in c("fig3_CART_L2_cluster_donut.png", "fig3_CART_L2_cluster_donut.pdf")) {
  ggsave(file.path(FIG_DIR, FIG_SUBDIR, nm), p3_donut, width = 100, height = 100, units = "mm", dpi = 300, device = tools::file_ext(nm), bg = 'transparent')
}

# ------------------------- Fig 4: UMAP split by group ---------------------
p4 <- DimPlot(seu, reduction = "umap", split.by = "group", group.by = "CellType_l2", label = FALSE) +
  coord_fixed(ratio = 1) + scale_color_manual(values = colorpanel)

ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig4_CART_UMAP_split.png"), p4, width = 250, height = 150, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig4_CART_UMAP_split.pdf"), p4, width = 250, height = 150, units = "mm", device = "pdf", bg = 'transparent')

# ------------------------- Fig 5: density on UMAP -------------------------
# Black theme for density raster
theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), legend.background = element_rect(color = NA, fill = "black"),
          legend.key = element_rect(color = "white", fill = "black"), legend.position = "none",
          panel.background = element_rect(fill = "black", color = NA), panel.border = element_blank(),
          panel.grid = element_blank(), plot.background = element_rect(color = "black", fill = "black"),
          strip.background = element_rect(fill = "grey30", color = "grey10"),
          strip.text = element_text(size = base_size*0.8, color = "white"), plot.margin = unit(rep(0,4), "lines"))
}

coord <- Embeddings(seu, reduction = "umap")[, 1:2] %>% as.data.frame() %>% tibble::rownames_to_column("ID") %>% dplyr::rename(UMAP_1 = 2, UMAP_2 = 3)
meta  <- seu@meta.data %>% tibble::rownames_to_column("ID")
meta  <- dplyr::left_join(meta, coord, by = "ID")

p5 <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_viridis(option = "magma") + coord_fixed(ratio = 1) + facet_wrap(~ group, nrow = 1) + theme_black()

ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig5_CART_L2_cluster_density.png"), p5, width = 200, height = 100, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig5_CART_L2_cluster_density.pdf"), p5, width = 200, height = 100, units = "mm", device = "pdf", bg = 'transparent')

# ------------------------- Fig 6: signature scores ------------------------
DefaultAssay(seu) <- "SCT"
# Define gene signatures
sig_list <- list(
  NaiveMemory   = c("SELL","CCR7","IL7R","FOXO1","TCF7"),
  Polyfunction  = c("GZMB","IFNG","FASLG","TNF","CCL1","CCL3","CCL4","CXCL11"),
  Proliferation = c("MKI67","PCNA","TOP2A","BIRC5","MCM5"),
  Exhaustion    = c("TOX","TOX2","NR4A2","NR4A3","PDCD1","TIGIT","LAG3","HAVCR2")
)

# Score (UCell). Columns appended as *_UCell
seu <- AddModuleScore_UCell(seu, features = sig_list, ncores = NWORKERS, maxRank = 10000)

uCell_cols <- grep("_UCell$", colnames(seu@meta.data), value = TRUE)
# Density map per signature
for (sig in uCell_cols) {
  p_den <- Nebulosa::plot_density(seu, features = sig, method = c("wkde"), adjust = 2, size = 0.3, joint = FALSE) +
           coord_fixed(ratio = 1) + ggtitle(sig) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(FIG_DIR, FIG_SUBDIR, sprintf("fig6_signature_density_%s_CART.png", sig)), p_den, width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR, FIG_SUBDIR, sprintf("fig6_signature_density_%s_CART.pdf", sig)), p_den, width = 100, height = 100, units = "mm", device = "pdf", bg = 'transparent')

  # Ridge by group
  if ("group" %in% colnames(seu@meta.data)) {
    grp_fill <- c(PostCAR = "#d5d5e7", PostCC5 = "#efd2d2")
    p_ridge <- RidgePlot(seu, features = sig, group.by = "group", stack = FALSE, same.y.lims = TRUE, combine = TRUE) +
               scale_fill_manual(values = grp_fill) + scale_y_discrete(limits = c("PostCC5", "PostCAR"))
    ggsave(file.path(FIG_DIR, FIG_SUBDIR, sprintf("fig6_signature_ridge_%s_CART.png", sig)), p_ridge, width = 200, height = 100, units = "mm", dpi = 300)
    ggsave(file.path(FIG_DIR, FIG_SUBDIR, sprintf("fig6_signature_ridge_%s_CART.pdf", sig)), p_ridge, width = 200, height = 100, units = "mm", device = "pdf", bg = 'transparent')
  }
}

# ------------------------- Fig 7: PostCC5 vs PostCAR DE (CD8T) ------------
if ("CellType_l1" %in% colnames(seu@meta.data)) {
  cells.use <- WhichCells(seu, expression = CellType_l1 == "CD8T")
  seu_sub   <- subset(seu, cells = cells.use)
  if (ncol(seu_sub) > 0) {
    seu_sub$group <- factor(seu_sub$group)

    # Exclude confounders from DE
    rn <- if ("RNA" %in% names(seu_sub@assays)) rownames(seu_sub@assays$RNA@counts) else rownames(seu_sub)
    hist_genes    <- grep("^Hist", rn, ignore.case = TRUE, value = TRUE)
    hb_genes      <- grep("^Hb[ab]-|^HB[^(P)]", rn, value = TRUE)
    mt_genes      <- grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l", rn, ignore.case = TRUE, value = TRUE)
    rps_genes     <- grep("^Rp[sl]|^RP[SL]", rn, ignore.case = TRUE, value = TRUE)
    rik_genes     <- grep("^Rik", rn, ignore.case = TRUE, value = TRUE)
    alu_genes     <- grep("^AL", rn, ignore.case = TRUE, value = TRUE)
    pseudo_genes  <- grep("-rs|-ps", rn, ignore.case = TRUE, value = TRUE)
    mir_genes     <- grep("^Mir", rn, ignore.case = TRUE, value = TRUE)
    gencode_genes <- grep("^Gm", rn, ignore.case = TRUE, value = TRUE)

    bad_features <- unique(c(hist_genes, hb_genes, mt_genes, rps_genes, rik_genes, alu_genes, pseudo_genes, mir_genes, gencode_genes))
    features_de  <- setdiff(rn, bad_features)

    DefaultAssay(seu_sub) <- "RNA"
    seu_sub <- NormalizeData(seu_sub)
    Idents(seu_sub) <- "group"

    deg <- FindMarkers(seu_sub, assay = "RNA", features = features_de, ident.1 = "PostCC5", ident.2 = "PostCAR",
                       test.use = "wilcox", min.pct = 0.10, logfc.threshold = 0.05, densify = TRUE)

    volcano_data <- deg %>% tibble::rownames_to_column(var = "symbol") %>% arrange(p_val_adj)

    wb <- wb_add_worksheet(wb, sheet = "CC5vsCAR_CD8T_DEG")
    wb <- wb_add_data(wb, sheet = "CC5vsCAR_CD8T_DEG", x = as.data.frame(volcano_data))

    label_features <- c("GZMA","GZMB","CCL3","IFNG","CCL1","GNLY","FOXO1","MKI67","PCNA","MCM7",
                        "CCL4","TNF","TOP2A","CCL2","KLRB1","KLRD1","TNFRSF9","PDCD1","TOX","TIGIT","EOMES")

    p7 <- EnhancedVolcano(volcano_data, lab = volcano_data$symbol, x = "avg_log2FC", y = "p_val_adj",
                          xlim = c(-2,2), ylim = c(0,200), selectLab = label_features, pCutoff = 0.05, FCcutoff = 0.25,
                          xlab = bquote(~ Log[2] ~ "fold change"), pointSize = 1.0, labSize = 6.0, labCol = "black",
                          labFace = "bold", boxedLabels = TRUE, colAlpha = 0.8, legendPosition = "bottom",
                          legendLabSize = 14, legendIconSize = 4.0, drawConnectors = TRUE, widthConnectors = 0.5,
                          lengthConnectors = 2, arrowheads = FALSE, colConnectors = "black")

    ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig7_CC5vsCAR_Volcano_CART.png"), p7, width = 200, height = 250, units = "mm", dpi = 300)
    ggsave(file.path(FIG_DIR, FIG_SUBDIR, "fig7_CC5vsCAR_Volcano_CART.pdf"), p7, width = 200, height = 250, units = "mm", device = "pdf", bg = 'transparent')
  }
}

# ------------------------- Save outputs -----------------------------------
wb_save(wb, file.path(RES_DIR, "Step3_results_CART.xlsx"), overwrite = TRUE)

res_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("03_Figures_CART_%s.rds", format(Sys.Date(), "%Y%m%d")) else "03_Figures_CART.rds"
saveRDS(seu, file = file.path(RES_DIR, res_file))

sess_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("session_info_step3_%s.txt", format(Sys.Date(), "%Y%m%d")) else "session_info_step3.txt"
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, sess_file))
