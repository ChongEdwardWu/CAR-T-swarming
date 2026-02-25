#!/usr/bin/env Rscript
# PBMC scRNA-seq — Step 3: Figures & signature analyses

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
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
})

set.seed(123)
# ------------------------- Parameters -------------------------------------
NWORKERS              <- as.integer(Sys.getenv("NWORKERS", unset = "8"))
PLAN_STRATEGY         <- Sys.getenv("PLAN", unset = "sequential")  # or "multicore"
PROJECT_BASE          <- Sys.getenv("PROJECT_BASE", unset = ".")
WORKDIR               <- Sys.getenv("WORKDIR", unset = file.path(PROJECT_BASE, "03_R"))
RES_DIR               <- file.path(WORKDIR, "results")
FIG_DIR               <- file.path(WORKDIR, "figures")
FIG_DIR_PBMC          <- file.path(FIG_DIR, "PBMC")
ANNOTATED_RDS         <- Sys.getenv("ANNOTATED_RDS", unset = file.path(RES_DIR, "02_Annotation_PBMC.rds"))
GROUP_LEVELS_ENV      <- Sys.getenv("GROUP_LEVELS", unset = "PBMC_CAR,PBMC_CC5")
INCLUDE_DATE_IN_NAMES <- as.logical(Sys.getenv("INCLUDE_DATE_IN_NAMES", unset = "FALSE"))
GENE_BLACKLIST_RDS    <- Sys.getenv("GENE_BLACKLIST_RDS", unset = "")  # optional

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR_PBMC)) dir.create(FIG_DIR_PBMC, recursive = TRUE)
setwd(WORKDIR)
plan(PLAN_STRATEGY, workers = NWORKERS)

# ------------------------- Load object ------------------------------------
seu <- readRDS(ANNOTATED_RDS)

# ------------------------- Factors & palettes -----------------------------
# l1 overview
if ("CellType_l1" %in% colnames(seu@meta.data)) print(table(seu$CellType_l1))
# group ordering
grp_levels <- trimws(strsplit(GROUP_LEVELS_ENV, ",")[[1]])
if ("group" %in% colnames(seu@meta.data)) {
  seu$group <- factor(seu$group, levels = grp_levels)
  print(table(seu$group))
}
# l2 ordering (kept from Step2)
if ("CellType_l2" %in% colnames(seu@meta.data)) {
  Idents(seu) <- seu$CellType_l2
  print(levels(seu$CellType_l2))
}

colorpanel <- c(
  # CD8 T (naïve → tex)
  "CD8T_Tn_CCR7"              = "#8CCBF8",
  "CD8T_TCM_IL7R"             = "#4FA5F6",
  "CD8T_TCM_Cycling"          = "#F4C84B",
  "CD8T_PreDys_GZMK"          = "#F79B3E",
  "CD8T_PreDys_CSF2"          = "#F05A2C",
  "CD8T_Tex_TNFRSF4"          = "#D32626",
  "CD8T_Tex_TNFRSF9"          = "#A1122E",
  # MAIT
  "CD8T_MAIT_CCL20"           = "#A56BFF",
  "CD8T_MAIT_Exh_CTLA4"       = "#7435BF",
  # CD4 T
  "CD4T_Tn_LEF1"              = "#C4DA5D",
  "CD4T_TCM_KLF2"             = "#4FBF8C",
  "CD4T_EarlyAct_TNFRSF4"     = "#27966B",
  "CD4T_Cycling"              = "#B4E4A0",
  # NK
  "NK_Cycling"                = "#F6A2C9",
  "NK_EarlyAct_NR4A"          = "#E053B2",
  "NK_IFNprimed_LDB2"         = "#B82485",
  "NK_Conventional_Cytotoxic" = "#79144F",
  # B
  "B_naive_FCER2"             = "#BBA78B",
  "B_Inflammatory_STAT1"      = "#8C7A61"
)

# ------------------------- Figure 1: UMAP (l2) ----------------------------
p_umap_l2 <- DimPlot(seu, reduction = "umap", group.by = "CellType_l2", label = FALSE) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = colorpanel)

ggsave(file.path(FIG_DIR_PBMC, "fig1_PBMC_UMAP_L2.png"), p_umap_l2, width = 250, height = 150, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig1_PBMC_UMAP_L2.pdf"), p_umap_l2, width = 250, height = 150, units = "mm", device = "pdf", bg = "transparent")

# ------------------------- Figure 2: l2 lineage markers -------------------
feature_genes <- list(
  T_cells  = c("CD3D","CD3E","TRAC","CD8A","CD8B","CD4"),
  NK_cells = c("NKG7","KLRD1","NCR1"),
  B_cells  = c("MS4A1","CD79A")
)

DefaultAssay(seu) <- "SCT"
l2_gene_panel <- intersect(unique(unlist(feature_genes)), rownames(seu))

p_dot_l2 <- DotPlot(seu, group.by = "CellType_l2", assay = "SCT", features = rev(l2_gene_panel)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        legend.direction = "horizontal",
        legend.position  = "bottom") +
  scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))

ggsave(file.path(FIG_DIR_PBMC, "fig2_PBMC_L2_CellTypeFeat_dotplot.png"), p_dot_l2, width = 200, height = 150, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig2_PBMC_L2_CellTypeFeat_dotplot.pdf"), p_dot_l2, width = 200, height = 150, units = "mm", device = "pdf", bg = "transparent")

# ------------------------- Figure 3: Cluster markers per lineage ----------
# Helper to set identities safely
.set_idents <- function(obj, id_col) {
  if (is.character(id_col)) Idents(obj) <- id_col else Idents(obj) <- id_col
  obj
}

## --- CD8 T cells
seu_CD8T <- seu[, seu$CellType_l1 == "CD8T"]
seu_CD8T$CellType_l2 <- droplevels(seu_CD8T$CellType_l2)
seu_CD8T <- .set_idents(seu_CD8T, "CellType_l2")
CD8T_genes <- list(
  `CD8T_Tn_CCR7`       = c("CCR7","TCF7"),
  `CD8T_TCM_IL7R`      = c("IL7R","LEF1","KLF2"),
  `CD8T_TCM_Cycling`   = c("MKI67","TOP2A","BIRC5"),
  `CD8T_PreDys_GZMK`   = c("GZMK","GZMH","HAVCR2"),
  `CD8T_PreDys_CSF2`   = c("CSF2"),
  `CD8T_Tex_TNFRSF4`   = c("TOX","TNFRSF4","BATF","LAG3"),
  `CD8T_Tex_TNFRSF9`   = c("TNFRSF9","TNFRSF18","TIGIT"),
  `CD8T_MAIT_CCL20`    = c("CCL20"),
  `CD8T_MAIT_Exh_CTLA4`= c("SLC4A10","IL23R","KLRB1","CTLA4")
)
CD8T_panel <- intersect(unique(unlist(CD8T_genes)), rownames(seu_CD8T))
CD8T_dotplot <- DotPlot(seu_CD8T, group.by = "CellType_l2", assay = "SCT", features = rev(CD8T_panel)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90"), legend.direction = "horizontal", legend.position = "bottom") +
  scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))

## --- CD4 T cells
seu_CD4T <- seu[, seu$CellType_l1 == "CD4T"]
seu_CD4T$CellType_l2 <- droplevels(seu_CD4T$CellType_l2)
seu_CD4T <- .set_idents(seu_CD4T, "CellType_l2")
CD4T_genes <- list(
  `CD4T_Tn_LEF1`      = c("LEF1","SELL","MAL","TCF7","CCR7","NME4","IFNGR2","IRF8","IGF1R"),
  `CD4T_TCM_KLF2`     = c("KLF2","IL7R","LTB"),
  `CD4T_EarlyAct_TNFRSF4` = c("TNFRSF4","TIGIT","TNFRSF9","TNFRSF18"),
  `CD4T_Cycling`      = c("MKI67","TOP2A","BIRC5","STMN1","PCNA")
)
CD4T_panel <- intersect(unique(unlist(CD4T_genes)), rownames(seu_CD4T))
CD4T_dotplot <- DotPlot(seu_CD4T, group.by = "CellType_l2", assay = "SCT", features = rev(CD4T_panel)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90"), legend.direction = "horizontal", legend.position = "bottom") +
  scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))

## --- NK cells
seu_NK <- seu[, seu$CellType_l1 == "NK"]
seu_NK$CellType_l2 <- droplevels(seu_NK$CellType_l2)
seu_NK <- .set_idents(seu_NK, "CellType_l2")
NK_genes <- list(
  `NK_Cycling`                = c("STMN1","H2AFX","MKI67","TOP2A","BIRC5"),
  `NK_EarlyAct_NR4A`          = c("LDB2","PECAM1","DPP4","ITGA1","KIR2DL4"),
  `NK_IFNprimed_LDB2`         = c("NR4A3","DUSP2","KLRF1","FGFBP2","RGS2","AREG"),
  `NK_Conventional_Cytotoxic` = c("FCER1G","TYROBP","GZMB","GZMM","HCST")
)
NK_panel <- intersect(unique(unlist(NK_genes)), rownames(seu_NK))
NK_dotplot <- DotPlot(seu_NK, group.by = "CellType_l2", assay = "SCT", features = rev(NK_panel)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90"), legend.direction = "horizontal", legend.position = "bottom") +
  scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))

## --- B cells
seu_B <- seu[, seu$CellType_l1 == "B"]
seu_B$CellType_l2 <- droplevels(seu_B$CellType_l2)
seu_B <- .set_idents(seu_B, "CellType_l2")
B_genes <- list(
  `B_naive_FCER2`           = c("FCER2","TCL1A","SELL","BACH2","BTLA","TNFRSF13C"),
  `B_Inflammatory_STAT1`    = c("STAT1","IRF1","CXCL10","CCL22","CCL17","CXCL9","TNFRSF13B")
)
B_panel <- intersect(unique(unlist(B_genes)), rownames(seu_B))
B_dotplot <- DotPlot(seu_B, group.by = "CellType_l2", assay = "SCT", features = rev(B_panel)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90"), legend.direction = "horizontal", legend.position = "bottom") +
  scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))

# Combine 4 dotplots vertically with proportional heights
get_n_ybreaks <- function(p){
  gb <- ggplot_build(p)
  pp <- gb$layout$panel_params[[1]]
  # ggplot2 pre/post 3.5 compatibility
  if (!is.null(pp$y$get_breaks())) return(length(pp$y$get_breaks()))
  if (!is.null(pp$y$breaks))      return(length(pp$y$breaks))
  return(1)
}

plot_list <- list(CD8T = CD8T_dotplot, CD4T = CD4T_dotplot, NK = NK_dotplot, B = B_dotplot)
heights_raw <- vapply(plot_list, get_n_ybreaks, numeric(1))
heights_rel <- if (sum(heights_raw) > 0) heights_raw / sum(heights_raw) else rep(1/length(plot_list), length(plot_list))

panel_fig <- wrap_plots(plot_list, ncol = 1, heights = heights_rel, guides = "keep") +
  plot_annotation(title = "Figure 3 │ Cluster Marker Dot-Plots") &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold", hjust = 0.02))

ggsave(file.path(FIG_DIR_PBMC, "fig3_PBMC_L2_ClusFeat_dotplot.png"), panel_fig, width = 150, height = 650, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig3_PBMC_L2_ClusFeat_dotplot.pdf"), panel_fig, width = 150, height = 650, units = "mm", device = "pdf", bg = "transparent")

# ------------------------- Figure 4: Abundance bar charts ------------------
wb <- wb_workbook()

## -- T (CD4 + CD8)
seu_T <- seu[, seu$CellType_l1 %in% c("CD8T","CD4T")]
seu_T$CellType_l2 <- droplevels(seu_T$CellType_l2)
source_cluster_T <- tibble(Cluster = seu_T$CellType_l2, Group = seu_T$group) %>%
  group_by(Cluster, Group) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no)
wb <- wb_add_worksheet(wb, sheet = "group_cluster_T"); wb <- wb_add_data(wb, sheet = "group_cluster_T", x = as.data.frame(source_cluster_T))

p_T <- ggplot(source_cluster_T, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") + scale_fill_manual(values = colorpanel) + coord_fixed(ratio = 1/10) +
  theme_bw() + xlab("group") + ylab("%")

ggsave(file.path(FIG_DIR_PBMC, "fig4_PBMC_L2_cluster_abundance_T.png"), p_T, width = 100, height = 100, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig4_PBMC_L2_cluster_abundance_T.pdf"), p_T, width = 100, height = 100, units = "mm", device = "pdf", bg = "transparent")

## -- NK
seu_NK2 <- seu[, seu$CellType_l1 == "NK"]
seu_NK2$CellType_l2 <- droplevels(seu_NK2$CellType_l2)
source_cluster_NK <- tibble(Cluster = seu_NK2$CellType_l2, Group = seu_NK2$group) %>%
  group_by(Cluster, Group) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no)
wb <- wb_add_worksheet(wb, sheet = "group_cluster_NK"); wb <- wb_add_data(wb, sheet = "group_cluster_NK", x = as.data.frame(source_cluster_NK))

p_NK <- ggplot(source_cluster_NK, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") + scale_fill_manual(values = colorpanel) + coord_fixed(ratio = 1/10) +
  theme_bw() + xlab("group") + ylab("%")

ggsave(file.path(FIG_DIR_PBMC, "fig4_PBMC_L2_cluster_abundance_NK.png"), p_NK, width = 100, height = 100, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig4_PBMC_L2_cluster_abundance_NK.pdf"), p_NK, width = 100, height = 100, units = "mm", device = "pdf", bg = "transparent")

## -- B
seu_B2 <- seu[, seu$CellType_l1 == "B"]
seu_B2$CellType_l2 <- droplevels(seu_B2$CellType_l2)
source_cluster_B <- tibble(Cluster = seu_B2$CellType_l2, Group = seu_B2$group) %>%
  group_by(Cluster, Group) %>% summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>% mutate(total.no = sum(no.cell), perc = 100 * no.cell / total.no)
wb <- wb_add_worksheet(wb, sheet = "group_cluster_B"); wb <- wb_add_data(wb, sheet = "group_cluster_B", x = as.data.frame(source_cluster_B))

p_B <- ggplot(source_cluster_B, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") + scale_fill_manual(values = colorpanel) + coord_fixed(ratio = 1/10) +
  theme_bw() + xlab("group") + ylab("%")

ggsave(file.path(FIG_DIR_PBMC, "fig4_PBMC_L2_cluster_abundance_B.png"), p_B, width = 100, height = 100, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig4_PBMC_L2_cluster_abundance_B.pdf"), p_B, width = 100, height = 100, units = "mm", device = "pdf", bg = "transparent")

wb_save(wb, file.path(RES_DIR, "Step3_results_PBMC.xlsx"), overwrite = TRUE)

# ------------------------- Figure 5: Split UMAP (l2) ----------------------
p_umap_split <- DimPlot(seu, reduction = "umap", split.by = "group", group.by = "CellType_l2", label = FALSE) +
  coord_fixed(ratio = 1) + scale_color_manual(values = colorpanel)

ggsave(file.path(FIG_DIR_PBMC, "fig5_PBMC_UMAP_split.png"), p_umap_split, width = 250, height = 150, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig5_PBMC_UMAP_split.pdf"), p_umap_split, width = 250, height = 150, units = "mm", device = "pdf", bg = "transparent")

# ------------------------- Figure 6: Density (by group) -------------------
# dark theme for density background tiles
theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
          legend.background=element_rect(color=NA, fill="black"), legend.key=element_rect(color="white", fill="black"),
          legend.position="none", panel.background=element_rect(fill="black"), panel.border=element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.spacing=unit(0, "lines"),
          strip.background=element_rect(fill="grey30", color="grey10"), strip.text=element_text(color="white"),
          plot.background=element_rect(color="black", fill="black"), plot.title=element_text(color="white"),
          plot.margin=unit(rep(0,4), "lines"))
}

coord <- Embeddings(seu, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1","UMAP_2")
coord <- tibble(ID = rownames(coord), coord)
meta  <- tibble(ID = rownames(seu@meta.data), seu@meta.data)
meta  <- left_join(meta, coord, by = "ID")

p_density <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_viridis(option = "magma") + coord_fixed(ratio = 1) +
  facet_wrap(~ group, nrow = 1) + theme_black()

ggsave(file.path(FIG_DIR_PBMC, "fig6_PBMC_L2_cluster_density.png"), p_density, width = 400, height = 300, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig6_PBMC_L2_cluster_density.pdf"), p_density, width = 400, height = 300, units = "mm", device = "pdf", bg = "transparent")

# ------------------------- Figure 7: Signature scores ---------------------
DefaultAssay(seu) <- "SCT"

## -- T cells (CD4+CD8)
seu_T <- seu[, seu$CellType_l1 %in% c("CD8T","CD4T")]

naivememory <- c("SELL","CCR7","TCF7","LEF1")
proliferation <- c("STMN1","MKI67","TOP2A","BIRC5")
exhaustion <- c("LAG3","HAVCR2","TOX","PDCD1")

sig_T <- list(Naive_Memory = naivememory, Proliferation = proliferation, Exhaustion = exhaustion)
seu_T <- AddModuleScore_UCell(seu_T, features = sig_T, ncores = min(16, NWORKERS), maxRank = 5000)

u_T_cols <- grep("_UCell$", colnames(seu_T@meta.data), value = TRUE)
for (sig in u_T_cols) {
  p <- Nebulosa::plot_density(seu_T, features = sig, method = "wkde", adjust = 2, size = 0.3, joint = FALSE) +
    coord_fixed(ratio = 1) + ggtitle(sig) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_", sig, "_PBMC_T.png")), p, width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_", sig, "_PBMC_T.pdf")), p, width = 100, height = 100, units = "mm", device = "pdf", bg = "transparent")

  grp_fill <- c(PBMC_CAR = "#d5d5e7", PBMC_CC5 = "#efd2d2")
  p_ridge <- RidgePlot(seu_T, features = sig, group.by = "group", stack = FALSE, same.y.lims = TRUE, combine = TRUE) +
    scale_fill_manual(values = grp_fill) + scale_y_discrete(limits = c("PBMC_CC5","PBMC_CAR")) + ggtitle(sig)
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_ridgeplot_", sig, "_PBMC_T.png")), p_ridge, width = 200, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_ridgeplot_", sig, "_PBMC_T.pdf")), p_ridge, width = 200, height = 100, units = "mm", device = "pdf", bg = "transparent")
}

T_gene_panel <- unique(unlist(sig_T))
p_dot_T <- DotPlot(seu_T, group.by = "group", assay = "SCT", features = T_gene_panel) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90"), legend.direction = "horizontal", legend.position = "bottom") +
  scale_color_gradientn(colours = c('#1f618d','#2E86C1','white','#ec7063','#AD272B'))

ggsave(file.path(FIG_DIR_PBMC, "fig7_signature_dotplot_PBMC_T.png"), p_dot_T, width = 150, height = 100, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig7_signature_dotplot_PBMC_T.pdf"), p_dot_T, width = 150, height = 100, units = "mm", device = "pdf", bg = "transparent")

## -- NK cells
seu_NKsig <- seu[, seu$CellType_l1 == "NK"]
cytotoxity <- c("GZMA","GZMB","IFNG","GNLY","FASLG")
prolif_NK  <- c("STMN1","MKI67","TOP2A","PCNA","BIRC5")
sig_NK <- list(Cytotoxity = cytotoxity, Proliferation = prolif_NK)
seu_NKsig <- AddModuleScore_UCell(seu_NKsig, features = sig_NK, ncores = min(16, NWORKERS), maxRank = 5000)

u_NK_cols <- grep("_UCell$", colnames(seu_NKsig@meta.data), value = TRUE)
for (sig in u_NK_cols) {
  p <- Nebulosa::plot_density(seu_NKsig, features = sig, method = "wkde", adjust = 2, size = 0.3, joint = FALSE) +
    coord_fixed(ratio = 1) + ggtitle(sig) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_", sig, "_PBMC_NK.png")), p, width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_", sig, "_PBMC_NK.pdf")), p, width = 100, height = 100, units = "mm", device = "pdf", bg = "transparent")

  grp_fill <- c(PBMC_CAR = "#d5d5e7", PBMC_CC5 = "#efd2d2")
  p_ridge <- RidgePlot(seu_NKsig, features = sig, group.by = "group", stack = FALSE, same.y.lims = TRUE, combine = TRUE) +
    scale_fill_manual(values = grp_fill) + scale_y_discrete(limits = c("PBMC_CC5","PBMC_CAR")) + ggtitle(sig)
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_ridgeplot_", sig, "_PBMC_NK.png")), p_ridge, width = 200, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR_PBMC, paste0("fig7_signature_score_ridgeplot_", sig, "_PBMC_NK.pdf")), p_ridge, width = 200, height = 100, units = "mm", device = "pdf", bg = "transparent")
}

p_dot_NK <- DotPlot(seu_NKsig, group.by = "group", assay = "SCT", features = unique(unlist(sig_NK))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90"), legend.direction = "horizontal", legend.position = "bottom") +
  scale_color_gradientn(colours = c('#1f618d','#2E86C1','white','#ec7063','#AD272B'))

ggsave(file.path(FIG_DIR_PBMC, "fig7_signature_dotplot_PBMC_NK.png"), p_dot_NK, width = 150, height = 100, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR_PBMC, "fig7_signature_dotplot_PBMC_NK.pdf"), p_dot_NK, width = 150, height = 100, units = "mm", device = "pdf", bg = "transparent")

# ----------------- Differential expression, enrichment & volcano plots -----------------
# This section reproduces Supplemental Fig. S4G and exports DEG/GO/GSEA results.
# Comparison is PBMC_CC5 (CAR.CCL5 + PBMC) vs PBMC_CAR (CAR.con + PBMC) within each CellType_l1.

# Output paths
DEG_XLSX <- Sys.getenv("DEG_XLSX", unset = file.path(RES_DIR, "Step3_CellType_DEG_GO_GSEA_results.xlsx"))
VOLCANO_DIR <- Sys.getenv("VOLCANO_DIR", unset = file.path(FIG_DIR_PBMC, "nonCAR_volcano"))
dir.create(VOLCANO_DIR, showWarnings = FALSE, recursive = TRUE)

# Parameters
DEG_IDENT1 <- Sys.getenv("DEG_IDENT1", unset = "PBMC_CC5")  # numerator
DEG_IDENT2 <- Sys.getenv("DEG_IDENT2", unset = "PBMC_CAR")  # denominator
MIN_CELLS_PER_GROUP <- as.numeric(Sys.getenv("MIN_CELLS_PER_GROUP", unset = "25"))
MIN_CELLS_PER_CELLTYPE <- as.numeric(Sys.getenv("MIN_CELLS_PER_CELLTYPE", unset = "100"))
P_ADJ_CUTOFF <- as.numeric(Sys.getenv("P_ADJ_CUTOFF", unset = "0.05"))
LOGFC_THRESH <- as.numeric(Sys.getenv("LOGFC_THRESH", unset = "0.05"))
N_TOP_LABEL <- as.integer(Sys.getenv("N_TOP_LABEL", unset = "10"))

# Confounder feature blacklist (can be extended via an optional RDS file)
BAD_FEATURES <- c(
  grep("^MT-", rownames(seu), value = TRUE),
  grep("^RPL", rownames(seu), value = TRUE),
  grep("^RPS", rownames(seu), value = TRUE),
  grep("^TR[ABDG][VDJC]", rownames(seu), value = TRUE),
  grep("^IG[HKL][VDJC]", rownames(seu), value = TRUE),
  grep("^HLA-", rownames(seu), value = TRUE),
  grep("^HIST", rownames(seu), value = TRUE),
  "MKI67", "TOP2A", "PCNA", "BIRC5", "STMN1", "TYMS", "HMGB2", "TUBA1B",
  "HBB", "HBA1", "HBA2"
)
BLACKLIST_RDS <- Sys.getenv("BLACKLIST_RDS", unset = "")
if (nzchar(BLACKLIST_RDS) && file.exists(BLACKLIST_RDS)) {
  message("Loading additional blacklist features from: ", BLACKLIST_RDS)
  extra_blacklist <- readRDS(BLACKLIST_RDS)
  BAD_FEATURES <- unique(c(BAD_FEATURES, extra_blacklist))
}
BAD_FEATURES <- unique(BAD_FEATURES)

# Run DE within each CellType_l1 subset
DefaultAssay(seu) <- "SCT"
Idents(seu) <- "CellType_l1"

celltypes_l1 <- levels(factor(seu$CellType_l1))
degs_list <- list()
go_up_list <- list()
go_dn_list <- list()
gsea_list <- list()

# MSigDB collections for GSEA
msig_collections <- c("H", "C2", "C5")

for (ct in celltypes_l1) {
  message("DEG: processing CellType_l1 = ", ct)
  sub <- subset(seu, idents = ct)

  if (ncol(sub) < MIN_CELLS_PER_CELLTYPE) {
    message("  -> skipped (", ncol(sub), " cells < MIN_CELLS_PER_CELLTYPE)")
    next
  }

  # Require both groups
  grp_tab <- table(sub$group)
  if (!(DEG_IDENT1 %in% names(grp_tab) && DEG_IDENT2 %in% names(grp_tab))) {
    message("  -> skipped (missing one of the groups: ", DEG_IDENT1, ", ", DEG_IDENT2, ")")
    next
  }
  if (any(grp_tab[c(DEG_IDENT1, DEG_IDENT2)] < MIN_CELLS_PER_GROUP)) {
    message("  -> skipped (insufficient cells per group: ", paste(names(grp_tab), grp_tab, collapse = ", "), ")")
    next
  }

  # Re-run SCT within subset for DE (mirrors raw analysis)
  DefaultAssay(sub) <- "RNA"
  if ("SCT" %in% Assays(sub)) sub[["SCT"]] <- NULL
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  sub <- SCTransform(sub, vst.flavor = "v2", variable.features.n = 3000, verbose = FALSE)
  DefaultAssay(sub) <- "SCT"
  sub <- PrepSCTFindMarkers(sub)

  feats <- setdiff(rownames(sub), BAD_FEATURES)

  deg <- FindMarkers(
    sub,
    group.by = "group",
    ident.1 = DEG_IDENT1,
    ident.2 = DEG_IDENT2,
    features = feats,
    test.use = "wilcox",
    min.pct = 0.10,
    logfc.threshold = LOGFC_THRESH,
    verbose = FALSE
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("symbol")

  degs_list[[ct]] <- deg

  # GO enrichment (BP) for up/down genes
  universe <- feats
  up_genes <- deg %>% filter(p_val_adj < P_ADJ_CUTOFF, avg_log2FC > 0.25) %>% pull(symbol)
  dn_genes <- deg %>% filter(p_val_adj < P_ADJ_CUTOFF, avg_log2FC < -0.25) %>% pull(symbol)

  if (length(up_genes) > 5) {
    go_up <- enrichGO(
      gene = up_genes,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    ) %>% as.data.frame()
    go_up_list[[ct]] <- go_up
  }

  if (length(dn_genes) > 5) {
    go_dn <- enrichGO(
      gene = dn_genes,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    ) %>% as.data.frame()
    go_dn_list[[ct]] <- go_dn
  }

  # GSEA with MSigDB collections
  ranks <- deg$avg_log2FC
  names(ranks) <- deg$symbol
  ranks <- sort(ranks, decreasing = TRUE)
  ranks <- ranks[!is.na(ranks)]

  gsea_list[[ct]] <- list()
  for (coll in msig_collections) {
    msig_df <- msigdbr(species = "Homo sapiens", category = coll)
    pathways <- split(msig_df$gene_symbol, msig_df$gs_name)

    fg <- fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = 10,
      maxSize = 500,
      nperm = 1000
    ) %>%
      as.data.frame() %>%
      arrange(padj)

    gsea_list[[ct]][[coll]] <- fg
  }
}

# Export DEG/GO/GSEA workbook
wb_deg <- wb_workbook()

for (ct in names(degs_list)) {
  # DEG
  sheet_deg <- paste0("DEG_", ct)
  wb_add_worksheet(wb_deg, sheet_deg)
  wb_add_data(wb_deg, sheet_deg, degs_list[[ct]])

  # GO UP/DN
  if (!is.null(go_up_list[[ct]])) {
    sheet_go_up <- paste0("GO_UP_", ct)
    wb_add_worksheet(wb_deg, sheet_go_up)
    wb_add_data(wb_deg, sheet_go_up, go_up_list[[ct]])
  }
  if (!is.null(go_dn_list[[ct]])) {
    sheet_go_dn <- paste0("GO_DN_", ct)
    wb_add_worksheet(wb_deg, sheet_go_dn)
    wb_add_data(wb_deg, sheet_go_dn, go_dn_list[[ct]])
  }

  # GSEA collections
  if (!is.null(gsea_list[[ct]])) {
    for (coll in names(gsea_list[[ct]])) {
      sheet_gsea <- paste0("GSEA_", coll, "_", ct)
      wb_add_worksheet(wb_deg, sheet_gsea)
      wb_add_data(wb_deg, sheet_gsea, gsea_list[[ct]][[coll]])
    }
  }
}

wb_save(wb_deg, file = DEG_XLSX, overwrite = TRUE)
message("DEG/GO/GSEA workbook written: ", DEG_XLSX)

# Volcano plots for major lymphoid subsets
volcano_celltypes <- c("CD8T", "CD4T", "NK")
for (ct in volcano_celltypes) {
  deg <- degs_list[[ct]]
  if (is.null(deg)) {
    message("Volcano: missing DEG table for ", ct, "; skipping.")
    next
  }

  top_genes <- deg %>%
    filter(p_val_adj < 1e-6, abs(avg_log2FC) > 0.5) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(N_TOP_LABEL) %>%
    pull(symbol)

  p <- EnhancedVolcano(
    deg,
    lab = deg$symbol,
    x = "avg_log2FC",
    y = "p_val_adj",
    pCutoff = P_ADJ_CUTOFF,
    FCcutoff = 0.25,
    pointSize = 2.0,
    labSize = 4.0,
    selectLab = top_genes,
    title = paste0(ct, ": ", DEG_IDENT1, " vs ", DEG_IDENT2),
    subtitle = NULL,
    caption = NULL
  )

  ggsave(file.path(VOLCANO_DIR, paste0("Volcano_", ct, ".png")), p, width = 6, height = 6, dpi = 300)
  ggsave(file.path(VOLCANO_DIR, paste0("Volcano_", ct, ".pdf")), p, width = 6, height = 6)
}

# ------------------------- Save object & session info ---------------------
fig_rds <- if (INCLUDE_DATE_IN_NAMES) sprintf("03_Figures_PBMC_%s.rds", format(Sys.Date(), "%Y%m%d")) else "03_Figures_PBMC.rds"
saveRDS(seu, file = file.path(RES_DIR, fig_rds))
sess_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("session_info_step3_%s.txt", format(Sys.Date(), "%Y%m%d")) else "session_info_step3.txt"
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, sess_file))
