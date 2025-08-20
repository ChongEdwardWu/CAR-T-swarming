#!/usr/bin/env Rscript
# CAR T scRNA-seq — Step 2: Clustering, marker discovery, and annotation

# ------------------------- Setup & parameters ------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(future)
  library(stringr)
  library(SeuratObject)
  library(tibble)
  library(openxlsx2)
})

set.seed(123)
NWORKERS              <- as.integer(Sys.getenv("NWORKERS", unset = "20"))
PLAN_STRATEGY         <- Sys.getenv("PLAN", unset = "sequential")  # or "multicore"
PROJECT_BASE          <- Sys.getenv("PROJECT_BASE", unset = ".")
WORKDIR               <- Sys.getenv("WORKDIR", unset = file.path(PROJECT_BASE, "03_R"))
RES_DIR               <- file.path(WORKDIR, "results")
FIG_DIR               <- file.path(WORKDIR, "figures")
INTEGRATED_RDS        <- Sys.getenv("INTEGRATED_RDS", unset = file.path(RES_DIR, "01_Integration_CART.rds"))
INCLUDE_DATE_IN_NAMES <- as.logical(Sys.getenv("INCLUDE_DATE_IN_NAMES", unset = "FALSE"))
GENE_BLACKLIST_RDS    <- Sys.getenv("GENE_BLACKLIST_RDS", unset = "")  # optional RDS path for TCR/BCR/Ig/etc.

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
setwd(WORKDIR)

plan(PLAN_STRATEGY, workers = NWORKERS)

# ------------------------- Load integrated object -------------------------
seu <- readRDS(INTEGRATED_RDS)

# ------------------------- Section 1: l1-type markers ---------------------
# Minimal lineage panel for quick sanity check via DotPlot.
marker_sets_l1 <- list(
  lineage_genes = c(
    # TCR / CD3
    "CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2",
    # NK receptors / granules
    "NCR1","KLRD1","KLRF1","NKG7",
    # Canonical CD markers
    "CD4","CD8A","CD8B",
    # MAIT
    "SLC4A10","RORC","IL23R",
    # B cell
    "CD79A","MS4A1"
  )
)

gene_panel <- unique(unlist(marker_sets_l1))
gene_panel <- gene_panel[gene_panel %in% rownames(seu)]

try({
  p_l1 <- DotPlot(seu, assay = "SCT", features = gene_panel) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_blank(),
          legend.direction = "vertical",
          legend.position = "bottom") +
    scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))
  ggsave(file.path(FIG_DIR, "Step2_dotplot_l1.png"), p_l1, width = 180, height = 160, units = "mm", dpi = 150, bg = "white")
}, silent = TRUE)

# Coarse l1 annotation (cluster-level). Keep the original mapping.
seu$CellType_l1 <- factor(dplyr::recode(
  seu$seurat_clusters,
  "1"="CD8T","2"="CD8T","3"="CD8T","4"="CD8T","5"="CD8T","6"="CD8T",
  "7"="CD8T","8"="CD4T","9"="CD8T"
))

# ------------------------- Section 2: l2 markers --------------------------
# UMAP by clusters for reference
p_umap <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + coord_fixed(ratio = 1)
ggsave(file.path(FIG_DIR, "Step2_umap_clusters.png"), p_umap, width = 160, height = 140, units = "mm", dpi = 150, bg = "white")

# Build blacklist for DE (remove confounders like histone/mito/ribo etc.)
rn <- if ("RNA" %in% names(seu@assays)) rownames(seu@assays$RNA@counts) else rownames(seu)
cc_genes      <- unique(unlist(Seurat::cc.genes))
hist_genes    <- grep("^Hist", rn, ignore.case = TRUE, value = TRUE)
hb_genes      <- grep("^Hb[ab]-|^HB[^(P)]", rn, value = TRUE)
mt_genes      <- grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l", rn, ignore.case = TRUE, value = TRUE)
rps_genes     <- grep("^Rp[sl]|^RP[SL]", rn, ignore.case = TRUE, value = TRUE)
rik_genes     <- grep("^Rik", rn, ignore.case = TRUE, value = TRUE)
alu_genes     <- grep("^AL", rn, ignore.case = TRUE, value = TRUE)
pseudo_genes  <- grep("-rs|-ps", rn, ignore.case = TRUE, value = TRUE)
mir_genes     <- grep("^Mir", rn, ignore.case = TRUE, value = TRUE)
gencode_genes <- grep("^Gm", rn, ignore.case = TRUE, value = TRUE)

# Optional external blacklist (e.g., TCR/BCR/Ig lists). If not found, skip gracefully.
extra_blacklist <- if (nzchar(GENE_BLACKLIST_RDS) && file.exists(GENE_BLACKLIST_RDS)) {
  tryCatch(readRDS(GENE_BLACKLIST_RDS), error = function(e) character())
} else character()

bad_features <- unique(c(
  extra_blacklist,
  # cc_genes,   # keep commented per the original script
  hist_genes, hb_genes, mt_genes, rps_genes, rik_genes,
  alu_genes, pseudo_genes, mir_genes, gencode_genes
))

bad_features <- intersect(bad_features, rn)
features <- setdiff(rn, bad_features)

# Conserved markers per cluster across groups (per-sample test)
seu$seurat_clusters <- factor(seu$seurat_clusters)
Idents(seu) <- "seurat_clusters"

# Ensure grouping var exists; if missing, create a dummy
if (!"group" %in% colnames(seu@meta.data)) seu$group <- "all"

DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)
cluster_ids <- levels(seu$seurat_clusters)

conserved_markers <- lapply(cluster_ids, function(cl) {
  FindConservedMarkers(
    seu,
    ident.1       = cl,
    features      = features,
    grouping.var  = "group",
    only.pos      = TRUE,
    min.pct       = 0.3,
    min.cells.per.ident = 10,
    logfc.threshold = 0.25,
    test.use      = "wilcox",
    combine       = 'replace',  # fill missing groups with zero instead of dropping
    densify       = TRUE
  ) %>% dplyr::mutate(cluster = cl, gene = rownames(.))
})

conserved_markers_df <- bind_rows(conserved_markers) %>%
  dplyr::group_by(cluster) %>% dplyr::arrange(max_pval, .by_group = TRUE)

# Export markers (cluster-wise)
wb <- wb_workbook()
wb <- wb_add_worksheet(wb, sheet = "all")
wb <- wb_add_data(wb, sheet = "all", x = as.data.frame(conserved_markers_df))
wb_save(wb, file.path(RES_DIR, "Step2_cluster_l2_CART.xlsx"), overwrite = TRUE)

# ------------------------- Section 3: l2 annotation -----------------------
# Optional cluster order for readability in figures
desired_order <- c("1","4","6","2","3","5","7","9","8")
seu$seurat_clusters <- factor(seu$seurat_clusters, levels = desired_order)

# Final l2 labels (keep the original names)
seu$CellType_l2 <- dplyr::recode(
  seu$seurat_clusters,
  "1"="CD8T_Cycling_MKI67",
  "4"="CD8T_PolyCTL_IFNG",
  "6"="CD8T_TCM_CCR7",
  "2"="CD8T_HyperAct_MHCII",
  "3"="CD8T_Tex_TOX",
  "5"="CD8T_Tex_TOX",
  "7"="CD8T_Tex_TOX",
  "9"="CD8T_Tex_TOX",
  "8"="CD4T_Helper_TNFRSF4"
)

# Order l2 labels for downstream plots
l2_order <- c("CD8T_Cycling_MKI67","CD8T_PolyCTL_IFNG","CD8T_TCM_CCR7",
              "CD8T_HyperAct_MHCII","CD8T_Tex_TOX","CD4T_Helper_TNFRSF4")
seu$CellType_l2 <- factor(seu$CellType_l2, levels = l2_order)

# ------------------------- Section 4: re-run conserved markers at l2 ------
Idents(seu) <- "CellType_l2"
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)

cluster_ids2 <- levels(seu$CellType_l2)
conserved_markers2 <- lapply(cluster_ids2, function(cl) {
  FindConservedMarkers(
    seu,
    ident.1       = cl,
    features      = features,
    grouping.var  = "group",
    only.pos      = TRUE,
    min.pct       = 0.3,
    min.cells.per.ident = 10,
    logfc.threshold = 0.25,
    test.use      = "wilcox",
    combine       = 'replace',
    densify       = TRUE
  ) %>% dplyr::mutate(cluster = cl, gene = rownames(.))
})

conserved_markers_df2 <- bind_rows(conserved_markers2) %>%
  dplyr::group_by(cluster) %>% dplyr::arrange(max_pval, .by_group = TRUE)

wb2 <- wb_workbook()
wb2 <- wb_add_worksheet(wb2, sheet = "all")
wb2 <- wb_add_data(wb2, sheet = "all", x = as.data.frame(conserved_markers_df2))
wb_save(wb2, file.path(RES_DIR, "Step2_cluster_l2_CART.xlsx"), overwrite = TRUE)

# Summary DotPlot for l2 classes
try({
  p_l2_final <- DotPlot(seu, group.by = "CellType_l2", assay = "SCT", features = rev(gene_panel2)) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title   = element_blank(),
          panel.grid.major = element_line(color = "grey90"),
          legend.direction = "horizontal",
          legend.position  = "bottom") +
    scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))
  ggsave(file.path(FIG_DIR, "Step2_dotplot_l2_final.png"), p_l2_final, width = 200, height = 160, units = "mm", dpi = 150, bg = "white")
}, silent = TRUE)

# ------------------------- Save annotated object & session info ------------
res_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("02_Annotation_CART_%s.rds", format(Sys.Date(), "%Y%m%d")) else "02_Annotation_CART.rds"
saveRDS(seu, file = file.path(RES_DIR, res_file))

sess_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("session_info_step2_%s.txt", format(Sys.Date(), "%Y%m%d")) else "session_info_step2.txt"
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, sess_file))
