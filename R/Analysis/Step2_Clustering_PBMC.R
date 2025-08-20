#!/usr/bin/env Rscript
# PBMC scRNA-seq — Step 2: Clustering, marker discovery, and annotation

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
# ------------------------- Config -----------------------------------------
NWORKERS              <- as.integer(Sys.getenv("NWORKERS", unset = "20"))
PLAN_STRATEGY         <- Sys.getenv("PLAN", unset = "sequential")  # or "multicore"
PROJECT_BASE          <- Sys.getenv("PROJECT_BASE", unset = ".")
WORKDIR               <- Sys.getenv("WORKDIR", unset = file.path(PROJECT_BASE, "03_R"))
RES_DIR               <- file.path(WORKDIR, "results")
FIG_DIR               <- file.path(WORKDIR, "figures")
INTEGRATED_RDS        <- Sys.getenv("INTEGRATED_RDS", unset = file.path(RES_DIR, "01_Integration_PBMC.rds"))
GENE_BLACKLIST_RDS    <- Sys.getenv("GENE_BLACKLIST_RDS", unset = "")  # optional; if not provided, skip
INCLUDE_DATE_IN_NAMES <- as.logical(Sys.getenv("INCLUDE_DATE_IN_NAMES", unset = "FALSE"))

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
setwd(WORKDIR)
plan(PLAN_STRATEGY, workers = NWORKERS)

# ------------------------- Load object ------------------------------------
seu <- readRDS(INTEGRATED_RDS)

# ------------------------- Section 1: l1-type markers ---------------------
# Quick lineage sanity panel (DotPlot saved is optional; here we only compute)
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
          axis.title = element_blank(), legend.direction = "vertical", legend.position = "bottom") +
    scale_color_gradientn(colours = c('white','#fde3d8','#ed684e','#d21e20','#981917'))
  ggsave(file.path(FIG_DIR, "PBMC_Step2_dotplot_l1.png"), p_l1, width = 180, height = 160, units = "mm", dpi = 150, bg = "white")
}, silent = TRUE)

# Coarse l1 annotation (cluster-level) — keep your mapping
seu$CellType_l1 <- factor(dplyr::recode(
  seu$seurat_clusters,
  "1"="CD4T","2"="CD4T","3"="CD8T","4"="CD8T","5"="CD8T","6"="NK","7"="CD8T",
  "8"="CD8T","9"="NK","10"="CD4T","11"="CD4T","12"="CD8T","13"="NK","14"="CD8T",
  "15"="B","16"="CD8T","17"="B","18"="NK","19"="CD8T","20"="NK"
))

# ------------------------- Section 2: conserved markers (per l1) ----------
# UMAP by clusters for reference
p_umap <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + coord_fixed(ratio = 1)
ggsave(file.path(FIG_DIR, "PBMC_Step2_umap_clusters.png"), p_umap, width = 160, height = 140, units = "mm", dpi = 150, bg = "white")

# Ensure RNA normalized for DE
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)

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

extra_blacklist <- if (nzchar(GENE_BLACKLIST_RDS) && file.exists(GENE_BLACKLIST_RDS)) {
  tryCatch(readRDS(GENE_BLACKLIST_RDS), error = function(e) character())
} else character()

bad_features <- unique(c(
  extra_blacklist,
  # cc_genes,  # comment out if you prefer to keep cc genes
  hist_genes, hb_genes, mt_genes, rps_genes, rik_genes,
  alu_genes, pseudo_genes, mir_genes, gencode_genes
))

bad_features <- intersect(bad_features, rn)
features <- setdiff(rn, bad_features)

# Per-CellType_l1 conserved markers across groups
seu$seurat_clusters <- factor(seu$seurat_clusters)
if (!"group" %in% colnames(seu@meta.data)) seu$group <- "all"

celltypes_l1 <- levels(seu$CellType_l1)
conserved_markers_all <- list()
wb <- wb_workbook()

plan(PLAN_STRATEGY, workers = NWORKERS)
for (ct in celltypes_l1) {
  message("=== Processing CellType L1: ", ct, " ===")
  seu_sub <- subset(seu, subset = CellType_l1 == ct)
  seu_sub$seurat_clusters <- droplevels(seu_sub$seurat_clusters)
  cl_ids <- levels(seu_sub$seurat_clusters)
  Idents(seu_sub) <- "seurat_clusters"

  conserved_list <- lapply(cl_ids, function(cl) {
    FindConservedMarkers(
      object = seu_sub,
      ident.1 = cl,
      features = features,
      grouping.var = "group",
      only.pos = TRUE,
      min.pct = 0.30,
      min.cells.per.ident = 10,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      combine = "replace",
      densify = TRUE
    ) %>% tibble::rownames_to_column("gene") %>% dplyr::mutate(cluster = cl, CellType_l1 = ct)
  })

  ct_df <- dplyr::bind_rows(conserved_list) %>% dplyr::group_by(cluster) %>% dplyr::arrange(max_pval, .by_group = TRUE) %>% dplyr::ungroup()
  wb <- wb_add_worksheet(wb, sheet = ct)
  wb <- wb_add_data(wb, sheet = ct, x = as.data.frame(ct_df))
  conserved_markers_all[[ct]] <- ct_df
}
plan("sequential")

all_df <- dplyr::bind_rows(conserved_markers_all)
wb <- wb_add_worksheet(wb, sheet = "ALL")
wb <- wb_add_data(wb, sheet = "ALL", x = as.data.frame(all_df))
wb_save(wb, file.path(RES_DIR, "Step2_cluster_l2_byCellType.xlsx"), overwrite = TRUE)
message("✓ Done → results/Step2_cluster_l2_byCellType.xlsx")

# ------------------------- Section 3: l2 annotation -----------------------
# Map clusters to detailed l2 labels (keep your mapping) and order them
seu$CellType_l2 <- dplyr::recode(
  seu$seurat_clusters,
  # CD8 T
  "4"="CD8T_Tn_CCR7","3"="CD8T_TCM_IL7R","5"="CD8T_TCM_Cycling","7"="CD8T_PreDys_GZMK",
  "19"="CD8T_PreDys_CSF2","16"="CD8T_Tex_TNFRSF4","8"="CD8T_Tex_TNFRSF9",
  "14"="CD8T_MAIT_CCL20","12"="CD8T_MAIT_Exh_CTLA4",
  # CD4/NK/B
  "1"="CD4T_Tn_LEF1","11"="CD4T_TCM_KLF2","2"="CD4T_EarlyAct_TNFRSF4","10"="CD4T_Cycling",
  "6"="NK_Cycling","20"="NK_EarlyAct_NR4A","9"="NK_IFNprimed_LDB2","18"="NK_Conventional_Cytotoxic",
  "13"="B_naive_FCER2","17"="B_InnateLike_CD5","15"="B_Inflammatory_STAT1"
)

l2_order <- c(
  # CD8 T (naïve → tex)
  "CD8T_Tn_CCR7","CD8T_TCM_IL7R","CD8T_TCM_Cycling","CD8T_PreDys_GZMK","CD8T_PreDys_CSF2",
  "CD8T_Tex_TNFRSF4","CD8T_Tex_TNFRSF9",
  # MAIT
  "CD8T_MAIT_CCL20","CD8T_MAIT_Exh_CTLA4",
  # CD4 T
  "CD4T_Tn_LEF1","CD4T_TCM_KLF2","CD4T_EarlyAct_TNFRSF4","CD4T_Cycling",
  # NK
  "NK_Cycling","NK_EarlyAct_NR4A","NK_IFNprimed_LDB2","NK_Conventional_Cytotoxic",
  # B
  "B_naive_FCER2","B_Inflammatory_STAT1"
)
seu$CellType_l2 <- factor(seu$CellType_l2, levels = l2_order)

# ------------------------- Section 4: conserved markers (l2) --------------
celltypes_l1 <- levels(seu$CellType_l1)
conserved_markers_all <- list()
wb2 <- wb_workbook()

plan(PLAN_STRATEGY, workers = NWORKERS)
for (ct in celltypes_l1) {
  message("=== Processing CellType L1 (l2 stage): ", ct, " ===")
  seu_sub <- subset(seu, subset = CellType_l1 == ct)
  seu_sub$CellType_l2 <- droplevels(seu_sub$CellType_l2)
  cl_ids <- levels(seu_sub$CellType_l2)
  Idents(seu_sub) <- "CellType_l2"

  conserved_list <- lapply(cl_ids, function(cl) {
    FindConservedMarkers(
      object = seu_sub, ident.1 = cl, features = features, grouping.var = "group",
      only.pos = TRUE, min.pct = 0.30, min.cells.per.ident = 10, logfc.threshold = 0.25,
      test.use = "wilcox", combine = "replace", densify = TRUE
    ) %>% tibble::rownames_to_column("gene") %>% dplyr::mutate(cluster = cl, CellType_l1 = ct)
  })

  ct_df <- dplyr::bind_rows(conserved_list) %>% dplyr::group_by(cluster) %>% dplyr::arrange(max_pval, .by_group = TRUE) %>% dplyr::ungroup()
  wb2 <- wb_add_worksheet(wb2, sheet = ct)
  wb2 <- wb_add_data(wb2, sheet = ct, x = as.data.frame(ct_df))
  conserved_markers_all[[ct]] <- ct_df
}
plan("sequential")

all_df2 <- dplyr::bind_rows(conserved_markers_all)
wb2 <- wb_add_worksheet(wb2, sheet = "ALL")
wb2 <- wb_add_data(wb2, sheet = "ALL", x = as.data.frame(all_df2))
wb_save(wb2, file.path(RES_DIR, "Step2_cluster_l2_byCellType.xlsx"), overwrite = TRUE)
message("✓ Done (l2) → results/Step2_cluster_l2_byCellType.xlsx")

# ------------------------- Save annotated object & session info ------------
res_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("02_Annotation_PBMC_%s.rds", format(Sys.Date(), "%Y%m%d")) else "02_Annotation_PBMC.rds"
saveRDS(seu, file = file.path(RES_DIR, res_file))

sess_file <- if (INCLUDE_DATE_IN_NAMES) sprintf("session_info_step2_%s.txt", format(Sys.Date(), "%Y%m%d")) else "session_info_step2.txt"
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, sess_file))
