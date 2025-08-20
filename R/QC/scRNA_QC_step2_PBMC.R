#!/usr/bin/env Rscript
# scRNA-seq QC Step 2 — per-sample script

suppressPackageStartupMessages({
  library(scater)
  library(edgeR)
  library(scran)
  library(BiocParallel)
  library(SingleR)
  library(scMCA)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(patchwork)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(Azimuth)
  library(matrixStats)
  library(Matrix)
  library(gridExtra)
})

## ---- Compatibility patch: Matrix ≥1.6 with older compiled binaries -------
if (!isGeneric("colSums")) setGeneric("colSums")
if (!isGeneric("rowMeans")) setGeneric("rowMeans")
if (!isGeneric("rowSums")) setGeneric("rowSums")

setMethod("colSums", signature(x = "dgCMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
            base::colSums(as.matrix(x), na.rm = na.rm))
setMethod("colSums", signature(x = "CsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
            base::colSums(as.matrix(x), na.rm = na.rm))
setMethod("rowMeans", signature(x = "dgCMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
            base::rowMeans(as.matrix(x), na.rm = na.rm))
setMethod("rowSums", signature(x = "dgCMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
            base::rowSums(as.matrix(x), na.rm = na.rm))
setMethod("rowSums", signature(x = "CsparseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
            base::rowSums(as.matrix(x), na.rm = na.rm))
## -------------------------------------------------------------------------

# -------------------------- Configurable parameters -----------------------
# You may override these via environment variables; defaults are generic.
#   export SAMPLE="PBMC_CAR" # or "PBMC_CC5"
#   export SPECIES="hs"              # or "mm"
#   export PROJECT_BASE="."
#   export WORK_BASE="03_R/QC"
#   export CR_PATH="/path/to/CR/outs/filtered_feature_bc_matrix"
#   export LOOM_PATH=""               # optional velocyto .loom
#   export NWORKERS=16
#   export INCLUDE_DATE_IN_LOG=FALSE
#   export EXCLUDE_CLUSTERS=""        # e.g. "4,15,22,24"
#   export IMMGEN_REF="/path/to/ImmGen_reference-Heng_2008.RData" # optional
#   export AZIMUTH_REF="/path/to/Azimuth/Human_PBMC"              # optional

set.seed(123)
NWORKERS            <- as.integer(Sys.getenv("NWORKERS", unset = "8"))
bp                  <- MulticoreParam(workers = NWORKERS)
SAMPLE              <- Sys.getenv("SAMPLE", unset = "PBMC_CAR")
SPECIES             <- Sys.getenv("SPECIES", unset = "hs")
PROJECT_BASE        <- Sys.getenv("PROJECT_BASE", unset = ".")
WORK_BASE           <- Sys.getenv("WORK_BASE", unset = file.path("03_R", "QC"))
CR_PATH             <- Sys.getenv("CR_PATH", unset = file.path(PROJECT_BASE, "02_CRcount", SAMPLE, "outs", "filtered_feature_bc_matrix"))
LOOM_PATH           <- Sys.getenv("LOOM_PATH", unset = "")
INCLUDE_DATE_IN_LOG <- as.logical(Sys.getenv("INCLUDE_DATE_IN_LOG", unset = "FALSE"))
EXCLUDE_CLUSTERS    <- Sys.getenv("EXCLUDE_CLUSTERS", unset = "")
IMMGEN_REF          <- Sys.getenv("IMMGEN_REF", unset = "")
AZIMUTH_REF         <- Sys.getenv("AZIMUTH_REF", unset = "")

# Working/figure directories (per-sample)
workdir <- file.path(PROJECT_BASE, WORK_BASE, SAMPLE)
figdir  <- file.path(workdir, "figures")
if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
if (!dir.exists(figdir))  dir.create(figdir)
setwd(workdir)

# -------------------------- Load Step 1 output ----------------------------
# Expect an RDS saved in Step 1: 01_QC_step1_<sample>.rds
norm <- readRDS(file = file.path(workdir, sprintf("01_QC_step1_%s.rds", SAMPLE)))

# ----------------- Identify unqualified clusters & re-threshold -----------
# Optionally tag pre-determined clusters for discard (e.g., doublets/RBC/platelets).
# Provide comma-separated cluster labels via EXCLUDE_CLUSTERS (empty = none).
colData(norm)$discard.cl <- 0L
if (nzchar(EXCLUDE_CLUSTERS)) {
  to_drop <- strsplit(EXCLUDE_CLUSTERS, ",")[[1]]
  to_drop <- trimws(to_drop)
  colData(norm)$discard.cl[ norm$label %in% to_drop ] <- 1L
}

# Dataset-specific QC thresholds (tune as needed for your data/species)
qc.nexprs  <- norm$detected < 1000         # low features
qc.mito    <- norm$subsets_Mito_percent > 10
qc.rbc     <- norm$subsets_RBC_percent > 1
qc.doublet <- norm$DoubletClass == "doublet"
qc.cluster <- norm$discard.cl == 1

discard.new <- qc.nexprs | qc.mito | qc.rbc | qc.doublet | qc.cluster
norm$discard <- discard.new

# ------------------------------- Logging ----------------------------------
log_stub <- sprintf("%s_scRNA_QC_step2", SAMPLE)
log_file <- if (isTRUE(INCLUDE_DATE_IN_LOG)) {
  file.path(workdir, sprintf("%s_%s.log", log_stub, format(Sys.time(), "%Y%m%d")))
} else {
  file.path(workdir, sprintf("%s.log", log_stub))
}
con <- file(log_file, open = "wt"); sink(con); sink(con, type = "message")
cat(sprintf("[%s] QC step 2 started for %s\n", as.character(Sys.time()), SAMPLE))

# Summary table of new discard flags
print("Summary of filtering"); print(table(discard.new))

# UMAP view of kept vs discarded cells
sceumap <- gridExtra::grid.arrange(
  plotReducedDim(norm[, norm$discard == 0], "UMAP", colour_by = "label", text_by = "label") + ggtitle("Cells remained"),
  plotReducedDim(norm[, norm$discard == 1], "UMAP", colour_by = "label", text_by = "label") + ggtitle("Cells discarded"),
  ncol = 2
)
ggsave(filename = file.path(figdir, "08_Final_Filtering_UMAP.png"), plot = sceumap,
       width = 250, height = 250, units = "mm", dpi = 150, device = "png", bg = "white")

# --------------------- Diagnose gene-level loss (lost vs kept) ------------
# Compare average expression between discarded and retained cells.
rownames(rowData(norm)) <- rowData(norm)$Symbol
lost <- calculateAverage(counts(norm)[,  norm$discard])
kept <- calculateAverage(counts(norm)[, !norm$discard])
logged    <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC     <- logged[, 1] - logged[, 2]
abundance <- rowMeans(logged)

png(file = file.path(figdir, "09_filtered_cells_DEGs.png"), width = 250, height = 250, units = "mm", res = 150)
plot(abundance, logFC, xlab = "Average count (log CPM)", ylab = "Log-FC (lost/kept)", pch = 16)
# (No overlay of cell-level mask here; points are genes.)
dev.off()

DEGs <- as.data.frame(lost[logFC > 1])
print("Top genes enriched in lost cells (logFC>1):")
print(head(rowData(norm)$Symbol[which(rowData(norm)$ID %in% rownames(DEGs))], 20))

# ------------------------- Final cell filtering ---------------------------
set.seed(123)
sce <- norm[, (!norm$discard)]
sce$group <- SAMPLE

# Cell cycle scoring (scran::cyclone) with species-specific markers
if (tolower(SPECIES) == "mm") {
  CCpairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package = "scran"))
} else {
  CCpairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))
}
assignments <- cyclone(sce, CCpairs, gene.names = rowData(sce)$ID, BPPARAM = bp)
sce$G1.score   <- assignments$normalized.scores$G1
sce$S.score    <- assignments$normalized.scores$S
sce$G2M.score  <- assignments$normalized.scores$G2M
sce$CCphase    <- assignments$phases
sce$CC.Diff    <- sce$S.score - sce$G2M.score  # recommended regressor in downstream models

# Stable, prefixed cell IDs (match loom/10x short barcodes)
sce$CellID   <- str_c(paste0(SAMPLE, ":"), str_sub(string = colnames(sce), start = 1, end = 16), "x")
colnames(sce) <- sce$CellID

# ------------------ Optional: mouse-only reference annotations ------------
if (tolower(SPECIES) == "mm") {
  rownames(sce) <- rowData(sce)$Symbol

  # scMCA coarse labels (operates on logcounts matrix)
  mca_result <- scMCA(scdata = assays(sce)@listData[["logcounts"]], numbers_plot = 3)
  sce$CellType_mca <- mca_result[["scMCA"]]

  # ImmGen fine labels via SingleR (if reference file is available)
  if (nzchar(IMMGEN_REF) && file.exists(IMMGEN_REF)) {
    load(IMMGEN_REF)  # should load an object named `immgen`
    pred <- SingleR(test = sce, ref = immgen, labels = immgen$label.fine, assay.type.test = 1, BPPARAM = bp)
    sce$CellType_immgen <- pred$pruned.labels
  } else {
    message("IMMGEN_REF not provided or not found; skipping SingleR ImmGen annotation.")
  }
}

# ---------------------- Export per-cell metadata --------------------------
colData_tbl <- as_tibble(sce@colData)
saveRDS(colData_tbl, file = sprintf("%s_colData.rds", SAMPLE))  # no date in filename

# --------------------- Build Seurat object for this sample ----------------
library(future)
plan("multicore", workers = NWORKERS)
options(Seurat.object.assay.version = "v3")
set.seed(123)

# Load raw counts from 10x filtered matrix
counts <- Read10X(data.dir = CR_PATH)
seu    <- CreateSeuratObject(counts = counts)

# Harmonize cell barcodes with SCE (prefix + truncated barcode + 'x')
newnames <- str_c(paste0(SAMPLE, ":"), str_sub(string = Cells(seu), start = 1, end = 16), "x")
seu      <- RenameCells(seu, new.names = newnames)

# Integrity check: expect all metadata CellIDs present among Seurat cells
message("Check CellID integrity between SCE colData and Seurat cells:")
print(table(colData_tbl$CellID %in% Cells(seu)))

# Filter Seurat to QC-passing cells only
seu <- seu[, Cells(seu) %in% colData_tbl$CellID]

# Optional: add RNA velocity assays from .loom, if available
if (nzchar(LOOM_PATH) && file.exists(LOOM_PATH)) {
  loom <- ReadVelocity(file = LOOM_PATH)
  seu_loom <- as.Seurat(loom)
  seu_loom <- subset(seu_loom, cells = Cells(seu))
  spliced_assay   <- CreateAssayObject(counts = seu_loom@assays$spliced$counts[, Cells(seu)])
  unspliced_assay <- CreateAssayObject(counts = seu_loom@assays$unspliced$counts[, Cells(seu)])
  stopifnot(all(Cells(seu) == colnames(spliced_assay)), all(Cells(seu) == colnames(unspliced_assay)))
  seu[["spliced"]]   <- spliced_assay
  seu[["unspliced"]] <- unspliced_assay
}

# Join per-cell metadata from SCE into Seurat
seu$CellID <- Cells(seu)
seu@meta.data <- dplyr::left_join(as_tibble(seu@meta.data), colData_tbl, by = "CellID") %>%
  column_to_rownames("CellID")

# Normalization (SCTransform v2) + PCA
vars_to_regress <- c("subsets_Mito_percent", "S.score", "G2M.score")
seu <- SCTransform(seu, vst.flavor = "v2", verbose = TRUE, return.only.var.genes = FALSE,
                   vars.to.regress = vars_to_regress) %>%
       RunPCA(npcs = 30, verbose = TRUE)

# Optional: human PBMC annotation via Azimuth (local reference directory)
options(timeout = 10000)
if (tolower(SPECIES) == "hs" && nzchar(AZIMUTH_REF) && dir.exists(AZIMUTH_REF)) {
  seu.annot <- RunAzimuth(seu, reference = AZIMUTH_REF)
  seu@meta.data <- seu.annot@meta.data
} else if (tolower(SPECIES) == "hs") {
  message("AZIMUTH_REF not provided or not found; skipping Azimuth annotation.")
}

# Save per-sample Seurat object (no date in filename)
saveRDS(seu, file = sprintf("%s_seu.rds", SAMPLE))
message(sprintf("Job %s successfully completed at %s.", SAMPLE, as.character(Sys.time())))

# Close log and clean up
sink(); sink(type = "message")
rm(list = ls()); graphics.off(); gc()
