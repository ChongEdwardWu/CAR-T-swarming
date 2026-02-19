# CAR‑T‑swarming
Code for the single-cell RNA-seq analyses used in the manuscript “Activation-gated, self-reinforcing CCL5–CCR5 relay drives CAR T cell swarming and immune remodeling in solid tumors”.

This repository contains R scripts to perform quality control, doublet detection, SCTransform v2 normalization, multi-sample integration, clustering, marker discovery, differential expression, gene signature scoring (UCell), and figure generation for two scRNA-seq datasets

- intratumoral human CAR T cells (tNGFR+)
- intratumoral human non-CAR leukocytes (hCD45+ tNGFR−)

---
**Manuscript title**  
**Activation-gated, self-reinforcing CCL5–CCR5 relay drives CAR T cell swarming and immune remodeling in solid tumors**

**Abstract**  
Limited trafficking and dysfunction constrain CAR T cell therapy in solid tumors. Here, we developed an activation-gated, self-reinforcing strategy for tumor-localized CAR T cell swarming. CCR5 was enriched on CD8+ CAR T cells and CCL5 was poised for activation-coupled release, yet ex vivo activation and expansion depleted preformed CCL5 stores and reduced CCL5 transcripts in CAR T products. Constitutive CCL5 expression restored CAR-gated chemokine pulses with minimal basal secretion. Recruited CCR5+ CAR T cells amplified these pulses upon antigen encounter, establishing a self-reinforcing relay that promoted CAR T cell swarming into collagen-embedded tumor spheroids, orthotopic hepatocellular carcinoma, and subcutaneous tumors. In xenografts, swarming CAR T cells retained anti-tumor activity and, with PBMC co-infusion, recruited bystander lymphocytes; in immunocompetent hosts without lymphodepletion, they remodeled the tumor immune microenvironment, improved antigen-positive tumor control, and restrained contralateral antigen-negative lesions. Activation-gated chemokine relays thus offer a modular strategy to coordinate engineered and host immunity against solid tumors.

---

## **Repository overview**

This repository contains lightly annotated, de‑personalized code to reproduce the single‑cell RNA‑seq analyses used in the study. The code is organized to run per‑sample QC, cohort‑level integration, clustering/annotation, and figure generation in a reproducible sequence.

**Note**: All absolute paths in scripts were anonymized to placeholders (e.g., path_to_data, path_to_.radian_profile). Date suffixes in file names were removed for portability. Replace these with your local paths before running.

---

### Directory layout

```
CAR-T-swarming/
├─ R/
│  ├─ Analysis/
│  │  ├─ Step1_Integration_CART.R       # Integration (Seurat SCT/RPCA) + DR/UMAP for CAR‑T
│  │  ├─ Step1_Integration_PBMC.R       # Integration (Seurat SCT/RPCA) + DR/UMAP for PBMC
│  │  ├─ Step2_Clustering_CART.R        # Clustering, automated labels, markers (CAR‑T)
│  │  ├─ Step2_Clustering_PBMC.R        # Clustering, automated labels, markers (PBMC)
│  │  ├─ Step3_Figures_CART.R           # Final figure panels (UMAP, dot plots, abundances) for CAR‑T
│  │  └─ Step3_Figures_PBMC.R           # Final figure panels for PBMC
│  └─ QC/
│     ├─ scRNA_QC_loop_step1.R          # Per-sample pre‑QC: filtering metrics, doublets, DR snapshot
│     ├─ scRNA_QC_step2_CART.R          # Post‑QC + SCTransform + metadata merge (CAR‑T)
│     └─ scRNA_QC_step2_PBMC.R          # Post‑QC + SCTransform + metadata merge (PBMC)
└─ README.md

```

---

### Pipeline overview

```
**CellRanger → CellBender → QC step1 → QC step2 → Integration
      │          │             │           │           │
      └─ run_cellbender_batch.sh
                 │
      └─ scRNA_QC_loop_CB_step1.R
                 │
      └─ scRNA_QC_loop_CB_step2.R        (per-sample Seurat RDS)
                 │
      └─ 01_hepaLSK_Data_Integration.R   (integrated Seurat)
                 │
      └─ 01.5_hepaLSK_pySCENIC.R         (export LOOM → run pySCENIC → import AUC/Bin)
                 │
      └─ 02_hepaLSK_Clustering.r         (annotations & marker tables)
                 │
      └─ 03_hepaLSK_Group_Comparisons.r  (RNA/AUC/Bin DE; intersections)
                 │
      └─ 04_hepaLSK_Trajectory.r         (Slingshot trajectory)
                 │
      └─ 05_hepaLSK_final.r              (final figures & UCell trends)**
```
---

## Input data

This repository does not include raw sequencing data.

The scripts expect Cell Ranger `count` outputs for each sample, and specifically use the `filtered_feature_bc_matrix/` directory produced by Cell Ranger (10x Genomics 3′ v3 chemistry in the manuscript). Please refer to the manuscript for the experimental design, cell sorting strategy, and sequencing details.

---

## Software requirements

The analysis was developed and tested in a recent R 4.x environment. Key packages include, but are not limited to

- `Seurat` and `SeuratObject`
- `SingleCellExperiment`, `scater`, `scran`
- `scDblFinder`
- `sctransform`
- `UCell`, `Nebulosa`
- `ggplot2`, `patchwork`, `cowplot`
- `dplyr`, `tibble`, `stringr`, `tidyr`
- `openxlsx`

Please refer to the `library()` calls at the top of each script for the full set of dependencies.

---

## Citation

If you use these scripts, please cite the corresponding manuscript above.  
Please also cite relevant software tools (Seurat, scran, pySCENIC, Slingshot, UCell, etc.) as appropriate.

---

## License

This repository is provided for academic and non-commercial use.  
For redistribution or commercial use, please open an issue to discuss licensing.

---

## Contact

For any questions or reproducibility requests, please email to wuchong5@mail(dot)sysu(dot)edu(dot)cn.

---

**Acknowledgements**  
We thank developers of Seurat, Bioconductor’s single-cell packages, pySCENIC, Slingshot, and UCell for their foundational tools enabling this analysis.
