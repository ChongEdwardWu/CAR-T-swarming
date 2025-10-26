# CAR‑T‑swarming

**Manuscript title**  
**Self‑Amplifying CCL5–CCR5 Circuit Drives CAR T Swarming and Tumor Immunity Normalization**

**Abstract**  
Poor trafficking and limited persistence hinder chimeric antigen receptor T cell (CAR-T) therapy in solid tumors. Here we convert CAR-T from solitary effectors into cooperative swarms by hard-wiring a stimulus-gated CCL5–CCR5 circuit. Because CAR/T-cell-receptor engagement rapidly releases pre-formed CCL5 yet represses its transcription, constitutive CCL5 co-expression restored a burst-on-demand chemokine pulse upon antigen engagement with minimal leakage. CCL5-armored CAR-T markedly increased peer recruitment, streamed from tumor-spheroid rims to tumor cores, and achieved remarkable infiltration in orthotopic tumor exografts in a CCL5-CCR5-dependent manner. Swarming boosted IFN-γ production while reducing PD-1/TIM-3 co-expression in CAR-T, recruited endogenous CCR5⁺ lymphoid cytotoxic effectors, and re-programmed suppressive myeloid populations toward MHC-II–rich phenotypes, collectively normalizing the tumor microenvironment and culminating in durable tumor control. Therefore, embedding a CCL5–CCR5 feedback loop rewires CAR-T into self-organizing swarms and highlights programmable chemokine circuits as plug-and-play tools for coordinated trafficking and durable immune normalization in solid-tumor immunotherapy.


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
