# Self-Amplifying CCL5–CCR5 Circuit Drives CAR T Swarming and Tumor Immunity Normalization

**Abstract.**  
Poor trafficking and limited persistence hinder chimeric antigen receptor T cell (CAR-T) therapy in solid tumors. Here we convert CAR-T from solitary effectors into cooperative swarms by hard-wiring a stimulus-gated CCL5–CCR5 circuit. Because CAR/T-cell-receptor engagement rapidly releases pre-formed CCL5 yet represses its transcription, constitutive CCL5 co-expression restored a burst-on-demand chemokine pulse upon antigen engagement with minimal leakage. CCL5-armored CAR-T markedly increased peer recruitment, streamed from tumor-spheroid rims to tumor cores, and achieved remarkable infiltration in orthotopic tumor exografts in a CCL5-CCR5-dependent manner. Swarming boosted IFN-γ production while reducing PD-1/TIM-3 co-expression in CAR-T, recruited endogenous CCR5⁺ lymphoid cytotoxic effectors, and re-programmed suppressive myeloid populations toward MHC-II–rich phenotypes, collectively normalizing the tumor microenvironment and culminating in durable tumor control. Therefore, embedding a CCL5–CCR5 feedback loop rewires CAR-T into self-organizing swarms and highlights programmable chemokine circuits as plug-and-play tools for coordinated trafficking and durable immune normalization in solid-tumor immunotherapy.

---

## Repository overview

This repo hosts lightly annotated, de-personalized R scripts to reproduce the single-cell RNA-seq (scRNA-seq) analysis used in the study. The pipeline spans **droplet QC**, **sample-level filtering**, **integration**, **annotation**, and **figure generation** for two cohorts: **CAR-T** and **PBMC**.


- **Inputs:** 10x Genomics `filtered_feature_bc_matrix` directories per sample.  
- **Outputs:** RDS objects, PNG/PDF figures, and XLSX tables are written under `R/results/` and `R/figures/`.

---


