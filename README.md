# mola-fsel-nf — Module 1 (MOFA → LASSO)

**Goal:** explain **why specific genes are selected (mutated) in a given tissue** by integrating multi-omics and identifying features that separate *gene-mut* vs *gene-WT* samples.

---

## What this pipeline does

- **Download & stage TCGA-like data** (Expression, CNV, Mutation).
- **Preprocess**
  - Expression: log2(TPM + 1)
  - CNV: thresholded/continuous as configured
  - Mutation: per-gene event counts (e.g., MISSENSE / TRUNC)
  - Align sample IDs and take the intersection across views.
- **Integrate omics with MOFA** (per cohort)
  - Train an initial (raw) model and inspect variance explained.
  - Retrain a final model with **selected factors** and export variance plots, factors and weights.
- **Per-gene LASSO** (per cohort × gene)
  - Build labels for a target gene (e.g., **TP53_mut** vs **TP53_wt**).
  - Use a **no-leak** mutation design (drop target-gene mutation columns) and run L1-logistic regression to rank **driver features**.
- **Predefined genes**
  - `--predef` forces genes into MOFA inputs / top-weights (expanded per view; e.g., `_MIS`, `_TRUNC` for Mutation).

> **Tested on HPC (SLURM)**; also runnable locally for small tests.

---

## Example run (SLURM)


```bash
nextflow run wf_main.nf -profile slurm  \
  --projects TCGA-BRCA  \
  --genes TP53,PIK3CA  \
  --predef "BRCA1,BRCA2,PIK3CA,KEAP1,EGFR,KRAS,TP53"  \
  --results_base /scratch/yazbecka/mola-fsel-nf/results_BRCA_release44 \
  --data_base    /scratch/yazbecka/mola-fsel-nf/data/raw/  \
  --nfolds 5 --seed 1 --mofa_topn 100  \
  --topn_expr 20 --topn_mut 20 --topn_cnv 20  \
  -with-dag dagRelease44.svg -with-report reportRelease44.html
```

```bash
data/<PROJECT>/ # from --data_base (cache/prepared inputs)
  ├─ raw/
      ├─ cnv_thresholded_by_genes.rds
      ├─ cnv_thresholded_by_genes.tsv
      ├─ expression_se.rds
      ├─ mutations_df.rds;
      └─ mutations_maf.rds
results/<PROJECT>/
  ├─ mofa/
  │   ├─ model_raw.rds / model_final.rds / model_final.hdf5
  │   ├─ Model_raw_variance_explained.pdf
  │   ├─ Model_final_variance_explained.pdf
  │   ├─ factors.tsv
  │   └─ weights_list.rds
  ├─ weights/<GENE>/                # top features per view (for LASSO)
  └─ lasso/<GENE>/                  # coefficients, CV metrics, (plots will be added later)
```
---

## To keep in mind until all modules are finished and containerised: some hard-coded paths & parameters to update

- **R libraries:** `scripts/00_globals.R` may set `.libPaths("/scratch/.../R/x86_64.../4.x")`. Edit to your path or remove and rely on the conda/container env.
- **MOFA Python env:** `scripts/legacy/PrepareMOFA.R` uses `reticulate::use_condaenv("/scratch/.../CondaMofa/conda_envs/mofa", required=TRUE)`. Point this to your env or switch to `use_basilisk = TRUE`.
- **MOFA options:** In the same file adjust  
  `model_opts$num_factors`,  
  `model_opts$likelihoods = c(Expression="gaussian", Mutation="poisson", CNV="gaussian")`,  
  `train_opts$drop_factor_threshold`,  
  `data_opts$scale_views`.
- **Top-N feature selection:** Controlled by CLI flags `--mofa_topn`, `--topn_expr`, `--topn_mut`, `--topn_cnv` (Nextflow) or `--topn` in legacy R scripts.
- **Mutation events:** `mut_mat_mis_trunc()` in `scripts/00_globals.R` defines how MIS/TRUNC counts are built—tweak here if needed.
- **CNV thresholds:** Thresholding/continuous handling lives in `scripts/00_globals.R` (CNV helpers) — edit cutoffs there.
- **Reproducibility knobs:** `--seed` and `--nfolds` (LASSO CV); set explicitly for stable outputs.
- **Predefined genes:** Use `--predef "KRAS,TP53,..."` to force inclusion (view-normalized names are auto-expanded, e.g., `_MIS/_TRUNC` for Mutation).

## Why this matters (short)

- Reveals **tissue-specific selective pressures** behind why a gene is mutated in one tissue but not another.
- Produces **interpretable features** across Expression / CNV / Mutation (not a black box).
- Uses **no-leak** modeling around the target gene to keep results valid.
- **Reproducible & HPC-ready** (Nextflow/SLURM), easy to extend (e.g., add other modules like xCell/PROGENy; in process and coming soon).

## Requirements (minimal)

- **Nextflow ≥ 21.04.3**
- **R ≥ 4.1**: `optparse`, `fs`, `dplyr`, `readr`, `SummarizedExperiment`, `MOFA2`, `reticulate`, `ggplot2`, `maftools`
- **Python** (for MOFA via reticulate): `mofapy2==0.7.0`  
  *(Alternatively, MOFA2 can use `use_basilisk = TRUE`.)*

> Scripts resolve paths relative to the repo (`workflow.projectDir`); no cluster-specific hardcoding is required.

## Notes

- Labels can be auto-generated per target gene (`*_mut` vs `*_wt`).
- LASSO uses a **no-leak** mutation design (drops target-gene mutation features).
- `--predef` ensures specific genes are included in the MOFA/weights stage.


## Citation & License

- Free to use, copy, change, and share.
- No guarantees. It may have bugs. Use at your own risk.
- I’m not liable for any problems or losses.
- Please keep my name and this license when you share it.
- Please cite **MOFA2** and your data source (e.g., TCGA).
- Full legal text: see the **LICENSE** file (MIT).


