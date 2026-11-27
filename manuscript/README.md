# Manuscript analyses

These scripts organize the available study analyses. The frozen model is evaluated from saved predictions without refitting. Some panels require additional source inputs; the table below identifies them.

## Setup

Run from the repository root:

```bash
Rscript scripts/test_scoring.R
cd manuscript
unzip data.zip
export OMWI_ROOT="$PWD"
```

The data archive contains the modeling matrices at six taxonomic ranks, metadata restricted to fields required for these analyses, frozen training and independent external predictions, repeated 10-fold OOF predictions, and rank-comparison summary. Raw reads are available via the study accession information in the supplementary data. Generated files are not included.

R packages used by the analysis scripts include `data.table`, `glmnet`, `pROC`, `dplyr`, `tidyr`, `tibble`, `stringr`, `ggplot2`, `ggpubr`, `patchwork`, `scales`, `vegan`, `sf`, `rnaturalearth`, and `rnaturalearthdata`. Install additional packages required by each script as indicated in its imports. The original analysis used R 4.4.3. Run each script from `manuscript/`.

## Analysis scripts

| Script | Purpose / inputs |
| --- | --- |
| `prepare_cohort_split.R`, `prepare_taxonomic_profiles.R` | Recreate study-level split and rank matrices from the original `data/all_oral/taxonomy.tsv` and `data/all_oral/all_oral_metadata.csv`; those original all-site source files are not included in this snapshot. The fixed modeling split is included. |
| `taxonomic_rank_comparison.R` | Six-rank nested CV, uncompressed species model included in this comparison. |
| `model_development.R` | Stability ranking and compact elastic-net refit with alpha=0.1, prevalence=0.10 for ranking, 40 resamples and 0.02 AUC loss budget. |
| `loocv.R` | Nested leave-one-sample-out refitting with alpha=0.1, ranking prevalence 0.10, and re-selection of K and gamma inside each training partition; corrected predictions are not supplied. |
| `Fig1_geographic_distribution.R`, `Fig1_phenotype_distribution.R` | Selected cohort geography and phenotype counts. |
| `Fig1_oral_niche_analysis.R` | All-site ordination and PERMANOVA; requires actual all-site metadata and `data/all_oral/level/taxonomy_species.tsv`. |
| `Fig2_model_performance.R`, `Fig2_score_intervals.R` | Available panels a, b, c, f and g; b plots available curves, adding corrected LOOCV only after it is supplied. ROC uses pooled repeated OOF predictions; cutoff uses per-sample mean OOF probabilities. Current study/phenotype exclusion inputs for panels d/e are not included. |
| `Fig3_phenotype_associations.R` | OMWI and diversity comparisons across studies, phenotypes and systems. |
| `Fig4_disease_severity.R` | Cross-sectional gastric lesion, colorectal neoplasia and rheumatoid arthritis activity analyses. |
| `Supplementary_LOSO_coefficients.R` | Refit coefficients from LOSO fold selections; regenerate model folds and final fit locally first. |
| `Supplementary_diversity_analysis.R` | Diversity and PCoA. |
| `Supplementary_subgroup_analysis.R` | Age, sex, geographic and oral-site analyses. |
| `Supplementary_covariate_adjustment.R` | Train-only age/sex adjustment and independent validation. |
| `Supplementary_model_comparison.R` | ROC comparison of conventional ML models; requires matching final prediction files in `differmodel/final/`. |
| `Supplementary_dental_validation.R` | Frozen dental and paired saliva/dental evaluation; requires `data/dental_validation/B/frozen_predictions_Dental_only.csv` and `data/dental_validation/C/Saliva_Dental_subject_pairs.csv`. |

The current snapshot does not include study/phenotype exclusion inputs for Fig. 2d/e, saliva-only model results, or MMUPHin sensitivity analysis inputs/scripts. No previous-version results are substituted.

## Model commands

```bash
Rscript model_development.R preflight
# Example outer CV task (100 tasks across repeats 1–10 and folds 1–10):
Rscript model_development.R cv 1 1
# Example LOSO task (one per training study, indices 1–17):
Rscript model_development.R loso 1
# After all CV and LOSO tasks finish:
Rscript model_development.R aggregate
Rscript model_development.R final-fit
Rscript model_development.R predict-external
# Rank analysis (repeat for each rank):
OMWI_RANK=species Rscript taxonomic_rank_comparison.R preflight
# Separate corrected LOOCV, one task for each sample index 1–1077:
Rscript loocv.R preflight
Rscript loocv.R loocv 1
# After every task finishes:
Rscript loocv.R aggregate-loocv
```

Refitting is computationally intensive. The frozen scoring model does not depend on these refits. The included plotting scripts preserve original analysis choices where applicable; the complete R workflow requires execution in an R environment before release.

From the repository root after extracting data, `Rscript scripts/test_real_sample_parity.R` checks scoring against the frozen scores.
