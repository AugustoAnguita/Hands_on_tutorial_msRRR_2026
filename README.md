# Unveiling the Exposome through Outcome-Wide Analysis: The msRRR Pipeline

**An end-to-end workflow for high-dimensional exposome data using
mixed-response sparse reduced-rank regression, bootstrap assessment, and
optional leave-one-cohort-out sensitivity analyses.**

Workflow and tutorial by:

- Dr Augusto Anguita-Ruiz, Junior Research Leader at the Barcelona Institute
  for Global Health (ISGlobal)
- María Arteaga Jover, Bioinformatics Technician at ISGlobal and PhD Student
  in Information and Communication Technologies (Bioinformatics) at the
  University of Granada

[YouTube tutorial companion](https://youtu.be/YarXHQBGgLA?si=F4EX3K43nnFaqjHz&t=1)

[![DOI](https://zenodo.org/badge/1155530150.svg)](https://doi.org/10.5281/zenodo.19615003)

## Overview

Exposome studies often involve many correlated environmental predictors and
several interconnected health outcomes, which may follow different
distributions. This repository provides a reproducible outcome-wide workflow
based on **mixed-response sparse reduced-rank regression (msRRR)**.

The workflow covers data preparation, leakage-aware cross-validation, model
diagnostics, whole-sample refitting, bootstrap analyses, publication-ready
summaries, and optional leave-one-cohort-out (LOCO) sensitivity analyses. It
illustrates two distinct post-selection approaches. They answer different
questions and should not be combined into a single inferential filter; an
applied analysis should pre-specify its primary approach.

## Main features

1. **Mixed-response reduced-rank modelling.** Gaussian and non-Gaussian
   outcomes can be modelled jointly through a lower-dimensional latent
   structure.

2. **Sparse exposure selection.** A group-lasso penalty controlled by lambda
   selects exposures across the multivariate outcome model.

3. **Leakage-aware model selection.** The example uses a stratified 80/20
   development-test split. Within cross-validation, transformations of
   continuous predictors, covariates, and Gaussian outcomes are learned from
   each training fold and then applied unchanged to its validation fold.
   Binary 0/1 variables remain on their original scale.

4. **Cross-validation diagnostics.** The revised `msrrr_cv()` interface
   reports potential issues involving the candidate rank range, lambda
   boundaries, the cross-validation criterion, missing outcomes, weak
   relationships among Gaussian outcomes, and rare dummy predictors.
   Diagnostic summaries can also be written to CSV files.

5. **Warm-start support.** During model selection, `warm = TRUE` reuses the
   preceding lambda solution only within the same rank and cross-validation
   fold. It does not transfer fitted information between folds. Setting
   `warm = FALSE` gives each lambda an independent initialisation but can be
   substantially more computationally expensive.

6. **Whole-sample refitting.** After rank and lambda selection, the penalised
   model is refitted on the whole analytical sample using the selected tuning
   parameters and the stored initialisation.

7. **Two separate post-selection approaches.**

   - **Approach A: penalised bootstrap stability.** Rank and lambda remain
     fixed. Each bootstrap fit uses an independently recalculated
     initialisation. `sel_prob` is the proportion of resamples in which an
     exposure-outcome coefficient is non-zero, and `sign_consistency` is the
     proportion of its selected coefficients whose sign agrees with the
     original whole-sample coefficient. These quantities are stability
     measures, not p-values. The tutorial illustrates thresholds of
     `sel_prob >= 0.90` and `sign_consistency >= 0.80`.
   - **Approach B: conditional inference after unpenalised refitting.** The
     exposure subset selected by the penalised whole-sample model is fixed and
     refitted with `lambda = 0`. Bootstrap estimates are used to calculate
     approximate standard errors, two-sided p-values, and 95% confidence
     intervals. These summaries are conditional on the preceding data-driven
     selection and should be interpreted as post-selection rather than
     selection-independent inference. The tutorial uses `p_value < 0.05` to
     highlight associations in the corresponding heatmaps.

8. **Optional LOCO sensitivity analyses.** Penalised and unpenalised LOCO
   analyses can be enabled independently. Each analysis removes one cohort at
   a time and recalculates preprocessing using only the retained participants.
   Separate switches control whether the corresponding LOCO bootstrap is run.

9. **Reproducible outputs.** The workflow saves model objects, diagnostic
   tables, coefficient summaries, heatmaps, and other figures. Heatmaps are
   exported to PDF and SVG, while detailed numerical results remain available
   in CSV files.

## Tutorial settings and definitive analyses

The current tutorial uses:

- 25 lambda values
- Candidate ranks 1 to 4
- 100 bootstrap resamples
- Penalised stability thresholds of 0.90 for selection probability and 0.80
  for sign consistency in Approach A
- An approximate two-sided p-value threshold of 0.05 for highlighting results
  from the conditional unpenalised refit in Approach B

These settings keep the complete demonstration computationally manageable. At
least 500 bootstrap resamples are recommended for a definitive analysis.
Candidate ranks, the lambda grid, bootstrap size, thresholds, and LOCO switches
should be reviewed for each application rather than treated as universal
defaults. The cross-validation diagnostic plots and messages should be checked
before interpreting the selected model.

For a definitive analysis, approximately 50 lambda values provide a reasonable
general starting point, together with at least 500 bootstrap resamples.
Consider increasing the lambda grid to 100 values when the CV curve changes
sharply near its minimum or when finer discrimination between neighbouring
lambda values is needed. This is not always necessary: routinely using 100
lambdas may add substantial computation without a meaningful gain when the
curve around the optimum is already smooth or flat.

## Repository structure

```text
├── data/
│   ├── codebook.csv
│   ├── covariates.csv
│   ├── exposome.csv
│   └── phenotype_NA.csv
├── functions/
│   ├── msrrr_v4.R
│   ├── bootstrapping_functions_v_31_07_26.r
│   ├── msrrr_plotting_functions.R
│   └── msrrr_tutorial_workflow_functions.R
├── results/
├── main_script_msrrr_exposome_v_03_08_2026.R
├── rmarkdown_main_script_msrrr_exposome_v_02_08_2026.html
└── quick_guide_msrrr_v4.html
```

The rendered HTML files provide the tutorial and practical guide without
requiring users to knit the R Markdown sources. The main R script contains the
executable full workflow. Output paths are based on the single `results_dir`
setting near the beginning of the script, so users can redirect all generated
results by changing it once.

## Running the workflow

1. Clone or download the repository and open R in the repository root.
2. Review the input files and variable definitions in `data/`.
3. Open `main_script_msrrr_exposome_v_03_08_2026.R` and review the settings
   near the beginning, particularly `results_dir`, `n_lambda`,
   `rank_sequence`, `n_boot_global`, the stability thresholds, and the LOCO
   switches.
4. Run the script from the repository root so that the relative paths to
   `data/` and `functions/` resolve correctly.
5. Inspect the cross-validation diagnostics before interpreting or reporting
   the whole-sample and bootstrap results.

The complete workflow, especially bootstrap and LOCO analyses, is
computationally intensive. Disable optional analyses that are not required and
use reduced settings only for workflow checks, not as a substitute for a
definitive run.

## Interpretation

The pipeline identifies multivariate exposure-outcome patterns and evaluates
their robustness under resampling and cohort exclusion. Stable selection or
sign consistency does not, by itself, establish causality. Results should be
interpreted alongside study design, data quality, subject-matter knowledge,
sensitivity analyses, and the limitations of post-selection inference.
