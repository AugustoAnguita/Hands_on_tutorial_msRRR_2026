# 🌟 Unveiling the Exposome through Outcome-Wide Analysis: The msRRR Pipeline

**Robust Inference for High-Dimensional Exposome Data using Mixed-Response Sparse Reduced-Rank Regression and Stability Selection.**

By Dr. Augusto Anguita-Ruiz
---

## 📖 Overview

Traditional environmental epidemiology often relies on univariate regression models. However, modern exposome research faces the **"Curse of Dimensionality" ($P \gg N$)**, extreme collinearity among environmental variables, and the need to evaluate multiple, interconnected health outcomes simultaneously.

This repository provides a complete, end-to-end analytical pipeline to tackle these challenges. It implements an **Outcome-Wide** approach using the **msRRR** (Mixed-Response Sparse Reduced-Rank Regression) statistical engine, coupled with a rigorous **Stability Selection** framework to filter out noise and discover true, pleiotropic environmental drivers.

### 🔑 Key Methodological Innovations
1. **The msRRR Engine:** Simultaneously models mixed-type responses (e.g., Gaussian for BMI z-score, Binomial for Asthma status) using Latent Factors (Rank $r$) to capture shared biological pathways.
2. **Strict Sparsity:** Applies a Group-Lasso penalty ($\lambda$) to force the coefficients of irrelevant exposures to absolute zero, preventing false-positive explosions.
3. **"Leakage-Proof" Harmonization:** A strict training-test split (80/20) where Z-score standardization and [0-1] scaling parameters are learned *exclusively* from the training set, preserving mathematical fairness.
4. **Advanced Balanced Bootstrapping:** Executes 500 resampling iterations (without global replacement) maintaining the optimal penalty ($\lambda_{opt}$) to account for variable selection uncertainty.
5. **The "Three-Lens" Inference System:** Replaces fragile standard p-values with robust stability metrics:
   * **Structural Stability:** Selection probability (`sel_prob` > 0.80).
   * **Directional Consistency:** Sign agreement across bootstraps (`pval_ratio` < 0.05).
   * **Empirical Significance:** Non-parametric robust p-value (`pval_corrected` < 0.05).

---

## 📂 Repository Structure

```text
├── data/
│   ├── codebook.csv            # Variable metadata and family groupings
│   ├── covariates.csv          # Adjustment variables (Z matrix)
│   ├── exposome.csv            # Environmental exposures (X matrix)
│   └── phenotype_NA.csv        # Health outcomes (Y matrix)
├── functions/
│   ├── msrrr_v3.R              # Core msRRR algorithm implementation
│   └── bootstrapping_functions.R # Custom parallel bootstrap logic
├── results/                    # Output directory for plots, heatmaps, and CSVs
├── main_script_msrrr.R         # Raw R script containing the full pipeline
└── Tutorial_msRRR_Exposome.Rmd # Interactive RMarkdown Tutorial with step-by-step theory
