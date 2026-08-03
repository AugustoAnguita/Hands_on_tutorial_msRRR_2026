# ==============================================================================
# PROJECT: msRRR Analysis Pipeline
# DESCRIPTION: Unified script for Model Fitting, Bootstrapping, and Sparsity.
# ORIGINAL AUTHOR: Augusto Anguita-Ruiz
# MAINTENANCE AND LATEST UPDATES: María Arteaga
# ==============================================================================
#
# CURRENT REVISION:
# IMPORTANT INPUT REQUIREMENT FOR msrrr_cv():
#    Pass analysis-ready numeric X, Z and Y matrices after variable coding and
#    dummy creation, but before standardising continuous variables or scaling
#    Gaussian outcomes. msrrr_cv() performs these transformations internally:
#    in each CV iteration it learns them from the training folds only and applies
#    them unchanged to the held-out fold. After model selection, it learns and
#    stores a final transformation from the complete development sample for
#    refitting and prediction. Binary 0/1 columns and non-Gaussian outcomes
#    remain on their original scales. This input requirement is specific to the
#    recommended msrrr_cv() interface; do not pre-transform its input matrices.
#
# 1) CV diagnostic audit. The following checks are displayed when each
#    msrrr_cv() model is fitted, stored in the fitted object and saved in the
#    scenario-specific file supplied through diagnostics_file:
#    (Set diagnostics = FALSE in msrrr_cv() to suppress checks A-F)
#    A) Rank range: warns when a requested rank exceeds min(ncol(X), ncol(Y)), 
#       CV selects a tested boundary that could reasonably be extended, or the
#       selected model contains fewer active exposures than its rank, making
#       some latent components redundant. Infeasible ranks are not removed
#       automatically.
#    B) Lambda range: issues a warning when CV selects a lambda at either 
#       endpoint or when the largest tested lambda still selects exposures, 
#       reporting the tested endpoint and an example of the next bound to try.
#    C) CV criterion: for all-Gaussian outcomes, notes that pMSE may be
#       preferable to deviance. If any outcome is non-Gaussian, warns when
#       pMSE is requested and recommends deviance.
#    D) Missing outcomes: reports outcomes with at least 20% missing values and
#       warns when any outcome has at least 50% missing values. msRRR allows 
#       missing outcome values, but high levels of missingness reduce the 
#       information available and may weaken model estimation.
#    E) Outcome relationships: for all-Gaussian outcomes, notes when more than
#       half of the evaluable outcome pairs have |r| < 0.25.
#    F) Rare dummy predictors: for 0/1 columns, reports when the less frequent
#       value represents 1% to <5% of observed rows, warns below 1%, and warns
#       when a CV training fold contains only one observed level.
#
# 2) Demonstrates two DISTINCT post-selection approaches. They answer different
#    questions and should not be combined by requiring an association to pass
#    both. A real analysis should pre-specify one approach as its primary method.
#
#    APPROACH A -- PENALISED BOOTSTRAP STABILITY
#    Cross-validation selects rank and lambda. The whole sample is refitted with
#    those values and the stored initialisation, retaining the penalty.
#    Bootstrap refits keep rank and lambda fixed and independently recalculate
#    their initialisation with init = NULL. sel_prob is the proportion of
#    resamples in which an exposure-outcome coefficient is selected.
#    sign_consistency is the proportion of its selected bootstrap coefficients
#    whose sign agrees with the original whole-sample coefficient. These are
#    stability measures, not
#    p-values. Because penalised bootstrap distributions combine exact zeros
#    from non-selection with non-zero estimates, classical confidence intervals
#    and p-values are not reported. CSVs retain the original beta, bootstrap
#    summaries, both stability measures, their thresholds and TRUE/FALSE
#    decisions. Heatmaps use grey for associations that fail the stability
#    filter without removing them from the result tables.
#
#    APPROACH B -- CONDITIONAL INFERENCE AFTER UNPENALISED REFITTING
#    The exposures selected by the penalised whole-sample model are fixed. A new
#    model is fitted using only that subset, lambda = 0, the selected rank and
#    init = NULL. If the selected subset contains fewer exposures than the
#    selected rank, the rank is reduced to the largest mathematically feasible
#    value. Every bootstrap refits that same unpenalised model and independently
#    recalculates its initialisation. For each exposure-outcome coefficient, the
#    bootstrap estimates provide a standard error used to calculate an
#    approximate two-sided p-value and 95% confidence interval. These summaries
#    are conditional on the exposure subset having first been selected using the
#    same data; they do not fully account for that selection step and should be
#    interpreted as post-selection rather than selection-independent inference.
#    CSVs retain all coefficients, standard errors, p-values and intervals.
#    Heatmaps colour associations with p_value < 0.05, while non-significant
#    results remain available in the tables.
#
#    Practical orientation: Approach A targets a stable sparse signature.
#    Approach B targets coefficient estimation and approximate uncertainty after
#    fixing the selected exposure set.
#
# 3) Uses warm starts during model selection. The whole-sample penalised refit
#    is performed with refit_msrrr(), reusing the exact initial values stored by
#    msrrr_cv() so that the selected solution is reproduced consistently.
# 4) Produces separate pre/post-bootstrap heatmaps for the penalised and
#    unpenalised approaches. An additional penalised-versus-unpenalised figure is
#    included only to illustrate their differences in this tutorial; a real
#    analysis would normally pre-specify and report one primary approach.
#    Plot filtering never deletes coefficients from the result tables.
# 5) Comparative heatmaps use fixed cell dimensions, preserve predictor order
#    and share a colour scale. They may optionally be split into separate pages
#    by exposure family for display only; this does not change the fitted models.
# 6) Standardises continuous X and Z columns while leaving 0/1 dummy variables
#    on their original scale. The same rule is used in training/test,
#    whole-sample, bootstrap and LOCO analyses. Rare dummy variables are still
#    reported by diagnostic F because limited representation can make their
#    estimates unstable.
# 7) Adds OPTIONAL penalised and unpenalised leave-one-cohort-out (LOCO)
#    sensitivity analyses. LOCO removes one cohort at a time to assess whether
#    results depend strongly on a particular cohort. For every exclusion,
#    standardisation and continuous-outcome scaling are recalculated using only
#    retained participants. Bootstrap samples are generated independently
#    within each cohort-excluded dataset.
#
#    PENALISED LOCO (activate with run_penalized_loco <- TRUE):
#    Rank and lambda remain fixed at their whole-sample values, but coefficients
#    may be selected differently after each cohort is excluded. The ordinary
#    LOCO heatmap counts the exclusions in which each coefficient is selected.
#    When LOCO bootstrap is enabled, the post-bootstrap heatmap instead counts
#    the exclusions in which that coefficient passes both sel_prob and
#    sign_consistency. Colour indicates the predominant coefficient sign.
#
#    UNPENALISED LOCO (activate with run_unpenalized_loco <- TRUE):
#    The exposure subset selected by the corresponding penalised whole-sample
#    model is fixed across cohort exclusions. Within each exclusion, a new model
#    is fitted with that subset, lambda = 0 and init = NULL, followed by its own
#    bootstrap. The heatmap counts the exclusions with p_value < 0.05; colour
#    indicates the predominant sign among those results.
#
# 8) Exports fixed-layout heatmaps as multi-page PDFs and editable SVGs, with
#    one scenario per page or file. LOCO comparison reports place the summaries
#    before and after bootstrap side by side.
#
#

# ------------------------------------------------------------------------------
# BLOCK 1: SET UP & LIBRARIES
# ------------------------------------------------------------------------------

# Clean environment
rm(list = ls())

# Create the common output directory before any plot, checkpoint or table is
# written. recursive = TRUE is harmless when the directory already exists.
results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(results_dir)) {
  stop("The output directory could not be created: ", results_dir)
}

# ------------------------------------------------------------------------------
# CORE MODEL AND BOOTSTRAP SETTINGS
# ------------------------------------------------------------------------------

# Number of lambda values evaluated along the penalisation path.
# Around 50 values usually provide sufficient resolution for model selection.
# Increase this value if the CV curve changes sharply near its minimum or if
# finer discrimination between neighbouring lambda values is required. This
# tutorial uses 25 because to substantially reduce computation time.
n_lambda <- 25

# Candidate reduced ranks evaluated by cross-validation.
# Rank 1 should normally be included. The upper value should reflect the
# expected latent structure and cannot meaningfully exceed
# min(ncol(X), ncol(Y)). Review the rank-boundary diagnostics after fitting.
rank_sequence <- 1:4

# Number of bootstrap resamples used by both post-selection approaches.
# At least 500 resamples are recommended for a definitive analysis. This
# tutorial uses 100 because it demonstrates both bootstrap approaches and
# optional LOCO analyses, which substantially increases computation time.
n_boot_global <- 100

# Close all open graphical devices from previous/interrupted runs
graphics.off()
if (!is.null(dev.list())) {dev.off()}

# Sanity check
print(dev.list())

# Define all required packages
required_packages <- c(
  "tidyverse", "rrpack", "foreach", "doParallel", "pheatmap", "glmnet", "fastDummies", 
  "gridExtra", "corrplot","caret","ggrepel","reshape2","RColorBrewer","grid",
  "doSNOW", "matrixStats", "RMTL" 
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(rrpack)
  library(foreach)
  library(doParallel)
  library(doSNOW)
  library(matrixStats)
  library(pheatmap)
  library(glmnet)
  library(fastDummies)
  library(gridExtra)
  library(corrplot)
  library(caret)
  library(ggrepel) 
  library(reshape2)
  library(RColorBrewer)
  library(grid) 
})

# SOURCE EXTERNAL FUNCTIONS
# Ensure these files are in a 'functions' subfolder
if(file.exists("functions/msrrr_v4.R")) {
  source("functions/msrrr_v4.R") 
  source("functions/bootstrapping_functions_v_31_07_26.r") 
  source("functions/msrrr_plotting_functions.R")
  source("functions/msrrr_tutorial_workflow_functions.R")
  selection_tol <- 1e-12
  message(">>> Setup Complete. Libraries and Functions loaded.")
} else {
  stop("ERROR: Functions not found. Please check your 'functions' folder.")
}


# ------------------------------------------------------------------------------

# BLOCK 2: LOAD DATA
# ------------------------------------------------------------------------------

# Helper to read European CSVs (semicolon separated, comma decimals)
# Workflow helper functions used below are defined in
# functions/msrrr_tutorial_workflow_functions.R.

message(">>> Loading data files...")
covs_raw <- read_euro_csv("data/covariates.csv")
exp_raw  <- read_euro_csv("data/exposome.csv")
# Loading the version with NAs as requested (Standard CSV read)
phen_raw <- read.csv("data/phenotype_NA.csv", row.names = 1, na.strings = "NA")


# ------------------------------------------------------------------------------

# BLOCK 3: VARIABLE SELECTION & PROCESSING
# ------------------------------------------------------------------------------
message(">>> [Block 3] Starting Variable Selection with fastDummies...")

# 3.1 Load Codebook
codebook <- read_euro_csv("data/codebook.csv")

# A. Target Families (Numeric variables)
target_families <- c(
  "Air Pollution", "Built environment", "Indoor air", "Lifestyle", 
  "Metals", "Meteorological", "Natural Spaces", "Organochlorines", 
  "Per- and polyfluoroalkyl substances (PFAS)", "Social and economic capital", 
  "Tobacco Smoke", "Traffic", "Water DBPs"
)

# B. Forced Factors (Variables to convert to Dummies)
forced_factors <- c(
  "e3_alcpreg_yn_None", "h_pamod_t3_None", "h_pavig_t3_None", 
  "FAS_cat_None", "hs_globalexp2_None"
)

# 3.2 Filter Variables from Codebook
vars_info <- codebook %>%
  dplyr::filter(
    (family %in% target_families & var_type == "numeric") |
      (variable_name %in% forced_factors)
  )

# 3.3 Split Variables by Period
vars_preg_names  <- vars_info %>% dplyr::filter(period == "Pregnancy") %>% dplyr::pull(variable_name)
vars_child_names <- vars_info %>% dplyr::filter(period == "Postnatal") %>% dplyr::pull(variable_name)

# Intersect with available data (Safety check)
vars_preg_names  <- intersect(vars_preg_names, colnames(exp_raw))
vars_child_names <- intersect(vars_child_names, colnames(exp_raw))


# 3.4 Process Exposome Matrices (X)
# Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.

message("    Processing Pregnancy Matrix...")
X_preg <- process_exposome_matrix(exp_raw, vars_preg_names, forced_factors)

message("    Processing Childhood Matrix...")
X_child <- process_exposome_matrix(exp_raw, vars_child_names, forced_factors)

message(paste("    -> X_preg dimensions:", nrow(X_preg), "x", ncol(X_preg)))
message(paste("    -> X_child dimensions:", nrow(X_child), "x", ncol(X_child)))


# 3.5 Process Covariates (Z)
covs_list <- c("h_mbmi_None", "e3_gac_None", "e3_sex_None", 
               "h_age_None", "h_cohort", "h_edumc_None")

Z_df <- covs_raw %>%
  dplyr::select(all_of(covs_list)) %>%
  dplyr::mutate(
    e3_sex_None = dplyr::case_when(
      e3_sex_None == "female" ~ 0,
      e3_sex_None == "male" ~ 1,
      TRUE ~ NA_real_
    ),
    h_cohort = as.factor(h_cohort),
    h_edumc_None = as.factor(h_edumc_None)
  )

# With K cohort categories, use K - 1 dummy columns because the model includes
# an intercept. remove_first_dummy = TRUE makes the first factor level the
# reference category (for six cohorts, five cohort dummies are therefore
# created). Every cohort-dummy coefficient is interpreted relative to that
# omitted reference cohort.
Z_df_dummies <- dummy_cols(Z_df, 
                           select_columns = c("h_cohort", "h_edumc_None"),
                           remove_first_dummy = TRUE, 
                           remove_selected_columns = TRUE)

Z <- as.matrix(Z_df_dummies)
message(paste("    -> Covariates dimensions:", nrow(Z), "vars.", ncol(Z)))

# 3.6 Process Outcomes (Y)
# UPDATED: Added "hs_bmi_c_cat" to the list
outcomes_list <- c(
  "e3_bw", "hs_asthma", "hs_zbmi_who", 
  "hs_correct_raven", "hs_Gen_Tot", "hs_bmi_c_cat"
)

Y_df <- phen_raw %>%
  dplyr::select(all_of(outcomes_list)) %>%
  # Force numeric on all columns (this handles hs_bmi_c_cat categorical -> numeric)
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

Y <- as.matrix(Y_df)
message(paste("    -> Outcomes (Y) Ready:", ncol(Y), "vars."))


# ------------------------------------------------------------------------------

# BLOCK 3.7: QUALITY CONTROL - NA COUNT PER VARIABLE
# ------------------------------------------------------------------------------
message(">>> [Block 3.7] Running NA Quality Checks (Counts)...")

# Function to report NAs per variable
# Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.

# Run Checks
report_nas(X_preg,  "Exposome Pregnancy (X_preg)")
report_nas(X_child, "Exposome Childhood (X_child)")
report_nas(Z,       "Covariates (Z)")
report_nas(Y,       "Phenotypes/Outcomes (Y)")


# ------------------------------------------------------------------------------
# BLOCK 3.8: MISSING VALUES VISUALIZATION (Phenotypes)
# ------------------------------------------------------------------------------
message("\n>>> [Block 3.8] Generating Missing Values Plots...")

# Use Y_df (dataframe) for plotting
df_ <- Y_df 

# 1. Prepare Data for Percentage Plot
missing.values <- df_ %>%
  tidyr::gather(key = "key", value = "val") %>%
  dplyr::mutate(isna = is.na(val)) %>%
  dplyr::group_by(key) %>%
  dplyr::mutate(total = dplyr::n()) %>%
  dplyr::group_by(key, total, isna) %>%
  dplyr::summarise(num.isna = dplyr::n(), .groups = 'drop') %>%
  dplyr::mutate(pct = num.isna / total * 100)

levels <- (
  missing.values %>%
    dplyr::filter(isna == TRUE) %>%
    dplyr::arrange(dplyr::desc(pct))
)$key
if(length(levels) == 0) levels <- unique(missing.values$key)

# 2. Percentage Bar Plot
percentage.plot <- missing.values %>%
  ggplot() +
  geom_bar(aes(x = reorder(key, desc(pct)), y = pct, fill=isna), 
           stat = 'identity', alpha=0.8) +
  scale_x_discrete(limits = levels) +
  scale_fill_manual(name = "", values = c('steelblue', 'tomato3'), 
                    labels = c("Present", "Missing")) +
  coord_flip() +
  labs(title = "Percentage of missing values", x = 'Variable', y = "% of missing values")

# 3. Row Missingness Plot
row.plot <- df_ %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  tidyr::gather(-id, key = "key", value = "val") %>%
  dplyr::mutate(isna = is.na(val)) %>%
  ggplot(aes(key, id, fill = isna)) +
  geom_raster(alpha=0.8) +
  scale_fill_manual(name = "", values = c('steelblue', 'tomato3'),
                    labels = c("Present", "Missing")) +
  scale_x_discrete(limits = levels) +
  labs(x = "Variable", y = "Row Number", title = "Missing values in rows") +
  coord_flip()

# 4. Combine, save and display plots
missingness_plot <- gridExtra::arrangeGrob(
  percentage.plot, row.plot, ncol = 2
)
missingness_file <- file.path(results_dir, "MISSING_VALUES_PHENOTYPES.pdf")
grDevices::pdf(missingness_file, width = 14, height = 8)
grid::grid.newpage()
grid::grid.draw(missingness_plot)
grDevices::dev.off()
grid::grid.newpage()
grid::grid.draw(missingness_plot)
message(">>> Missing-values plots saved to: ", missingness_file)

# 5. Complete Cases Analysis Report
cc_matrix <- na.omit(df_)
message(paste0("\n>>> COMPLETE CASES ANALYSIS for PHENOTYPES:"))
message(paste0("    Original Dimensions: ", paste(dim(df_), collapse = " x ")))
message(paste0("    Complete Dimensions: ", paste(dim(cc_matrix), collapse = " x ")))
message(paste0("    Rows dropped: ", nrow(df_) - nrow(cc_matrix)))


# ------------------------------------------------------------------------------

# BLOCK 4: CORRELATION ANALYSIS (Corrplots)
# ------------------------------------------------------------------------------
message("\n>>> [Block 4] Generating Correlation Plots...")

# Plotting helper functions used below are defined in
# functions/msrrr_plotting_functions.R.

# Generate plots for all 3 datasets
generate_corrplot(
  X_preg, "Correlation Matrix: Pregnancy Exposome",
  file.path(results_dir, "CORRPLOT_predictors_Pregnancy.pdf"),
  label_cex = 0.55, number_cex = 0.35
)
generate_corrplot(
  X_child, "Correlation Matrix: Childhood Exposome",
  file.path(results_dir, "CORRPLOT_predictors_Childhood.pdf"),
  label_cex = 0.50, number_cex = 0.35
)
generate_corrplot(
  Y, "Correlation Matrix: Phenotypes",
  file.path(results_dir, "CORRPLOT_Phenotypes.pdf"),
  label_cex = 0.95, number_cex = 0.80
)

# Save only the necessary processed matrices and the codebook info
save(X_preg, X_child, Z, Y, Y_df, vars_info, 
     file = file.path(results_dir, "checkpoint_preproc.RData"))



# ------------------------------------------------------------------------------

# BLOCK 5: DATA SPLITTING & STRATIFICATION (Training/Test & Folds)
# ------------------------------------------------------------------------------
message(">>> [Block 5] Starting Data Splitting and Stratification...")

# 5.1 Configuration
# ------------------------------------------------------------------------------
set.seed(12345)
train_prop <- 0.8
num_folds  <- 5

# 5.2 Create Stratification Variable (Combined Asthma & BMI Category)
# ------------------------------------------------------------------------------
# We combine both variables to ensure balance in both. 
# We handle NAs as a separate level so we don't lose any individuals.

# Temporary factors for stratification
s_asthma  <- ifelse(is.na(Y[, "hs_asthma"]), "NA", as.character(Y[, "hs_asthma"]))
s_bmi_cat <- ifelse(is.na(Y[, "hs_bmi_c_cat"]), "NA", as.character(Y[, "hs_bmi_c_cat"]))

# Create the interaction strata
combined_strata <- as.factor(paste(s_asthma, s_bmi_cat, sep = "_"))

message("    Strata distribution (Asthma_BMIcat):")
print(table(combined_strata))

# 5.3 Split Population (Training 80% / Test 20%)
# ------------------------------------------------------------------------------
# Using caret::createDataPartition to ensure the 80/20 split is balanced
train_index <- caret::createDataPartition(combined_strata, p = train_prop, list = FALSE)
train_index <- as.vector(train_index)

# Split Exposome
X_preg_train  <- X_preg[train_index, ]
X_preg_test   <- X_preg[-train_index, ]
X_child_train <- X_child[train_index, ]
X_child_test  <- X_child[-train_index, ]

# Split Covariates and Outcomes
Z_train  <- Z[train_index, ]
Z_test   <- Z[-train_index, ]
Y_train  <- Y[train_index, ]
Y_test   <- Y[-train_index, ]

message(paste("    Split complete. Training N:", nrow(Y_train), "| Test N:", nrow(Y_test)))

# 5.4 Generate Stratified Fold IDs for Cross-Validation
# ------------------------------------------------------------------------------
# We repeat the stratification logic inside the Training set for the 5 folds.
train_strata <- combined_strata[train_index]

# Create folds
folds_list <- caret::createFolds(train_strata, k = num_folds, list = TRUE)

# Initialize foldid vector
foldid <- numeric(nrow(Y_train))
for(i in seq_along(folds_list)) {
  foldid[folds_list[[i]]] <- i
}
# ------------------------------------------------------------------------------
# BLOCK 5.5: COMPACT BALANCE VERIFICATION (Fixed Alignment)
# ------------------------------------------------------------------------------
message("\n>>> [Block 5.5] Verifying Balance across Folds (Fixed Table)...")

# Variables to check
check_vars <- c("hs_asthma", "hs_bmi_c_cat")

for (var_ in check_vars) {
  message(paste0("\n--- Proportion Table for: ", var_, " ---"))
  
  # 1. Identify all possible levels (including NA) across the whole dataset
  # This ensures all tables have the same rows even if a fold is missing a category
  all_levels <- names(table(Y[, var_], useNA = "always"))
  
  # 2. Function to get consistent proportions per fold
  # Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.
  
  # 3. Build the comparison matrix
  comp_matrix <- sapply(1:num_folds, function(f) {
    get_prop(Y_train[foldid == f, var_])
  })
  
  # 4. Add the Test Set column
  test_col <- get_prop(Y_test[, var_])
  final_tab <- cbind(comp_matrix, test_col)
  
  # 5. Formatting
  colnames(final_tab) <- c(paste0("Fold_", 1:num_folds), "Test_Set")
  rownames(final_tab) <- ifelse(is.na(all_levels), "Missing/NA", all_levels)
  
  print(final_tab)
}

# ------------------------------------------------------------------------------
# BLOCK 5.6: FINAL DIMENSIONS & SPLIT SUMMARY
# ------------------------------------------------------------------------------
message("\n>>> Final sample distribution:")
print(table(foldid, dnn = "Fold ID"))

message(paste0("\n    Exposome Matrix (Train): ", nrow(X_child_train), " x ", ncol(X_child_train)))
message(paste0("    Exposome Matrix (Test):  ", nrow(X_child_test), " x ", ncol(X_child_test)))

message("\n>>> [Block 5] Step completed. The proportions are now perfectly aligned.")


# ------------------------------------------------------------------------------
# BLOCK 5.7: TEST SET PREPARATION FOR EVALUATION
# ------------------------------------------------------------------------------
# In a tutorial context, we want the Test Set to have ground truth values
# to calculate prediction metrics (RMSE, R2, etc.) later.

# We check if there are any samples in Y_test that present any NA
full_na_test <- rowSums(is.na(Y_test)) >= 1

if(any(full_na_test)) {
  message(paste(">>> Removing", sum(full_na_test), "subjects from Test set with NA outcome data."))
  
  X_preg_test  <- X_preg_test[!full_na_test, ]
  X_child_test <- X_child_test[!full_na_test, ]
  Z_test       <- Z_test[!full_na_test, ]
  Y_test       <- Y_test[!full_na_test, ]
}

# Note: Individual NAs within Y_test are fine, but when calculating metrics later,
# we will need to use functions that handle NAs (e.g., cor(..., use = "complete.obs")).

message(">>> Block 5 fully optimized for Model Evaluation.")


# ------------------------------------------------------------------------------

# BLOCK 5.8: FINAL DATA INTEGRITY CHECK (Dimensions & Names)
# ------------------------------------------------------------------------------
message("\n>>> [Block 5.8] Final Data Integrity Check before Modeling...")

# 1. List of main datasets to inspect
# ------------------------------------------------------------------------------
datasets_to_check <- list(
  X_preg_train  = X_preg_train, 
  X_preg_test   = X_preg_test,
  X_child_train = X_child_train, 
  X_child_test  = X_child_test,
  Z_train       = Z_train,       
  Z_test        = Z_test,
  Y_train       = Y_train,       
  Y_test        = Y_test
)

# 2. Automated Inspection Loop
# ------------------------------------------------------------------------------
for (name in names(datasets_to_check)) {
  mat <- datasets_to_check[[name]]
  
  # Basic Info
  cat(paste0("\n--- Dataset: ", name, " ---\n"))
  cat(paste0("  [Dim]: ", nrow(mat), " rows x ", ncol(mat), " columns\n"))
  
  # Column check (First 10)
  all_cols <- colnames(mat)
  n_to_show <- min(10, length(all_cols))
  cat(paste0("  [First ", n_to_show, " vars]: ", 
             paste(all_cols[1:n_to_show], collapse = ", "), "...\n"))
  
  # NA Check (Should be 0 for X and Z, but can have some in Y)
  n_nas <- sum(is.na(mat))
  if(n_nas > 0) {
    cat(paste0("  [Alert]: Found ", n_nas, " missing values (Expected in Y, check if in X/Z)\n"))
  }
}

# 3. Validation Fold (foldid) Check
# ------------------------------------------------------------------------------
cat("\n--- CV Folds (foldid) ---\n")
cat(paste0("  Length: ", length(foldid), " (Matches Y_train: ", length(foldid) == nrow(Y_train), ")\n"))
cat("  Frequency per fold:\n")
print(table(foldid))

message("\n>>> [Integrity Check] Done. All datasets are aligned and ready.")


# ------------------------------------------------------------------------------

# BLOCK 5.9: RAW MATRICES FOR LEAKAGE-FREE INTERNAL CV
# ------------------------------------------------------------------------------
message(
  ">>> [Block 5.9] Keeping X, Z and Y untransformed for msrrr_cv()..."
)

# msrrr_cv() learns X/Z standardisation and Gaussian Y min-max scaling using
# only the four training folds in each CV iteration, then applies those
# parameters unchanged to the held-out validation fold. Its final training fit
# learns a new transformation from the complete 80% development sample.
# Binary X/Z columns remain 0/1.
#
# hs_bmi_c_cat is stored as four ordered categories (1-4), not as an existing
# dummy. Convert its scientifically defined lower/higher grouping explicitly
# before CV; this is recoding, not data-derived scaling.
# 0 (Under+Normal) vs 1 (Over+Obese)
Y_train[, "hs_bmi_c_cat"] <- as.numeric(Y_train[, "hs_bmi_c_cat"] >= 3)
Y_test[, "hs_bmi_c_cat"]  <- as.numeric(Y_test[, "hs_bmi_c_cat"] >= 3)

message("\n--- BMI recoding check (Train/Test) ---")
print(table(Y_train[, "hs_bmi_c_cat"], dnn = "Train_BMI_Binary"))
print(table(Y_test[, "hs_bmi_c_cat"], dnn = "BMI_Binary"))



# ------------------------------------------------------------------------------

# BLOCK 6: msRRR MODEL FITTING (Sequential Scenarios)
# ------------------------------------------------------------------------------
message(">>> [Block 6] Fitting Models for Pregnancy, Childhood, and Combined...")

# 6.0 Shared Model Configuration
# ------------------------------------------------------------------------------
# Outcome Family Mapping:
# 1: Birthweight (Gaussian), 2: Asthma (Binomial), 3: zBMI (Gaussian), 
# 4: Raven (Gaussian), 5: Gen_Health (Gaussian), 6: Binary BMI Cat (Binomial)
outcome_families <- list(gaussian(), binomial())
family_mapping   <- c(1, 2, 1, 1, 1, 2) 

# Lambda and rank options:
# - nrankseq contains the candidate ranks compared by CV.
# - Supply lamseq explicitly to test a user-defined lambda sequence.
# - With lamseq = NULL, msrrr_cv() automatically creates nlam lambda values and
#   reuses that same sequence in every fold. nlam is ignored when lamseq is
#   supplied. With warm = TRUE, lambdas are evaluated from largest to smallest.
#
# CV and diagnostics:
# - foldid fixes the fold assignment; nfold must match its number of levels.
# - cv.criteria selects "deviance" or "pMSE".
# - diagnostics = TRUE runs and stores checks in fit$diagnostics.
# - print_diagnostics_at_end = TRUE prints their compact end-of-fit summary.
# - diagnostics_file is optional; when a path is supplied, the same audit is
#   also written to CSV. Use NULL if no separate file is required.
# - progress = TRUE displays one overall bar based on completed model fits.
#   control$trace = TRUE is much more verbose and is intended for debugging.

# TECHNICAL NOTE ON MISSING Y VALUES DURING AUTOMATIC LAMBDA CALCULATION:
# crossprod() and svd() cannot operate with NA values. When lamseq = NULL,
# msrrr_cv() creates a temporary copy of Y and replaces missing entries with
# their observed outcome-column mean solely to calculate the initial lambda
# range. This temporary replacement is not used to fit the models and does not
# modify the original outcome matrix or its missing values.

# ------------------------------------------------------------------------------
# 6.1 SCENARIO A: PREGNANCY EXPOSOME ONLY
# ------------------------------------------------------------------------------
message("\n>>> Scenario A: Fitting Pregnancy Model...")

cv.crit_preg <- "deviance"
# cv.crit_preg <- "pMSE"

fit_preg <- msrrr_cv(
  Y = Y_train, X = X_preg_train, Z = Z_train,
  family = outcome_families, familygroup = family_mapping,
  nrankseq = rank_sequence, lamseq = NULL, nlam = n_lambda,
  foldid = foldid, nfold = num_folds,
  method = "CV", cv.criteria = cv.crit_preg, warm = TRUE,
  diagnostics = TRUE, print_diagnostics_at_end = TRUE, progress = TRUE,
  diagnostics_file = file.path(results_dir, "MSRRR_CV_DIAGNOSTICS_PREGNANCY.csv")
)


# ------------------------------------------------------------------------------
# 6.2 SCENARIO B: CHILDHOOD EXPOSOME ONLY
# ------------------------------------------------------------------------------
message("\n>>> Scenario B: Fitting Childhood Model...")

cv.crit_child <- "deviance"
# cv.crit_child <- "pMSE"

fit_child <- msrrr_cv(
  Y = Y_train, X = X_child_train, Z = Z_train,
  family = outcome_families, familygroup = family_mapping,
  nrankseq = rank_sequence, lamseq = NULL, nlam = n_lambda,
  foldid = foldid, nfold = num_folds,
  method = "CV", cv.criteria = cv.crit_child, warm = TRUE,
  diagnostics = TRUE, print_diagnostics_at_end = TRUE, progress = TRUE,
  diagnostics_file = file.path(results_dir, "MSRRR_CV_DIAGNOSTICS_CHILDHOOD.csv")
)

# ------------------------------------------------------------------------------
# 6.3 SCENARIO C: COMBINED EXPOSOME (Pregnancy + Childhood)
# ------------------------------------------------------------------------------
message("\n>>> Scenario C: Fitting Combined Model...")

# Concatenating sources to evaluate joint effects
X_combined_train <- cbind(X_preg_train, X_child_train)

cv.crit_comb <- "deviance"
# cv.crit_comb <- "pMSE"

fit_comb <- msrrr_cv(
  Y = Y_train, X = X_combined_train, Z = Z_train,
  family = outcome_families, familygroup = family_mapping,
  nrankseq = rank_sequence, lamseq = NULL, nlam = n_lambda,
  foldid = foldid, nfold = num_folds,
  method = "CV", cv.criteria = cv.crit_comb, warm = TRUE,
  diagnostics = TRUE, print_diagnostics_at_end = TRUE, progress = TRUE,
  diagnostics_file = file.path(results_dir, "MSRRR_CV_DIAGNOSTICS_COMBINED.csv")
)


# ------------------------------------------------------------------------------
# BLOCK 7.1: INITIALIZE RESULTS SUMMARY TABLE
# ------------------------------------------------------------------------------
message(">>> [Block 7.1] Summarizing Model Results...")

# 1. Prepare Joint Test Matrix for Scenario C
X_combined_test <- cbind(X_preg_test, X_child_test)

# 2. Initialize the results dataframe
results_df <- data.frame(
  Scenario      = c("Pregnancy", "Childhood", "Combined"),
  Seed          = rep(12345, 3),
  Opt_Rank      = rep(NA, 3),
  Opt_Lambda    = rep(NA, 3),
  N_Exposures   = rep(NA, 3),
  CV_Crit_Name   = c(cv.crit_preg, cv.crit_child, cv.crit_comb),
  CV_Crit_value   = rep(NA, 3),
  pMSE_training = rep(NA, 3),
  pMSE_test     = rep(NA, 3)
)

# 3. Helper function to fill the table
# Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.

fill_results(1, fit_preg)
fill_results(2, fit_child)
fill_results(3, fit_comb)



# ------------------------------------------------------------------------------
# BLOCK 7.2: PREDICTION PERFORMANCE (pMSE)
# ------------------------------------------------------------------------------
message(">>> [Block 7.2] Calculating pMSE for Training and Test sets...")

# Scenario A: Pregnancy
fit_train_preg <- predict(
  fit_preg,
  Y.new = Y_train, X.new = X_preg_train, Z.new = Z_train,
  family = outcome_families, familygroup = family_mapping,
  cv.criteria = "pMSE"
)
fit_test_preg <- predict(
  fit_preg,
  Y.new = Y_test, X.new = X_preg_test, Z.new = Z_test,
  family = outcome_families, familygroup = family_mapping,
  cv.criteria = "pMSE"
)
results_df$pMSE_training[1] <- fit_train_preg$pred.perf
results_df$pMSE_test[1]     <- fit_test_preg$pred.perf

# Scenario B: Childhood
fit_train_child <- predict(
  fit_child,
  Y.new = Y_train, X.new = X_child_train, Z.new = Z_train,
  family = outcome_families, familygroup = family_mapping,
  cv.criteria = "pMSE"
)
fit_test_child <- predict(
  fit_child,
  Y.new = Y_test, X.new = X_child_test, Z.new = Z_test,
  family = outcome_families, familygroup = family_mapping,
  cv.criteria = "pMSE"
)
results_df$pMSE_training[2] <- fit_train_child$pred.perf
results_df$pMSE_test[2]     <- fit_test_child$pred.perf

# Scenario C: Combined
fit_train_comb <- predict(
  fit_comb,
  Y.new = Y_train, X.new = X_combined_train, Z.new = Z_train,
  family = outcome_families, familygroup = family_mapping,
  cv.criteria = "pMSE"
)
fit_test_comb <- predict(
  fit_comb,
  Y.new = Y_test, X.new = X_combined_test, Z.new = Z_test,
  family = outcome_families, familygroup = family_mapping,
  cv.criteria = "pMSE"
)
results_df$pMSE_training[3] <- fit_train_comb$pred.perf
results_df$pMSE_test[3]     <- fit_test_comb$pred.perf

results_df$Perf_Dif <- abs(results_df$pMSE_training - results_df$pMSE_test)
print(results_df)

# Save the primary model objects and the final results table
save(fit_preg, fit_child, fit_comb, results_df, 
     file = file.path(results_dir, "msrrr_objects_model_metrics_results_training_CV.RData"))


write.csv2(results_df, file.path(results_dir, "model_metrics_results_training_CV.csv"), row.names = FALSE)


# ------------------------------------------------------------------------------
# BLOCK 7.3: DIAGNOSTIC PLOTS (Generalized Function)
# ------------------------------------------------------------------------------

message(">>> [Block 7.3] Defining and generating Diagnostic Plots...")

# Function to generate Rank and Lambda diagnostics for any scenario
# Plotting helper functions used below are defined in functions/msrrr_plotting_functions.R.

# Example: Generate plots for each scenario
diag_preg  <- generate_diagnostics(fit_preg, "Pregnancy", cv.crit = cv.crit_preg)
diag_child <- generate_diagnostics(fit_child, "Childhood", cv.crit = cv.crit_child)
diag_comb  <- generate_diagnostics(fit_comb, "Combined", cv.crit = cv.crit_comb)

pdf(file.path(results_dir, "Model_Diagnostics_training_CV.pdf"), width = 8.5, height = 11)
plot(diag_preg)
plot(diag_child)
plot(diag_comb)
dev.off()


# ------------------------------------------------------------------------------
# BLOCK 7.4: COEFFICIENT LOLLIPOP PLOTS (Top Exposures)
# ------------------------------------------------------------------------------
message(">>> [Block 7.4] Plotting Coefficients for Key Outcomes...")

# Plotting helper functions used below are defined in functions/msrrr_plotting_functions.R.

pdf(file.path(results_dir, "Coefficient_Lollipop_Plots_training_CV.pdf"), width = 12, height = 8)
coefficient_plot_scenarios <- list(
  Pregnancy = list(fit = fit_preg, X = X_preg_train),
  Childhood = list(fit = fit_child, X = X_child_train),
  Combined = list(fit = fit_comb, X = X_combined_train)
)
for (scenario_label in names(coefficient_plot_scenarios)) {
  scenario <- coefficient_plot_scenarios[[scenario_label]]
  for (outcome_idx in seq_len(ncol(Y_train))) {
    coefficient_plot <- plot_coeffs(
      scenario$fit, outcome_idx, scenario$X, scenario_label
    )
    if (!is.null(coefficient_plot)) print(coefficient_plot)
  }
}
dev.off()





# ------------------------------------------------------------------------------
# BLOCK 8: PROFESSIONAL COEFFICIENT HEATMAPS (With Interpretation Guide)
# ------------------------------------------------------------------------------

message(">>> [Block 8] Generating custom heatmaps and documentation...")

# 8.0 TECHNICAL NOTE: INTERPRETING msRRR COEFFICIENTS (Beta)
# ------------------------------------------------------------------------------
# 1. CONTINUOUS OUTCOMES (e.g., Birthweight, zBMI):
#    Since outcomes were scaled [0-1] and predictors were standardized (Z-score),
#    Beta represents the change in the outcome range (0 to 1) for every 
#    1 standard deviation (SD) increase in the exposure, adjusted for covariates.
#    - Beta = 0.1 means: "1 SD increase in exposure is associated with a 10% 
#      increase in the total range of the outcome."
#
# 2. BINARY OUTCOMES (e.g., Asthma, BMI Category):
#    These use a Logit link. Beta represents the change in the LOG-ODDS of 
#    the outcome for every 1 SD increase in exposure. 
#    Exp(Beta) would provide the Odds Ratio (OR).
#    Beta = 0.4 means: OR = exp(0.4) ≈ 1.50 (a 50% increase in risk).

# 3. ARE EFFECT SIZES COMPARABLE?
#    Technically, NO. Beta_gaussian and Beta_binomial are on different scales 
#    (Linear vs. Logit). However, within this msRRR framework, they are 
#    comparable in terms of "Selection Strength": the penalty treats them 
#    similarly to find the most influential predictors across all domains. 
#    A strong color in both indicates the exposure is a powerful driver for
#    both health domains, regardless of the unit.
# ------------------------------------------------------------------------------


# 8.1 CUSTOM "SHARP" COLOR PALETTE (GREENS AND REDS)
# ------------------------------------------------------------------------------
# Define two gradients that do NOT pass through white or yellow.
# Negative (green): dark green to light green
neg_palette <- colorRampPalette(c("#1A9850", "#A6DBA0"))(50) 
# Neutral: grey (only for exact zero)
zero_color  <- "grey" 
# Positive (red): light red to dark red
pos_palette <- colorRampPalette(c("#F4A582", "#D73027"))(50) 

custom_colors_pro <- c(neg_palette, zero_color, pos_palette)

# 8.2 CUSTOM BREAKS (KEEPING NON-ZERO VALUES COLOURED)
# ------------------------------------------------------------------------------
# Plotting helper functions used below are defined in functions/msrrr_plotting_functions.R.

# 8.3 HEATMAP GENERATION FUNCTION
# ------------------------------------------------------------------------------
# Plotting helper functions used below are defined in functions/msrrr_plotting_functions.R.

# 8.4 Export
# ------------------------------------------------------------------------------
graphics.off()

cv_heatmap_grobs <- list(
  plot_final_heatmap(fit_preg, X_preg_train, "Pregnancy", draw = FALSE)$gtable,
  plot_final_heatmap(fit_child, X_child_train, "Childhood", draw = FALSE)$gtable,
  plot_final_heatmap(fit_comb, X_combined_train, "Combined", draw = FALSE)$gtable
)
save_grobs_pdf(
  cv_heatmap_grobs,
  file.path(results_dir, "Exposome_Signatures_training_CV.pdf"),
  width = 18, height = 30
)
graphics.off()
message(">>> Block 8 Complete. Heatmaps updated with non-linear scale.")

load(file.path(results_dir, "msrrr_objects_model_metrics_results_training_CV.RData"))
results_df <- read.csv2(file.path(results_dir, "model_metrics_results_training_CV.csv"))



# ------------------------------------------------------------------------------
# BLOCK 9: SHARED PENALISED WHOLE-SAMPLE MODELS
# ------------------------------------------------------------------------------
# These are the main whole-sample models reported in Approach A. Their selected
# exposure subsets also define the candidate exposures that are subsequently
# refitted with lambda = 0 in Approach B. Block 9 is therefore required for
# BOTH approaches.
message(">>> [Block 9] Fitting shared penalised whole-sample models...")
# ------------------------------------------------------------------------------
set.seed(12345)

# 9.1 Prepare Whole Sample Matrices
# ------------------------------------------------------------------------------
outcome_families <- list(gaussian(), binomial())
family_mapping   <- c(1, 2, 1, 1, 1, 2) 

# Supply raw matrices to refit_msrrr(). Because the CV objects contain
# preprocessing metadata, the refit learns new parameters from all subjects.
Y_whole_raw <- Y
Y_whole_raw[, "hs_bmi_c_cat"] <- as.numeric(
  Y_whole_raw[, "hs_bmi_c_cat"] >= 3
)
X_preg_whole_raw <- X_preg
X_child_whole_raw <- X_child
X_comb_whole_raw <- cbind(X_preg, X_child)
Z_whole_raw <- Z

# 9.2 Refit Final Models
# ------------------------------------------------------------------------------
# Reuse the exact initial values stored by msrrr_cv(). This keeps the current
# warm-start approach reproducible while applying the selected rank and lambda
# to the whole sample.

# Scenario A: Pregnancy
fit_preg_whole <- refit_msrrr(
  object = fit_preg, Y = Y_whole_raw, X = X_preg_whole_raw, Z = Z_whole_raw
)

# Scenario B: Childhood
fit_child_whole <- refit_msrrr(
  object = fit_child, Y = Y_whole_raw, X = X_child_whole_raw, Z = Z_whole_raw
)

# Scenario C: Combined
fit_comb_whole <- refit_msrrr(
  object = fit_comb, Y = Y_whole_raw, X = X_comb_whole_raw, Z = Z_whole_raw
)

# The bootstrap keeps the fixed whole-sample transformation. Obtain the
# processed matrices explicitly from each final fit; do not re-estimate them
# inside individual bootstrap samples.
preg_whole_processed <- transform_msrrr_data(
  fit_preg_whole, X = X_preg_whole_raw, Z = Z_whole_raw, Y = Y_whole_raw
  )

child_whole_processed <- transform_msrrr_data(
  fit_child_whole, X = X_child_whole_raw, Z = Z_whole_raw, Y = Y_whole_raw
  )

comb_whole_processed <- transform_msrrr_data(
  fit_comb_whole, X = X_comb_whole_raw, Z = Z_whole_raw, Y = Y_whole_raw
  )

X_preg_whole <- preg_whole_processed$X
X_child_whole <- child_whole_processed$X
X_comb_whole <- comb_whole_processed$X
Z_whole <- preg_whole_processed$Z
Y_whole <- preg_whole_processed$Y

stopifnot(
  isTRUE(all.equal(Z_whole, child_whole_processed$Z)),
  isTRUE(all.equal(Z_whole, comb_whole_processed$Z)),
  isTRUE(all.equal(Y_whole, child_whole_processed$Y)),
  isTRUE(all.equal(Y_whole, comb_whole_processed$Y))
)

# 9.3 Generate Final Heatmaps
# ------------------------------------------------------------------------------
# Plotting helper functions used below are defined in functions/msrrr_plotting_functions.R.
# ------------------------------------------------------------------------------
# BLOCK 9.3: EXPORT FINAL PDF
# ------------------------------------------------------------------------------

X_comb_whole <- cbind(X_preg_whole, X_child_whole)
whole_heatmap_grobs <- list(
  plot_final_heatmap(
    fit_preg_whole, X_preg_whole, "Pregnancy (Whole Sample)", draw = FALSE
  )$gtable,
  plot_final_heatmap(
    fit_child_whole, X_child_whole, "Childhood (Whole Sample)", draw = FALSE
  )$gtable,
  plot_final_heatmap(
    fit_comb_whole, X_comb_whole, "Combined (Whole Sample)", draw = FALSE
  )$gtable
)
save_grobs_pdf(
  whole_heatmap_grobs,
  file.path(results_dir, "Exposome_Signatures_WHOLE_SAMPLE_final_prebootstrap.pdf"),
  width = 18, height = 30
)
graphics.off()
if (!is.null(dev.list())) {dev.off()}

save(fit_preg_whole, fit_child_whole, fit_comb_whole, 
     file = file.path(results_dir, "msrrr_WHOLE_SAMPLE_models.RData"))


# ------------------------------------------------------------------------------
# BLOCK 9.4: SUMMARY TABLE FOR WHOLE SAMPLE MODELS
# ------------------------------------------------------------------------------
message(">>> [Block 9.4] Creating summary table for Whole Sample Models...")

# Initialize the dataframe
results_whole_df <- data.frame(
  Scenario      = c("Pregnancy", "Childhood", "Combined"),
  Opt_Rank      = results_df$Opt_Rank,
  Opt_Lambda    = results_df$Opt_Lambda,
  N_Exposures   = rep(NA, 3)
  #,In_Sample_pMSE = rep(NA, 3)
)

# refit_msrrr() returns the whole-sample refit directly, so its exposure
# coefficient matrix is stored in $B rather than $fit$B. An exposure is counted
# as selected when at least one outcome coefficient exceeds selection_tol.

results_whole_df$N_Exposures[1] <- sum(
  rowSums(abs(fit_preg_whole$B) > selection_tol) > 0
)
results_whole_df$N_Exposures[2] <- sum(
  rowSums(abs(fit_child_whole$B) > selection_tol) > 0
)
results_whole_df$N_Exposures[3] <- sum(
  rowSums(abs(fit_comb_whole$B) > selection_tol) > 0
)
# Print and Save
results_whole_df
results_df

write.csv2(results_whole_df, file.path(results_dir, "model_metrics_results_whole_sample_summary.csv"), row.names = FALSE)

message(">>> Whole Sample results table saved to: ", file.path(results_dir, "model_metrics_results_whole_sample_summary.csv"))

# Whole-sample heatmaps are generated by canonical Block 9.3 above.

# ==============================================================================
# BLOCK 10: APPROACH A - PENALISED BOOTSTRAP
# ==============================================================================

# APPROACH A: PENALISED BOOTSTRAP STABILITY AFTER MODEL SELECTION
# ------------------------------------------------------------------------------
# The final whole-sample models are bootstrapped while keeping the selected
# optimal rank and lambda fixed. In this penalized sparse setting, the bootstrap
# output is used to assess the stability of the selected exposure-outcome
# coefficients, rather than as classical post-selection inference.
#
# The two main post-bootstrap criteria used below are:
#
# 1. Selection probability (sel_prob)
#    Proportion of bootstrap samples in which each exposure-outcome coefficient
#    is selected, i.e. non-zero after applying the same penalized model.
#
# 2. Directional consistency (sign_consistency)
#    Among the bootstrap samples in which the coefficient is selected, the
#    proportion that keep the same sign as the original whole-sample coefficient.
#
# Two filters are visualised in Block 11:
#
# 1) Selection-stability filter:
#    sel_prob >= sel_prob_threshold
#
# 2) Stricter stability filter:
#    sel_prob >= sel_prob_threshold and
#    sign_consistency >= sign_consistency_threshold
#
# There is no universal best threshold. Higher thresholds are more conservative
# and may retain fewer associations; lower thresholds may retain more associations
# but can also include less stable ones. The default values used here are
# intended as practical starting points, not fixed rules.
# ==============================================================================

# Shared definition of a selected coefficient and of the two stability filters.
# These values are passed to the reusable bootstrap summary so that the CSV
# records both the continuous metrics and the TRUE/FALSE decisions.
sel_prob_threshold <- 0.90
sign_consistency_threshold <- 0.80

# ==============================================================================
# BLOCK 10A: APPROACH A - PENALISED BOOTSTRAP, PREGNANCY
# ==============================================================================
message(">>> [Block 10A] Starting penalised bootstrap: Pregnancy...")

# 10.1 Data Preparation
# ------------------------------------------------------------------------------
X_boot <- X_preg_whole
Z_boot <- Z_whole
Y_boot <- Y_whole
n_boot <- n_boot_global

# 10.2 Generate Resampling Indices
# ------------------------------------------------------------------------------
boot_sampling <- generate_resampling_msrrr(
  blocks             = list(X_boot, Z_boot, Y_boot),
  n_boot             = n_boot,
  balanced           = TRUE,
  keep_all_variables = FALSE,
  verbose            = TRUE
)

# 10.3 Pre-processing: Filter Covariates (Z) for Variance Stability
# ------------------------------------------------------------------------------
sd_null <- boot_sampling$sd_null
if (!is.null(sd_null)) {
  message("    -> Removing covariates with zero variance in some partitions...")
  Z_boot <- Z_boot[, -which(colnames(Z_boot) %in% names(sd_null[[2]]))]
} else {
  message("    -> OK: No covariate filtering needed.")
}

# 10.4 Parallel Computing Setup (PSOCK)
# ------------------------------------------------------------------------------
require(doSNOW)
num_cores <- detectCores() - 1
cl        <- parallel::makeCluster(num_cores, type = "PSOCK")

registerDoParallel(cl)
registerDoSNOW(cl)

# Progress Bar Setup
pb       <- txtProgressBar(max = n_boot, style = 3)
progress <- make_progress_callback(pb)
opts     <- list(progress = progress)

# 10.5 Core Bootstrap Loop (foreach)
# ------------------------------------------------------------------------------
message(paste("    -> Running", n_boot, "iterations on", num_cores, "cores..."))

boot_output <- foreach(
  i           = 1:n_boot,
  .combine    = "cbind",
  .packages   = c("glmnet", "rrpack"),
  .options.snow = opts
) %dopar% {
  
  # Core Worker Function
  unlist(msrrr_bootstrap_k(
    Y           = Y_boot,
    X           = X_boot,
    Z           = Z_boot,
    inds        = boot_sampling$full_idx[i],
    family      = outcome_families,
    familygroup = family_mapping,
    nrank       = results_whole_df$Opt_Rank[1], # Pregnancy Rank
    lambda      = results_whole_df$Opt_Lambda[1]      # Pregnancy Lambda
  ))
}

stopCluster(cl)
close(pb)

# 10.6 Advanced Statistical Extraction
# ------------------------------------------------------------------------------
message("    -> Calculating stability metrics...")

original_betas <- fit_preg_whole$B
colnames(original_betas) <- colnames(Y_boot)
rownames(original_betas) <- colnames(X_boot)

res_preg <- summarize_msrrr_bootstrap(
  boot_output = boot_output,
  original_B = original_betas,
  lambda = results_whole_df$Opt_Lambda[1],
  selection_tol = selection_tol,
  sel_prob_threshold = sel_prob_threshold,
  sign_consistency_threshold = sign_consistency_threshold
)

View(res_preg)

# 10.7 Save Scenario Results
# ------------------------------------------------------------------------------
write.csv2(res_preg, file.path(results_dir, "BOOTSTRAP_STATISTICS_PREGNANCY_wholesample.csv"), row.names = FALSE)
write.csv2(
  res_preg[res_preg$passes_stability_filter, , drop = FALSE],
  file.path(results_dir, "STABLE_ASSOCIATIONS_PREGNANCY_wholesample.csv"),
  row.names = FALSE
)


# ==============================================================================
# BLOCK 10B: APPROACH A - PENALISED BOOTSTRAP, CHILDHOOD
# ==============================================================================
message(">>> [Block 10B] Starting penalised bootstrap: Childhood...")

# 10.1 Data Preparation
# ------------------------------------------------------------------------------
X_boot <- X_child_whole
Z_boot <- Z_whole
Y_boot <- Y_whole
n_boot <- n_boot_global

# 10.2 Generate Resampling Indices
# ------------------------------------------------------------------------------
boot_sampling <- generate_resampling_msrrr(
  blocks             = list(X_boot, Z_boot, Y_boot),
  n_boot             = n_boot,
  balanced           = TRUE,
  keep_all_variables = FALSE,
  verbose            = TRUE
)

# 10.3 Pre-processing: Filter Covariates (Z) for Variance Stability
# ------------------------------------------------------------------------------
sd_null <- boot_sampling$sd_null
if (!is.null(sd_null)) {
  message("    -> Removing covariates with zero variance in some partitions...")
  Z_boot <- Z_boot[, -which(colnames(Z_boot) %in% names(sd_null[[2]]))]
} else {
  message("    -> OK: No covariate filtering needed.")
}

# 10.4 Parallel Computing Setup (PSOCK)
# ------------------------------------------------------------------------------
require(doSNOW)
num_cores <- detectCores() - 1
cl        <- parallel::makeCluster(num_cores, type = "PSOCK")

registerDoParallel(cl)
registerDoSNOW(cl)

# Progress Bar Setup
pb       <- txtProgressBar(max = n_boot, style = 3)
progress <- make_progress_callback(pb)
opts     <- list(progress = progress)

# 10.5 Core Bootstrap Loop (foreach)
# ------------------------------------------------------------------------------
message(paste("    -> Running", n_boot, "iterations on", num_cores, "cores..."))

boot_output <- foreach(
  i           = 1:n_boot,
  .combine    = "cbind",
  .packages   = c("glmnet", "rrpack"),
  .options.snow = opts
) %dopar% {
  
  # Core Worker Function
  unlist(msrrr_bootstrap_k(
    Y           = Y_boot,
    X           = X_boot,
    Z           = Z_boot,
    inds        = boot_sampling$full_idx[i],
    family      = outcome_families,
    familygroup = family_mapping,
    nrank       = results_whole_df$Opt_Rank[2], # Childhood Rank
    lambda      = results_whole_df$Opt_Lambda[2]      # Childhood Lambda
  ))
}

stopCluster(cl)
close(pb)

# 10.6 Advanced Statistical Extraction
# ------------------------------------------------------------------------------
message("    -> Calculating stability metrics...")

original_betas <- fit_child_whole$B
colnames(original_betas) <- colnames(Y_boot)
rownames(original_betas) <- colnames(X_boot)

res_child <- summarize_msrrr_bootstrap(
  boot_output = boot_output,
  original_B = original_betas,
  lambda = results_whole_df$Opt_Lambda[2],
  selection_tol = selection_tol,
  sel_prob_threshold = sel_prob_threshold,
  sign_consistency_threshold = sign_consistency_threshold
)

View(res_child)

# 10.7 Save Scenario Results
# ------------------------------------------------------------------------------
write.csv2(res_child, file.path(results_dir, "BOOTSTRAP_STATISTICS_CHILDHOOD_wholesample.csv"), row.names = FALSE)
write.csv2(
  res_child[res_child$passes_stability_filter, , drop = FALSE],
  file.path(results_dir, "STABLE_ASSOCIATIONS_CHILDHOOD_wholesample.csv"),
  row.names = FALSE
)


# ==============================================================================
# BLOCK 10C: APPROACH A - PENALISED BOOTSTRAP, COMBINED
# ==============================================================================
message(">>> [Block 10C] Starting penalised bootstrap: Combined...")

# 10.1 Data Preparation
# ------------------------------------------------------------------------------
X_boot <- X_comb_whole
Z_boot <- Z_whole
Y_boot <- Y_whole
n_boot <- n_boot_global

# 10.2 Generate Resampling Indices
# ------------------------------------------------------------------------------
boot_sampling <- generate_resampling_msrrr(
  blocks             = list(X_boot, Z_boot, Y_boot),
  n_boot             = n_boot,
  balanced           = TRUE,
  keep_all_variables = FALSE,
  verbose            = TRUE
)

# 10.3 Pre-processing: Filter Covariates (Z) for Variance Stability
# ------------------------------------------------------------------------------
sd_null <- boot_sampling$sd_null
if (!is.null(sd_null)) {
  message("    -> Removing covariates with zero variance in some partitions...")
  Z_boot <- Z_boot[, -which(colnames(Z_boot) %in% names(sd_null[[2]]))]
} else {
  message("    -> OK: No covariate filtering needed.")
}

# 10.4 Parallel Computing Setup (PSOCK)
# ------------------------------------------------------------------------------
require(doSNOW)
num_cores <- detectCores() - 1
cl        <- parallel::makeCluster(num_cores, type = "PSOCK")

registerDoParallel(cl)
registerDoSNOW(cl)

# Progress Bar Setup
pb       <- txtProgressBar(max = n_boot, style = 3)
progress <- make_progress_callback(pb)
opts     <- list(progress = progress)

# 10.5 Core Bootstrap Loop (foreach)
# ------------------------------------------------------------------------------
message(paste("    -> Running", n_boot, "iterations on", num_cores, "cores..."))

boot_output <- foreach(
  i           = 1:n_boot,
  .combine    = "cbind",
  .packages   = c("glmnet", "rrpack"),
  .options.snow = opts
) %dopar% {
  
  # Core Worker Function
  unlist(msrrr_bootstrap_k(
    Y           = Y_boot,
    X           = X_boot,
    Z           = Z_boot,
    inds        = boot_sampling$full_idx[i],
    family      = outcome_families,
    familygroup = family_mapping,
    nrank       = results_whole_df$Opt_Rank[3],       # Combined Rank
    lambda      = results_whole_df$Opt_Lambda[3]      # Combined Lambda
  ))
}

stopCluster(cl)
close(pb)

# 10.6 Advanced Statistical Extraction
# ------------------------------------------------------------------------------
message("    -> Calculating stability metrics...")

original_betas <- fit_comb_whole$B
colnames(original_betas) <- colnames(Y_boot)
rownames(original_betas) <- colnames(X_boot)

res_comb_ <- summarize_msrrr_bootstrap(
  boot_output = boot_output,
  original_B = original_betas,
  lambda = results_whole_df$Opt_Lambda[3],
  selection_tol = selection_tol,
  sel_prob_threshold = sel_prob_threshold,
  sign_consistency_threshold = sign_consistency_threshold
)

View(res_comb_)

# 10.7 Save Scenario Results
# ------------------------------------------------------------------------------
write.csv2(res_comb_, file.path(results_dir, "BOOTSTRAP_STATISTICS_COMB_wholesample.csv"), row.names = FALSE)
write.csv2(
  res_comb_[res_comb_$passes_stability_filter, , drop = FALSE],
  file.path(results_dir, "STABLE_ASSOCIATIONS_COMB_wholesample.csv"),
  row.names = FALSE
)


# ==============================================================================
# BLOCK 11: APPROACH A PENALISED PRE- VS POST-BOOTSTRAP HEATMAPS
# ==============================================================================
# Each page compares the original whole-sample coefficients (left) with the
# same coefficient matrix after the strict bootstrap stability filter (right).
# Both panels use identical rows, row order, fixed cell sizes and colour scale.
#
# Set TRUE if a scenario contains too many rows for one readable page. The PDF
# will then contain one page per exposure family/domain within each scenario.
split_heatmaps_by_domain <- FALSE
heatmap_cellwidth <- 32
heatmap_cellheight <- 13

# Plotting helper functions used below are defined in functions/msrrr_plotting_functions.R.

output_pdf <- file.path(
  results_dir,
  "PENALISED_PRE_vs_POST_BOOTSTRAP_EXPOSOME_SIGNATURES_wholesample.pdf"
)
X_comb_whole <- cbind(X_preg_whole, X_child_whole)
heatmap_scenarios <- list(
  list(fit = fit_preg_whole, boot = res_preg, X = X_preg_whole,
       name = "Pregnancy"),
  list(fit = fit_child_whole, boot = res_child, X = X_child_whole,
       name = "Childhood"),
  list(fit = fit_comb_whole, boot = res_comb_, X = X_comb_whole,
       name = "Combined")
)

penalised_pre_post_grobs <- list()
for (scenario in heatmap_scenarios) {
  family_pages <- character(0)
  if (split_heatmaps_by_domain) {
    clean_names <- gsub(
      "_None|_Ter_2|_Ter_3", "", colnames(scenario$X)
    )
    row_info <- vars_info %>%
      dplyr::filter(variable_name %in% clean_names) %>%
      dplyr::distinct(variable_name, .keep_all = TRUE)
    family_pages <- unique(na.omit(
      row_info$family[match(clean_names, row_info$variable_name)]
    ))
  }
  
  if (length(family_pages) == 0L) {
    penalised_pre_post_grobs[[length(penalised_pre_post_grobs) + 1L]] <-
      plot_original_vs_stable(
      fit_whole = scenario$fit,
      boot_res = scenario$boot,
      X_mat = scenario$X,
      scenario_name = scenario$name,
      cellwidth = heatmap_cellwidth,
      cellheight = heatmap_cellheight,
      draw = FALSE
    )
  } else {
    for (family_page in family_pages) {
      penalised_pre_post_grobs[[length(penalised_pre_post_grobs) + 1L]] <-
        plot_original_vs_stable(
        fit_whole = scenario$fit,
        boot_res = scenario$boot,
        X_mat = scenario$X,
        scenario_name = scenario$name,
        family_filter = family_page,
        cellwidth = heatmap_cellwidth,
        cellheight = heatmap_cellheight,
        draw = FALSE
      )
    }
  }
}
save_grobs_pdf(
  penalised_pre_post_grobs, output_pdf, width = 20, height = 26
)
graphics.off()
message(">>> Comparative heatmaps exported to: ", output_pdf)


# ==============================================================================
# BLOCK 12: OPTIONAL APPROACH A PENALISED LOCO SENSITIVITY ANALYSIS
# ==============================================================================
# Set to TRUE to run this block. It is deliberately independent from the
# ordinary bootstrap above and from the unpenalised approach below.
#
# For each h_cohort level, this block:
#   1) excludes every participant from that cohort
#   2) recalculates X/Z standardisation and continuous-outcome [0,1] scaling
#      using only the retained participants
#   3) removes Z columns that become constant; when the original reference
#      cohort is excluded, it also removes one remaining cohort dummy to avoid
#      perfect collinearity with the intercept
#   4) refits the penalised Pregnancy, Childhood and Combined models with the
#      whole-sample selected rank and lambda, but init = NULL
#   5) allows the selected exposures to change between LOCO fits;
#   6) generates independent bootstrap samples within each retained LOCO sample
#      and applies the same sel_prob/sign_consistency criteria; and
#   7) summarises both raw fitted selection and bootstrap-stable selection
#      frequency/sign across the six cohort exclusions; and
#   8) exports the LOCO summary plots and pre- versus post-bootstrap comparisons
#
# REMINDER: change run_penalized_loco to TRUE to execute this entire block.
# Computational cost: 6 cohorts x 3 scenarios x n_boot_global bootstrap fits
# (9,000 fits when n_boot_global = 500).
# ==============================================================================
run_penalized_loco <- TRUE
run_penalized_loco_bootstrap <- TRUE

if (run_penalized_loco) {
  message(">>> [Optional LOCO] Starting penalised leave-one-cohort-out analysis...")
  if (run_penalized_loco_bootstrap) {
    message(
      ">>> [Optional LOCO] Bootstrap is enabled: ", n_boot_global,
      " independent resamples will be generated inside each LOCO scenario."
    )
  }
  
  loco_output_dir <- file.path(results_dir, "LOCO_PENALISED")
  dir.create(loco_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cohort_id <- as.character(covs_raw$h_cohort)
  if (length(cohort_id) != nrow(Y) ||
      nrow(X_preg) != nrow(Y) ||
      nrow(X_child) != nrow(Y) ||
      nrow(Z) != nrow(Y)) {
    stop(
      paste(
        "LOCO requires cohort_id, X_preg, X_child, Z and Y to contain",
        "the same participants in the same row order."
      )
    )
  }
  if (anyNA(cohort_id) || any(cohort_id == "")) {
    stop("h_cohort contains missing/empty values; LOCO membership is undefined.")
  }
  
  cohort_levels <- sort(unique(cohort_id))
  if (length(cohort_levels) < 2L) {
    stop("LOCO requires at least two observed cohorts.")
  }
  reference_cohort <- levels(as.factor(covs_raw$h_cohort))[1L]
  cohort_dummy_columns <- grep(
    "^h_cohort_", colnames(Z), value = TRUE
  )
  
  # Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.
  
  loco_scenario_specifications <- list(
    Pregnancy = list(
      X_source = X_preg,
      rank = results_whole_df$Opt_Rank[1],
      lambda = results_whole_df$Opt_Lambda[1],
      original_fit = fit_preg_whole
    ),
    Childhood = list(
      X_source = X_child,
      rank = results_whole_df$Opt_Rank[2],
      lambda = results_whole_df$Opt_Lambda[2],
      original_fit = fit_child_whole
    ),
    Combined = list(
      X_source = cbind(X_preg, X_child),
      rank = results_whole_df$Opt_Rank[3],
      lambda = results_whole_df$Opt_Lambda[3],
      original_fit = fit_comb_whole
    )
  )
  
  # Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.
  
  loco_fits <- setNames(vector("list", length(cohort_levels)), cohort_levels)
  loco_bootstrap_results <- setNames(
    vector("list", length(cohort_levels)), cohort_levels
  )
  loco_run_information <- vector("list", length(cohort_levels))
  
  for (cohort_index in seq_along(cohort_levels)) {
    omitted_cohort <- cohort_levels[cohort_index]
    keep_subject <- cohort_id != omitted_cohort
    
    message(
      ">>> [Optional LOCO] Excluding cohort ", omitted_cohort,
      " (", sum(!keep_subject), " participant(s)); retaining ",
      sum(keep_subject), "."
    )
    
    Z_loco_raw <- as.matrix(Z[keep_subject, , drop = FALSE])
    
    # If the original reference cohort is absent, the remaining cohort dummies
    # span all retained cohort categories and one must become the new reference.
    if (identical(omitted_cohort, reference_cohort)) {
      remaining_cohort_dummies <- intersect(
        cohort_dummy_columns, colnames(Z_loco_raw)
      )
      if (length(remaining_cohort_dummies) > 0L) {
        new_reference_dummy <- remaining_cohort_dummies[1L]
        message(
          "    -> Original reference cohort excluded; removing ",
          new_reference_dummy, " to define a new reference category."
        )
        Z_loco_raw <- Z_loco_raw[
          , colnames(Z_loco_raw) != new_reference_dummy, drop = FALSE
        ]
      }
    }
    
    Z_loco <- scale_loco_matrix(
      Z_loco_raw, matrix_name = "Z", drop_constant = TRUE
    )
    Y_loco <- scale_loco_outcomes(Y[keep_subject, , drop = FALSE])
    
    loco_fits[[omitted_cohort]] <- list()
    loco_bootstrap_results[[omitted_cohort]] <- list()
    for (scenario_name in names(loco_scenario_specifications)) {
      scenario <- loco_scenario_specifications[[scenario_name]]
      X_loco <- scale_loco_matrix(
        scenario$X_source[keep_subject, , drop = FALSE],
        matrix_name = paste0("X (", scenario_name, ")"),
        drop_constant = FALSE
      )
      
      fit_loco <- msrrr.fit(
        Y = Y_loco,
        X = X_loco,
        Z = Z_loco,
        family = outcome_families,
        familygroup = family_mapping,
        nrank = scenario$rank,
        lambda = scenario$lambda,
        init = NULL,
        ensure_intercept = TRUE
      )
      rownames(fit_loco$B) <- colnames(X_loco)
      colnames(fit_loco$B) <- colnames(Y_loco)
      loco_fits[[omitted_cohort]][[scenario_name]] <- fit_loco
      
      saveRDS(
        fit_loco,
        file.path(
          loco_output_dir,
          paste0(
            "PENALISED_LOCO_", scenario_name,
            "_without_cohort_", omitted_cohort, ".rds"
          )
        )
      )
      
      if (run_penalized_loco_bootstrap) {
        loco_bootstrap_result <- run_one_penalized_loco_bootstrap(
          fit_loco = fit_loco,
          X_loco = X_loco,
          Y_loco = Y_loco,
          Z_loco = Z_loco,
          nrank = scenario$rank,
          lambda = scenario$lambda,
          family = outcome_families,
          familygroup = family_mapping,
          scenario_name = scenario_name,
          omitted_cohort = omitted_cohort,
          n_boot = n_boot_global
        )
        loco_bootstrap_results[[omitted_cohort]][[scenario_name]] <-
          loco_bootstrap_result
        utils::write.csv2(
          loco_bootstrap_result,
          file.path(
            loco_output_dir,
            paste0(
              "PENALISED_LOCO_BOOTSTRAP_", toupper(scenario_name),
              "_WITHOUT_COHORT_", omitted_cohort, ".csv"
            )
          ),
          row.names = FALSE
        )
      }
    }
    
    loco_run_information[[cohort_index]] <- data.frame(
      omitted_cohort = omitted_cohort,
      n_omitted = sum(!keep_subject),
      n_retained = sum(keep_subject),
      n_covariates_retained = ncol(Z_loco),
      stringsAsFactors = FALSE
    )
  }
  
  loco_run_information <- do.call(rbind, loco_run_information)
  utils::write.csv2(
    loco_run_information,
    file.path(loco_output_dir, "PENALISED_LOCO_RUN_INFORMATION.csv"),
    row.names = FALSE
  )
  
  # Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.
  
  loco_summaries <- lapply(
    names(loco_scenario_specifications),
    function(scenario_name) {
      summarize_penalized_loco(
        scenario_name = scenario_name,
        scenario = loco_scenario_specifications[[scenario_name]],
        loco_fits = loco_fits,
        cohort_levels = cohort_levels
      )
    }
  )
  names(loco_summaries) <- names(loco_scenario_specifications)
  
  for (scenario_name in names(loco_summaries)) {
    utils::write.csv2(
      loco_summaries[[scenario_name]]$summary,
      file.path(
        loco_output_dir,
        paste0("PENALISED_LOCO_SUMMARY_", toupper(scenario_name), ".csv")
      ),
      row.names = FALSE
    )
  }
  saveRDS(
    loco_summaries,
    file.path(loco_output_dir, "PENALISED_LOCO_ALL_SUMMARIES.rds")
  )
  
  # Build a second set of LOCO coefficient matrices in which a coefficient is
  # retained only when it passes the bootstrap sel_prob + sign_consistency
  # filter within that specific cohort-excluded sample.
  loco_bootstrap_stable_summaries <- NULL
  if (run_penalized_loco_bootstrap) {
    loco_bootstrap_stable_fits <- loco_fits
    
    for (omitted_cohort in cohort_levels) {
      for (scenario_name in names(loco_scenario_specifications)) {
        fit_loco <- loco_fits[[omitted_cohort]][[scenario_name]]
        bootstrap_result <-
          loco_bootstrap_results[[omitted_cohort]][[scenario_name]]
        
        stable_mask <- matrix(
          FALSE,
          nrow = nrow(fit_loco$B),
          ncol = ncol(fit_loco$B),
          dimnames = dimnames(fit_loco$B)
        )
        exposure_index <- match(
          bootstrap_result$exposure, rownames(stable_mask)
        )
        outcome_index <- match(
          bootstrap_result$outcome, colnames(stable_mask)
        )
        if (anyNA(exposure_index) || anyNA(outcome_index)) {
          stop(
            "A LOCO bootstrap result could not be aligned with its fitted B ",
            "matrix for ", scenario_name, " without cohort ", omitted_cohort, "."
          )
        }
        stable_mask[cbind(exposure_index, outcome_index)] <-
          bootstrap_result$passes_stability_filter
        stable_mask[is.na(stable_mask)] <- FALSE
        
        stable_fit <- fit_loco
        stable_fit$B[!stable_mask] <- 0
        loco_bootstrap_stable_fits[[omitted_cohort]][[scenario_name]] <-
          stable_fit
      }
    }
    
    loco_bootstrap_stable_summaries <- lapply(
      names(loco_scenario_specifications),
      function(scenario_name) {
        summarize_penalized_loco(
          scenario_name = scenario_name,
          scenario = loco_scenario_specifications[[scenario_name]],
          loco_fits = loco_bootstrap_stable_fits,
          cohort_levels = cohort_levels
        )
      }
    )
    names(loco_bootstrap_stable_summaries) <-
      names(loco_scenario_specifications)
    
    for (scenario_name in names(loco_bootstrap_stable_summaries)) {
      utils::write.csv2(
        loco_bootstrap_stable_summaries[[scenario_name]]$summary,
        file.path(
          loco_output_dir,
          paste0(
            "PENALISED_LOCO_BOOTSTRAP_STABLE_SUMMARY_",
            toupper(scenario_name), ".csv"
          )
        ),
        row.names = FALSE
      )
    }
    saveRDS(
      loco_bootstrap_stable_summaries,
      file.path(
        loco_output_dir,
        "PENALISED_LOCO_BOOTSTRAP_STABLE_ALL_SUMMARIES.rds"
      )
    )
  }
  
  # LOCO summary and SVG-export plotting helpers used below are defined in
  # functions/msrrr_plotting_functions.R.

  loco_pdf <- file.path(
    loco_output_dir, "PENALISED_LOCO_SELECTION_AND_SIGN_SUMMARY.pdf"
  )
  loco_summary_grobs <- list()
  for (scenario_name in names(loco_summaries)) {
    loco_summary_grobs[[scenario_name]] <- plot_penalized_loco_summary(
      loco_summaries[[scenario_name]], scenario_name, draw = FALSE
    )
  }
  save_grobs_pdf(loco_summary_grobs, loco_pdf, 14, 30)
  for (scenario_name in names(loco_summaries)) {
    save_grob_svg(
      loco_summary_grobs[[scenario_name]],
      file.path(
        loco_output_dir,
        paste0("PENALISED_LOCO_SELECTION_", toupper(scenario_name), ".svg")
      ),
      width = 14, height = 30
    )
  }
  
  if (run_penalized_loco_bootstrap) {
    loco_bootstrap_pdf <- file.path(
      loco_output_dir,
      "PENALISED_LOCO_BOOTSTRAP_STABLE_SELECTION_AND_SIGN_SUMMARY.pdf"
    )
    stable_loco_grobs <- list()
    for (scenario_name in names(loco_bootstrap_stable_summaries)) {
      stable_definition <- paste0(
        "stable: sel_prob >= ", sel_prob_threshold,
        " and sign_consistency >= ", sign_consistency_threshold
      )
      stable_loco_grobs[[scenario_name]] <- plot_penalized_loco_summary(
        loco_bootstrap_stable_summaries[[scenario_name]],
        scenario_name,
        selection_definition = stable_definition,
        draw = FALSE
      )
    }
    save_grobs_pdf(stable_loco_grobs, loco_bootstrap_pdf, 14, 30)
    for (scenario_name in names(loco_bootstrap_stable_summaries)) {
      save_grob_svg(
        stable_loco_grobs[[scenario_name]],
        file.path(
          loco_output_dir,
          paste0(
            "PENALISED_LOCO_BOOTSTRAP_STABLE_",
            toupper(scenario_name), ".svg"
          )
        ),
        width = 14, height = 30
      )
    }
    save_grob_svg(
      gridExtra::arrangeGrob(grobs = stable_loco_grobs, ncol = 1),
      sub("\\.pdf$", ".svg", loco_bootstrap_pdf), 14, 90
    )

    # Direct pre- versus post-bootstrap LOCO comparison.
    loco_pre_post_pdf <- file.path(
      loco_output_dir,
      "PENALISED_LOCO_PRE_vs_POST_BOOTSTRAP_SUMMARY.pdf"
    )
    loco_pre_post_grobs <- list()
    for (scenario_name in names(loco_summaries)) {
      common_rows <- rownames(loco_summaries[[scenario_name]]$n_selected)[
        rowSums(
          loco_summaries[[scenario_name]]$n_selected > 0L |
            loco_bootstrap_stable_summaries[[scenario_name]]$n_selected > 0L
        ) > 0L
      ]
      pre_grob <- plot_penalized_loco_summary(
        loco_summaries[[scenario_name]], scenario_name,
        selection_definition = "selected in penalised fit",
        method_label = "Pre-bootstrap", rows_to_plot = common_rows,
        cellwidth = 32, cellheight = 13, draw = FALSE
      )
      post_grob <- plot_penalized_loco_summary(
        loco_bootstrap_stable_summaries[[scenario_name]], scenario_name,
        selection_definition = stable_definition,
        method_label = "Post-bootstrap", rows_to_plot = common_rows,
        cellwidth = 32, cellheight = 13, draw = FALSE
      )
      comparison_grob <- gridExtra::arrangeGrob(
        pre_grob, grid::nullGrob(), post_grob, ncol = 3,
        widths = grid::unit.c(
          grid::unit(1, "null"), grid::unit(0.25, "in"),
          grid::unit(1, "null")
        ),
        top = grid::textGrob(
          paste0(scenario_name, " — penalised LOCO sensitivity"),
          gp = grid::gpar(fontface = "bold", fontsize = 14)
        )
      )
      loco_pre_post_grobs[[scenario_name]] <- comparison_grob
      svg_name <- file.path(
        loco_output_dir,
        paste0(
          "PENALISED_LOCO_PRE_vs_POST_BOOTSTRAP_",
          toupper(scenario_name), ".svg"
        )
      )
      save_grob_svg(comparison_grob, svg_name, 22, 30)
    }
    save_grobs_pdf(loco_pre_post_grobs, loco_pre_post_pdf, 22, 30)
  }
  
  save(
    loco_fits,
    loco_bootstrap_results,
    loco_bootstrap_stable_summaries,
    loco_run_information,
    file = file.path(loco_output_dir, "PENALISED_LOCO_FITS.RData")
  )
  message(">>> [Optional LOCO] Results saved under: ", loco_output_dir)
}


# ==============================================================================
# BLOCK 13: APPROACH B - UNPENALISED WHOLE-SAMPLE MODELS
# ==============================================================================
# After completing the shared models in Block 9, users choosing Approach B may
# continue directly here. Blocks 10-12 belong only to Approach A and are not
# required for the unpenalised model, bootstrap, heatmaps or LOCO analysis.
#
# This is an alternative to the penalised stability approach above. Exposures
# are selected once from each final penalised WHOLE-SAMPLE model. The selected
# subset and rank are then fixed, the model is refitted with lambda = 0 and
# init = NULL, and the same unpenalised model is fitted in every bootstrap.
# P-values and 95% intervals are approximate and conditional on the preceding
# exposure-selection step. No multiple-testing correction is applied.
# ==============================================================================
message(">>> [Approach B] Refitting selected whole-sample models with lambda = 0...")

# BLOCK 13A: PREGNANCY UNPENALISED WHOLE-SAMPLE MODEL
fit_preg_unpenalized <- refit_msrrr_unpenalized(
  object = fit_preg_whole,
  Y = Y_whole_raw,
  X = X_preg_whole_raw,
  Z = Z_whole_raw,
  init = NULL
)

# BLOCK 13B: CHILDHOOD UNPENALISED WHOLE-SAMPLE MODEL
fit_child_unpenalized <- refit_msrrr_unpenalized(
  object = fit_child_whole,
  Y = Y_whole_raw,
  X = X_child_whole_raw,
  Z = Z_whole_raw,
  init = NULL
)

# BLOCK 13C: COMBINED UNPENALISED WHOLE-SAMPLE MODEL
fit_comb_unpenalized <- refit_msrrr_unpenalized(
  object = fit_comb_whole,
  Y = Y_whole_raw,
  X = X_comb_whole_raw,
  Z = Z_whole_raw,
  init = NULL
)

save(
  fit_preg_unpenalized,
  fit_child_unpenalized,
  fit_comb_unpenalized,
  file = file.path(results_dir, "msrrr_UNPENALISED_WHOLE_SAMPLE_models.RData")
)

# ==============================================================================
# BLOCK 14: APPROACH B - CONDITIONAL UNPENALISED BOOTSTRAP
# ==============================================================================
# The exposure subsets and ranks come from the shared penalised whole-sample
# models in Block 9. Each bootstrap fit keeps that candidate set fixed, uses
# lambda = 0 and starts independently with init = NULL.
message(">>> [Block 14] Running conditional unpenalised bootstrap...")

# Run the fixed-subset, lambda = 0 bootstrap for one scenario.
# Workflow helper functions used below are defined in functions/msrrr_tutorial_workflow_functions.R.

# BLOCK 14A: PREGNANCY CONDITIONAL UNPENALISED BOOTSTRAP
res_preg_unpenalized <- run_unpenalized_bootstrap(
  fit_preg_unpenalized, Y_whole, Z_whole,
  outcome_families, family_mapping, n_boot_global, "Pregnancy"
)

# BLOCK 14B: CHILDHOOD CONDITIONAL UNPENALISED BOOTSTRAP
res_child_unpenalized <- run_unpenalized_bootstrap(
  fit_child_unpenalized, Y_whole, Z_whole,
  outcome_families, family_mapping, n_boot_global, "Childhood"
)

# BLOCK 14C: COMBINED CONDITIONAL UNPENALISED BOOTSTRAP
res_comb_unpenalized <- run_unpenalized_bootstrap(
  fit_comb_unpenalized, Y_whole, Z_whole,
  outcome_families, family_mapping, n_boot_global, "Combined"
)

write.csv2(
  res_preg_unpenalized,
  file.path(results_dir, "UNPENALISED_BOOTSTRAP_INFERENCE_PREGNANCY.csv"),
  row.names = FALSE
)
write.csv2(
  res_child_unpenalized,
  file.path(results_dir, "UNPENALISED_BOOTSTRAP_INFERENCE_CHILDHOOD.csv"),
  row.names = FALSE
)
write.csv2(
  res_comb_unpenalized,
  file.path(results_dir, "UNPENALISED_BOOTSTRAP_INFERENCE_COMBINED.csv"),
  row.names = FALSE
)

# ==============================================================================
# BLOCK 15: APPROACH B UNPENALISED PRE- VS POST-BOOTSTRAP HEATMAPS
# ==============================================================================
# The right panel colours only coefficients with p_value < 0.05.
# Non-significant betas remain unchanged in model objects and CSV files; NA is
# used only in the temporary plotting matrix so those cells are not coloured.

# Plotting helper functions used throughout this block are defined in
# functions/msrrr_plotting_functions.R.

split_heatmaps_by_domain <- FALSE
unpenalized_heatmap_scenarios <- list(
  list(
    fit = fit_preg_unpenalized, results = res_preg_unpenalized,
    name = "Pregnancy"
  ),
  list(
    fit = fit_child_unpenalized, results = res_child_unpenalized,
    name = "Childhood"
  ),
  list(
    fit = fit_comb_unpenalized, results = res_comb_unpenalized,
    name = "Combined"
  )
)

unpenalized_pdf <- file.path(
  results_dir,
  "UNPENALISED_PRE_vs_POST_BOOTSTRAP_EXPOSOME_SIGNATURES_wholesample.pdf"
)
unpenalized_pre_post_grobs <- list()
for (scenario in unpenalized_heatmap_scenarios) {
  family_pages <- character(0)
  if (split_heatmaps_by_domain) {
    family_pages <- unique(na.omit(
      get_heatmap_row_family(colnames(scenario$fit$X_selected))
    ))
  }
  if (length(family_pages) == 0L) {
    unpenalized_pre_post_grobs[[length(unpenalized_pre_post_grobs) + 1L]] <-
      plot_unpenalized_pre_post(
        scenario$fit, scenario$results, scenario$name, draw = FALSE
      )
  } else {
    for (family_page in family_pages) {
      unpenalized_pre_post_grobs[[length(unpenalized_pre_post_grobs) + 1L]] <-
        plot_unpenalized_pre_post(
          scenario$fit, scenario$results, scenario$name, family_page,
          draw = FALSE
        )
    }
  }
}
save_grobs_pdf(
  unpenalized_pre_post_grobs, unpenalized_pdf, width = 20, height = 26
)
graphics.off()
message(">>> Unpenalised pre/post heatmaps exported to: ", unpenalized_pdf)

# SVG is a single-page format: export one editable file per scenario.
for (scenario_index in seq_along(unpenalized_heatmap_scenarios)) {
  scenario <- unpenalized_heatmap_scenarios[[scenario_index]]
  svg_file <- file.path(
    results_dir,
    paste0("UNPENALISED_PRE_vs_POST_BOOTSTRAP_", toupper(scenario$name), ".svg")
  )
  scenario_grob <- plot_unpenalized_pre_post(
    scenario$fit, scenario$results, scenario$name, draw = FALSE
  )
  save_grob_svg(scenario_grob, svg_file, width = 20, height = 26)
}
unpenalized_combined_svg <- sub("\\.pdf$", ".svg", unpenalized_pdf)
save_grob_svg(
  gridExtra::arrangeGrob(grobs = unpenalized_pre_post_grobs, ncol = 1),
  unpenalized_combined_svg, width = 20, height = 72
)


# BLOCK 16: OPTIONAL APPROACH B UNPENALISED LOCO SENSITIVITY ANALYSIS
# ==============================================================================
# This block deliberately fixes the whole-sample penalised exposure subset for
# every cohort exclusion and every bootstrap sample. It does not reselect
# exposures inside a LOCO subsample. Consequently, the resulting p-values and
# intervals remain conditional on the original whole-sample selection.
# This block creates its own cohort definitions and does not require the
# penalised LOCO block. The preprocessing helper functions themselves are
# sourced from functions/msrrr_tutorial_workflow_functions.R.
# It also exports its own LOCO significance/sign summaries and pre- versus
# post-bootstrap comparison plots.
# Computational cost: another 6 cohorts x 3 scenarios x n_boot_global fits.
# ==============================================================================
run_unpenalized_loco <- TRUE
run_unpenalized_loco_bootstrap <- TRUE

if (run_unpenalized_loco) {
  cohort_id <- as.character(covs_raw$h_cohort)
  if (length(cohort_id) != nrow(Y) ||
      nrow(X_preg) != nrow(Y) ||
      nrow(X_child) != nrow(Y) ||
      nrow(Z) != nrow(Y)) {
    stop(
      paste(
        "LOCO requires covs_raw, X_preg, X_child, Z and Y to contain",
        "the same participants in the same row order."
      )
    )
  }
  if (anyNA(cohort_id) || any(cohort_id == "")) {
    stop("h_cohort contains missing/empty values; LOCO membership is undefined.")
  }
  cohort_levels <- sort(unique(cohort_id))
  if (length(cohort_levels) < 2L) {
    stop("LOCO requires at least two observed cohorts.")
  }
  reference_cohort <- levels(as.factor(covs_raw$h_cohort))[1L]
  cohort_dummy_columns <- grep(
    "^h_cohort_", colnames(Z), value = TRUE
  )
  
  message(
    ">>> [Optional LOCO] Starting unpenalised LOCO with fixed ",
    "whole-sample-selected exposure subsets..."
  )
  unpenalized_loco_dir <- file.path(results_dir, "LOCO_UNPENALISED")
  dir.create(unpenalized_loco_dir, recursive = TRUE, showWarnings = FALSE)
  
  unpenalized_loco_specifications <- list(
    Pregnancy = list(
      X_source = X_preg,
      whole_object = fit_preg_unpenalized
    ),
    Childhood = list(
      X_source = X_child,
      whole_object = fit_child_unpenalized
    ),
    Combined = list(
      X_source = cbind(X_preg, X_child),
      whole_object = fit_comb_unpenalized
    )
  )
  unpenalized_loco_fits <- setNames(
    vector("list", length(cohort_levels)), cohort_levels
  )
  unpenalized_loco_results <- setNames(
    vector("list", length(cohort_levels)), cohort_levels
  )
  
  for (omitted_cohort in cohort_levels) {
    keep_subject <- cohort_id != omitted_cohort
    message(
      ">>> [Optional unpenalised LOCO] Excluding cohort ", omitted_cohort,
      "; retaining ", sum(keep_subject), " participant(s)."
    )
    
    Z_loco_raw <- as.matrix(Z[keep_subject, , drop = FALSE])
    if (identical(omitted_cohort, reference_cohort)) {
      remaining_cohort_dummies <- intersect(
        cohort_dummy_columns, colnames(Z_loco_raw)
      )
      if (length(remaining_cohort_dummies) > 0L) {
        Z_loco_raw <- Z_loco_raw[
          , colnames(Z_loco_raw) != remaining_cohort_dummies[1L],
          drop = FALSE
        ]
      }
    }
    Z_loco <- scale_loco_matrix(
      Z_loco_raw, matrix_name = "Z", drop_constant = TRUE
    )
    Y_loco <- scale_loco_outcomes(Y[keep_subject, , drop = FALSE])
    unpenalized_loco_fits[[omitted_cohort]] <- list()
    unpenalized_loco_results[[omitted_cohort]] <- list()
    
    for (scenario_name in names(unpenalized_loco_specifications)) {
      scenario <- unpenalized_loco_specifications[[scenario_name]]
      selected_names <- colnames(scenario$whole_object$X_selected)
      X_loco_all <- scale_loco_matrix(
        scenario$X_source[keep_subject, , drop = FALSE],
        matrix_name = paste0("X (", scenario_name, ")"),
        drop_constant = FALSE
      )
      missing_selected <- setdiff(selected_names, colnames(X_loco_all))
      if (length(missing_selected) > 0L) {
        stop(
          "Selected exposure(s) missing from ", scenario_name, " LOCO data: ",
          paste(missing_selected, collapse = ", ")
        )
      }
      X_loco_selected <- X_loco_all[, selected_names, drop = FALSE]
      loco_rank <- min(
        scenario$whole_object$nrank,
        ncol(X_loco_selected),
        ncol(Y_loco)
      )
      
      fit_loco_unpenalized <- msrrr.fit(
        Y = Y_loco,
        X = X_loco_selected,
        Z = Z_loco,
        family = outcome_families,
        familygroup = family_mapping,
        nrank = loco_rank,
        lambda = 0,
        init = NULL,
        ensure_intercept = TRUE
      )
      rownames(fit_loco_unpenalized$B) <- selected_names
      colnames(fit_loco_unpenalized$B) <- colnames(Y_loco)
      unpenalized_object_loco <- list(
        fit = fit_loco_unpenalized,
        X_selected = X_loco_selected,
        nrank = loco_rank,
        lambda = 0
      )
      unpenalized_loco_fits[[omitted_cohort]][[scenario_name]] <-
        unpenalized_object_loco
      
      if (run_unpenalized_loco_bootstrap) {
        inference_loco <- run_unpenalized_bootstrap(
          unpenalized_object = unpenalized_object_loco,
          Y = Y_loco,
          Z = Z_loco,
          family = outcome_families,
          familygroup = family_mapping,
          n_boot = n_boot_global,
          scenario_name = paste0(
            scenario_name, " without cohort ", omitted_cohort
          )
        )
        inference_loco$omitted_cohort <- omitted_cohort
        inference_loco$scenario <- scenario_name
        unpenalized_loco_results[[omitted_cohort]][[scenario_name]] <-
          inference_loco
        write.csv2(
          inference_loco,
          file.path(
            unpenalized_loco_dir,
            paste0(
              "UNPENALISED_LOCO_INFERENCE_", toupper(scenario_name),
              "_WITHOUT_COHORT_", omitted_cohort, ".csv"
            )
          ),
          row.names = FALSE
        )
      }
    }
  }
  
  unpenalized_loco_summaries <- list()
  unpenalized_loco_pre_summaries <- list()
  for (scenario_name in names(unpenalized_loco_specifications)) {
    selected_names <- colnames(
      unpenalized_loco_specifications[[scenario_name]]$whole_object$X_selected
    )
    outcome_names <- colnames(Y_whole)
    n_significant <- n_positive <- n_negative <- matrix(
      0L, nrow = length(selected_names), ncol = length(outcome_names),
      dimnames = list(selected_names, outcome_names)
    )
    n_fitted <- n_fitted_positive <- n_fitted_negative <- matrix(
      0L, nrow = length(selected_names), ncol = length(outcome_names),
      dimnames = list(selected_names, outcome_names)
    )
    for (omitted_cohort in cohort_levels) {
      if (run_unpenalized_loco_bootstrap) {
        result <- unpenalized_loco_results[[omitted_cohort]][[scenario_name]]
        row_index <- match(result$exposure, selected_names)
        col_index <- match(result$outcome, outcome_names)
        valid <- !is.na(row_index) & !is.na(col_index) & result$significant
        for (i in which(valid)) {
          n_significant[row_index[i], col_index[i]] <-
            n_significant[row_index[i], col_index[i]] + 1L
          if (result$beta_unpenalized[i] > 0) {
            n_positive[row_index[i], col_index[i]] <-
              n_positive[row_index[i], col_index[i]] + 1L
          } else if (result$beta_unpenalized[i] < 0) {
            n_negative[row_index[i], col_index[i]] <-
              n_negative[row_index[i], col_index[i]] + 1L
          }
        }
      }
      fitted_B <- as.matrix(
        unpenalized_loco_fits[[omitted_cohort]][[scenario_name]]$fit$B
      )
      rownames(fitted_B) <- selected_names
      colnames(fitted_B) <- outcome_names
      fitted_selected <- abs(fitted_B) > selection_tol
      n_fitted <- n_fitted + fitted_selected
      n_fitted_positive <- n_fitted_positive +
        (fitted_selected & fitted_B > 0)
      n_fitted_negative <- n_fitted_negative +
        (fitted_selected & fitted_B < 0)
    }
    unpenalized_loco_pre_summaries[[scenario_name]] <- list(
      n_selected = n_fitted,
      n_positive = n_fitted_positive,
      n_negative = n_fitted_negative
    )
    if (run_unpenalized_loco_bootstrap) {
      summary_long <- expand.grid(
        exposure = selected_names,
        outcome = outcome_names,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      summary_long$n_significant_loco <- as.vector(n_significant)
      summary_long$n_positive_significant <- as.vector(n_positive)
      summary_long$n_negative_significant <- as.vector(n_negative)
      unpenalized_loco_summaries[[scenario_name]] <- list(
        summary = summary_long,
        n_selected = n_significant,
        n_positive = n_positive,
        n_negative = n_negative
      )
      write.csv2(
        summary_long,
        file.path(
          unpenalized_loco_dir,
          paste0(
            "UNPENALISED_LOCO_SIGNIFICANCE_SUMMARY_",
            toupper(scenario_name), ".csv"
          )
        ),
        row.names = FALSE
      )
    }
  }

  # Always export the unpenalised LOCO fits before bootstrap. This remains
  # available when run_unpenalized_loco_bootstrap = FALSE.
  unpenalized_loco_pre_pdf <- file.path(
    unpenalized_loco_dir,
    "UNPENALISED_LOCO_SELECTION_AND_SIGN_SUMMARY.pdf"
  )
  unpenalized_loco_pre_grobs <- list()
  for (scenario_name in names(unpenalized_loco_pre_summaries)) {
    unpenalized_loco_pre_grobs[[scenario_name]] <-
      plot_penalized_loco_summary(
        unpenalized_loco_pre_summaries[[scenario_name]],
        scenario_name,
        selection_definition = paste0("abs(beta) > ", selection_tol),
        method_label = "unpenalised LOCO - Pre-bootstrap",
        draw = FALSE
      )
  }
  save_grobs_pdf(
    unpenalized_loco_pre_grobs, unpenalized_loco_pre_pdf, 14, 30
  )
  for (scenario_name in names(unpenalized_loco_pre_summaries)) {
    save_grob_svg(
      unpenalized_loco_pre_grobs[[scenario_name]],
      file.path(
        unpenalized_loco_dir,
        paste0(
          "UNPENALISED_LOCO_SELECTION_", toupper(scenario_name), ".svg"
        )
      ),
      width = 14, height = 30
    )
  }
  
  # Significance and pre/post outputs require bootstrap results.
  if (run_unpenalized_loco_bootstrap) {
  unpenalized_loco_pdf <- file.path(
    unpenalized_loco_dir,
    "UNPENALISED_LOCO_SIGNIFICANCE_AND_SIGN_SUMMARY.pdf"
  )
  unpenalized_loco_grobs <- list()
  for (scenario_name in names(unpenalized_loco_summaries)) {
    unpenalized_loco_grobs[[scenario_name]] <- plot_penalized_loco_summary(
      unpenalized_loco_summaries[[scenario_name]],
      scenario_name,
      selection_definition = "p_value < 0.05",
      method_label = "unpenalised LOCO",
      draw = FALSE
    )
  }
  save_grobs_pdf(unpenalized_loco_grobs, unpenalized_loco_pdf, 14, 30)
  for (scenario_name in names(unpenalized_loco_summaries)) {
    save_grob_svg(
      unpenalized_loco_grobs[[scenario_name]],
      file.path(
        unpenalized_loco_dir,
        paste0(
          "UNPENALISED_LOCO_SIGNIFICANCE_", toupper(scenario_name), ".svg"
        )
      ),
      width = 14, height = 30
    )
  }
  save_grob_svg(
    gridExtra::arrangeGrob(grobs = unpenalized_loco_grobs, ncol = 1),
    sub("\\.pdf$", ".svg", unpenalized_loco_pdf), 14, 90
  )

  unpenalized_loco_pre_post_pdf <- file.path(
    unpenalized_loco_dir,
    "UNPENALISED_LOCO_PRE_vs_POST_BOOTSTRAP_SUMMARY.pdf"
  )
  unpenalized_loco_pre_post_grobs <- list()
  for (scenario_name in names(unpenalized_loco_summaries)) {
    common_rows <- rownames(
      unpenalized_loco_pre_summaries[[scenario_name]]$n_selected
    )[
      rowSums(
        unpenalized_loco_pre_summaries[[scenario_name]]$n_selected > 0L |
          unpenalized_loco_summaries[[scenario_name]]$n_selected > 0L
      ) > 0L
    ]
    pre_grob <- plot_penalized_loco_summary(
      unpenalized_loco_pre_summaries[[scenario_name]], scenario_name,
      selection_definition = "coefficient differs from zero",
      method_label = "Pre-bootstrap", rows_to_plot = common_rows,
      cellwidth = 32, cellheight = 13, draw = FALSE
    )
    post_grob <- plot_penalized_loco_summary(
      unpenalized_loco_summaries[[scenario_name]], scenario_name,
      selection_definition = "p_value < 0.05",
      method_label = "Post-bootstrap", rows_to_plot = common_rows,
      cellwidth = 32, cellheight = 13, draw = FALSE
    )
    comparison_grob <- gridExtra::arrangeGrob(
      pre_grob, grid::nullGrob(), post_grob, ncol = 3,
      widths = grid::unit.c(
        grid::unit(1, "null"), grid::unit(0.25, "in"),
        grid::unit(1, "null")
      ),
      top = grid::textGrob(
        paste0(scenario_name, " — unpenalised LOCO sensitivity"),
        gp = grid::gpar(fontface = "bold", fontsize = 14)
      )
    )
    unpenalized_loco_pre_post_grobs[[scenario_name]] <- comparison_grob
    svg_name <- file.path(
      unpenalized_loco_dir,
      paste0(
        "UNPENALISED_LOCO_PRE_vs_POST_BOOTSTRAP_",
        toupper(scenario_name), ".svg"
      )
    )
    save_grob_svg(comparison_grob, svg_name, 22, 30)
  }
  save_grobs_pdf(
    unpenalized_loco_pre_post_grobs,
    unpenalized_loco_pre_post_pdf,
    22, 30
  )
  }
  
  save(
    unpenalized_loco_fits,
    unpenalized_loco_results,
    unpenalized_loco_summaries,
    unpenalized_loco_pre_summaries,
    file = file.path(
      unpenalized_loco_dir, "UNPENALISED_LOCO_RESULTS.RData"
    )
  )
  message(
    ">>> [Optional LOCO] Unpenalised results saved under: ",
    unpenalized_loco_dir
  )
}


# ==============================================================================
