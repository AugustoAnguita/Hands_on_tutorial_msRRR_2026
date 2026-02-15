# ==============================================================================
# PROJECT: msRRR Analysis Pipeline (Local Version)
# DESCRIPTION: Unified script for Model Fitting, Bootstrapping, and Sparsity.
# AUTHOR: Augusto Anguita-Ruiz
# ==============================================================================

# ------------------------------------------------------------------------------
# BLOCK 1: SET UP & LIBRARIES
# ------------------------------------------------------------------------------

# Clean environment
rm(list = ls())

# Define all required packages
required_packages <- c(
  "tidyverse", "rrpack", "foreach", "doParallel", 
  "pheatmap", "glmnet", "fastDummies", 
  "gridExtra", "corrplot","caret","ggrepel","reshape2","RColorBrewer","grid"
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
if(file.exists("functions/msrrr_v3.R")) {
  source("functions/msrrr_v3.R") 
  source("functions/bootstrapping_functions_aug_25_03_24.r") 
  message(">>> Setup Complete. Libraries and Functions loaded.")
} else {
  stop("ERROR: Functions not found. Please check your 'functions' folder.")
}


# ------------------------------------------------------------------------------
# BLOCK 2: LOAD DATA
# ------------------------------------------------------------------------------

# Helper to read European CSVs (semicolon separated, comma decimals)
read_euro_csv <- function(filename) {
  read_delim(filename, delim = ";", 
             locale = locale(decimal_mark = ","), 
             show_col_types = FALSE)
}

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
  filter(
    (family %in% target_families & var_type == "numeric") |
      (variable_name %in% forced_factors)
  )

# 3.3 Split Variables by Period
vars_preg_names  <- vars_info %>% filter(period == "Pregnancy") %>% pull(variable_name)
vars_child_names <- vars_info %>% filter(period == "Postnatal") %>% pull(variable_name)

# Intersect with available data (Safety check)
vars_preg_names  <- intersect(vars_preg_names, colnames(exp_raw))
vars_child_names <- intersect(vars_child_names, colnames(exp_raw))


# 3.4 Process Exposome Matrices (X)
process_exposome_matrix <- function(df_raw, var_names, factors_list) {
  
  # 1. Select specific variables
  df_subset <- df_raw %>% select(all_of(var_names))
  
  # 2. Identify which 'forced_factors' are present
  present_factors <- intersect(factors_list, colnames(df_subset))
  
  if(length(present_factors) > 0) {
    # Ensure they are factors
    df_subset <- df_subset %>%
      mutate(across(all_of(present_factors), as.factor))
    
    # 3. Apply dummy_cols (remove first dummy to avoid collinearity)
    df_subset <- dummy_cols(df_subset, 
                            select_columns = present_factors,
                            remove_first_dummy = TRUE, 
                            remove_selected_columns = TRUE)
  }
  
  # 4. Convert to Matrix
  mat <- as.matrix(df_subset)
  return(mat)
}

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
  select(all_of(covs_list)) %>%
  mutate(
    e3_sex_None = case_when(
      e3_sex_None == "female" ~ 0,
      e3_sex_None == "male" ~ 1,
      TRUE ~ NA_real_
    ),
    h_cohort = as.factor(h_cohort),
    h_edumc_None = as.factor(h_edumc_None)
  )

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
  select(all_of(outcomes_list)) %>%
  # Force numeric on all columns (this handles hs_bmi_c_cat categorical -> numeric)
  mutate(across(everything(), as.numeric))

Y <- as.matrix(Y_df)
message(paste("    -> Outcomes (Y) Ready:", ncol(Y), "vars."))


# ------------------------------------------------------------------------------
# BLOCK 3.7: QUALITY CONTROL - NA COUNT PER VARIABLE
# ------------------------------------------------------------------------------
message(">>> [Block 3.7] Running NA Quality Checks (Counts)...")

# Function to report NAs per variable
report_nas <- function(data_obj, name) {
  if(is.matrix(data_obj)) data_obj <- as.data.frame(data_obj)
  
  message(paste0("\n--- Checking Dataset: ", name, " ---"))
  
  na_counts <- colSums(is.na(data_obj))
  vars_with_na <- na_counts[na_counts > 0]
  
  if(length(vars_with_na) > 0) {
    message(paste("    Found", length(vars_with_na), "variables with missing values."))
    print(vars_with_na)
  } else {
    message("    OK: No missing values found (Complete).")
  }
  
  n_complete <- sum(complete.cases(data_obj))
  pct_complete <- round((n_complete / nrow(data_obj)) * 100, 1)
  message(paste0("    Complete Cases (Rows): ", n_complete, "/", nrow(data_obj), 
                 " (", pct_complete, "%)"))
}

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
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, isna) %>%
  summarise(num.isna = n(), .groups = 'drop') %>%
  mutate(pct = num.isna / total * 100)

levels <- (missing.values %>% filter(isna == TRUE) %>% arrange(desc(pct)))$key
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
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  ggplot(aes(key, id, fill = isna)) +
  geom_raster(alpha=0.8) +
  scale_fill_manual(name = "", values = c('steelblue', 'tomato3'),
                    labels = c("Present", "Missing")) +
  scale_x_discrete(limits = levels) +
  labs(x = "Variable", y = "Row Number", title = "Missing values in rows") +
  coord_flip()

# 4. Display Plots
grid.arrange(percentage.plot, row.plot, ncol = 2)

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

generate_corrplot <- function(data_mat, plot_title, file_name, tl) {
  
  message(paste("    -> Processing:", plot_title))
  
  # Calculate Correlation (Pairwise Complete)
  M <- cor(data_mat, use = "pairwise.complete.obs")
  
  col_palette <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  
  # Save PDF
  pdf(file = file_name, width = 12, height = 10)
  corrplot(M, 
           method = "color", 
           col = col_palette(200), 
           type = "upper", 
           order = "hclust", 
           addCoef.col = "black", 
           tl.col = "black", 
           tl.cex = 0.6, 
           number.cex = tl, 
           diag = FALSE, 
           title = plot_title,
           mar = c(0,0,2,0))
  dev.off()
  
  message(paste("       Saved to:", file_name))
}

# Generate plots for all 3 datasets
generate_corrplot(X_preg, "Correlation Matrix: Pregnancy Exposome", "results/CORRPLOT_predictors_Pregnancy.pdf",tl=0.4)
generate_corrplot(X_child, "Correlation Matrix: Childhood Exposome", "results/CORRPLOT_predictors_Childhood.pdf",tl=0.2)
generate_corrplot(Y, "Correlation Matrix: Phenotypes", "results/CORRPLOT_Phenotypes.pdf",tl=1)

# Save only the necessary processed matrices and the codebook info
save(X_preg, X_child, Z, Y, Y_df, vars_info, 
     file = "results/checkpoint_preproc.RData")



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
  get_prop <- function(data_vec) {
    tab <- table(factor(data_vec, levels = all_levels, exclude = NULL))
    return(round(prop.table(tab), 3))
  }
  
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
# BLOCK 5.9: STANDARDIZATION & SCALING (Corrected for Data Leakage)
# ------------------------------------------------------------------------------
message(">>> [Block 5.9] Standardizing predictors (Z-score) and scaling outcomes [0-1]...")

# A. Helper function for Z-score (Predictors X and Covariates Z)
# ------------------------------------------------------------------------------
# We calculate Mean and SD on Train, then apply to both Train and Test.
standardize_pair <- function(train_mat, test_mat) {
  # Standardize Train and save the attributes (center and scale)
  train_scaled <- scale(train_mat)
  center_val   <- attr(train_scaled, "scaled:center")
  scale_val    <- attr(train_scaled, "scaled:scale")
  
  # Standardize Test using Train's parameters
  test_scaled  <- scale(test_mat, center = center_val, scale = scale_val)
  
  return(list(train = train_scaled, test = test_scaled))
}

# Apply to Exposome and Covariates
res_preg  <- standardize_pair(X_preg_train, X_preg_test)
X_preg_train <- res_preg$train; X_preg_test <- res_preg$test

res_child <- standardize_pair(X_child_train, X_child_test)
X_child_train <- res_child$train; X_child_test <- res_child$test

res_z     <- standardize_pair(Z_train, Z_test)
Z_train <- res_z$train; Z_test <- res_z$test


# B. Helper function for [0-1] Scaling (Outcomes Y)
# ------------------------------------------------------------------------------
# We use Min and Max from Train to scale both sets.
scale_01_pair <- function(train_mat, test_mat) {
  # Copy structures
  sc_train <- train_mat
  sc_test  <- test_mat
  
  for (col in colnames(train_mat)) {
    # Only scale if variable is not binary (more than 2 unique values)
    if (length(unique(na.omit(train_mat[, col]))) > 2) {
      min_val <- min(train_mat[, col], na.rm = TRUE)
      max_val <- max(train_mat[, col], na.rm = TRUE)
      range_val <- max_val - min_val
      
      if (range_val > 0) {
        sc_train[, col] <- (train_mat[, col] - min_val) / range_val
        sc_test[, col]  <- (test_mat[, col] - min_val) / range_val
      }
    }
  }
  return(list(train = sc_train, test = sc_test))
}

# Apply to Outcomes
res_y <- scale_01_pair(Y_train, Y_test)
Y_train <- res_y$train; Y_test <- res_y$test


# C. Specific Binarization: BMI Category
# ------------------------------------------------------------------------------
# 0 (Under+Normal) vs 1 (Over+Obese)
# After [0-1] scaling, the threshold 0.5 separates the groups correctly.
Y_train[, "hs_bmi_c_cat"] <- as.numeric(Y_train[, "hs_bmi_c_cat"] > 0.5)
Y_test[, "hs_bmi_c_cat"]  <- as.numeric(Y_test[, "hs_bmi_c_cat"] > 0.5)


# ------------------------------------------------------------------------------
# FINAL CHECKS
# ------------------------------------------------------------------------------
message("\n--- Scaling Verification (Train Set) ---")
print(apply(Y_train, 2, range, na.rm = TRUE))

message("\n--- BMI Binarization Check (Test Set) ---")
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
n_obs            <- nrow(Y_train)

# TECHNICAL NOTE ON MEAN IMPUTATION FOR LAMBDA CALCULATION:
# The svd() function and crossprod() cannot handle missing values (NA). 
# To calculate the gradient at zero (lambda_max), we temporarily replace NAs 
# in Y with the column means. This represents the most neutral statistical 
# assumption, effectively zeroing out the influence of missing observations 
# in the cross-product matrix.
Y_imputed <- Y_train
for(i in 1:ncol(Y_imputed)) {
  col_mean <- mean(Y_imputed[, i], na.rm = TRUE)
  Y_imputed[is.na(Y_imputed[, i]), i] <- col_mean
}

# ------------------------------------------------------------------------------
# 6.1 SCENARIO A: PREGNANCY EXPOSOME ONLY
# ------------------------------------------------------------------------------
message("\n>>> Scenario A: Fitting Pregnancy Model...")

# TECHNICAL RATIONALE FOR LAMBDA_MAX:
# Based on Karush-Kuhn-Tucker (KKT) conditions, lambda_max is the point where 
# all coefficients shrink to zero. We calculate this as the spectral norm 
# (largest singular value) of the cross-product (X'Y), normalized by n.
C_preg          <- crossprod(X_preg_train, Y_imputed)
spec_norm_preg  <- svd(C_preg, nu = 0, nv = 0)$d[1]
lam_max_preg    <- (spec_norm_preg / n_obs) * 1.1 #The 1.1 factor provides a numerical safety margin to ensure the regularization path begins strictly at the Null Model. 
# It accounts for potential floating-point inaccuracies, guaranteeing that all exposure coefficients start at exactly zero.

# TECHNICAL RATIONALE FOR LOGARITHMIC SEQUENCE:
# We use a 10^seq(log10...) approach to generate a geometric progression. 
# Penalized models are highly sensitive at low lambda values (where variables 
# enter the model) and less sensitive at high values. A logarithmic sequence 
# provides higher resolution for the critical "elbow" of the regularization path.
# Depth of 1e-4 will define how close to the lambda max we finish the sequence. It allows for the inclusion of denser models if the signal is weak. Decreasing it to 1e-6 forces the model to include more variables (higher sensitivity for weak signals), 
# while increasing it to 1e-2 restricts selection to only the strongest predictors (higher specificity).
# For epidemiologists, use a larger depth if you expect many small, pleiotropic effects, 
# and a smaller depth if the goal is to identify a few dominant environmental drivers.
lam_seq_preg <- 10^(seq(log10(lam_max_preg), log10(lam_max_preg * 1e-4), length.out = 100))

fit_preg <- msrrr(
  Y = Y_train, X = X_preg_train, Z = Z_train,
  family = outcome_families, familygroup = family_mapping,
  nrankseq = 2:4, lamseq = lam_seq_preg,
  foldid = foldid, method = 'CV', cv.criteria = "deviance", warm = TRUE
)


# ------------------------------------------------------------------------------
# 6.2 SCENARIO B: CHILDHOOD EXPOSOME ONLY
# ------------------------------------------------------------------------------
message("\n>>> Scenario B: Fitting Childhood Model...")

# Recalculating lambda_max for the Childhood covariate structure
C_child         <- crossprod(X_child_train, Y_imputed)
spec_norm_child <- svd(C_child, nu = 0, nv = 0)$d[1]
lam_max_child   <- (spec_norm_child / n_obs) * 1.1
lam_seq_child   <- 10^(seq(log10(lam_max_child), log10(lam_max_child * 1e-4), length.out = 100))

fit_child <- msrrr(
  Y = Y_train, X = X_child_train, Z = Z_train,
  family = outcome_families, familygroup = family_mapping,
  nrankseq = 2:4, lamseq = lam_seq_child,
  foldid = foldid, method = 'CV', cv.criteria = "deviance", warm = TRUE
)

# ------------------------------------------------------------------------------
# 6.3 SCENARIO C: COMBINED EXPOSOME (Pregnancy + Childhood)
# ------------------------------------------------------------------------------
message("\n>>> Scenario C: Fitting Combined Model...")

# Concatenating sources to evaluate joint effects
X_combined_train <- cbind(X_preg_train, X_child_train)

# Recalculating lambda_max for the joint dimension
C_comb          <- crossprod(X_combined_train, Y_imputed)
spec_norm_comb  <- svd(C_comb, nu = 0, nv = 0)$d[1]
lam_max_comb    <- (spec_norm_comb / n_obs) * 1.1
lam_seq_comb    <- 10^(seq(log10(lam_max_comb), log10(lam_max_comb * 1e-4), length.out = 100))

fit_comb <- msrrr(
  Y = Y_train, X = X_combined_train, Z = Z_train,
  family = outcome_families, familygroup = family_mapping,
  nrankseq = 2:4, lamseq = lam_seq_comb,
  foldid = foldid, method = 'CV', cv.criteria = "deviance", warm = TRUE
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
  CV_Deviance   = rep(NA, 3),
  pMSE_training = rep(NA, 3),
  pMSE_test     = rep(NA, 3)
)

# 3. Helper function to fill the table
fill_results <- function(row_idx, model_obj) {
  results_df$Opt_Rank[row_idx]    <<- model_obj$nrank
  results_df$Opt_Lambda[row_idx]  <<- model_obj$lam.opt
  # Count selected predictors (non-zero coefficients in B matrix)
  results_df$N_Exposures[row_idx] <<- sum(rowSums(model_obj$fit$B != 0) > 0)
  # Extract the minimum deviance from the tuning path
  results_df$CV_Deviance[row_idx] <<- unlist(model_obj$tunepath.opt)[which.min(unlist(model_obj$tunepath.opt))]
}

fill_results(1, fit_preg)
fill_results(2, fit_child)
fill_results(3, fit_comb)



# ------------------------------------------------------------------------------
# BLOCK 7.2: PREDICTION PERFORMANCE (pMSE)
# ------------------------------------------------------------------------------
message(">>> [Block 7.2] Calculating pMSE for Training and Test sets...")

# Scenario A: Pregnancy
fit_train_preg <- predict.msrrr(fit_preg, Y.new = Y_train, X.new = X_preg_train, Z = Z_train, family=outcome_families, familygroup=family_mapping, cv.criteria="pMSE")
fit_test_preg  <- predict.msrrr(fit_preg, Y.new = Y_test,  X.new = X_preg_test,  Z = Z_test,  family=outcome_families, familygroup=family_mapping, cv.criteria="pMSE")
results_df$pMSE_training[1] <- fit_train_preg$pred.perf
results_df$pMSE_test[1]     <- fit_test_preg$pred.perf

# Scenario B: Childhood
fit_train_child <- predict.msrrr(fit_child, Y.new = Y_train, X.new = X_child_train, Z = Z_train, family=outcome_families, familygroup=family_mapping, cv.criteria="pMSE")
fit_test_child  <- predict.msrrr(fit_child, Y.new = Y_test,  X.new = X_child_test,  Z = Z_test,  family=outcome_families, familygroup=family_mapping, cv.criteria="pMSE")
results_df$pMSE_training[2] <- fit_train_child$pred.perf
results_df$pMSE_test[2]     <- fit_test_child$pred.perf

# Scenario C: Combined
fit_train_comb <- predict.msrrr(fit_comb, Y.new = Y_train, X.new = X_combined_train, Z = Z_train, family=outcome_families, familygroup=family_mapping, cv.criteria="pMSE")
fit_test_comb  <- predict.msrrr(fit_comb, Y.new = Y_test,  X.new = X_combined_test,  Z = Z_test,  family=outcome_families, familygroup=family_mapping, cv.criteria="pMSE")
results_df$pMSE_training[3] <- fit_train_comb$pred.perf
results_df$pMSE_test[3]     <- fit_test_comb$pred.perf

results_df$Perf_Dif <- abs(results_df$pMSE_training - results_df$pMSE_test)
print(results_df)

# Save the primary model objects and the final results table
save(fit_preg, fit_child, fit_comb, results_df, 
     file = "results/msrrr_objects_model_metrics_results_training_CV.RData")


write.csv2(results_df,"results/model_metrics_results_training_CV.csv")


# ------------------------------------------------------------------------------
# BLOCK 7.3: DIAGNOSTIC PLOTS (Generalized Function)
# ------------------------------------------------------------------------------

message(">>> [Block 7.3] Defining and generating Diagnostic Plots...")

# Function to generate Rank and Lambda diagnostics for any scenario
generate_diagnostics <- function(model_obj, scenario_name) {
  
  # 1. Rank vs Deviance Plot
  rank_data <- data.frame(
    Rank = model_obj$nrankseq, 
    Deviance = unlist(model_obj$tunepath.opt)
  )
  
  p_rank <- ggplot(rank_data, aes(x = Rank, y = Deviance)) +
    geom_line(color = "grey") + 
    geom_point(size = 3, color = "#08737f") +
    theme_minimal() + 
    labs(title = paste(scenario_name, "Model: Rank vs Deviance"), 
         subtitle = "Optimization of latent factors")
  
  # 2. Lambda Boxplot (for the optimal Rank)
  opt_rank_name <- paste0("nrank_", model_obj$nrank)
  # Extract fold data (removing the summary cv.mean column)
  dtplot_raw <- model_obj$Tunepath[[opt_rank_name]]
  dtplot_folds <- t(dtplot_raw[, -ncol(dtplot_raw)]) 
  
  colnames(dtplot_folds) <- sprintf("%.1e", model_obj$lamseq)
  dtplot_long <- reshape2::melt(dtplot_folds)
  dtplot_long$Var2 <- factor(dtplot_long$Var2, levels = unique(dtplot_long$Var2))
  
  my_colors <- colorRampPalette(c("#08737f", "#f95559"))(100)
  
  p_lambda <- ggplot(dtplot_long, aes(x = Var2, y = value, fill = Var2)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = my_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1), 
          legend.position = "none") +
    labs(x = "Lambda (Log Scale)", y = "Deviance / Log-Likelihood", 
         title = "Regularization Path: Error across Folds")
  
  return(grid.arrange(p_rank, p_lambda, ncol = 1))
}

# Example: Generate plots for each scenario
diag_preg  <- generate_diagnostics(fit_preg, "Pregnancy")
diag_child <- generate_diagnostics(fit_child, "Childhood")
diag_comb  <- generate_diagnostics(fit_comb, "Combined")

pdf("results/Model_Diagnostics_training_CV.pdf", width = 8.5, height = 11)
plot(diag_preg)
plot(diag_child)
plot(diag_comb)
dev.off()



# ------------------------------------------------------------------------------
# BLOCK 7.4: COEFFICIENT LOLLIPOP PLOTS (Top Exposures)
# ------------------------------------------------------------------------------
message(">>> [Block 7.4] Plotting Coefficients for Key Outcomes...")

plot_coeffs <- function(model_obj, outcome_idx, X_matrix, scenario_label) {
  toplot <- data.frame(
    Coefficient = model_obj$fit$B[, outcome_idx],
    Exposure    = colnames(X_matrix)
  )
  
  toplot$vargroup <- ifelse(toplot$Coefficient > 0, "Positive", 
                            ifelse(toplot$Coefficient < 0, "Negative", "Null"))
  
  toplot_active <- subset(toplot, vargroup != "Null")
  
  # Ensure there is something to plot
  if(nrow(toplot_active) == 0) {
    message(paste("No variables selected for outcome", outcome_idx, "in", scenario_label))
    return(NULL)
  }
  
  ggplot(toplot, aes(x = Exposure, y = Coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = vargroup), size = 2) +
    geom_label_repel(data = toplot_active, 
                     aes(label = Exposure, color = vargroup), 
                     size = 2.5, max.overlaps = 20, segment.size = 0.2) +
    scale_color_manual(values = c("Negative" = "#CD534CFF", "Null" = "black", "Positive" = "#227CAD")) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), legend.position = "bottom") +
    labs(title = paste(colnames(Y_train)[outcome_idx], "-", scenario_label), 
         x = "Exposome Variables", y = "Estimated Beta")
}

# Generate specific plots for the tutorial
p_bw_comb <- plot_coeffs(fit_comb, 1, X_combined_train, "Combined")
p_as_comb <- plot_coeffs(fit_comb, 2, X_combined_train, "Combined")







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
#    similarly to find the most influential predictors across all domains. A strong color in both indicates the exposure is a powerful 
#      driver for both health domains, regardless of the unit.
# ------------------------------------------------------------------------------


# 8.1 CUSTOM COLOR PALETTE (Red-Yellow-Green)
# ------------------------------------------------------------------------------
# We create a palette where:
# - Deep Red: Strong Negative
# - Light Red/Yellow: Approaching Zero (Negative)
# - Pure Yellow: EXACTLY ZERO
# - Light Green/Yellow: Approaching Zero (Positive)
# - Deep Green: Strong Positive
# This ensures that "Yellow" is only seen at the sparsity threshold.

# 8.1 CUSTOM "SHARP" COLOR PALETTE (Verdes y Rojos)
# ------------------------------------------------------------------------------
# Definimos dos rampas que NO pasan por el blanco/amarillo hasta el final.
# Negativo (Verde): De Verde oscuro a Verde muy claro
neg_palette <- colorRampPalette(c("#1A9850", "#A6DBA0"))(50) 
# Neutral: Amarillo (Solo para el CERO absoluto)
zero_color  <- "grey" 
# Positivo (Rojo): De Rojo muy claro a Rojo oscuro
pos_palette <- colorRampPalette(c("#F4A582", "#D73027"))(50) 

custom_colors_pro <- c(neg_palette, zero_color, pos_palette)

# 8.2 CUSTOM BREAKS (La clave para evitar los blancos)
# ------------------------------------------------------------------------------
get_pro_breaks <- function(mat) {
  limit <- max(abs(mat), na.rm = TRUE)
  if(limit == 0) limit <- 0.1
  
  # Creamos una secuencia que va desde el límite hasta un valor ínfimo (1e-5)
  # Esto deja al color central (grey) un espacio casi invisible,
  # obligando a que cualquier valor pequeño tenga color.
  eps <- 1e-99
  breaks <- c(
    seq(-limit, -eps, length.out = 51), # 50 intervalos negativos
    0,                                  # El punto cero exacto
    seq(eps, limit, length.out = 51)    # 50 intervalos positivos
  )
  return(unique(sort(breaks)))
}

# 8.3 HEATMAP GENERATION FUNCTION
# ------------------------------------------------------------------------------
plot_final_heatmap <- function(model_obj, X_matrix, scenario_name) {
  
  B_mat <- model_obj$fit$B
  rownames(B_mat) <- colnames(X_matrix)
  colnames(B_mat) <- colnames(Y_train)
  
  # Filtramos filas donde TODAS las celdas sean exactamente 0
  active_rows <- rowSums(B_mat != 0) > 0
  coef_subset <- B_mat[active_rows, , drop = FALSE]
  
  if(nrow(coef_subset) == 0) return(NULL)
  
  # Anotaciones
  clean_names <- gsub("_None|_Ter_2|_Ter_3", "", rownames(coef_subset))
  row_info <- vars_info %>% 
    filter(variable_name %in% clean_names) %>%
    distinct(variable_name, .keep_all = TRUE)
  
  df_row_ann <- data.frame(
    Family = row_info$family[match(clean_names, row_info$variable_name)],
    row.names = rownames(coef_subset)
  )
  
  if(scenario_name == "Combined") {
    df_row_ann$Period <- ifelse(rownames(coef_subset) %in% colnames(X_preg_train), 
                                "Pregnancy", "Childhood")
  }
  
  pheatmap(
    coef_subset,
    main = paste("Environmental Signature:", scenario_name),
    color = custom_colors_pro,
    breaks = get_pro_breaks(coef_subset),
    annotation_row = df_row_ann,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    border_color = "white",
    fontsize_row = 8,
    cellwidth = 30
  )
}

# 8.4 Export
# ------------------------------------------------------------------------------
pdf("results/Exposome_Signatures_training_CV.pdf", width = 10, height = 12)
plot_final_heatmap(fit_preg, X_preg_train, "Pregnancy")
plot_final_heatmap(fit_child, X_child_train, "Childhood")

X_comb_train <- cbind(X_preg_train, X_child_train)
plot_final_heatmap(fit_comb, X_comb_train, "Combined")
dev.off()

message(">>> Block 8 Complete. Heatmaps updated with non-linear scale.")



load("results/msrrr_objects_model_metrics_results_training_CV.RData")
results_df <- read.csv2("results/model_metrics_results_training_CV.csv")








# ------------------------------------------------------------------------------
# BLOCK 9: WHOLE SAMPLE REFITTING
# ------------------------------------------------------------------------------
message(">>> [Block 9] Refitting models using the WHOLE SAMPLE...")
# ------------------------------------------------------------------------------
set.seed(12345)
# 9.1 Prepare Whole Sample Matrices
# ------------------------------------------------------------------------------
# Standardize the full set (internally consistent for the final model)
X_preg_whole  <- scale(X_preg)
X_child_whole <- scale(X_child)
Z_whole       <- scale(Z)
outcome_families <- list(gaussian(), binomial())
family_mapping   <- c(1, 2, 1, 1, 1, 2) 

# Scale Y [0-1] for the whole sample
Y_whole <- Y
for (col in colnames(Y_whole)) {
  if (length(unique(na.omit(Y_whole[, col]))) > 2) {
    min_val <- min(Y_whole[, col], na.rm = TRUE)
    max_val <- max(Y_whole[, col], na.rm = TRUE)
    Y_whole[, col] <- (Y_whole[, col] - min_val) / (max_val - min_val)
  }
}
Y_whole[, "hs_bmi_c_cat"] <- as.numeric(Y_whole[, "hs_bmi_c_cat"] > 0.5)

# 9.2 Refit Final Models (using msrrr.fit for fixed parameters)
# ------------------------------------------------------------------------------
# Scenario A: Pregnancy
fit_preg_whole <- msrrr.fit(
  Y = Y_whole, X = X_preg_whole, Z = Z_whole,
  family = outcome_families, familygroup = family_mapping,
  nrank = results_df$Opt_Rank[1], lambda = results_df$Opt_Lambda[1]
)

# Scenario B: Childhood
fit_child_whole <- msrrr.fit(
  Y = Y_whole, X = X_child_whole, Z = Z_whole,
  family = outcome_families, familygroup = family_mapping,
  nrank = results_df$Opt_Rank[2], lambda = results_df$Opt_Lambda[2]
)

# Scenario C: Combined
X_comb_whole <- cbind(X_preg_whole, X_child_whole)
fit_comb_whole <- msrrr.fit(
  Y = Y_whole, X = X_comb_whole, Z = Z_whole,
  family = outcome_families, familygroup = family_mapping,
  nrank = results_df$Opt_Rank[3], lambda = results_df$Opt_Lambda[3]
)

# 9.3 Generate Final Heatmaps
# ------------------------------------------------------------------------------
# Note: We create a small wrapper list to reuse our 'plot_final_heatmap' function


# ------------------------------------------------------------------------------
# FUNCIÓN DE HEATMAP ACTUALIZADA (Para PDF y Whole Sample)
# ------------------------------------------------------------------------------
plot_final_heatmap <- function(model_obj, X_matrix, scenario_name) {
  
  # 1. Extraer Beta (Manejamos tanto estructura msrrr como msrrr.fit)
  if(!is.null(model_obj$fit)) {
    B_mat <- model_obj$fit$B
  } else {
    B_mat <- model_obj$B
  }
  
  rownames(B_mat) <- colnames(X_matrix)
  colnames(B_mat) <- colnames(Y_train)
  
  # 2. Filtrar filas activas
  active_rows <- rowSums(B_mat != 0) > 0
  coef_subset <- B_mat[active_rows, , drop = FALSE]
  
  if(nrow(coef_subset) == 0) {
    message(paste("    ! No hay coeficientes para el heatmap en:", scenario_name))
    return(NULL)
  }
  
  # 3. Anotaciones de fila
  clean_names <- gsub("_None|_Ter_2|_Ter_3", "", rownames(coef_subset))
  row_info <- vars_info %>% 
    filter(variable_name %in% clean_names) %>%
    distinct(variable_name, .keep_all = TRUE)
  
  df_row_ann <- data.frame(
    Family = row_info$family[match(clean_names, row_info$variable_name)],
    row.names = rownames(coef_subset)
  )
  
  if(scenario_name == "Combined (Whole Sample)") {
    df_row_ann$Period <- ifelse(rownames(coef_subset) %in% colnames(X_preg_whole), 
                                "Pregnancy", "Childhood")
  }
  
  # 4. CREAR EL GRÁFICO
  p <- pheatmap(
    coef_subset,
    main = paste("Environmental Signature:", scenario_name),
    color = custom_colors_pro,
    breaks = get_pro_breaks(coef_subset),
    annotation_row = df_row_ann,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    border_color = "white",
    fontsize_row = 8,
    cellwidth = 30,
    silent = TRUE # Evita que se dibuje antes de tiempo
  )
  
  # 5. IMPRESIÓN EXPLÍCITA (Crucial para el PDF)
  print(p) 
}
# ------------------------------------------------------------------------------
# Bloque 9.3: Generar PDF Final
# ------------------------------------------------------------------------------

pdf("results/Exposome_Signatures_WHOLE_SAMPLE_final_prebootstrap.pdf", width = 10, height = 12)

# Escenario A
grid.newpage() # Limpia el lienzo para la página 1
plot_final_heatmap(fit_preg_whole, X_preg_whole, "Pregnancy (Whole Sample)")

# Escenario B
grid.newpage() # Salto a la página 2
plot_final_heatmap(fit_child_whole, X_child_whole, "Childhood (Whole Sample)")

# Escenario C
grid.newpage() # Salto a la página 3
X_comb_whole <- cbind(X_preg_whole, X_child_whole)
plot_final_heatmap(fit_comb_whole, X_comb_whole, "Combined (Whole Sample)")

dev.off()


save(fit_preg_whole, fit_child_whole, fit_comb_whole, 
     file = "results/msrrr_WHOLE_SAMPLE_models.RData")


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

# Function to extract counts and performance
# Note: Since msrrr.fit returns a list directly (not a nested $fit), 
# we access model_obj$B directly.
# summarize_whole <- function(row_idx, model_obj, X_mat,Y_whole,Z_whole) {
#   # 1. Count selected exposures (at least one non-zero coefficient across outcomes)
#   results_whole_df$N_Exposures[row_idx] <<- sum(rowSums(model_obj$B != 0) > 0)
#   
#   # 2. Calculate In-sample pMSE
#   # We wrap the model in a list $fit to satisfy predict.msrrr's expectation if needed,
#   # but our previous custom predict calls handle the object structure.
#   # Let's use a temporary object that mimics the msrrr class structure for the predictor
#   temp_obj <- list(B = model_obj$B, phi = model_obj$phi, family = outcome_families, familygroup = family_mapping)
#   
#   perf <- predict.msrrr(temp_obj, 
#                         Y.new = Y_whole, 
#                         X.new = X_mat, 
#                         Z.new = Z_whole, 
#                         family = outcome_families, 
#                         familygroup = family_mapping)
#   
#   results_whole_df$In_Sample_pMSE[row_idx] <<- perf$pred.perf
# }
# Fill the table
# summarize_whole(1, fit_preg_whole,  X_preg_whole, Y_whole, Z_whole)
# summarize_whole(2, fit_child_whole, X_child_whole, Y_whole, Z_whole)
# summarize_whole(3, fit_comb_whole,  X_comb_whole, Y_whole, Z_whole)


results_whole_df$N_Exposures[1] <- sum(rowSums(fit_preg_whole$B != 0) > 0)
results_whole_df$N_Exposures[2] <- sum(rowSums(fit_child_whole$B != 0) > 0)
results_whole_df$N_Exposures[3] <- sum(rowSums(fit_comb_whole$B != 0) > 0)
# Print and Save
results_whole_df
results_df


write.csv2(results_whole_df, "results/model_metrics_results_whole_sample_summary.csv", row.names = FALSE)

message(">>> Whole Sample results table saved to 'results/results_whole_sample_summary.csv'")








# ==============================================================================
# SECTION: ADVANCED BOOTSTRAP INFERENCE & STABILITY SELECTION
# ==============================================================================
#
# METHODOLOGICAL CHOICE: Stability Selection vs. Post-Selection Inference
# ------------------------------------------------------------------------------
# We implement the 'Stability Selection' approach by maintaining the optimal 
# penalty (lambda_opt) during bootstrapping. 
#
# Why avoid lambda = 0? 
# Refitting only selected variables without a penalty leads to 'p-value hacking' 
# or over-optimistic inference. By keeping the penalty, we account for 
# 'Selection Uncertainty'—the probability that a different variable might 
# have been chosen in a different sub-sample.
#
# THREE COMPLEMENTARY INFERENTIAL LENSES:
#
# 1. Sign Consistency (pval_ratio):
#    Formula: min(count_pos, count_neg) / max(count_pos, count_neg)
#    Interpretation: Measures 'Directional Stability'. A value < 0.05 
#    means >95% of bootstraps agree on the effect sign (e.g., always harmful).
#
# 2. Empirical Significance (pval_corrected):
#    Formula: 2 * min(mean(beta <= 0), mean(beta >= 0))
#    Interpretation: Non-parametric p-value. It represents the probability 
#    of the effect magnitude 'crossing zero'. It is robust against non-normal 
#    data distributions (Statistical Robustness).
#
# 3. Parametric Significance (th_pval):
#    Formula: 2 * pnorm(-abs(Beta_orig / SD_boot))
#    Interpretation: Classical p-value based on the 'Bootstrap Ratio'. 
#    It assumes the sampling distribution is Gaussian (Standard Inference).
#
# 4. Selection Probability (sel_prob):
#    Interpretation: The most critical metric for High-Dimensional data. 
#    It measures the frequency of a variable being selected across 500 
#    different versions of the dataset.
# ==============================================================================

# Final decision rule for "Stable Environmental Signatures":
# ------------------------------------------------------------------------------
# A finding is considered a "Priority Environmental Driver" if:
# - Stability: sel_prob > 0.80  (Selected in >80% of samples)
# - Consistency: pval_ratio < 0.05 (Sign is stable)
# - Robustness: pval_corrected < 0.05 (Empirically significant)




# ==============================================================================
# BLOCK 10: BOOTSTRAPPING WHOLE SAMPLE - PREGNANCY
# ==============================================================================
message(">>> [Block 10] Starting Bootstrap: Pregnancy Scenario...")

# 10.1 Data Preparation
# ------------------------------------------------------------------------------
X_boot <- X_preg_whole
Z_boot <- Z_whole
Y_boot <- Y_whole
n_boot <- 500  

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
progress <- function(n) setTxtProgressBar(pb, n)
opts     <- list(progress = progress)

# 10.5 Core Bootstrap Loop (foreach)
# ------------------------------------------------------------------------------
message(paste("    -> Running", n_boot, "iterations on", num_cores, "cores..."))

boot_output <- foreach(
  i           = 1:n_boot,
  .combine    = "cbind",
  .packages   = c(
    "RMTL", "data.table", "rrpack", "MCMCpack", "GIGrvg", "utils", "MASS", 
    "nortest", "MBSP", "foreach", "Matrix", "parallel", "glmnet", "spls", 
    "MXM", "dlnm", "splines", "mgcv", "doParallel", "ranger", 
    "palmerpenguins", "tidyverse", "kableExtra", "haven", "corrplot", "pheatmap"
  ),
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

# Prepare grid for mapping results
grid_df       <- expand.grid(exposures = colnames(X_boot), 
                             outcomes  = colnames(Y_boot), 
                             bootstraps = 1:n_boot)
grid_df$value <- as.vector(boot_output)

# Matrix representation for matrixStats performance
val_matrix <- matrix(grid_df$value, ncol = n_boot)
count_pos  <- matrixStats::rowSums2(val_matrix > 0)
count_neg  <- matrixStats::rowSums2(val_matrix < 0)
z_counts   <- cbind(count_pos, count_neg)

# Build Results Dataframe
res_preg <- cbind(
  matrixStats::rowMeans2(val_matrix),                    # 1. Mean Beta
  matrixStats::rowSds(val_matrix),                       # 2. SD
  matrixStats::rowMins(z_counts) / matrixStats::rowMaxs(z_counts), # 3. pval_ratio
  
  # Corrected P-value (Empirical density crossing zero)
  pmin(2 * matrixStats::rowMins(cbind(
    matrixStats::rowMeans2(val_matrix <= 0),
    matrixStats::rowMeans2(val_matrix >= 0)
  )), 1),                                                # 4. pval_corrected
  
  matrixStats::rowQuantiles(val_matrix, probs = c(0.025, 0.975)), # 5. 95% CI
  matrixStats::rowMeans2(val_matrix != 0)                # 6. Selection Probability
  
  
)

colnames(res_preg) <- c("mean", "sd", "pval_ratio", "pval_corrected", 
                        "lower_bound", "upper_bound", "sel_prob")

res_preg          <- data.frame(res_preg)
res_preg$exposure <- rep(colnames(X_boot), times = ncol(Y_boot))
res_preg$outcome  <- rep(colnames(Y_boot), each = ncol(X_boot))

# --- ADDITIONAL STATISTICS (Theoretical/Parametric) ---
#-----------Compute additional statistics using aggregated statistics
original_betas <- fit_preg_whole$B
colnames(original_betas) <- colnames(Y_boot)
rownames(original_betas) <- colnames(X_boot)
original_betas <- as.data.frame(original_betas)
original_betas$Predictor <- rownames(original_betas)
original_betas <- melt(original_betas, id.vars = "Predictor", variable.name = "Outcome", value.name = "Value")
original_betas$Predictor == res_preg$exposure
original_betas$Outcome == res_preg$outcome
tail <- qnorm(1 - .05 / 2)
res_preg$bootstrap_ratio <- original_betas$Value / res_preg$sd
res_preg$th_pval <- 2 * pnorm(abs(res_preg$bootstrap_ratio), lower.tail = FALSE)
res_preg$th_lower_bound <- original_betas$Value - res_preg$sd * tail
res_preg$th_upper_bound <- original_betas$Value + res_preg$sd * tail

View(res_preg)
# 10.7 Save Scenario Results
# ------------------------------------------------------------------------------
write.csv2(res_preg, "results/BOOTSTRAP_STATISTICS_PREGNANCY_wholesample.csv", row.names = FALSE)








# ==============================================================================
# BLOCK 10: BOOTSTRAPPING WHOLE SAMPLE - CHILDHOOD
# ==============================================================================
message(">>> [Block 10] Starting Bootstrap: Childhood Scenario...")

# 10.1 Data Preparation
# ------------------------------------------------------------------------------
X_boot <- X_child_whole
Z_boot <- Z_whole
Y_boot <- Y_whole
n_boot <- 500  

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
progress <- function(n) setTxtProgressBar(pb, n)
opts     <- list(progress = progress)

# 10.5 Core Bootstrap Loop (foreach)
# ------------------------------------------------------------------------------
message(paste("    -> Running", n_boot, "iterations on", num_cores, "cores..."))

boot_output <- foreach(
  i           = 1:n_boot,
  .combine    = "cbind",
  .packages   = c(
    "RMTL", "data.table", "rrpack", "MCMCpack", "GIGrvg", "utils", "MASS", 
    "nortest", "MBSP", "foreach", "Matrix", "parallel", "glmnet", "spls", 
    "MXM", "dlnm", "splines", "mgcv", "doParallel", "ranger", 
    "palmerpenguins", "tidyverse", "kableExtra", "haven", "corrplot", "pheatmap"
  ),
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

# Prepare grid for mapping results
grid_df       <- expand.grid(exposures = colnames(X_boot), 
                             outcomes  = colnames(Y_boot), 
                             bootstraps = 1:n_boot)
grid_df$value <- as.vector(boot_output)

# Matrix representation for matrixStats performance
val_matrix <- matrix(grid_df$value, ncol = n_boot)
count_pos  <- matrixStats::rowSums2(val_matrix > 0)
count_neg  <- matrixStats::rowSums2(val_matrix < 0)
z_counts   <- cbind(count_pos, count_neg)

# Build Results Dataframe
res_child <- cbind(
  matrixStats::rowMeans2(val_matrix),                    # 1. Mean Beta
  matrixStats::rowSds(val_matrix),                       # 2. SD
  matrixStats::rowMins(z_counts) / matrixStats::rowMaxs(z_counts), # 3. pval_ratio
  
  # Corrected P-value (Empirical density crossing zero)
  pmin(2 * matrixStats::rowMins(cbind(
    matrixStats::rowMeans2(val_matrix <= 0),
    matrixStats::rowMeans2(val_matrix >= 0)
  )), 1),                                                # 4. pval_corrected
  
  matrixStats::rowQuantiles(val_matrix, probs = c(0.025, 0.975)), # 5. 95% CI
  matrixStats::rowMeans2(val_matrix != 0)                # 6. Selection Probability
  
  
)

colnames(res_child) <- c("mean", "sd", "pval_ratio", "pval_corrected", 
                        "lower_bound", "upper_bound", "sel_prob")

res_child          <- data.frame(res_child)
res_child$exposure <- rep(colnames(X_boot), times = ncol(Y_boot))
res_child$outcome  <- rep(colnames(Y_boot), each = ncol(X_boot))

# --- ADDITIONAL STATISTICS (Theoretical/Parametric) ---
#-----------Compute additional statistics using aggregated statistics
original_betas <- fit_child_whole$B
colnames(original_betas) <- colnames(Y_boot)
rownames(original_betas) <- colnames(X_boot)
original_betas <- as.data.frame(original_betas)
original_betas$Predictor <- rownames(original_betas)
original_betas <- melt(original_betas, id.vars = "Predictor", variable.name = "Outcome", value.name = "Value")
original_betas$Predictor == res_child$exposure
original_betas$Outcome == res_child$outcome
tail <- qnorm(1 - .05 / 2)
res_child$bootstrap_ratio <- original_betas$Value / res_child$sd
res_child$th_pval <- 2 * pnorm(abs(res_child$bootstrap_ratio), lower.tail = FALSE)
res_child$th_lower_bound <- original_betas$Value - res_child$sd * tail
res_child$th_upper_bound <- original_betas$Value + res_child$sd * tail

View(res_child)
# 10.7 Save Scenario Results
# ------------------------------------------------------------------------------
write.csv2(res_child, "results/BOOTSTRAP_STATISTICS_CHILDHOOD_wholesample.csv", row.names = FALSE)






# ==============================================================================
# BLOCK 10: BOOTSTRAPPING WHOLE SAMPLE - COMBINED
# ==============================================================================
message(">>> [Block 10] Starting Bootstrap: Combined Scenario...")

# 10.1 Data Preparation
# ------------------------------------------------------------------------------
X_boot <- X_comb_whole
Z_boot <- Z_whole
Y_boot <- Y_whole
n_boot <- 500  

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
progress <- function(n) setTxtProgressBar(pb, n)
opts     <- list(progress = progress)

# 10.5 Core Bootstrap Loop (foreach)
# ------------------------------------------------------------------------------
message(paste("    -> Running", n_boot, "iterations on", num_cores, "cores..."))

boot_output <- foreach(
  i           = 1:n_boot,
  .combine    = "cbind",
  .packages   = c(
    "RMTL", "data.table", "rrpack", "MCMCpack", "GIGrvg", "utils", "MASS", 
    "nortest", "MBSP", "foreach", "Matrix", "parallel", "glmnet", "spls", 
    "MXM", "dlnm", "splines", "mgcv", "doParallel", "ranger", 
    "palmerpenguins", "tidyverse", "kableExtra", "haven", "corrplot", "pheatmap"
  ),
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
    nrank       = results_whole_df$Opt_Rank[3], # Combined Rank
    lambda      = results_whole_df$Opt_Lambda[3]      # Combined Lambda
  ))
}

stopCluster(cl)
close(pb)

# 10.6 Advanced Statistical Extraction
# ------------------------------------------------------------------------------
message("    -> Calculating stability metrics...")

# Prepare grid for mapping results
grid_df       <- expand.grid(exposures = colnames(X_boot), 
                             outcomes  = colnames(Y_boot), 
                             bootstraps = 1:n_boot)
grid_df$value <- as.vector(boot_output)

# Matrix representation for matrixStats performance
val_matrix <- matrix(grid_df$value, ncol = n_boot)
count_pos  <- matrixStats::rowSums2(val_matrix > 0)
count_neg  <- matrixStats::rowSums2(val_matrix < 0)
z_counts   <- cbind(count_pos, count_neg)

# Build Results Dataframe
res_comb_ <- cbind(
  matrixStats::rowMeans2(val_matrix),                    # 1. Mean Beta
  matrixStats::rowSds(val_matrix),                       # 2. SD
  matrixStats::rowMins(z_counts) / matrixStats::rowMaxs(z_counts), # 3. pval_ratio
  
  # Corrected P-value (Empirical density crossing zero)
  pmin(2 * matrixStats::rowMins(cbind(
    matrixStats::rowMeans2(val_matrix <= 0),
    matrixStats::rowMeans2(val_matrix >= 0)
  )), 1),                                                # 4. pval_corrected
  
  matrixStats::rowQuantiles(val_matrix, probs = c(0.025, 0.975)), # 5. 95% CI
  matrixStats::rowMeans2(val_matrix != 0)                # 6. Selection Probability
  
  
)

colnames(res_comb_) <- c("mean", "sd", "pval_ratio", "pval_corrected", 
                         "lower_bound", "upper_bound", "sel_prob")

res_comb_          <- data.frame(res_comb_)
res_comb_$exposure <- rep(colnames(X_boot), times = ncol(Y_boot))
res_comb_$outcome  <- rep(colnames(Y_boot), each = ncol(X_boot))

# --- ADDITIONAL STATISTICS (Theoretical/Parametric) ---
#-----------Compute additional statistics using aggregated statistics
original_betas <- fit_comb_whole$B
colnames(original_betas) <- colnames(Y_boot)
rownames(original_betas) <- colnames(X_boot)
original_betas <- as.data.frame(original_betas)
original_betas$Predictor <- rownames(original_betas)
original_betas <- melt(original_betas, id.vars = "Predictor", variable.name = "Outcome", value.name = "Value")
original_betas$Predictor == res_comb_$exposure
original_betas$Outcome == res_comb_$outcome
tail <- qnorm(1 - .05 / 2)
res_comb_$bootstrap_ratio <- original_betas$Value / res_comb_$sd
res_comb_$th_pval <- 2 * pnorm(abs(res_comb_$bootstrap_ratio), lower.tail = FALSE)
res_comb_$th_lower_bound <- original_betas$Value - res_comb_$sd * tail
res_comb_$th_upper_bound <- original_betas$Value + res_comb_$sd * tail

View(res_comb_)
# 10.7 Save Scenario Results
# ------------------------------------------------------------------------------
write.csv2(res_comb_, "results/BOOTSTRAP_STATISTICS_COMB_wholesample.csv", row.names = FALSE)





# ==============================================================================
# BLOCK 13: UNIVERSAL STABLE HEATMAP GENERATOR
# ==============================================================================
# DESCRIPTION: Generates a cleaned heatmap by masking non-significant 
#              associations based on chosen bootstrap metrics.
# ==============================================================================

# 13.1 Define the Universal Plotting Function
# ------------------------------------------------------------------------------
plot_stable_signature <- function(fit_whole, boot_res, X_mat, 
                                  scenario_name, 
                                  filter_by = "pval_ratio", 
                                  alpha = 0.05) {
  
  message(paste0(">>> Generating Stable Heatmap for: ", scenario_name))
  message(paste0("    -> Filter: ", filter_by, " | Threshold: < ", alpha))
  
  # A. EXTRACT ORIGINAL BETAS
  # ----------------------------------------------------------------------------
  # Get the B matrix from the whole sample fit
  B_orig <- fit_whole$B
  rownames(B_orig) <- colnames(X_mat)
  colnames(B_orig) <- colnames(Y_whole)
  
  # B. PREPARE THE STABILITY MASK
  # ----------------------------------------------------------------------------
  # Pivot the chosen bootstrap metric into a matrix matching B
  stability_mask <- boot_res %>%
    select(exposure, outcome, !!sym(filter_by)) %>%
    tidyr::pivot_wider(names_from = outcome, values_from = !!sym(filter_by)) %>%
    tibble::column_to_rownames("exposure")
  
  # Align rows and columns with B matrix
  stability_mask <- stability_mask[rownames(B_orig), colnames(B_orig)]
  
  # C. APPLY THE FILTER (MASKING)
  # ----------------------------------------------------------------------------
  # Create a copy for filtering
  B_stable <- B_orig
  
  # Logic: If metric >= alpha, it is NOT significant -> Set to 0
  B_stable[as.matrix(stability_mask) >= alpha] <- 0
  
  # D. POST-FILTER CLEANUP
  # ----------------------------------------------------------------------------
  # Remove rows that are entirely zero after masking (cleaner visualization)
  active_rows <- rowSums(B_stable != 0) > 0
  if(sum(active_rows) == 0) {
    message("    ! WARNING: No associations survived the stability filter.")
    return(NULL)
  }
  coef_subset <- B_stable[active_rows, , drop = FALSE]
  
  # E. ROW ANNOTATIONS (Families)
  # ----------------------------------------------------------------------------
  clean_names <- gsub("_None|_Ter_2|_Ter_3", "", rownames(coef_subset))
  row_info <- vars_info %>% 
    filter(variable_name %in% clean_names) %>%
    distinct(variable_name, .keep_all = TRUE)
  
  df_row_ann <- data.frame(
    Family = row_info$family[match(clean_names, row_info$variable_name)],
    row.names = rownames(coef_subset)
  )
  
  # F. RENDER PHEATMAP
  # ----------------------------------------------------------------------------
  p <- pheatmap(
    coef_subset,
    main              = paste0(scenario_name, " (Stable Signature: ", filter_by, " < ", alpha, ")"),
    color             = custom_colors_pro,    # Defined in Block 8
    breaks            = get_pro_breaks(coef_subset), # Defined in Block 8
    annotation_row    = df_row_ann,
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    border_color      = "white",
    fontsize_row      = 8,
    cellwidth         = 30,
    silent            = TRUE
  )
  
  print(p)
}
    
    
    # Define output PDF path
    output_pdf <- "results/STABLE_EXPOSOME_SIGNATURES_wholesample_afterbootstrapping_preg.pdf"
    
    pdf(output_pdf, width = 12, height = 14)
    
    # Page 1: Pregnancy (Filtered by Sign Consistency)
    grid.newpage()
    plot_stable_signature(fit_preg_whole, res_preg, X_preg_whole, 
                          "Pregnancy", filter_by = "pval_ratio", alpha = 0.05)
    
    # Page 2: Pregnancy (Filtered by Parametric P-value)
    grid.newpage()
    plot_stable_signature(fit_preg_whole, res_preg, X_preg_whole, 
                          "Pregnancy", filter_by = "th_pval", alpha = 0.05)
    
    # Page 3: Pregnancy (Strict Filter: Sign Consistency (corrected) < 0.1)
    grid.newpage()
    plot_stable_signature(fit_preg_whole, res_preg, X_preg_whole, 
                          "Pregnancy", filter_by = "pval_corrected", alpha = 0.1)
    
    dev.off()
    
    message(paste0(">>> Final Report exported to: ", output_pdf))
    
    
    
    
    
    # Define output PDF path
    output_pdf <- "results/STABLE_EXPOSOME_SIGNATURES_wholesample_afterbootstrapping.pdf"
    
    pdf(output_pdf, width = 12, height = 14)
    
    # Page 1: Pregnancy (Filtered by Sign Consistency corrected < 0.1)
    grid.newpage()
    plot_stable_signature(fit_preg_whole, res_preg, X_preg_whole, 
                          "Pregnancy", filter_by = "pval_corrected", alpha = 0.1)
    
    # Page 2: Childhood (Filtered by Sign Consistency corrected < 0.1)
    grid.newpage()
    plot_stable_signature(fit_child_whole, res_child, X_child_whole, 
                          "Childhood", filter_by = "pval_corrected", alpha = 0.1)
    
    # Page 3: Combined (Strict Filter: Sign Consistency corrected < 0.1)
    grid.newpage()
    X_comb_whole <- cbind(X_preg_whole, X_child_whole)
    plot_stable_signature(fit_comb_whole, res_comb_, X_comb_whole, 
                          "Combined", filter_by = "pval_corrected", alpha = 0.1)
    
    dev.off()
