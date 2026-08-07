## =============================================================================
## msRRR tutorial workflow functions
## VERSION: v2 (7 August 2026)
## =============================================================================
##
## Reusable data-preparation, reporting, bootstrap orchestration and LOCO
## helpers used by both the standalone R script and the R Markdown tutorial.
## Statistical model functions remain in msrrr_v4.R; graphical functions
## remain in msrrr_plotting_functions.R.
## =============================================================================

make_progress_callback <- function(progress_bar) {
  force(progress_bar)
  function(n) utils::setTxtProgressBar(progress_bar, n)
}

read_euro_csv <- function(filename) {
  read_delim(filename, delim = ";", 
             locale = locale(decimal_mark = ","), 
             show_col_types = FALSE)
}

process_exposome_matrix <- function(df_raw, var_names, factors_list) {

  # 1. Select specific variables
  df_subset <- df_raw %>% dplyr::select(all_of(var_names))

  # 2. Identify which 'forced_factors' are present
  present_factors <- intersect(factors_list, colnames(df_subset))

  if(length(present_factors) > 0) {
    # Ensure they are factors
    df_subset <- df_subset %>%
      dplyr::mutate(dplyr::across(all_of(present_factors), as.factor))

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

get_prop <- function(
    data_vec,
    levels = get0("all_levels", ifnotfound = NULL, inherits = TRUE)) {
  if (is.null(levels)) {
    levels <- names(table(data_vec, useNA = "always"))
  }
  tab <- table(factor(data_vec, levels = levels, exclude = NULL))
  return(round(prop.table(tab), 3))
}

fill_results <- function(row_idx, model_obj) {
  updated_results <- results_df
  updated_results$Opt_Rank[row_idx] <- model_obj$nrank
  updated_results$Opt_Lambda[row_idx] <- model_obj$lam.opt
  # Count selected predictors (non-zero coefficients in B matrix)
  updated_results$N_Exposures[row_idx] <- count_selected_exposures(
    model_obj, tol = selection_tol
  )
  # Extract the best value for the selected CV criterion.
  criterion_column <- intersect(
    c("CV_Crit_Value", "CV_Crit_value", "CV_Deviance"),
    names(updated_results)
  )[1L]
  if (is.na(criterion_column)) {
    stop("results_df has no recognised CV-criterion column.")
  }
  updated_results[[criterion_column]][row_idx] <- min(
    unlist(model_obj$tunepath.opt), na.rm = TRUE
  )
  results_df <<- updated_results
}

# Return TRUE only when every observed value is coded as 0 or 1. Keeping this
# helper in the workflow file makes the LOCO preprocessing functions
# self-contained instead of depending on a similarly named private helper in
# msrrr_v4.R.
is_binary_01 <- function(x) {
  observed <- unique(stats::na.omit(x))
  length(observed) > 0L &&
    length(observed) <= 2L &&
    all(observed %in% c(0, 1))
}

scale_loco_matrix <- function(M, matrix_name, drop_constant = FALSE) {
  M <- as.matrix(M)
  column_sd <- apply(M, 2L, stats::sd, na.rm = TRUE)
  constant <- !is.finite(column_sd) | column_sd == 0

  if (any(constant)) {
    constant_names <- colnames(M)[constant]
    if (!drop_constant) {
      stop(
        matrix_name, " contains constant column(s) after cohort exclusion: ",
        paste(constant_names, collapse = ", "),
        ". Penalised X columns are not removed automatically."
      )
    }
    message(
      "    -> Removing constant ", matrix_name, " column(s): ",
      paste(constant_names, collapse = ", ")
    )
    M <- M[, !constant, drop = FALSE]
  }

  # Reuse the central msRRR rule after cohort exclusion: z-score continuous
  # columns and preserve binary 0/1 columns.
  scaling_parameters <- .fit_standardization(M)
  scaled <- .apply_standardization(M, scaling_parameters)
  if (any(!is.finite(scaled))) {
    stop(matrix_name, " produced non-finite values during LOCO scaling.")
  }
  scaled
}

scale_loco_outcomes <- function(Y_subset, family, familygroup) {
  Y_subset <- as.matrix(Y_subset)
  # Use the same family-based rule as CV and whole-sample refits. Parameters
  # are deliberately relearned from the cohorts retained in this LOCO fit.
  scaling_parameters <- .fit_y_scaling(Y_subset, family, familygroup)
  .apply_y_scaling(Y_subset, scaling_parameters)
}

run_one_penalized_loco_bootstrap <- function(
  fit_loco,
  X_loco,
  Y_loco,
  Z_loco,
  nrank,
  lambda,
  family,
  familygroup,
  scenario_name,
  omitted_cohort,
  n_boot = n_boot_global,
  show_progress = TRUE) {

  message(
    "    -> Bootstrap: ", scenario_name,
    " without cohort ", omitted_cohort,
    " (", n_boot, " resamples)."
  )

  boot_sampling <- generate_resampling_msrrr(
    blocks = list(X_loco, Z_loco, Y_loco),
    n_boot = n_boot,
    balanced = TRUE,
    keep_all_variables = FALSE,
    verbose = TRUE
  )

  Z_boot <- Z_loco
  sd_null <- boot_sampling$sd_null
  if (!is.null(sd_null)) {
    removed_covariates <- names(sd_null[[2]])
    if (length(removed_covariates) > 0L) {
      message(
        "       Removing Z columns with zero variance in at least one ",
        "bootstrap sample: ",
        paste(removed_covariates, collapse = ", ")
      )
      Z_boot <- Z_boot[
        , !colnames(Z_boot) %in% removed_covariates, drop = FALSE
      ]
    }
  }

  num_cores <- max(1L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(num_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Load sourced project functions explicitly on the fresh PSOCK workers.
  msrrr_worker_file <- normalizePath(
    "functions/msrrr_v4.R", winslash = "/", mustWork = TRUE
  )
  bootstrap_worker_file <- normalizePath(
    "functions/bootstrapping_functions_v_31_07_26.r",
    winslash = "/", mustWork = TRUE
  )
  parallel::clusterCall(
    cl,
    function(msrrr_file, bootstrap_file) {
      sys.source(msrrr_file, envir = .GlobalEnv)
      sys.source(bootstrap_file, envir = .GlobalEnv)
      invisible(TRUE)
    },
    msrrr_file = msrrr_worker_file,
    bootstrap_file = bootstrap_worker_file
  )

  doParallel::registerDoParallel(cl)
  doSNOW::registerDoSNOW(cl)

  opts <- list()
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(max = n_boot, style = 3)
    on.exit(close(pb), add = TRUE)
    progress <- make_progress_callback(pb)
    opts <- list(progress = progress)
  }

  boot_output <- foreach(
    i = seq_len(n_boot),
    .combine = "cbind",
    .packages = c("glmnet", "rrpack"),
    .options.snow = opts
  ) %dopar% {
    unlist(msrrr_bootstrap_k(
      Y = Y_loco,
      X = X_loco,
      Z = Z_boot,
      inds = boot_sampling$full_idx[i],
      family = family,
      familygroup = familygroup,
      nrank = nrank,
      lambda = lambda,
      init = NULL
    ))
  }

  original_B <- as.matrix(fit_loco$B)
  rownames(original_B) <- colnames(X_loco)
  colnames(original_B) <- colnames(Y_loco)

  summarize_msrrr_bootstrap(
    boot_output = boot_output,
    original_B = original_B,
    lambda = lambda,
    selection_tol = selection_tol,
    sel_prob_threshold = sel_prob_threshold,
    sign_consistency_threshold = sign_consistency_threshold
  )
}

summarize_penalized_loco <- function(
  scenario_name,
  scenario,
  loco_fits,
  cohort_levels,
  tol = selection_tol) {

  original_B <- as.matrix(scenario$original_fit$B)
  rownames(original_B) <- colnames(scenario$X_source)
  colnames(original_B) <- colnames(Y_whole)

  B_array <- array(
    NA_real_,
    dim = c(nrow(original_B), ncol(original_B), length(cohort_levels)),
    dimnames = list(
      rownames(original_B), colnames(original_B), cohort_levels
    )
  )
  for (cohort_index in seq_along(cohort_levels)) {
    omitted_cohort <- cohort_levels[cohort_index]
    B_array[, , cohort_index] <-
      loco_fits[[omitted_cohort]][[scenario_name]]$B
  }

  selected_array <- abs(B_array) > tol
  n_selected <- apply(selected_array, c(1L, 2L), sum)
  n_positive <- apply(B_array > tol, c(1L, 2L), sum)
  n_negative <- apply(B_array < -tol, c(1L, 2L), sum)
  dominant_count <- pmax(n_positive, n_negative)
  dominant_direction <- ifelse(
    n_selected == 0L, "Not selected",
    ifelse(
      n_positive > n_negative, "Positive",
      ifelse(n_negative > n_positive, "Negative", "Tie")
    )
  )
  dominant_sign_consistency <- ifelse(
    n_selected > 0L, dominant_count / n_selected, NA_real_
  )

  original_sign <- sign(original_B)
  same_original_sign <- array(FALSE, dim = dim(B_array))
  for (cohort_index in seq_along(cohort_levels)) {
    same_original_sign[, , cohort_index] <- selected_array[
      , , cohort_index
    ] & sign(B_array[, , cohort_index]) == original_sign
  }
  n_same_original_sign <- apply(
    same_original_sign, c(1L, 2L), sum
  )
  original_sign_consistency <- ifelse(
    n_selected > 0L & abs(original_B) > tol,
    n_same_original_sign / n_selected,
    NA_real_
  )

  summary_data <- data.frame(
    scenario = scenario_name,
    exposure = rep(rownames(original_B), times = ncol(original_B)),
    outcome = rep(colnames(original_B), each = nrow(original_B)),
    original_beta = as.vector(original_B),
    n_loco_fits = length(cohort_levels),
    n_selected = as.vector(n_selected),
    selection_frequency = as.vector(n_selected) / length(cohort_levels),
    n_positive = as.vector(n_positive),
    n_negative = as.vector(n_negative),
    dominant_direction = as.vector(dominant_direction),
    dominant_sign_consistency = as.vector(dominant_sign_consistency),
    original_sign_consistency = as.vector(original_sign_consistency),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_data,
    B_array = B_array,
    n_selected = n_selected,
    n_positive = n_positive,
    n_negative = n_negative
  )
}

is_unpenalized_refit_available <- function(object) {
  is.list(object) &&
    !identical(object$available, FALSE) &&
    !is.null(object$fit) &&
    !is.null(object$fit$B) &&
    !is.null(object$X_selected) &&
    ncol(object$X_selected) > 0L
}


empty_unpenalized_inference <- function(reason = NA_character_) {
  out <- data.frame(
    exposure = character(0),
    outcome = character(0),
    lambda = numeric(0),
    beta_unpenalized = numeric(0),
    bootstrap_mean = numeric(0),
    bootstrap_standard_error = numeric(0),
    ci95_lower = numeric(0),
    ci95_upper = numeric(0),
    bootstrap_q025 = numeric(0),
    bootstrap_q975 = numeric(0),
    z_value = numeric(0),
    p_value = numeric(0),
    significant = logical(0),
    n_valid_bootstrap = integer(0),
    stringsAsFactors = FALSE
  )
  attr(out, "status") <- "not_applicable_no_selected_exposures"
  attr(out, "reason") <- reason
  out
}


unpenalized_refit_status <- function(object, scenario_name) {
  available <- is_unpenalized_refit_available(object)
  data.frame(
    scenario = scenario_name,
    approach_b_available = available,
    n_selected_exposures = if (available) ncol(object$X_selected) else 0L,
    status = if (!is.null(object$status)) object$status else if (available) "fitted" else "not_applicable",
    reason = if (!is.null(object$reason)) object$reason else NA_character_,
    stringsAsFactors = FALSE
  )
}


run_unpenalized_bootstrap <- function(
    unpenalized_object,
    Y,
    Z,
    family,
    familygroup,
    n_boot,
    scenario_name,
    show_progress = TRUE) {

  if (!is_unpenalized_refit_available(unpenalized_object)) {
    reason <- if (!is.null(unpenalized_object$reason)) {
      unpenalized_object$reason
    } else {
      "No exposures were selected by the penalised whole-sample model."
    }
    message(
      ">>> [Approach B] Bootstrap skipped for ", scenario_name, ": ", reason
    )
    return(empty_unpenalized_inference(reason))
  }

  X_selected <- unpenalized_object$X_selected
  Z_boot <- Z

  message(
    ">>> [Approach B] Bootstrap for ", scenario_name,
    ": using ", ncol(X_selected),
    " candidate exposure(s) retained from the penalised whole-sample model; ",
    "each bootstrap sample is refitted with rank = ",
    unpenalized_object$nrank, ", lambda = 0 and init = NULL."
  )

  boot_sampling <- generate_resampling_msrrr(
    blocks = list(X_selected, Z_boot, Y),
    n_boot = n_boot,
    balanced = TRUE,
    keep_all_variables = FALSE,
    verbose = TRUE
  )

  sd_null <- boot_sampling$sd_null
  if (!is.null(sd_null)) {
    removed_covariates <- names(sd_null[[2]])
    message(
      "    -> Removing covariates with zero variance in some bootstrap ",
      "samples: ", paste(removed_covariates, collapse = " - ")
    )
    Z_boot <- Z_boot[
      , !colnames(Z_boot) %in% removed_covariates, drop = FALSE
    ]
  }

  num_cores <- max(1L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(num_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # PSOCK workers are fresh R sessions. Explicitly source the complete msRRR
  # function chain on every worker; foreach auto-export is not sufficient here
  # because msrrr_bootstrap_k() calls several indirectly referenced helpers.
  msrrr_worker_file <- normalizePath(
    "functions/msrrr_v4.R", winslash = "/", mustWork = TRUE
  )
  bootstrap_worker_file <- normalizePath(
    "functions/bootstrapping_functions_v_31_07_26.r",
    winslash = "/", mustWork = TRUE
  )
  parallel::clusterCall(
    cl,
    function(msrrr_file, bootstrap_file) {
      sys.source(msrrr_file, envir = .GlobalEnv)
      sys.source(bootstrap_file, envir = .GlobalEnv)
      invisible(TRUE)
    },
    msrrr_file = msrrr_worker_file,
    bootstrap_file = bootstrap_worker_file
  )

  doParallel::registerDoParallel(cl)
  doSNOW::registerDoSNOW(cl)

  opts <- list()
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(max = n_boot, style = 3)
    on.exit(close(pb), add = TRUE)
    progress <- make_progress_callback(pb)
    opts <- list(progress = progress)
  }

  boot_output <- foreach(
    i = seq_len(n_boot),
    .combine = "cbind",
    .packages = c("glmnet", "rrpack"),
    .options.snow = opts
  ) %dopar% {
    unlist(msrrr_bootstrap_k(
      Y = Y,
      X = X_selected,
      Z = Z_boot,
      inds = boot_sampling$full_idx[i],
      family = family,
      familygroup = familygroup,
      nrank = unpenalized_object$nrank,
      lambda = 0,
      init = NULL
    ))
  }

  original_B <- unpenalized_object$fit$B
  rownames(original_B) <- colnames(X_selected)
  colnames(original_B) <- colnames(Y)

  summarize_msrrr_unpenalized_bootstrap(
    boot_output = boot_output,
    original_B = original_B,
    lambda = 0
  )
}
