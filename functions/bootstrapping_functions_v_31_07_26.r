#' Generate `n_boot` bootstrap samples.
#' @details
#' The bootstrap samples are generated so that it is very unlikely that they
#' will be associated with a 0 variance variable. However such situation can
#' appear in very rare cases.
#'
#' The simplest procedure to avoid such situation would be to compute the
#' variance of each variable in each bootstrap sample. However, in
#' high-dimensional cases, this is time consuming. Thus, the function starts
#' by identifying variables that present a risk of being of null variance in
#' at least one bootstrap sample. A variable is detected as such if the
#' probability of sampling `N` times (where `N` is the number of observations)
#' its most frequent observed value is higher than a specific threshold
#' named `pval` (by default set to `1e-15`) corrected by the number of bootstrap
#' sample `n_boot` and the number of variables in the whole data-set
#' (\eqn{\sum_{j=1}^Jp_{j}}). In the end, a variable is defined as risky if the
#' proportion of its most frequent observed value is strictly higher than:
#' \deqn{risky\textunderscore threshold = \left(\frac{pval}{n_{boot}*
#' \left(\sum_{j=1}^J p_{j}\right)}\right)^{1/N}.}
#' This value is necessarily lower than `1` as `pval` is lower than `1`.
#' However, it could be lower than \eqn{1/N}, which means that all variables
#' are defined as risky. That is why the maximum value between
#' \eqn{risky\textunderscore threshold} and \eqn{1/N} is taken.
#'
#' Once risky variables have been identified, the function computes the
#' variance of each risky variable in each bootstrap sample. If there isn't any
#' 0 variance risky variable, the bootstrap samples are returned as is.
#' Otherwise, three cases are possible:
#' \itemize{
#' \item `keep_all_variables = F`: problematic covariates in `Z` are marked for
#' global removal from every bootstrap fit; problematic exposures in `X` are
#' retained with a warning; and problematic outcomes in `Y` stop execution.
#' \item `keep_all_variables = T`, two possibilities :
#' \itemize{
#' \item If `balanced = T`, the procedure is repeated at most `5` times until all
#' bootstrap samples do not present any 0 variance variable.
#' \item If `balanced = F`, an heuristic procedure is used to modify the sampling
#' probability of each observation in order to keep all variables.
#' }
#' }
#'
#' A short detail of this lastly mentioned heuristic consist in replacing each
#' observed value of the risky variables by \eqn{1-\frac{N_k}{N}}
#' (where \eqn{N_k} corresponds to the number of times this value is observed),
#' normalized so that the sum through all the observations (for each variable)
#' equals `1`. Then, if \eqn{prob_i} defines the sampling probability for the
#' observation associated with line \eqn{i}, it is defined as the maximum value
#' of the previous matrix through all risky variables (again normalized so that
#' \eqn{\sum_{i=1}^N prob_i = 1}). Thus the sampling probability of an
#' observation is more of less associated with `1 - its proportion in the
#' variable where this observation is in the lowest frequent group
#' (through all risky variables)`.
#' @inheritParams rgcca_bootstrap
#' @param pval For all the variables, a threshold for the proportion of the most
#' frequent value of this variable is computed. This threshold is evaluated so
#' that the probability to sample only this value is below `pval`. This
#' probability is corrected by the number of bootstrap samples and the number
#' of variables (default is `1e-15`).
#' @return \item{full_idx}{A list of size n_boot containing the observations
#' kept for each bootstrap sample.}
#' @return \item{sd_null}{When `keep_all_variables = FALSE`, a block-aligned
#' list containing only the `Z` covariates to remove globally because they lose
#' observed variation in at least one bootstrap sample. It is `NULL` when no
#' covariate needs removal. Exposures are never placed in `sd_null`, and an
#' affected outcome stops execution.}
#' @return \item{sd_null_detected}{A block-aligned diagnostic list retaining
#' every `X`, `Z`, and `Y` variable found to lose observed variation in at
#' least one bootstrap sample before the block-specific policy is applied.}
#' @title Generate bootstrap samples.
#' @noRd

## VERSION: v1 (31 July 2026)
##
## Changes:
## 1) Missing values are excluded when checking whether a variable has
##    observed variation. A bootstrap sample with no observed value, or with
##    only one unique observed value, is treated as having null observed
##    variance. Missingness is allowed.
## 2) Variables that are already constant in the original data are identified
##    before resampling: constant exposures and outcomes stop execution;
##    constant covariates produce a warning. A named intercept column is exempt.
## 3) Documents the one-element-list format expected by msrrr_bootstrap_k()
##    for bootstrap indices.
## 4) Replaces the unavailable RGCCA-specific error helper with a local,
##    self-contained error function based on base R stop().
## 5) Defines block-specific handling of variables that lose observed
##    variation in at least one bootstrap sample when keep_all_variables =
##    FALSE: retain X with a warning, mark Z for global removal with a warning,
##    and stop for Y. sd_null contains only Z variables to be removed by the
##    calling pipeline; sd_null_detected retains the complete X/Z/Y diagnostic.
## 6) Adds summarize_msrrr_bootstrap(), which reports the original penalised
##    coefficient, the bootstrap mean including zero coefficients, and the
##    conditional mean/median among selected bootstrap fits. It also reports
##    selection and sign stability using the common 1e-8 selection tolerance.
##    Classical p-values and confidence intervals are deliberately not
##    calculated for penalised fits.
## 7) Makes init explicit in msrrr_bootstrap_k(). It remains NULL by default,
##    so every bootstrap sample calculates its own ridge initialisation, but
##    callers can supply another initialisation when deliberately required.
## 8) Adds summarize_msrrr_unpenalized_bootstrap() as a separate conditional
##    inference summary for fixed-subset, lambda = 0 fits. It refuses penalised
##    input and reports bootstrap standard errors, normal-approximation 95%
##    intervals and p-values, plus percentile intervals. No
##    multiple-testing correction is applied automatically.
## 9) Retains a complete set of descriptive coefficient-distribution summaries:
##    the penalised summary reports mean, median, standard deviation and
##    quantiles across all fits and, separately, among non-zero selected fits;
##    the unpenalised summary reports mean, median, standard deviation (as the
##    bootstrap standard error) and percentile quantiles.

## Raise a bootstrap error without depending on RGCCA internals.
## exit_code is accepted for compatibility with the inherited checks.
stop_msrrr_bootstrap <- function(..., exit_code = NULL) {
  stop(paste0(...), call. = FALSE)
}

## Return the columns with at most one unique non-missing value.
.constant_observed_columns <- function(block) {
  if (is.null(block) || NCOL(block) == 0L) return(integer(0))

  which(vapply(
    seq_len(NCOL(block)),
    function(j) {
      observed <- block[!is.na(block[, j]), j]
      length(unique(observed)) <= 1L
    },
    logical(1)
  ))
}


## Use column names in diagnostics, with stable fallbacks when names are absent.
.bootstrap_variable_names <- function(block, idx, prefix) {
  variable_names <- colnames(block)
  if (is.null(variable_names)) {
    variable_names <- paste0(prefix, seq_len(NCOL(block)))
  }
  variable_names[idx]
}


generate_resampling_msrrr <- function(blocks, n_boot, balanced = TRUE,
                                keep_all_variables = FALSE, pval = 1e-15,
                                verbose = TRUE) {
                      
  if (verbose) {
    packageStartupMessage("Bootstrap samples sanity check...",
      appendLF = FALSE
    )
  }

  # Initialization
  pval <- min(pval, 1)
  NO_null_sd_var <- FALSE
  iter <- 0
  sd_null_detected <- NULL
  raw_blocks <- blocks
  N <- NROW(raw_blocks[[1]])
  prob <- rep(1 / N, N)

  # Check variables that have no observed variation in the original data.
  # blocks are expected in msRRR order: X (exposures), Z (covariates), Y
  # (outcomes). Missing values are allowed and excluded from the check.
  constant_original <- lapply(raw_blocks, .constant_observed_columns)
  constant_covariate_idx <- integer(0)

  if (length(constant_original) >= 1L &&
      length(constant_original[[1]]) > 0L) {
    constant_exposures <- .bootstrap_variable_names(
      raw_blocks[[1]], constant_original[[1]], "X"
    )
    stop(
      "Constant exposure variable(s) detected among observed values: ",
      paste(constant_exposures, collapse = " - "),
      ". Bootstrap fitting cannot estimate effects for constant exposures.",
      call. = FALSE
    )
  }

  if (length(constant_original) >= 3L &&
      length(constant_original[[3]]) > 0L) {
    constant_outcomes <- .bootstrap_variable_names(
      raw_blocks[[3]], constant_original[[3]], "Y"
    )
    stop(
      "Constant outcome variable(s) detected among observed values: ",
      paste(constant_outcomes, collapse = " - "),
      ". Bootstrap fitting requires observed variation in every outcome.",
      call. = FALSE
    )
  }

  if (length(constant_original) >= 2L &&
      length(constant_original[[2]]) > 0L) {
    constant_covariate_idx <- constant_original[[2]]
    covariate_names <- .bootstrap_variable_names(
      raw_blocks[[2]], constant_covariate_idx, "Z"
    )

    # A deliberately constant intercept is part of the model, not a
    # problematic covariate. Exempt it when it is clearly named as such.
    is_named_intercept <- grepl(
      "^\\(?intercept\\)?$",
      covariate_names,
      ignore.case = TRUE
    )
    constant_covariates <- covariate_names[!is_named_intercept]

    if (length(constant_covariates) > 0L) {
      warning(
        "Constant covariate variable(s) detected among observed values: ",
        paste(constant_covariates, collapse = " - "),
        ". These covariates contain no information for coefficient ",
        "estimation and may cause instability in bootstrap fitting.",
        call. = FALSE
      )
    }
  }

  # For any variable, threshold for the proportion of the most frequent value
  # of a variable. This threshold is computed so that the probability to sample
  # only this value is below `pval`. This probability is corrected by the number
  # of bootstrap samples and the number of variables.
  risky_threshold <- max(
    1 / N,
    (pval / (n_boot * sum(vapply(raw_blocks, NCOL, FUN.VALUE = 1L))))^(1 / N)
  )

  # Identify variables with value having an observed proportion higher than
  # risky_threshold.
  risky_var <- lapply(
    raw_blocks,
    function(block) {
      which(apply(
        block, 2,
        function(x) {
          observed <- x[!is.na(x)]
          observed_concentration <- if (length(observed) == 0L) {
            1
          } else {
            max(table(observed) / length(observed))
          }
          missing_concentration <- mean(is.na(x))

          # Missing values are not treated as a response category. The
          # separate missingness term only detects a material probability
          # that a bootstrap sample contains no observed value.
          observed_concentration > risky_threshold ||
            missing_concentration > risky_threshold
        }
      ))
    }
  )

  # Constant covariates have already produced their dedicated warning.
  # Do not subsequently convert that warning into a zero-variance resampling
  # error. They remain in the supplied Z block; the calling pipeline decides
  # whether to retain or remove them.
  if (length(risky_var) >= 2L && length(constant_covariate_idx) > 0L) {
    risky_var[[2]] <- setdiff(
      risky_var[[2]], constant_covariate_idx
    )
  }

  # Keep only risky variables for each block.
  raw_blocks_filtered <- Map(
    function(x, y) x[, y, drop = FALSE],
    raw_blocks, risky_var
  )
  # While there are variables with null variance among the risky variables.
  while (!NO_null_sd_var) {
    if (balanced) { # Balanced bootstrap sampling.
      full_idx <- rep(x = seq(N), each = n_boot)
      full_idx <- sample(x = full_idx, size = N * n_boot, replace = FALSE)
      full_idx <- split(full_idx, ceiling(seq_along(full_idx) / N))
    } else { # Unbalanced bootstrap sampling.
      full_idx <- lapply(
        seq(n_boot),
        function(x) {
          sample(x = seq(N), replace = TRUE, prob = prob)
        }
      )
    }
    # Compute blocks for each bootstrap sample (only with risky variables).
    boot_blocks_filtered <-
      lapply(
        full_idx,
        function(idx) {
          lapply(
            raw_blocks_filtered,
            function(x) {
              y <- x[idx, , drop = FALSE]
              rownames(y) <- paste("S", seq_along(idx))
              return(y)
            }
          )
        }
      )

    # For each sample, identify variables with a single unique value.

    boot_column_sd_null <-
      lapply(
        boot_blocks_filtered,
        function(boot) {
          lapply(boot, function(boot_bl) {
            which(apply(
              boot_bl, 2,
              function(x) {
                observed <- x[!is.na(x)]
                length(unique(observed)) <= 1L
              }
            ))
          })
        }
      )

    # Summarize through all the samples.
    eval_boot_sample <- vapply(
      boot_column_sd_null,
      function(x) sum(vapply(x, length, FUN.VALUE = 1L)),
      FUN.VALUE = 1L
    )
    NO_null_sd_var <- (sum(eval_boot_sample) == 0)

    if (NO_null_sd_var) {
      # through all samples, all variables have non null variances.
      sd_null <- NULL
    } else {
      # If at least one sample have been identified with a null
      # variance variable.
      if (!keep_all_variables) {
        # Extract every variable that loses observed variation in at least one
        # bootstrap sample. Keep the complete result for diagnostics.
        sd_null_detected <- Reduce("rbind", boot_column_sd_null)
        rownames(sd_null_detected) <- NULL
        sd_null_detected <- apply(
          sd_null_detected, 2,
          function(x) unique(names(Reduce("c", x)))
        )
        sd_null_detected <- Map(function(x, y) {
          z <- match(x, y)
          names(z) <- x
          return(z)
        }, sd_null_detected, lapply(raw_blocks, colnames))

        null_X <- if (length(sd_null_detected) >= 1L) {
          names(sd_null_detected[[1]])
        } else {
          character(0)
        }
        null_Z <- if (length(sd_null_detected) >= 2L) {
          names(sd_null_detected[[2]])
        } else {
          character(0)
        }
        null_Y <- if (length(sd_null_detected) >= 3L) {
          names(sd_null_detected[[3]])
        } else {
          character(0)
        }

        # Outcomes must remain estimable and identical in every bootstrap fit.
        if (length(null_Y) > 0L) {
          stop_msrrr_bootstrap(paste(
            "The following outcomes have no observed variation, or no observed",
            " values, in at least one bootstrap sample: ",
            paste(null_Y, collapse = " - "),
            ". Bootstrap fitting cannot continue because the same outcomes",
            " must be estimable in every replicate."
          ))
        }

        # Exposures are retained to preserve the coefficient structure.
        if (length(null_X) > 0L) {
          warning(
            paste(
              "The following exposures have no observed variation in at least",
              " one bootstrap sample: ",
              paste(null_X, collapse = " - "),
              ". They are retained because keep_all_variables = FALSE only",
              " removes problematic covariates. Bootstrap coefficients for",
              " these exposures may be unstable or unavailable in affected",
              " replicates."
            ),
            call. = FALSE
          )
        }

        # Only problematic Z variables are returned in sd_null for global
        # removal by the calling pipeline. If Z is unaffected, return NULL so
        # that a pipeline does not accidentally subset Z with integer(0).
        if (length(null_Z) > 0L) {
          warning(
            paste(
              "The following covariates have no observed variation in at least",
              " one bootstrap sample and must be excluded from all bootstrap",
              " fits: ",
              paste(null_Z, collapse = " - "),
              ". Bootstrap models will therefore use a reduced covariate",
              " adjustment set."
            ),
            call. = FALSE
          )
          sd_null <- vector("list", length(raw_blocks))
          sd_null[[2]] <- sd_null_detected[[2]]
        } else {
          sd_null <- NULL
        }

        NO_null_sd_var <- TRUE
      } else { # It is NOT allowed to remove variables.
        # Generate at most five different re-sampling until not a single
        # variable has a null variance.
        if (iter > 5) { # Otherwise STOP.
          # Extract the troublesome variables.
          sd_null <- Reduce("rbind", boot_column_sd_null)
          rownames(sd_null) <- NULL
          sd_null <- apply(
            sd_null, 2,
            function(x) unique(names(Reduce("c", x)))
          )
          sd_null <- Map(function(x, y) {
            z <- match(x, y)
            names(z) <- x
            return(z)
          }, sd_null, lapply(raw_blocks, colnames))

          error_message <- paste(
            "Impossible to define all bootstrap samples",
            "without variables with null variance. Please",
            "consider removing these variables: ",
            paste(names(Reduce("c", sd_null)), collapse = " - ")
          )
          # In the balanced case, you CANNOT play with the sampling probability
          # of the different observations as it is unbalanced.
          if (balanced) {
            error_message <- paste0(
              error_message,
              ". Please, consider unbalanced bootstrap by",
              " setting 'balanced' to FALSE."
            )
          }
          stop_msrrr_bootstrap(error_message)
        }
        # In the unbalanced case, you CAN play with the sampling probability
        # of the different observations.
        if (!balanced) {
          if (iter == 0) { # The first time, you define your unbalancedness.
            # Each observed value of the risky variables is replaced by
            # `1 - the proportion of this observed value`, normalized so that
            # the sum through all the observations (for each variable)
            # equals `1`.
            prob <- lapply(
              raw_blocks_filtered,
              function(block) {
                apply(block, 2, function(var) {
                  occurences <- table(var, useNA = "ifany") / length(var)
                  new_idx <- match(as.character(var), names(occurences))
                  new_var <- as.matrix(occurences[new_idx])
                  new_var <- (1 - new_var) / sum(1 - new_var)
                  return(new_var)
                })
              }
            )
            # The sampling probability for each observation is associated with
            # the maximum value of the previous matrix through all risky
            # variables (again normalize so that `sum(prob) = 1`). Thus
            # the sampling probability of an observation is more of less
            # associated with `1 - its proportion in the variable where this
            # observation is in the lowest frequent group (through all
            # risky variables)`.
            prob <- apply(Reduce("cbind", prob), 1, max) / sum(
              apply(Reduce("cbind", prob), 1, max)
            )
          }
        }
        iter <- iter + 1
      }
    }
  }
  if (verbose) packageStartupMessage("OK")
  return(list(
    full_idx = full_idx,
    sd_null = sd_null,
    sd_null_detected = sd_null_detected
  ))
}


#' Parallel lapply with progress bar
#' @inheritParams msrrr.fit
#' @param X a vector (atomic or list) or an expression object.
#' @param FUN the function to be applied to each element of X.
#' @param ... optional arguments to FUN.
#' @return The result of lapply(X, FUN, ...)
#' @noRd
par_pblapply <- function(X, FUN, ..., n_cores = 1, verbose = TRUE) {

  require(parallel)
  check_integer("n_cores", n_cores, min = 0)

  verbose <- verbose & interactive()

  if (!verbose) {
    pbapply::pboptions(type = "none")
  } else {
    pbapply::pboptions(type = "timer")
  }

  is_windows <- Sys.info()["sysname"] == "Windows"

  if (n_cores <= 1) {
    cl <- NULL
  } else if (is_windows) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, NULL, envir = environment())
  } else {
    cl <- parallel::makeCluster( # create the cluster
  n_cores, 
  type = "PSOCK"
)
clusterEvalQ(cl, {
  library(rrpack)  # Load required package
  source("./SCRIPTs/bootstrapping_functions_v_31_07_26.r")  # Load custom functions from script
})

# Export custom functions to each worker
clusterExport(cl, c("msrrr_bootstrap_k"))
  }

  W <- pbapply::pblapply(X, FUN, ..., cl = cl)

  if (is_windows && !is.null(cl)) parallel::stopCluster(cl)

  return(W)
}


#' Internal function for computing bootstrap of msRRR.
#' @blocks Y, X and Z datasets
#' @param inds A vector of integers defining the index of the observations
#' taken into account for this bootstrap sample.
#' @param type A character indicating the type of the second object to return.
#' @return \item{W}{A list of RGCCA bootstrap weights. Returned only if there
#' are no missing variables, otherwise the name of the missing variables are
#' returned.}
#' @return \item{L}{If type == "loadings", a list of RGCCA bootstrap loadings.
#' Returned only if there are no missing variables, otherwise the name of the
#' missing variables are returned.
#' If type == "AVE", the AVE of the fitted RGCCA model.}
#' @title Compute bootstrap (internal).
#' @noRd

msrrr_bootstrap_k <- function(Y,X,Z, inds = NULL,family=NULL, familygroup=NULL,
      nrank=NULL, lambda=NULL, init=NULL) {

# The calling pipeline must pass `inds` as a one-element list containing the
# complete bootstrap index vector, for example: boot_sampling$full_idx[i].
# Passing boot_sampling$full_idx[[i]] directly would not match this interface.
indsidx <- as.numeric(inds[[1]])

Z_bootstrap <- if (is.null(Z)) {
  NULL
} else {
  Z[indsidx, , drop = FALSE]
}

msrrr_res_boot <- msrrr.fit(
  Y = Y[indsidx, , drop = FALSE],
  X = X[indsidx, , drop = FALSE],
  Z = Z_bootstrap,
  family = family,
  familygroup = familygroup,
  nrank = nrank,
  lambda = lambda,
  init = init
)

colnames(msrrr_res_boot$B) <- colnames(Y)
rownames(msrrr_res_boot$B) <- colnames(X)

list_msrrr_res_boot <- lapply(names(as.data.frame(msrrr_res_boot$B)), function(col_name) {
  matrix(as.data.frame(msrrr_res_boot$B)[,col_name], ncol = 1, dimnames =list(rownames(msrrr_res_boot$B), col_name))
})

return(list_msrrr_res_boot)

}


## Summarise bootstrap coefficients as penalised coefficient and stability
## measures, rather than as classical post-selection inference.
summarize_msrrr_bootstrap <- function(
    boot_output,
    original_B,
    exposure_names = rownames(original_B),
    outcome_names = colnames(original_B),
    lambda,
    selection_tol = 1e-8,
    sel_prob_threshold = 0.90,
    sign_consistency_threshold = 0.80) {

  boot_output <- as.matrix(boot_output)
  original_B <- as.matrix(original_B)

  if (is.null(exposure_names) || length(exposure_names) != nrow(original_B)) {
    stop("exposure_names must contain one name per row of original_B.",
         call. = FALSE)
  }
  if (is.null(outcome_names) || length(outcome_names) != ncol(original_B)) {
    stop("outcome_names must contain one name per column of original_B.",
         call. = FALSE)
  }
  if (nrow(boot_output) != length(original_B)) {
    stop(
      "boot_output must have nrow(original_B) * ncol(original_B) rows.",
      call. = FALSE
    )
  }
  if (length(lambda) != 1L || is.na(lambda) || !is.finite(lambda) ||
      lambda < 0) {
    stop("lambda must be one finite non-negative value.", call. = FALSE)
  }
  if (length(selection_tol) != 1L || !is.finite(selection_tol) ||
      selection_tol < 0) {
    stop("selection_tol must be one finite non-negative value.",
         call. = FALSE)
  }

  selected <- abs(boot_output) > selection_tol
  positive <- boot_output > selection_tol
  negative <- boot_output < -selection_tol
  n_boot <- ncol(boot_output)
  n_selected <- rowSums(selected, na.rm = TRUE)

  selected_values <- boot_output
  selected_values[!selected] <- NA_real_

  row_quantile <- function(x, prob) {
    if (all(is.na(x))) return(NA_real_)
    as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
  }
  row_median <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    stats::median(x, na.rm = TRUE)
  }

  beta_original <- as.vector(original_B)
  original_sign <- ifelse(
    abs(beta_original) > selection_tol,
    sign(beta_original),
    0
  )
  same_sign <- (sign(boot_output) == original_sign) & selected
  sign_consistency <- ifelse(
    n_selected > 0L & original_sign != 0,
    rowSums(same_sign, na.rm = TRUE) / n_selected,
    NA_real_
  )

  sel_prob <- n_selected / n_boot
  passes_sel_prob <- sel_prob >= sel_prob_threshold
  passes_sign_consistency <- !is.na(sign_consistency) &
    sign_consistency >= sign_consistency_threshold

  out <- data.frame(
    exposure = rep(exposure_names, times = length(outcome_names)),
    outcome = rep(outcome_names, each = length(exposure_names)),
    lambda = lambda,
    selection_tol = selection_tol,
    beta_original = beta_original,
    beta_boot_mean_all = rowMeans(boot_output, na.rm = TRUE),
    beta_boot_median_all = apply(boot_output, 1L, row_median),
    beta_boot_mean_selected = rowMeans(selected_values, na.rm = TRUE),
    beta_boot_median_selected = apply(
      selected_values, 1L, row_median
    ),
    beta_boot_sd_all = apply(boot_output, 1L, stats::sd, na.rm = TRUE),
    beta_boot_sd_selected = apply(
      selected_values, 1L, stats::sd, na.rm = TRUE
    ),
    beta_boot_q025_all = apply(
      boot_output, 1L, row_quantile, prob = 0.025
    ),
    beta_boot_q975_all = apply(
      boot_output, 1L, row_quantile, prob = 0.975
    ),
    beta_boot_q025_selected = apply(
      selected_values, 1L, row_quantile, prob = 0.025
    ),
    beta_boot_q975_selected = apply(
      selected_values, 1L, row_quantile, prob = 0.975
    ),
    n_selected = n_selected,
    sel_prob = sel_prob,
    sel_prob_threshold = sel_prob_threshold,
    passes_sel_prob = passes_sel_prob,
    sign_consistency = sign_consistency,
    sign_consistency_threshold = sign_consistency_threshold,
    passes_sign_consistency = passes_sign_consistency,
    passes_stability_filter = passes_sel_prob &
      passes_sign_consistency,
    prop_positive = rowSums(positive, na.rm = TRUE) / n_boot,
    prop_negative = rowSums(negative, na.rm = TRUE) / n_boot,
    prop_zero = rowSums(!selected, na.rm = TRUE) / n_boot,
    stringsAsFactors = FALSE
  )

  ## Make the intended result explicit for coefficients never selected.
  out$beta_boot_mean_selected[n_selected == 0L] <- NA_real_
  out$beta_boot_median_selected[n_selected == 0L] <- NA_real_
  out$beta_boot_sd_selected[n_selected < 2L] <- NA_real_

  attr(out, "lambda") <- lambda
  attr(out, "selection_tol") <- selection_tol
  attr(out, "sel_prob_threshold") <- sel_prob_threshold
  attr(out, "sign_consistency_threshold") <- sign_consistency_threshold
  attr(out, "inference_note") <- if (lambda > 0) {
    paste(
      "Penalised fit: bootstrap summaries describe coefficient and selection",
      "stability. Classical p-values and confidence intervals are not",
      "calculated."
    )
  } else {
    paste(
      "Unpenalised fit: no p-values or confidence intervals are calculated",
      "automatically; their validity requires a separately justified",
      "inferential procedure."
    )
  }

  out
}


## Summarise fixed-subset, unpenalised (lambda = 0) bootstrap coefficients.
## Approximate two-sided p-values and 95% intervals use a normal approximation
## with the bootstrap standard error. Percentile intervals are also retained
## for evaluation in simulation studies. Inference is conditional on the
## exposure subset and rank having been fixed before this bootstrap.
summarize_msrrr_unpenalized_bootstrap <- function(
    boot_output,
    original_B,
    exposure_names = rownames(original_B),
    outcome_names = colnames(original_B),
    lambda = 0,
    alpha = 0.05) {

  boot_output <- as.matrix(boot_output)
  original_B <- as.matrix(original_B)

  if (length(lambda) != 1L || is.na(lambda) || !is.finite(lambda)) {
    stop("lambda must be one finite value.", call. = FALSE)
  }
  if (lambda > 0) {
    stop(
      paste(
        "summarize_msrrr_unpenalized_bootstrap() only accepts lambda = 0.",
        "Use summarize_msrrr_bootstrap() for penalised bootstrap fits."
      ),
      call. = FALSE
    )
  }
  if (lambda < 0) {
    stop("lambda must be zero; negative penalties are not supported.",
         call. = FALSE)
  }
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite value strictly between 0 and 1.",
         call. = FALSE)
  }
  if (is.null(exposure_names) || length(exposure_names) != nrow(original_B)) {
    stop("exposure_names must contain one name per row of original_B.",
         call. = FALSE)
  }
  if (is.null(outcome_names) || length(outcome_names) != ncol(original_B)) {
    stop("outcome_names must contain one name per column of original_B.",
         call. = FALSE)
  }
  if (nrow(boot_output) != length(original_B)) {
    stop(
      "boot_output must have nrow(original_B) * ncol(original_B) rows.",
      call. = FALSE
    )
  }
  if (ncol(boot_output) < 2L) {
    stop("At least two bootstrap samples are required.", call. = FALSE)
  }

  n_valid_bootstrap <- rowSums(is.finite(boot_output))
  bootstrap_values <- boot_output
  bootstrap_values[!is.finite(bootstrap_values)] <- NA_real_

  row_quantile <- function(x, prob) {
    if (sum(is.finite(x)) < 2L) return(NA_real_)
    as.numeric(stats::quantile(
      x, probs = prob, na.rm = TRUE, names = FALSE
    ))
  }

  beta_unpenalized <- as.vector(original_B)
  bootstrap_mean <- rowMeans(bootstrap_values, na.rm = TRUE)
  bootstrap_median <- apply(
    bootstrap_values, 1L, stats::median, na.rm = TRUE
  )
  bootstrap_standard_error <- apply(
    bootstrap_values, 1L, stats::sd, na.rm = TRUE
  )

  valid_standard_error <- is.finite(bootstrap_standard_error) &
    bootstrap_standard_error > 0 &
    is.finite(beta_unpenalized)

  z_value <- rep(NA_real_, length(beta_unpenalized))
  z_value[valid_standard_error] <- beta_unpenalized[valid_standard_error] /
    bootstrap_standard_error[valid_standard_error]

  # Approximate two-sided p-value using a normal approximation and the
  # bootstrap standard error.
  p_value <- rep(NA_real_, length(beta_unpenalized))
  p_value[valid_standard_error] <- 2 * stats::pnorm(
    abs(z_value[valid_standard_error]), lower.tail = FALSE
  )

  critical_value <- stats::qnorm(1 - alpha / 2)
  ci95_lower <- rep(NA_real_, length(beta_unpenalized))
  ci95_upper <- rep(NA_real_, length(beta_unpenalized))
  ci95_lower[valid_standard_error] <- beta_unpenalized[
    valid_standard_error
  ] - critical_value * bootstrap_standard_error[valid_standard_error]
  ci95_upper[valid_standard_error] <- beta_unpenalized[
    valid_standard_error
  ] + critical_value * bootstrap_standard_error[valid_standard_error]

  out <- data.frame(
    exposure = rep(exposure_names, times = length(outcome_names)),
    outcome = rep(outcome_names, each = length(exposure_names)),
    lambda = lambda,
    beta_unpenalized = beta_unpenalized,
    bootstrap_mean = bootstrap_mean,
    bootstrap_median = bootstrap_median,
    bootstrap_standard_error = bootstrap_standard_error,
    ci95_lower = ci95_lower,
    ci95_upper = ci95_upper,
    bootstrap_q025 = apply(
      bootstrap_values, 1L, row_quantile, prob = alpha / 2
    ),
    bootstrap_q975 = apply(
      bootstrap_values, 1L, row_quantile, prob = 1 - alpha / 2
    ),
    z_value = z_value,
    p_value = p_value,
    significant = !is.na(p_value) & p_value < alpha,
    n_valid_bootstrap = n_valid_bootstrap,
    stringsAsFactors = FALSE
  )

  invalid_rows <- n_valid_bootstrap < 2L
  if (any(invalid_rows)) {
    warning(
      sum(invalid_rows),
      " coefficient(s) had fewer than two finite bootstrap estimates; ",
      "their uncertainty summaries are unavailable.",
      call. = FALSE
    )
  }
  zero_se_rows <- is.finite(bootstrap_standard_error) &
    bootstrap_standard_error == 0
  if (any(zero_se_rows)) {
    warning(
      sum(zero_se_rows),
      " coefficient(s) had zero bootstrap standard error; their p-values ",
      "and normal-approximation intervals are unavailable.",
      call. = FALSE
    )
  }

  attr(out, "alpha") <- alpha
  attr(out, "inference_note") <- paste(
    "Conditional inference after an unpenalised refit with a fixed exposure",
    "subset and rank. Normal-approximation p-values and intervals use the",
    "bootstrap standard error. No multiple-testing correction is applied",
    "automatically."
  )

  out
}














############### CHECKS

check_integer <- function(x, y = x, type = "scalar", float = FALSE, min = 1,
                          max = Inf, max_message = NULL, exit_code = NULL,
                          min_message = NULL) {
  if (type %in% c("matrix", "data.frame")) {
    y_temp <- y
  }

  y <- tryCatch(
    as.double(as.matrix(y)),
    warning = function(w) {
      stop_msrrr_bootstrap(paste(x, "must be numeric."))
    }
  )

  if (any(is.na(y))) {
    stop_msrrr_bootstrap(paste(x, "must not be NA."))
  }

  if (type == "scalar" && length(y) != 1) {
    stop_msrrr_bootstrap(paste(x, "must be of length 1."))
  }

  if (!float) {
    if (any((y %% 1) != 0)) {
      stop_msrrr_bootstrap(paste(x, "must be an integer."))
    }
    y <- as.integer(y)
  }

  if (any(y < min)) {
    if (!is.null(min_message)) {
      stop_msrrr_bootstrap(min_message, exit_code = exit_code)
    } else {
      stop_msrrr_bootstrap(x, " must be higher than or equal to ", min, ".",
        exit_code = exit_code
      )
    }
  }

  if (any(y > max)) {
    if (!is.null(max_message)) {
      stop_msrrr_bootstrap(max_message, exit_code = exit_code)
    } else {
      stop_msrrr_bootstrap(x, " must be lower than or equal to ", max, ".",
        exit_code = exit_code
      )
    }
  }

  if (type %in% c("matrix", "data.frame")) {
    y <- matrix(
      y,
      NROW(y_temp),
      NCOL(y_temp),
      dimnames = dimnames(y_temp)
    )
  }

  if (type == "data.frame") {
    y <- as.data.frame(y)
  }

  return(y)
}
