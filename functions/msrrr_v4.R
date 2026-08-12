## engine: msRRR algorithm via MM + SRRR
## VERSION: v4 (12 August 2026)
## ORIGINAL METHODOLOGY AND v3 IMPLEMENTATION: Augusto Anguita-Ruiz
## v4 UPDATES, ROBUSTNESS AND MAINTENANCE: María Arteaga Jover
##
## Practical recommendations and automatic diagnostics:
## A) Rank range: warns when a requested rank exceeds min(ncol(X), ncol(Y)),
##    and when CV selects the smallest tested rank (if > 1) or the largest
##    tested rank (if below the maximum feasible rank). It also warns when the
##    final penalised model selects fewer exposures than its selected rank,
##    because its effective coefficient rank must then be lower. Infeasible
##    requested ranks are not removed automatically.
## B) Lambda range: checks every candidate rank, warning when its CV optimum is
##    at either endpoint or its largest lambda still selects exposures. It then
##    proposes (without running) a common grid for all ranks, extended tenfold
##    where needed and with an integer size that preserves the log resolution.
##    A sparse solution at the smallest lambda is not itself treated as an error.
## C) CV criterion: for all-Gaussian outcomes, notes that pMSE may be preferable
##    when outcomes have different scales. For any non-Gaussian outcome, warns
##    when pMSE is requested and recommends deviance.
## D) Missing outcomes: reports outcomes with at least 20% missing values and
##    warns when any outcome has at least 50% missing values.
## E) Outcome relationships: for all-Gaussian outcomes, notes when more than
##    half of the evaluable outcome pairs have |r| < 0.25.
## F) Rare dummy predictors: for 0/1 columns, reports when the least frequent
##    value represents 1% to <5% of observed rows, warns when it represents
##    <1%, and warns when a CV training fold contains only one observed level.
## These checks belong to msrrr_cv(). Set diagnostics = FALSE to suppress
## checks A-F. When diagnostics = TRUE, the emitted messages and warnings are
## always stored in object$diagnostics, are printed together at the end by
## default, and are additionally written to disk only when diagnostics_file is
## supplied.
##
## Compatibility and recommended interfaces:
## - msrrr() and msrrr.tuning() retain the historical v3 selection workflow.
## - msrrr_cv() and msrrr_cv.tuning() provide the recommended self-contained
##   CV workflow. IMPORTANT: supply analysis-ready numeric X, Z and Y matrices
##   after variable coding/dummy creation but before standardisation or outcome
##   scaling. Do not pass matrices that have already been standardised. Within
##   each CV iteration, preprocessing parameters are
##   estimated using only the training folds and then applied unchanged to
##   validation.
##   Continuous X/Z columns and Gaussian outcomes are z-score standardised,
##   0/1 X/Z columns remain unchanged, and non-Gaussian outcomes remain on
##   their original scales. The final selected
##   model is refitted after learning preprocessing from the complete dataset
##   passed to msrrr_cv(), and those parameters are stored for prediction.
##
## Primary corrections in v4:
## 1) Preserves A, V and C as matrices when rank = 1 and when values are
##    extracted from stored solution paths, preventing dropped dimensions
##    (Marco's change).
## 2) Corrects the Gaussian likelihood used by AIC, BIC, BICP and GIC.
##    The previous code passed phi (a variance) to dnorm() as if it were a
##    standard deviation; v4 passes sqrt(phi), which is the required SD.
## 3) Orders every lambda sequence from largest to smallest when warm starts
##    are enabled, including sequences supplied by the user.
## 4) Resets model initialisation at the beginning of every CV fold. Warm
##    starts are reused only across lambda values within the same training
##    fold, preventing information leakage between folds.
## 5) Preserves the exact initial values used for the final refit and provides
##    refit_msrrr() to reproduce it with the same rank, lambda, control,
##    intercept handling and initialisation. This is essential because msRRR
##    is non-convex and a different initialisation can select a different model.
##
## Secondary and robustness improvements in v4:
## 6) Handles the intercept consistently. The previous code always added one,
##    even if Z already contained it, and an intercept-only Z became a vector
##    inside CV. v4 adds an intercept only when needed and keeps fold subsets
##    as matrices. Z = NULL therefore fits an intercept-only model safely.
## 7) Checks the assignment of families to outcomes before fitting. With one
##    family and familygroup = NULL, every column of Y uses that family. With
##    multiple families, familygroup must contain one valid family index per
##    column of Y; otherwise the function stops with an explanatory error.
## 8) Associates each outcome group with family[[group]] explicitly, so the
##    result no longer depends on the order returned by unique(familygroup).
## 9) Labels fitted objects with class "msrrr". This makes the generic call
##    predict(object, ...) automatically use predict.msrrr(); without the
##    class, R does not know which model-specific prediction method to call.
## 10) Stops before fitting if a family, link or outcome value is unsupported.
##    The accepted combinations are Gaussian/identity with finite values,
##    binomial/logit with 0 or 1, and Poisson/log with non-negative integers.
## 11) Makes prediction usable with or without observed outcomes. type selects
##    response or link-scale predictions; pMSE or deviance is calculated only
##    when Y.new is supplied. Invalid type or cv.criteria values produce a
##    clear error.
## 12) Reports convergence status and the final relative objective change.
## 13) Provides count_selected_exposures(), using a numerical tolerance to
##    ignore tiny floating-point residuals. With the default threshold 1e-8, an
##    exposure is selected only if the L2 norm of its complete coefficient row
##    is greater than 1e-8. This changes only the reported count, not the fitted
##    coefficients or the selected model.
## 14) Stores the response-specific Kappa candidates and scalar kappa along
##    the lambda path and CV folds for later diagnosis. This does not change
##    how kappa is calculated or used.
## 15) Uses one numerical tolerance (1e-8) whenever the script distinguishes
##    zero from non-zero coefficients. Exposure selection uses the L2 norm of
##    the complete coefficient row, so all rank components are considered.
## 16) Provides refit_msrrr_unpenalized() as a separate optional post-selection
##    helper. It can take the final whole-sample penalised fit returned by
##    refit_msrrr(), restricts X to the exposures selected in that fit, keeps
##    its rank when feasible, and fits lambda = 0 with init = NULL by default.
##    It does not alter the penalised fit or its bootstrap.
## 17) Passes the family name ("gaussian", "binomial" or "poisson") rather
##    than the complete R family object to glmnet() during ridge
##    initialisation. Passing binomial() made glmnet use internally normalised
##    fractional weights and produced the misleading warning "non-integer
##    #successes in a binomial glm!" even when the response contained only
##    0/1 values.
## 18) Z-score standardises Gaussian outcomes instead of min-max scaling them.
##    The same family-based preprocessing is used by CV, whole-sample refits,
##    prediction and LOCO; each bootstrap keeps its parent fit's transformation
##    fixed. Legacy fitted objects storing min-max parameters remain usable for
##    transformation and prediction.
## 19) Uses warm = FALSE by default in the recommended msrrr_cv() workflow, so
##    every rank/lambda fit starts independently. The historical msrrr() and
##    msrrr.tuning() compatibility interfaces retain their original defaults.
##
## HELIX reproducibility extension (analysis-specific; not the default rule):
## - msrrr_cv() accepts optional X_no_scale and Z_no_scale character vectors.
## - NULL preserves the standard v4 behaviour: binary 0/1 columns retain their
##   original scale and all other eligible columns are standardised.
## - Supplied columns retain their original scale in every CV fold, the final
##   development fit, prediction and subsequent refits. Parameters for all
##   other columns remain leakage-free. This permits reproduction of the
##   historical HELIX <7-distinct-values rule and is not a new recommendation.

require(glmnet)
require(rrpack)

## 15) Common zero/non-zero tolerance used throughout the script.
.msrrr_selection_tol <- 1e-8

## Practical checks C), D), E) and F): inspect inputs without changing them.
.diagnose_msrrr_inputs <- function(Y, X, family, familygroup,
                                   method, cv.criteria) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  spec <- .validate_family_spec(Y, family, familygroup)
  family_names <- vapply(spec$fg, function(fam) fam$family, character(1))
  outcome_names <- colnames(Y)
  if (is.null(outcome_names)) {
    outcome_names <- paste0("Y", seq_len(ncol(Y)))
  }

  # Practical check C) Match the CV criterion to the declared families.
  if (identical(method, "CV")) {
    if (all(family_names == "gaussian") && identical(cv.criteria, "deviance")) {
      message(
        "All outcomes are Gaussian. Deviance is valid, but pMSE may be ",
        "preferable when outcomes have different scales because it ",
        "standardises prediction error by outcome variance."
      )
    }
    if (any(family_names != "gaussian") && identical(cv.criteria, "pMSE")) {
      warning(
        "pMSE was requested with non-Gaussian or mixed-family outcomes. ",
        "Deviance is recommended for binomial, Poisson and mixed-family ",
        "models; verify that pMSE is intentional.",
        call. = FALSE
      )
    }
  }

  # Practical check D) Report missingness separately for each outcome.
  missing_fraction <- colMeans(is.na(Y))
  moderate <- which(missing_fraction >= 0.20 & missing_fraction < 0.50)
  high <- which(missing_fraction >= 0.50)
  format_missing <- function(idx) {
    paste0(
      outcome_names[idx], " (",
      formatC(100 * missing_fraction[idx], format = "f", digits = 1),
      "%)",
      collapse = ", "
    )
  }
  if (length(moderate) > 0L) {
    message(
      "Some outcomes contain at least 20% missing values: ",
      format_missing(moderate),
      ". Greater outcome missingness reduces the information available for ",
      "fitting and CV. Review the missingness pattern, its possible causes and ",
      "the effective number of observed cases. Depending on the study design, ",
      "consider whether preprocessing, excluding selected outcomes or ",
      "participants, or sensitivity analyses would be appropriate."
    )
  }
  if (length(high) > 0L) {
    warning(
      "Some outcomes contain at least 50% missing values: ",
      format_missing(high),
      ". msRRR can accommodate missing outcomes, but estimates and CV ",
      "performance may be unstable when the observed information is limited. ",
      "Carefully review the missingness pattern and mechanism. Depending on ",
      "the study design, consider excluding highly incomplete outcomes, ",
      "excluding participants with extensive missingness, using alternative ",
      "missing-data strategies, or evaluating these choices in sensitivity ",
      "analyses. Any decision should be justified.",
      call. = FALSE
    )
  }

  # Practical check E) Quantify weak pairwise Gaussian-outcome correlations.
  if (ncol(Y) >= 2L && all(family_names == "gaussian")) {
    correlation_matrix <- suppressWarnings(
      stats::cor(Y, use = "pairwise.complete.obs")
    )
    diag(correlation_matrix) <- NA_real_
    abs_correlations <- abs(
      correlation_matrix[upper.tri(correlation_matrix)]
    )
    abs_correlations <- abs_correlations[is.finite(abs_correlations)]
    if (length(abs_correlations) > 0L) {
      weak_pair_fraction <- mean(abs_correlations < 0.25)
      if (weak_pair_fraction > 0.50) {
        message(
          "More than half of the evaluable Gaussian outcome pairs have ",
          "|r| < 0.25 (",
          formatC(100 * weak_pair_fraction, format = "f", digits = 1),
          "% of pairs). msRRR is most useful when outcomes are expected to ",
          "share a latent or exposure-related structure; review whether joint ",
          "reduced-rank modelling is scientifically justified."
        )
      }
    }
  }

  # Practical check F) Report rare 0/1 columns using percentage thresholds.
  predictor_names <- colnames(X)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(X)))
  }
  uncommon_dummy <- logical(ncol(X))
  rare_dummy <- logical(ncol(X))
  rare_details <- character(ncol(X))
  for (j in seq_len(ncol(X))) {
    observed <- X[!is.na(X[, j]), j]
    values <- sort(unique(observed))
    if (length(values) >= 1L && all(values %in% c(0, 1))) {
      counts <- table(factor(observed, levels = c(0, 1)))
      least_count <- min(counts)
      least_fraction <- if (length(observed) > 0L) {
        least_count / length(observed)
      } else {
        0
      }
      uncommon_dummy[j] <- least_fraction >= 0.01 &&
        least_fraction < 0.05
      rare_dummy[j] <- least_fraction < 0.01
      rare_details[j] <- paste0(
        predictor_names[j], " (least frequent value: ", least_count,
        "/", length(observed), ")"
      )
    }
  }
  if (any(uncommon_dummy)) {
    message(
      "Some binary/dummy predictors have a least frequent value representing ",
      "between 1% and <5% of observed rows: ",
      paste(rare_details[uncommon_dummy], collapse = ", "),
      ". Review whether these categories have sufficient representation for ",
      "stable estimation and CV."
    )
  }
  if (any(rare_dummy)) {
    warning(
      "Very rare binary/dummy predictors detected (<1% representation for ",
      "the least frequent value): ",
      paste(rare_details[rare_dummy], collapse = ", "),
      ". Coefficients and CV results may be unstable; reconsider whether these ",
      "categories can be estimated reliably and inspect their representation ",
      "within CV folds.",
      call. = FALSE
    )
  }
}


## Practical check F): detect dummy columns losing a level in CV training data.
.diagnose_dummy_folds <- function(X, foldid) {
  if (is.null(foldid)) return(invisible(NULL))
  X <- as.matrix(X)
  predictor_names <- colnames(X)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(X)))
  }

  binary_columns <- vapply(seq_len(ncol(X)), function(j) {
    values <- sort(unique(X[!is.na(X[, j]), j]))
    length(values) >= 1L && all(values %in% c(0, 1))
  }, logical(1))

  affected <- character()
  for (ifold in sort(unique(foldid))) {
    training_rows <- foldid != ifold
    for (j in which(binary_columns)) {
      values <- unique(X[training_rows & !is.na(X[, j]), j])
      if (length(values) < 2L) {
        affected <- c(
          affected,
          paste0(predictor_names[j], " (training fold ", ifold, ")")
        )
      }
    }
  }
  affected <- unique(affected)
  if (length(affected) > 0L) {
    warning(
      "Some binary/dummy predictors have only one observed level in a CV ",
      "training set: ", paste(affected, collapse = ", "),
      ". The corresponding coefficients and fold performance may be unstable.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


## Practical checks A) and B): inspect the selected rank and lambda boundaries.
.diagnose_requested_ranks <- function(nrankseq, p, q) {
  max_feasible_rank <- min(p, q)
  infeasible_ranks <- nrankseq[nrankseq > max_feasible_rank]
  if (length(infeasible_ranks) > 0L) {
    warning(
      "Requested rank value(s) ",
      paste(unique(infeasible_ranks), collapse = ", "),
      " exceed the maximum feasible rank min(ncol(X), ncol(Y)) = ",
      max_feasible_rank,
      ". These values are not removed automatically, but the corresponding ",
      "fits will generally fail during SVD-based initialisation because the ",
      "requested number of components does not exist.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


## Practical checks A) and B): inspect the selected rank and lambda boundaries.
.diagnose_msrrr_selection <- function(object, nrankseq, p, q,
                                      tol = .msrrr_selection_tol) {
  if (!identical(object$method, "CV")) return(invisible(NULL))

  # Practical check A) Warn when the selected rank lies on a testable boundary.
  min_rank <- min(nrankseq)
  max_rank <- max(nrankseq)
  max_feasible_rank <- min(p, q)
  if (object$nrank.opt == min_rank && min_rank > 1L) {
    warning(
      "The selected rank (", object$nrank.opt,
      ") is the smallest value tested. Consider extending nrankseq downward ",
      "and refitting with the same CV folds.",
      call. = FALSE
    )
  }
  if (object$nrank.opt == max_rank && max_rank < max_feasible_rank) {
    warning(
      "The selected rank (", object$nrank.opt,
      ") is the largest value tested. Consider extending nrankseq upward ",
      "and refitting with the same CV folds. The maximum feasible rank is ",
      max_feasible_rank, ".",
      call. = FALSE
    )
  }

  active_in_selected_model <- count_selected_exposures(
    object$fit$B, tol = tol
  )
  if (active_in_selected_model < object$nrank.opt) {
    warning(
      "The final penalised model selects ", active_in_selected_model,
      " exposure(s), fewer than its selected rank (", object$nrank.opt,
      "). Its coefficient matrix can therefore have effective rank at most ",
      active_in_selected_model,
      ", so one or more latent components are redundant. Consider including ",
      "lower ranks in nrankseq and reviewing the CV results.",
      call. = FALSE
    )
  }

  # Practical check B) Inspect the complete lambda path for every candidate
  # rank, not only for the rank selected by CV. A truncated path in another
  # rank could conceal a better rank/lambda combination.
  rank_objects <- object$out.allrank
  if (is.null(rank_objects)) rank_objects <- list(object)
  rank_values <- vapply(
    rank_objects,
    function(rank_object) as.integer(rank_object$nrank),
    integer(1)
  )
  common_lamseq <- object$lamseq
  min_lambda <- min(common_lamseq)
  max_lambda <- max(common_lamseq)
  format_lambda <- function(value) {
    formatC(value, format = "e", digits = 4)
  }

  optimal_at_min <- vapply(
    rank_objects,
    function(rank_object) rank_object$lam.opt == min(rank_object$lamseq),
    logical(1)
  )
  optimal_at_max <- vapply(
    rank_objects,
    function(rank_object) rank_object$lam.opt == max(rank_object$lamseq),
    logical(1)
  )
  active_at_max <- vapply(
    rank_objects,
    function(rank_object) {
      largest_index <- which.max(rank_object$lamseq)
      A_largest <- .as_matrix_dim(
        rank_object$Apath[largest_index, , ],
        nrow = p, ncol = rank_object$nrank
      )
      V_largest <- .as_matrix_dim(
        rank_object$Vpath[largest_index, , ],
        nrow = q, ncol = rank_object$nrank
      )
      count_selected_exposures(A_largest %*% t(V_largest), tol = tol)
    },
    integer(1)
  )

  extend_lower <- min_lambda > 0 && any(optimal_at_min)
  extend_upper <- any(optimal_at_max) || any(active_at_max > 0L)

  if (extend_lower || extend_upper) {
    proposed_min <- if (extend_lower) min_lambda / 10 else min_lambda
    proposed_max <- if (extend_upper) max_lambda * 10 else max_lambda

    # Preserve approximately the current number of positive lambda intervals
    # per log10 decade. ceiling() followed by as.integer() guarantees that the
    # suggested length.out is a whole number and never reduces path resolution.
    positive_lambdas <- sort(unique(common_lamseq[common_lamseq > 0]))
    if (length(positive_lambdas) >= 2L &&
        max(positive_lambdas) > min(positive_lambdas) &&
        proposed_min > 0) {
      current_decades <- log10(
        max(positive_lambdas) / min(positive_lambdas)
      )
      proposed_decades <- log10(proposed_max / proposed_min)
      intervals_per_decade <-
        (length(positive_lambdas) - 1L) / current_decades
      proposed_positive_n <- as.integer(
        ceiling(intervals_per_decade * proposed_decades) + 1L
      )
    } else {
      proposed_positive_n <- as.integer(max(2L, length(positive_lambdas)))
    }
    include_zero <- any(common_lamseq == 0)
    proposed_nlam <- as.integer(proposed_positive_n + as.integer(include_zero))

    details <- character()
    if (any(optimal_at_min)) {
      details <- c(
        details,
        paste0(
          "CV selects the smallest lambda for rank(s) ",
          paste(rank_values[optimal_at_min], collapse = ", "), "."
        )
      )
    }
    if (any(optimal_at_max)) {
      details <- c(
        details,
        paste0(
          "CV selects the largest lambda for rank(s) ",
          paste(rank_values[optimal_at_max], collapse = ", "), "."
        )
      )
    }
    if (any(active_at_max > 0L)) {
      affected <- which(active_at_max > 0L)
      details <- c(
        details,
        paste0(
          "The largest lambda still selects exposures for rank(s) ",
          paste0(
            rank_values[affected], " (n = ", active_at_max[affected], ")",
            collapse = ", "
          ), "."
        )
      )
    }

    sequence_code <- if (include_zero) {
      paste0(
        "lamseq <- sort(unique(c(0, 10^seq(log10(",
        format_lambda(proposed_max), "), log10(",
        format_lambda(proposed_min), "), length.out = ",
        proposed_positive_n, "L))), decreasing = TRUE)"
      )
    } else {
      paste0(
        "lamseq <- 10^seq(log10(", format_lambda(proposed_max),
        "), log10(", format_lambda(proposed_min),
        "), length.out = ", proposed_nlam, "L)"
      )
    }

    warning(
      paste(details, collapse = " "), " ",
      "The common lambda grid may not cover the optimum for every rank. ",
      "Rerun all ranks with the following common sequence (",
      proposed_nlam, " total lambda values, preserving approximately the ",
      "current log-scale resolution):\n", sequence_code,
      call. = FALSE
    )
  }
  invisible(NULL)
}


## 6) Add an intercept only when Z does not already contain one.
## This keeps msrrr.fit() compatible with both Z and cbind(1, Z).
.ensure_intercept <- function(Z, n, tol = sqrt(.Machine$double.eps)) {
  if (is.null(Z)) {
    return(matrix(1, nrow = n, ncol = 1,
                  dimnames = list(NULL, "(Intercept)")))
  }

  Z <- as.matrix(Z)
  if (nrow(Z) != n) stop("Z must have the same number of rows as Y/X.")

  has_intercept <- ncol(Z) > 0 && any(vapply(
    seq_len(ncol(Z)),
    function(j) all(is.finite(Z[, j])) && max(abs(Z[, j] - 1)) <= tol,
    logical(1)
  ))

  if (!has_intercept) {
    Z <- cbind("(Intercept)" = rep(1, n), Z)
  }

  Z
}

## 7) Validate and standardise familygroup before it is used.
.validate_family_spec <- function(Y, family, familygroup = NULL,
                                  check_response = TRUE) {
  Y <- as.matrix(Y)
  q <- ncol(Y)

  if (!is.list(family) || length(family) < 1L) {
    stop("family must be a non-empty list of R family objects.")
  }

  if (is.null(familygroup)) {
    if (length(family) != 1L) {
      stop("familygroup is required when more than one family is supplied.")
    }
    familygroup <- rep.int(1L, q)
  }

  if (length(familygroup) != q) {
    stop("familygroup must have one element per outcome (one per column of Y).")
  }
  if (anyNA(familygroup) ||
      any(familygroup != as.integer(familygroup)) ||
      any(!familygroup %in% seq_along(family))) {
    stop("familygroup must contain valid integer indices into family.")
  }
  familygroup <- as.integer(familygroup)

  used_groups <- sort(unique(familygroup))
  if (!identical(used_groups, seq_along(family))) {
    stop("Each supplied family must be used and numbered consecutively in familygroup.")
  }

  # 10) Restrict fitting to supported families, canonical links and valid data.
  supported_links <- c(
    gaussian = "identity",
    binomial = "logit",
    poisson = "log"
  )

  for (g in seq_along(family)) {
    fam <- family[[g]]
    if (is.null(fam$family) || is.null(fam$link) ||
        !fam$family %in% names(supported_links)) {
      stop(
        "Unsupported family in family[[", g,
        "]]. Supported families are gaussian, binomial and poisson."
      )
    }
    if (!identical(fam$link, unname(supported_links[fam$family]))) {
      stop(
        "family[[", g, "]] must use the canonical ",
        supported_links[fam$family], " link for ", fam$family, "."
      )
    }

    if (check_response) {
      yg <- Y[, familygroup == g, drop = FALSE]
      observed <- yg[!is.na(yg)]
      if (any(!is.finite(observed))) {
        stop("Observed outcomes assigned to family[[", g, "]] must be finite.")
      }
      if (fam$family == "binomial" && any(!observed %in% c(0, 1))) {
        stop("Binomial outcomes must contain only 0, 1 or NA.")
      }
      if (fam$family == "poisson" &&
          any(observed < 0 | observed != floor(observed))) {
        stop("Poisson outcomes must contain only non-negative integers or NA.")
      }
    }
  }

  list(
    familygroup = familygroup,
    fg = family[familygroup],
    groups = used_groups
  )
}


## 1) Preserve dimensions after array subsetting, especially at rank = 1.
.as_matrix_dim <- function(x, nrow, ncol) {
  if (length(x) != nrow * ncol) {
    stop("Initial-value dimensions are incompatible with the requested model.")
  }
  matrix(x, nrow = nrow, ncol = ncol)
}


## 13/15) Count active exposure rows using the common numerical tolerance.
count_selected_exposures <- function(model_or_B,
                                     tol = .msrrr_selection_tol) {
  B <- if (is.matrix(model_or_B)) {
    model_or_B
  } else if (!is.null(model_or_B$fit$B)) {
    model_or_B$fit$B
  } else if (!is.null(model_or_B$B)) {
    model_or_B$B
  } else {
    stop("Provide a coefficient matrix B or an msrrr/msrrr.fit object.")
  }

  sum(sqrt(rowSums(B^2)) > tol)
}


## matrix row-wise 2_1 norm for group sparsity 
norm21 <- function(A) {
  A <- as.matrix(A)
  sum(sqrt(rowSums(A^2)))
}

##' @importFrom stats dnorm dpois dbinom
logLikehood <- function(Y, MU, Sigma = 1, family){
  out <- switch(
    family,
    "gaussian" = dnorm(Y, MU, Sigma, log = TRUE), 
    "poisson"  = dpois(Y, MU, log = TRUE),
    "binomial" = dbinom(Y, 1, MU, log = TRUE),
    NULL
  )
  if (is.null(out)) stop("Unsupported family: ", family)
  out
}

## neg log-lik, to minimize
objFun <- function(Y, mu, Phi, family, familygroup = NULL) {
  spec <- .validate_family_spec(Y, family, familygroup)
  familygroup <- spec$familygroup
  n = nrow(Y)
  temp <- vector()
  # 8) Use the actual group identifier to select both outcomes and family.
  for (g in spec$groups) {
    idx <- familygroup == g
    if(family[[g]]$family == "gaussian"){
      Sig = t(matrix(Phi[idx], sum(idx), n))
    }else{Sig = 1}
    temp[g] <- -sum(logLikehood(Y[, idx, drop = FALSE],
                                mu[, idx, drop = FALSE],
                                sqrt(Sig), family[[g]]$family), na.rm = TRUE)
  }
  sum(temp)
}

# pseudo MSE...
pmse <- function(Y, mu, Phi=NULL, family, familygroup = NULL) {
  spec <- .validate_family_spec(Y, family, familygroup)
  familygroup <- spec$familygroup
  n = nrow(Y)
  m = sum(!is.na(Y))
  Yvar = apply(Y,2,var,na.rm=T)
  Yvar = t(matrix(Yvar, ncol(Y), n))
  temp <- vector()
  for (g in spec$groups) {
    idx <- familygroup == g
    if(family[[g]]$family == "gaussian"){
      # Sig = t(matrix(Phi[familygroup == cfamily[j]],
      #                sum(familygroup == cfamily[j]), n))
      Sig = Yvar[, idx, drop = FALSE]
      temp[g] <- sum((Y-mu)[, idx, drop = FALSE]^2 / Sig, na.rm = T)
    }else if(family[[g]]$family == "binomial"){
      # Sig = ((1-mu)*mu)[, familygroup == cfamily[j]]
      Sig = Yvar[, idx, drop = FALSE]
      temp[g] <- sum((Y-mu)[, idx, drop = FALSE]^2 / Sig, na.rm = T)
    }else if(family[[g]]$family == "poisson"){
      # Sig = mu[, familygroup == cfamily[j]]
      Sig = Yvar[, idx, drop = FALSE]
      temp[g] <- sum((Y-mu)[, idx, drop = FALSE]^2 / Sig, na.rm = T)
    } 
  }
  return(sum(temp)/m)
}


## msrrr fit with prespecified nrank and lambda
# Z may be supplied with or without an intercept column
msrrr.fit <- function(Y, X, Z, family, familygroup = NULL, nrank = 2, lambda, init = NULL,
                      control = list(epsilon = 1e-4, maxit = 200, trace = F),
                      ensure_intercept = TRUE) {
  
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  spec <- .validate_family_spec(Y, family, familygroup)
  familygroup <- spec$familygroup
  fg <- spec$fg

  # 6) Accept Z, cbind(1, Z), or Z = NULL with consistent intercept handling.
  if (ensure_intercept) {
    Z <- .ensure_intercept(Z, n)
  } else {
    Z <- as.matrix(Z)
  }

  pz = ncol(Z)
  Y_mis = Y
  id.mis = which(is.na(Y))
  r = nrank
  
  ## get init values via ridge
  if (is.null(init$A0)) {  # use ridge to get init working Y
    C0 = matrix(NA, pz, q)
    B0 = matrix(NA, p, q)
    
    for (iq in 1:q) {
      # fit0 <- glm(Y[,iq] ~ 0+Z+X, family = fg[[iq]])
      idx = !is.na(Y_mis[, iq])
      
      # ridge
      fit0 <- glmnet(
        cbind(Z, X)[idx, ],
        Y[idx, iq],
        # 17) Use glmnet's native family implementation. Supplying the full
        # binomial() object triggers a spurious fractional-success warning.
        family = fg[[iq]]$family,
        alpha = 0,
        lambda = 0.03,
        intercept = F,
        penalty.factor = c(rep(0, pz), rep(1, p))
      )
      
      C0[, iq] = as.numeric(fit0$beta)[1:pz]
      B0[, iq] = as.numeric(fit0$beta)[-c(1:pz)]
    }
    
    svdB0 = svd(B0)
    A0 = svdB0$u[, 1:r, drop = F] %*% diag(svdB0$d[1:r], r, r)
    V0 = svdB0$v[, 1:r, drop = F]
    
    # 1) Keep the rank dimension when nrank = 1.
    if (is.null(dim(A0))) A0 <- matrix(A0, ncol = 1)
    if (is.null(dim(V0))) V0 <- matrix(V0, ncol = 1)
    
    init$A0 = A0
    init$V0 = V0
    init$C0 = C0
    
  } else {
    # 1) Restore A/V/C matrix dimensions for supplied warm starts.
    A0 = .as_matrix_dim(init$A0, nrow = p, ncol = r)
    V0 = .as_matrix_dim(init$V0, nrow = q, ncol = r)
    C0 = .as_matrix_dim(init$C0, nrow = pz, ncol = q)
  }

  # 5) Retain the exact unscaled initial values used by this fit.
  init_used <- list(A0 = A0, V0 = V0, C0 = C0)
  # if (!is.null(init$kappa)) init_used$kappa <- init$kappa
  
  B0 = A0 %*% t(V0) # plot(Bt, B0)
  eta = (X %*% B0 + Z %*% C0) # nature parameter matrix
  mu = matrix(NA, n, q)        # mean matrix
  phi = rep(NA, q)             # dispersion par
  
  for (iq in 1:q) {
    mu[, iq] = fg[[iq]]$linkinv(eta[, iq])
  }
  
  Y[id.mis] = mu[id.mis]  # Y: imputed working outcome
  
  phi = ifelse(
    unlist(lapply(fg, function(a) a$family)) == "gaussian",
    colMeans((Y_mis - mu)^2, na.rm = T),
    1
  )
  
  # plot(Y_mis, mu)
  
  ## scaling predictor matrix
  Kappa <- rep(NA_real_, q)
  svdX0d1 <- svd(cbind(Z, X))$d[1]

  for (j in 1:q) {
    Kappa[j] <- switch(
      fg[[j]]$family,
      "gaussian" = svdX0d1 / phi[j],
      "binomial" = svdX0d1 / 2,
      "poisson" = svdX0d1 * quantile(Y_mis[, j], 0.9, na.rm = T)
    )
  }

  if (is.null(init$kappa)) {
    kappa <- max(Kappa)
    init$kappa = kappa
  } else {
    kappa = init$kappa
  }
  
  init_used$kappa <- kappa
  
  # kappa = 1
  X = X / kappa
  Z = Z / kappa
  A0 = A0 * kappa
  B0 = B0 * kappa
  C0 = C0 * kappa
  
  C = C0
  B = B0
  A = A0
  V = V0
  
  # mean((B/kappa - Bt)^2) # 0.00417
  dif = obj = rep(NA, control$maxit + 1)
  
  obj[1] = objFun(Y_mis, mu, phi, family, familygroup) +
    lambda * norm21(A0)
  
  for (iter in 1:control$maxit) {
    # if (control$trace) message(gettextf("iteration %d", iter), domain = NA)
    
    R = (t(X) %*% (Y - mu) %*% diag(1 / phi) + B0)
    
    solve.srrr = srrr(
      Y = R,
      X = diag(1, p),
      nrank,
      A0 = A0,
      V0 = V0,
      modstr = list(lamA = c(lambda), nlam = 1),
      control = list(epsilon = control$epsilon, maxit = 10)
    )
    
    A = solve.srrr$A.path[1, , ]
    V = solve.srrr$V.path[1, , ]
    
    # 1) Array subsetting can drop the rank dimension; restore it.
    if (is.null(dim(A))) A <- matrix(A, ncol = 1)
    if (is.null(dim(V))) V <- matrix(V, ncol = 1)
    
    B = A %*% t(V)
    
    # plot(A0, A)
    
    ## update C
    for (iq in 1:q) {
      # cat('.')
      # fine to use glm here as Z is low-d
      fit0 <- glm(
        Y_mis[, iq] ~ 0 + Z,
        offset = X %*% B[, iq],
        family = fg[[iq]]
      )
      
      C[, iq] = fit0$coef
    }
    
    ## update working Y, and phi
    eta <- X %*% B + Z %*% C
    
    for (iq in 1:q) {
      fm = fg[[iq]]
      mu[, iq] = fm$linkinv(eta[, iq])
    }
    
    Y[id.mis] = mu[id.mis]
    
    phi = ifelse(
      unlist(lapply(fg, function(a) a$family)) == "gaussian",
      colMeans((Y_mis - mu)^2, na.rm = T),
      1
    )
    
    ## conv
    # if (sum((eta - etaold)^2) < tol * sum(eta^2)) break
    # dif[iter] <- sum((B - B0)^2)/(sum(B0^2) + 1e-6)
    
    obj[iter + 1] = objFun(Y_mis, mu, phi, family, familygroup) +
      lambda * norm21(A)
    
    dif[iter] = obj[iter + 1] / obj[iter] - 1
    
    B0 = B
    A0 = A
    V0 = V
    C0 = C
    
    # if (dif[iter] < control$epsilon) break
    if (dif[iter] > 0 & control$trace == T) {
      warning("obj increased")
    }
    
    if (abs(dif[iter]) < control$epsilon) {
      break
    }
  }
  
  dif = dif[1:iter]
  obj = obj[1:(iter + 1)]
  
  A = A / kappa
  B = B / kappa
  C = C / kappa
  
  return(
    list(
      B = B,
      A = A,
      V = V,
      C = C,
      phi = phi,
      mu = mu,
      dif = dif,
      obj = obj,
      Y_imp = Y,
      Kappa = Kappa,
      kappa = kappa,
      iter = iter,
      # 12) Report whether the stopping criterion was reached.
      converged = isTRUE(abs(dif[iter]) < control$epsilon),
      final_relative_obj_change = dif[iter],
      init_used = init_used
    )
  )
}


## wrapper: select lambda with prespecified nrank, by IC or CV selection 
# Z is handled internally with an intercept column


## ---------------------------------------------------------------------------
## Historical v3 tuning and rank-selection wrappers
## ---------------------------------------------------------------------------
## These two public functions retain the original v3 selection workflow for
## backward compatibility. They use the corrected v4 msrrr.fit() engine, but
## they intentionally do not use the fold-specific preprocessing and
## diagnostics implemented in msrrr_cv().

msrrr.tuning = function(Y, X, Z, family, familygroup, nrank=2, init=NULL, 
                        lamseq=NULL, nlam=50, warm=T,  # lam.max, lam.min, 
                        method='CV', cv.criteria='pMSE', foldid=NULL, nfold=5, c.BIC=2, 
                        # c('CV', 'BIC', 'BICP', 'AIC', 'GIC')
                        control=list(epsilon=1e-4, maxit=200, trace=F, conv.obj=T)){
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  # Z = cbind(1,Z)
  pz = ncol(Z) 
  Y_mis = Y
  fg = family[familygroup]
  id.mis = which(is.na(Y))
  r = nrank  
  
  ## distribution families ##
  nfamily <- length(family)
  if (nfamily == 1 & is.null(familygroup)) familygroup <- rep(1, q)
  ## characters of families
  cfamily <- unique(familygroup)
  
  ## default lambda seq
  if (is.null(lamseq)) {
    cat("No lambda sequence provided, use default!\n") 
    yt <- rep(0, n)
    for (i in 1:nfamily) {   
      yt <- yt + switch(family[[i]]$family,
                        'gaussian' = apply(abs(2 * Y[, familygroup == i]), 1, function(a) sum(a, na.rm = TRUE)),
                        'binomial' = apply(abs(2 - Y[, familygroup == i]), 1, function(a) sum(a, na.rm = TRUE)),
                        'poisson'  = apply(abs(2 * Y[, familygroup == i] - 2), 1, function(a) sum(a, na.rm = TRUE)) )
    }
    lam.max <- sum(yt * apply(abs(X), 1, max)) / 1000 ## ???? 30-30000 too large?
    lamseq <- 10 ^ (seq(log10(lam.max * 0.001), log10(lam.max), len = nlam)) 
    if(control$trace==T) cat('lambda:', lamseq)
  }
  nlam = length(lamseq)
  
  ## sol path
  Apath = array(NA, c(nlam,p,nrank))
  Vpath = array(NA, c(nlam,q,nrank))
  Cpath = array(NA, c(nlam,pz,q))
  phipath = matrix(NA, nlam,q)
  mupath = array(NA, c(nlam,n,q))
  for(il in 1:nlam){
    if(control$trace==T) cat(il,'...')
    init.i = init 
    if(il!=1 & warm==T) init.i=list(A0=fit$A, V0=fit$V, C0=fit$C, kappa=fit$kappa)
    fit = msrrr.fit(Y, X, Z, family, familygroup, nrank, lamseq[il], init.i, control)
    Apath[il,,] = fit$A
    Vpath[il,,] = fit$V
    Cpath[il,,] = fit$C
    phipath[il,] = fit$phi
    mupath[il,,] = fit$mu
  }
  
  ## tuning
  if(method!='CV'){ 
    xrank = min(n, p) # sum(svd(XX)$d > 0.0001)
    dev = rep(0, nlam) # -2*logL
    df = nz = rep(NA, nlam)
    for(il in 1:nlam){ 
      for(iq in 1:q){   
        # tt= logLikehood(Y[,iq],mupath[il,,iq],phipath[il,iq],fg[[iq]]$family) 
        dev[il] <- dev[il] - 2*sum(logLikehood(Y[,iq],mupath[il,,iq],phipath[il,iq],fg[[iq]]$family), na.rm=T)
        # objFun(Y_test, mu.test, fit1$dispersion)
      }
      # d.f. 
      dfu0 = sum(Apath[il,,] != 0) 
      dfv0 = sum(Vpath[il,,] != 0) 
      df[il] = dfu0 * xrank/p + dfv0 - nrank*nrank + pz*q
      nz[il] = mean(Apath[il,,1] != 0)
    }
    ## IC
    n.obs = sum(!is.na(Y_mis))
    logqn = log(n.obs)
    BIC = (dev + df*logqn)/n.obs 
    BICP = (dev + df*logqn*c.BIC)/n.obs 
    AIC = (dev + 2*df)/n.obs 
    GIC = (dev + df*log(logqn)*log(p*q))/n.obs 
    # dfqn2 = (1 - df/(q*n))^2 
    # GCV = dev/(q*n*dfqn2) 
    ICpath = data.frame(nz, df, dev, BIC, BICP, AIC, GIC)
    lam.idx = switch(method,
                     'BIC'=which.min(BIC),
                     'BICP'=which.min(BICP),
                     'AIC'=which.min(AIC),
                     'GIC'=which.min(GIC)) 
    Tunepath = ICpath
    tunepath.opt = switch(method,
                          'BIC'=min(BIC),
                          'BICP'=min(BICP),
                          'AIC'=min(AIC),
                          'GIC'=min(GIC))
  }else{ 
    # store the deviance/pseudo MSE of the test data
    cv <- matrix(NA, nlam, nfold) 
    if (is.null(foldid)) {  ## CV by row  
      foldid <- rep(1:nfold, len = n)
      foldid <- sample(foldid, n, replace = FALSE)
    }  
    
    init.i=init 
    for (ifold in 1:nfold) { 
      Y_test <- Y[foldid==ifold,]
      X_test <- X[foldid==ifold,]
      Z_test <- Z[foldid==ifold,] 
      Y_train <- Y[foldid!=ifold,]
      X_train <- X[foldid!=ifold,]
      Z_train <- Z[foldid!=ifold,]  
      mu_test <- matrix(NA, sum(foldid==ifold), q)
      for (il in 1:nlam) {
        if(il!=1 & warm==T) init.i =list(A0=fit1$A, V0=fit1$V, C0=fit1$C, kappa=fit1$kappa)
        fit1 <- msrrr.fit(Y_train, X_train, Z_train, family, familygroup, nrank, lamseq[il], init.i, control)
        eta_test = X_test%*%fit1$B + Z_test%*%fit1$C
        for(iq in 1:q){
          fm = fg[[iq]]
          mu_test[,iq] = fm$linkinv(eta_test[,iq])      
        }  
        if(cv.criteria=='deviance') cv[il, ifold] <- 2 * objFun(Y_test, mu_test, fit1$phi, family, familygroup)
        if(cv.criteria=='pMSE') cv[il, ifold] <- pmse(Y_test, mu_test, fit1$phi, family, familygroup)
      }
    }
    cv.mean <- apply(cv, 1, mean)
    lam.idx <- which.min(cv.mean)

    Tunepath = data.frame(cv, cv.mean = cv.mean)
    tunepath.opt = min(cv.mean) 
  }
  
  lam.opt = lamseq[lam.idx] 
  ## refit
  init.i = NULL
  if(lam.idx>1 & warm==T) init.i = list(A0=Apath[lam.idx,,], V0=Vpath[lam.idx,,], C0=Cpath[lam.idx,,])
  fit = msrrr.fit(Y, X, Z, family, familygroup, nrank, lam.opt, init.i, control)
     
  out <- list(lamseq=lamseq, nrank=nrank, Apath=Apath, Vpath=Vpath, Cpath=Cpath,phipath=phipath,mupath=mupath,  
              method=method, cv.criteria=cv.criteria, foldid=foldid, nfold=nfold, c.BIC=c.BIC, # c('CV', 'BIC', 'BICP', 'AIC', 'GIC'), 
              Tunepath = Tunepath, lam.idx = lam.idx, lam.opt = lam.opt, tunepath.opt = tunepath.opt,
              fit = fit, A=fit$A, V=fit$V, B=fit$B, C=fit$C)
  return(out)
}

 
## final wrapper: select nrank  
# Z not include a intercept col
msrrr = function(Y, X, Z=NULL, family, familygroup, nrankseq=c(1:3), init=NULL, 
                  lamseq=NULL, nlam=50, warm=T,   
                  method='CV', cv.criteria='pMSE', foldid=NULL, nfold=5, c.BIC=2, 
                  # c('CV', 'BIC', 'BICP', 'AIC', 'GIC')
                  control=list(epsilon=1e-4, maxit=200, trace=F, conv.obj=T)){
  n = nrow(X)
  Z = cbind(rep(1,n), Z) # 
  # pz = ncol(Z)
  nr = length(nrankseq)
  out.allrank = list()
  
  if(method=='CV' & is.null(foldid)) {   
    foldid <- rep(1:nfold, len = n)
    foldid <- sample(foldid, n, replace = FALSE)
  }  
  for(ir in 1:nr){
    if(control$trace==T) cat(ir)
    out.allrank[[ir]] = msrrr.tuning(Y, X, Z, family, familygroup, nrankseq[ir], init, lamseq, 
                                     nlam, warm, method, cv.criteria, foldid, nfold, c.BIC, control)
  }
  names(out.allrank) = paste0('nrank_', nrankseq)
  Tunepath = lapply(out.allrank, function(a) a$Tunepath)
  tunepath.opt = lapply(out.allrank, function(a) a$tunepath.opt)
  idx = which.min(tunepath.opt)
  out = out.allrank[[idx]]
  
  out$nrankseq = nrankseq
  out$out.allrank = out.allrank
  out$Tunepath = Tunepath
  out$tunepath.opt = tunepath.opt
  out$nrank.opt = nrankseq[idx]
  out$family = family
  out$familygroup = familygroup
  return(out)
}


## ---------------------------------------------------------------------------
## Fold-specific preprocessing for the robust CV workflow
## ---------------------------------------------------------------------------

.msrrr_is_binary_01 <- function(x) {
  observed <- unique(stats::na.omit(x))
  length(observed) > 0L &&
    length(observed) <= 2L &&
    all(observed %in% c(0, 1))
}

.validate_no_scale_columns <- function(M, no_scale, argument_name) {
  M <- as.matrix(M)
  if (is.null(no_scale)) return(character(0))
  if (!is.character(no_scale) || anyNA(no_scale) || any(!nzchar(no_scale))) {
    stop(argument_name, " must be NULL or a character vector of column names.")
  }
  if (is.null(colnames(M))) {
    stop(argument_name, " requires the corresponding matrix to have column names.")
  }
  if (anyDuplicated(no_scale)) no_scale <- unique(no_scale)
  unknown <- setdiff(no_scale, colnames(M))
  if (length(unknown) > 0L) {
    stop(
      argument_name, " contains column(s) not found in the corresponding matrix: ",
      paste(unknown, collapse = ", "), "."
    )
  }
  no_scale
}

.fit_standardization <- function(M, no_scale = NULL) {
  M <- as.matrix(M)
  no_scale <- .validate_no_scale_columns(M, no_scale, "no_scale")
  is_binary <- apply(M, 2L, .msrrr_is_binary_01)
  excluded <- if (is.null(colnames(M))) {
    rep(FALSE, ncol(M))
  } else {
    colnames(M) %in% no_scale
  }
  center <- rep(0, ncol(M))
  scale <- rep(1, ncol(M))
  names(center) <- names(scale) <- colnames(M)

  continuous <- which(!is_binary & !excluded)
  if (length(continuous) > 0L) {
    center[continuous] <- colMeans(M[, continuous, drop = FALSE], na.rm = TRUE)
    scale[continuous] <- apply(
      M[, continuous, drop = FALSE], 2L, stats::sd, na.rm = TRUE
    )
    # Constant Z columns are allowed to reach the existing warning rather than
    # turning the complete transformed column into NA.
    invalid <- !is.finite(scale) | scale <= 0
    scale[invalid] <- 1
  }

  list(
    center = center, scale = scale, is_binary = is_binary,
    no_scale = excluded,
    no_scale_names = colnames(M)[excluded]
  )
}

.apply_standardization <- function(M, parameters) {
  M <- as.matrix(M)
  if (ncol(M) != length(parameters$center)) {
    stop("The matrix does not match the stored preprocessing parameters.")
  }
  out <- sweep(M, 2L, parameters$center, FUN = "-")
  out <- sweep(out, 2L, parameters$scale, FUN = "/")
  if (any(parameters$is_binary)) {
    out[, parameters$is_binary] <- M[, parameters$is_binary, drop = FALSE]
  }
  # `no_scale` was added for analysis-specific reproducibility. Older fitted
  # objects do not contain it and retain their previous behaviour.
  if (!is.null(parameters$no_scale) && any(parameters$no_scale)) {
    out[, parameters$no_scale] <- M[, parameters$no_scale, drop = FALSE]
  }
  dimnames(out) <- dimnames(M)
  out
}

.fit_y_scaling <- function(Y, family, familygroup) {
  Y <- as.matrix(Y)
  spec <- .validate_family_spec(Y, family, familygroup)
  is_gaussian <- vapply(
    spec$fg, function(fm) identical(fm$family, "gaussian"), logical(1)
  )
  center <- rep(0, ncol(Y))
  scale <- rep(1, ncol(Y))
  names(center) <- names(scale) <- colnames(Y)

  for (j in which(is_gaussian)) {
    center[j] <- mean(Y[, j], na.rm = TRUE)
    scale[j] <- stats::sd(Y[, j], na.rm = TRUE)
    if (!is.finite(center[j])) center[j] <- 0
    if (!is.finite(scale[j]) || scale[j] <= 0) scale[j] <- 1
  }

  list(
    center = center, scale = scale, is_gaussian = is_gaussian,
    method = "zscore"
  )
}

.apply_y_scaling <- function(Y, parameters) {
  Y <- as.matrix(Y)
  parameter_length <- if (!is.null(parameters$center)) {
    length(parameters$center)
  } else {
    length(parameters$lower)
  }
  if (ncol(Y) != parameter_length) {
    stop("Y does not match the stored preprocessing parameters.")
  }
  out <- Y
  gaussian <- which(parameters$is_gaussian)
  if (length(gaussian) > 0L) {
    if (!is.null(parameters$center) && !is.null(parameters$scale)) {
      out[, gaussian] <- sweep(
        Y[, gaussian, drop = FALSE],
        2L, parameters$center[gaussian], FUN = "-"
      )
      out[, gaussian] <- sweep(
        out[, gaussian, drop = FALSE],
        2L, parameters$scale[gaussian], FUN = "/"
      )
    } else if (!is.null(parameters$lower) && !is.null(parameters$range)) {
      # Backward compatibility for objects fitted before Gaussian outcomes
      # changed from min-max scaling to z-score standardisation.
      out[, gaussian] <- sweep(
        Y[, gaussian, drop = FALSE],
        2L, parameters$lower[gaussian], FUN = "-"
      )
      out[, gaussian] <- sweep(
        out[, gaussian, drop = FALSE],
        2L, parameters$range[gaussian], FUN = "/"
      )
    } else {
      stop("Stored Y preprocessing parameters are incomplete.")
    }
  }
  dimnames(out) <- dimnames(Y)
  out
}

.fit_msrrr_preprocessing <- function(
    Y, X, Z, family, familygroup,
    X_no_scale = NULL, Z_no_scale = NULL) {
  parameters <- list(
    X = .fit_standardization(X, no_scale = X_no_scale),
    Z = .fit_standardization(Z, no_scale = Z_no_scale),
    Y = .fit_y_scaling(Y, family, familygroup)
  )
  list(
    Y = .apply_y_scaling(Y, parameters$Y),
    X = .apply_standardization(X, parameters$X),
    Z = .apply_standardization(Z, parameters$Z),
    parameters = parameters
  )
}

.apply_msrrr_preprocessing <- function(
    preprocessing, X, Z, Y = NULL) {
  result <- list(
    X = .apply_standardization(X, preprocessing$X),
    Z = .apply_standardization(Z, preprocessing$Z)
  )
  if (!is.null(Y)) {
    result$Y <- .apply_y_scaling(Y, preprocessing$Y)
  }
  result
}

## Robust CV tuning implementation
msrrr_cv.tuning <- function(Y, X, Z, family, familygroup = NULL, nrank=2, init=NULL, 
                        lamseq=NULL, nlam=50, warm=FALSE,  # lam.max, lam.min, 
                        method='CV', cv.criteria='pMSE', foldid=NULL, nfold=5, c.BIC=2, 
                        # c('CV', 'BIC', 'BICP', 'AIC', 'GIC')
                        control=list(epsilon=1e-4, maxit=200, trace=F, conv.obj=T),
                        progress_callback = NULL,
                        X_no_scale = NULL, Z_no_scale = NULL){
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  spec <- .validate_family_spec(Y, family, familygroup)
  familygroup <- spec$familygroup
  fg <- spec$fg

  # msrrr_cv() receives untransformed matrices. Preserve them for the
  # fold-specific transformations, then create the common whole-development
  # transformation used for the solution path and final refit.
  Y_raw <- as.matrix(Y)
  X_raw <- as.matrix(X)
  Z_raw <- .ensure_intercept(Z, n)
  full_preprocessing <- .fit_msrrr_preprocessing(
    Y_raw, X_raw, Z_raw, family, familygroup,
    X_no_scale = X_no_scale, Z_no_scale = Z_no_scale
  )
  Y <- full_preprocessing$Y
  X <- full_preprocessing$X
  Z <- full_preprocessing$Z
  pz = ncol(Z) 
  Y_mis = Y_raw
  id.mis = which(is.na(Y))
  r = nrank  
  
  ## distribution families ##
  nfamily <- length(family)
  ## characters of families
  cfamily <- spec$groups
  
  ## default lambda seq
  if (is.null(lamseq)) {
    message(
      "No lambda sequence supplied. A common sequence is calculated from ",
      "the complete development sample after preprocessing; it is then reused ",
      "unchanged in every CV fold."
    )
    Y_for_lambda <- Y
    for (j in seq_len(ncol(Y_for_lambda))) {
      missing_j <- is.na(Y_for_lambda[, j])
      if (any(missing_j)) {
        Y_for_lambda[missing_j, j] <- mean(
          Y_for_lambda[, j], na.rm = TRUE
        )
      }
    }
    spectral_norm <- svd(
      crossprod(X, Y_for_lambda), nu = 0, nv = 0
    )$d[1]
    # This spectral value is a useful scale estimate, but it is not an exact
    # zero-solution bound for the non-convex mixed-response reduced-rank model.
    # Extend only the strongly penalised end by one decade. The original weak-
    # penalty endpoint and nlam are retained, so this broadens the range without
    # adding model fits or losing the low-lambda region already explored.
    lam.base <- (spectral_norm / n) * 1.1
    if (!is.finite(lam.base) || lam.base <= 0) {
      stop("A positive finite lambda maximum could not be calculated.")
    }
    lam.max <- lam.base * 10
    lam.min <- lam.base * 1e-4
    lamseq <- 10^seq(
      log10(lam.max), log10(lam.min), length.out = nlam
    )
    if(control$trace==T) cat('lambda:', lamseq)
  }
  if (anyNA(lamseq) || any(!is.finite(lamseq)) || any(lamseq < 0)) {
    stop("lamseq must contain finite, non-negative values.")
  }
  # 3) Warm starts always proceed from greater to smaller penalisation.
  if (warm) {
    lamseq <- sort(lamseq, decreasing = TRUE)
  }
  nlam = length(lamseq)
  
  ## sol path
  Apath = array(NA, c(nlam,p,nrank))
  Vpath = array(NA, c(nlam,q,nrank))
  Cpath = array(NA, c(nlam,pz,q))
  phipath = matrix(NA, nlam,q)
  mupath = array(NA, c(nlam,n,q))
  # 14) Store candidate and applied kappa values for later diagnosis.
  Kappapath = matrix(NA_real_, nrow = nlam, ncol = q)
  kappapath = rep(NA_real_, nlam)
  for(il in 1:nlam){
    if(control$trace==T) cat(il,'...')
    init.i = init 
    if(il!=1 & warm==T) init.i=list(A0=fit$A, V0=fit$V, C0=fit$C, kappa=fit$kappa)
    fit = msrrr.fit(Y, X, Z, family, familygroup, nrank, lamseq[il], init.i, control)
    if (!is.null(progress_callback)) progress_callback()
    Apath[il,,] = fit$A
    Vpath[il,,] = fit$V
    Cpath[il,,] = fit$C
    phipath[il,] = fit$phi
    mupath[il,,] = fit$mu
    Kappapath[il,] = fit$Kappa
    kappapath[il] = fit$kappa
  }
  
  ## tuning
  if(method!='CV'){ 
    xrank = min(n, p) # sum(svd(XX)$d > 0.0001)
    dev = rep(0, nlam) # -2*logL
    df = nz = rep(NA, nlam)
    for(il in 1:nlam){ 
      for(iq in 1:q){   
        # 2) dnorm() requires a standard deviation; stored phi is a variance.
        sigma_i <- if (fg[[iq]]$family == "gaussian") {
          sqrt(phipath[il, iq])
        } else {
          1
        }
        dev[il] <- dev[il] - 2*sum(
          logLikehood(
            Y[, iq], mupath[il, , iq], sigma_i, fg[[iq]]$family
          ),
          na.rm = TRUE
        )
        # objFun(Y_test, mu.test, fit1$dispersion)
      }
      # 15) Treat numerical values within the common tolerance as zero.
      A.il <- .as_matrix_dim(
        Apath[il, , ], nrow = p, ncol = nrank
      )
      V.il <- .as_matrix_dim(
        Vpath[il, , ], nrow = q, ncol = nrank
      )

      # d.f. (the formula is retained; only its zero test is harmonised).
      dfu0 = sum(abs(A.il) > .msrrr_selection_tol)
      dfv0 = sum(abs(V.il) > .msrrr_selection_tol)
      df[il] = dfu0 * xrank/p + dfv0 - nrank*nrank + pz*q
      # Exposure support uses every rank component, not only A[, 1].
      nz[il] = mean(sqrt(rowSums(A.il^2)) > .msrrr_selection_tol)
    }
    ## IC
    n.obs = sum(!is.na(Y_mis))
    logqn = log(n.obs)
    BIC = (dev + df*logqn)/n.obs 
    BICP = (dev + df*logqn*c.BIC)/n.obs 
    AIC = (dev + 2*df)/n.obs 
    GIC = (dev + df*log(logqn)*log(p*q))/n.obs 
    # dfqn2 = (1 - df/(q*n))^2 
    # GCV = dev/(q*n*dfqn2) 
    ICpath = data.frame(nz, df, dev, BIC, BICP, AIC, GIC)
    lam.idx = switch(method,
                     'BIC'=which.min(BIC),
                     'BICP'=which.min(BICP),
                     'AIC'=which.min(AIC),
                     'GIC'=which.min(GIC)) 
    Tunepath = ICpath
    tunepath.opt = switch(method,
                          'BIC'=min(BIC),
                          'BICP'=min(BICP),
                          'AIC'=min(AIC),
                          'GIC'=min(GIC))
  }else{ 
    # store the deviance/pseudo MSE of the test data
    cv <- matrix(NA, nlam, nfold) 
    cv.Kappapath <- array(
      NA_real_, dim = c(nlam, nfold, q),
      dimnames = list(lambda = NULL, fold = NULL, outcome = NULL)
    )
    cv.kappapath <- matrix(
      NA_real_, nrow = nlam, ncol = nfold,
      dimnames = list(lambda = NULL, fold = NULL)
    )
    if (is.null(foldid)) {  ## CV by row  
      foldid <- rep(1:nfold, len = n)
      foldid <- sample(foldid, n, replace = FALSE)
    }  
    
    for (ifold in 1:nfold) {
      # 4) Reset each CV fold to prevent leakage between training folds.
      # Warm starts are reused only across lambdas within the same fold.
      init.i <- init
      fit1 <- NULL

      # Fit every transformation using only the four training folds. Apply the
      # resulting parameters unchanged to the held-out validation fold.
      training_rows <- foldid != ifold
      validation_rows <- !training_rows
      fold_preprocessing <- .fit_msrrr_preprocessing(
        Y_raw[training_rows, , drop = FALSE],
        X_raw[training_rows, , drop = FALSE],
        Z_raw[training_rows, , drop = FALSE],
        family, familygroup,
        X_no_scale = X_no_scale, Z_no_scale = Z_no_scale
      )
      validation_data <- .apply_msrrr_preprocessing(
        fold_preprocessing$parameters,
        X = X_raw[validation_rows, , drop = FALSE],
        Z = Z_raw[validation_rows, , drop = FALSE],
        Y = Y_raw[validation_rows, , drop = FALSE]
      )
      Y_train <- fold_preprocessing$Y
      X_train <- fold_preprocessing$X
      Z_train <- fold_preprocessing$Z
      Y_test <- validation_data$Y
      X_test <- validation_data$X
      Z_test <- validation_data$Z
      mu_test <- matrix(NA, sum(foldid==ifold), q)
      for (il in 1:nlam) {
        if(il!=1 & warm==T) init.i =list(A0=fit1$A, V0=fit1$V, C0=fit1$C, kappa=fit1$kappa)
        fit1 <- msrrr.fit(Y_train, X_train, Z_train, family, familygroup, nrank, lamseq[il], init.i, control)
        if (!is.null(progress_callback)) progress_callback()
        cv.Kappapath[il, ifold, ] <- fit1$Kappa
        cv.kappapath[il, ifold] <- fit1$kappa
        eta_test = X_test%*%fit1$B + Z_test%*%fit1$C
        for(iq in 1:q){
          fm = fg[[iq]]
          mu_test[,iq] = fm$linkinv(eta_test[,iq])      
        }  
        if(cv.criteria=='deviance') cv[il, ifold] <- 2 * objFun(Y_test, mu_test, fit1$phi, family, familygroup)
        if(cv.criteria=='pMSE') cv[il, ifold] <- pmse(Y_test, mu_test, fit1$phi, family, familygroup)
      }
    }
    cv.mean <- apply(cv, 1, mean)
    lam.idx <- which.min(cv.mean)
    
    Tunepath = data.frame(cv, cv.mean = cv.mean)
    tunepath.opt = min(cv.mean) 
  }
  
  lam.opt = lamseq[lam.idx] 
  
  ## refit
  # 5) Preserve a supplied init when no path warm start is available.
  refit.init <- init
  
  if (lam.idx > 1 & warm == T) {
    # 1) Restore A, V and C matrix shapes before the final refit.
    A0 <- .as_matrix_dim(Apath[lam.idx, , ], nrow = p, ncol = nrank)
    V0 <- .as_matrix_dim(Vpath[lam.idx, , ], nrow = q, ncol = nrank)
    C0 <- .as_matrix_dim(Cpath[lam.idx, , ], nrow = pz, ncol = q)
    
    refit.init <- list(A0 = A0, V0 = V0, C0 = C0)
  }
  
  fit = msrrr.fit(
    Y, X, Z, family, familygroup,
    nrank, lam.opt, refit.init, control
  )
  if (!is.null(progress_callback)) progress_callback()

  # 5) Expose the exact matrix-valued initialisation actually used.
  refit.init <- fit$init_used
  
  out <- list(lamseq=lamseq, nrank=nrank, Apath=Apath, Vpath=Vpath, Cpath=Cpath,
              phipath=phipath, mupath=mupath, Kappapath=Kappapath,
              kappapath=kappapath,
              method=method, cv.criteria=cv.criteria, foldid=foldid, nfold=nfold, c.BIC=c.BIC, # c('CV', 'BIC', 'BICP', 'AIC', 'GIC'), 
              Tunepath = Tunepath, lam.idx = lam.idx, lam.opt = lam.opt, tunepath.opt = tunepath.opt,
              refit.init = refit.init, control = control, family = family, familygroup = familygroup,
              X_no_scale = X_no_scale, Z_no_scale = Z_no_scale,
              preprocessing = full_preprocessing$parameters,
              fit = fit, A=fit$A, V=fit$V, B=fit$B, C=fit$C)
  if (method == "CV") {
    out$cv.Kappapath <- cv.Kappapath
    out$cv.kappapath <- cv.kappapath
  }
  # 9) Enable S3 methods such as predict(object, ...).
  class(out) <- "msrrr"
  out
}


## final wrapper: select nrank  
# Z may be supplied with or without an intercept column
.msrrr_cv_impl <- function(Y, X, Z=NULL, family, familygroup = NULL, nrankseq=c(1:3), init=NULL, 
                 lamseq=NULL, nlam=50, warm=FALSE,   
                 method='CV', cv.criteria='pMSE', foldid=NULL, nfold=5, c.BIC=2, 
                 # c('CV', 'BIC', 'BICP', 'AIC', 'GIC')
                 control=list(epsilon=1e-4, maxit=200, trace=F, conv.obj=T),
                 diagnostics=TRUE, progress=TRUE,
                 X_no_scale = NULL, Z_no_scale = NULL){
  n = nrow(X)
  spec <- .validate_family_spec(Y, family, familygroup)
  familygroup <- spec$familygroup

  if (isTRUE(diagnostics)) {
    # Practical check A) Warn, but do not stop, for infeasible requested ranks.
    .diagnose_requested_ranks(
      nrankseq, p = ncol(X), q = ncol(Y)
    )
    .diagnose_msrrr_inputs(
      Y, X, family, familygroup,
      method = method, cv.criteria = cv.criteria
    )
  }

  # 6) Add an intercept only if it is not already present.
  Z <- .ensure_intercept(Z, n)
  X_no_scale <- .validate_no_scale_columns(X, X_no_scale, "X_no_scale")
  Z_no_scale <- .validate_no_scale_columns(Z, Z_no_scale, "Z_no_scale")

  # pz = ncol(Z)
  nr = length(nrankseq)
  out.allrank = list()

  # Display progress by completed model fits. This is intentionally separate
  # from control$trace, which prints verbose optimisation details for every
  # internal iteration and is mainly useful for low-level debugging.
  progress_callback <- NULL
  if (isTRUE(progress)) {
    path_length <- if (is.null(lamseq)) nlam else length(lamseq)
    fits_per_rank <- path_length + 1L
    if (method == "CV") {
      fits_per_rank <- path_length * (nfold + 1L) + 1L
    }
    progress_total <- nr * fits_per_rank
    progress_done <- 0L
    progress_bar <- utils::txtProgressBar(
      min = 0, max = progress_total, style = 3
    )
    on.exit(close(progress_bar), add = TRUE)
    progress_callback <- function() {
      progress_done <<- progress_done + 1L
      utils::setTxtProgressBar(
        progress_bar, min(progress_done, progress_total)
      )
    }
  }
  
  if(method=='CV' & is.null(foldid)) {   
    foldid <- rep(1:nfold, len = n)
    foldid <- sample(foldid, n, replace = FALSE)
  }
  if (method == "CV") {
    if (length(foldid) != n || anyNA(foldid)) {
      stop("foldid must contain one non-missing fold identifier per row.")
    }
    observed_folds <- sort(unique(foldid))
    if (length(observed_folds) != nfold) {
      stop(
        "nfold does not match the number of distinct identifiers in foldid."
      )
    }
  }
  if (isTRUE(diagnostics) && method == "CV") {
    .diagnose_dummy_folds(X, foldid)
  }
  for(ir in 1:nr){
    if(control$trace==T) cat(ir)
    out.allrank[[ir]] = msrrr_cv.tuning(
      Y, X, Z, family, familygroup, nrankseq[ir], init, lamseq,
      nlam, warm, method, cv.criteria, foldid, nfold, c.BIC, control,
      progress_callback = progress_callback,
      X_no_scale = X_no_scale, Z_no_scale = Z_no_scale
    )
  }
  names(out.allrank) = paste0('nrank_', nrankseq)
  Tunepath = lapply(out.allrank, function(a) a$Tunepath)
  tunepath.opt = lapply(out.allrank, function(a) a$tunepath.opt)
  idx = which.min(tunepath.opt)
  out = out.allrank[[idx]]
  
  out$nrankseq = nrankseq
  out$out.allrank = out.allrank
  out$Tunepath = Tunepath
  out$tunepath.opt = tunepath.opt
  out$nrank.opt = nrankseq[idx]
  out$family = family
  out$familygroup = familygroup
  out$call <- match.call()
  # 9) Preserve the class after selecting the optimal rank.
  class(out) <- "msrrr"
  if (isTRUE(diagnostics)) {
    .diagnose_msrrr_selection(
      out, nrankseq = nrankseq, p = ncol(X), q = ncol(Y)
    )
  }
  out
}

## Robust public wrapper: fold-specific preprocessing and self-contained
## diagnostics. X, Z and Y must be supplied on their original, untransformed
## scales. When diagnostics=TRUE, every captured message/warning is retained in
## object$diagnostics. A compact audit is printed once, after fitting, by
## default. Supplying diagnostics_file additionally writes that table to disk.
## X_no_scale and Z_no_scale are optional character vectors containing exact
## column names that must remain on their original scale. NULL preserves the
## standard v4 preprocessing rule. For historical HELIX reproduction, derive
## these vectors before fitting using the <7-distinct-values rule; this column
## classification is then held fixed across all CV folds.
msrrr_cv <- function(
    Y, X, Z = NULL, family, familygroup = NULL, nrankseq = c(1:3),
    init = NULL, lamseq = NULL, nlam = 50, warm = FALSE,
    method = "CV", cv.criteria = "pMSE", foldid = NULL, nfold = 5,
    c.BIC = 2,
    control = list(epsilon = 1e-4, maxit = 200, trace = FALSE,
                   conv.obj = TRUE),
    diagnostics = TRUE,
    print_diagnostics_at_end = TRUE,
    diagnostics_file = NULL,
    progress = TRUE,
    X_no_scale = NULL,
    Z_no_scale = NULL) {
  diagnostic_log <- data.frame(
    condition = character(),
    text = character(),
    stringsAsFactors = FALSE
  )

  fit_call <- function() {
    .msrrr_cv_impl(
      Y = Y, X = X, Z = Z,
      family = family, familygroup = familygroup,
      nrankseq = nrankseq, init = init,
      lamseq = lamseq, nlam = nlam, warm = warm,
      method = method, cv.criteria = cv.criteria,
      foldid = foldid, nfold = nfold, c.BIC = c.BIC,
      control = control, diagnostics = diagnostics, progress = progress,
      X_no_scale = X_no_scale, Z_no_scale = Z_no_scale
    )
  }

  if (isTRUE(diagnostics)) {
    object <- withCallingHandlers(
      fit_call(),
      message = function(m) {
        diagnostic_log <<- rbind(
          diagnostic_log,
          data.frame(
            condition = "message",
            text = conditionMessage(m),
            stringsAsFactors = FALSE
          )
        )
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        diagnostic_log <<- rbind(
          diagnostic_log,
          data.frame(
            condition = "warning",
            text = conditionMessage(w),
            stringsAsFactors = FALSE
          )
        )
        invokeRestart("muffleWarning")
      }
    )
    diagnostic_log <- unique(diagnostic_log)
  } else {
    object <- fit_call()
  }

  object$diagnostics <- diagnostic_log
  object$call <- match.call()

  if (isTRUE(diagnostics) && isTRUE(print_diagnostics_at_end)) {
    message("\n>>> END-OF-RUN msRRR CV DIAGNOSTIC AUDIT")
    if (nrow(diagnostic_log) == 0L) {
      message("No diagnostic messages or warnings were generated.")
    } else {
      for (i in seq_len(nrow(diagnostic_log))) {
        message(
          "[", toupper(diagnostic_log$condition[i]), "] ",
          trimws(diagnostic_log$text[i])
        )
      }
    }
  }

  if (!is.null(diagnostics_file)) {
    if (!isTRUE(diagnostics)) {
      warning(
        "diagnostics_file was supplied but diagnostics=FALSE; no file was written."
      )
    } else {
      if (!is.character(diagnostics_file) ||
          length(diagnostics_file) != 1L ||
          is.na(diagnostics_file) ||
          !nzchar(diagnostics_file)) {
        stop("diagnostics_file must be NULL or one non-empty file path.")
      }
      output_directory <- dirname(diagnostics_file)
      if (!dir.exists(output_directory)) {
        dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
      }
      if (!dir.exists(output_directory)) {
        stop("The diagnostics output directory could not be created.")
      }
      utils::write.csv2(
        diagnostic_log, file = diagnostics_file, row.names = FALSE
      )
    }
  }

  object
}

transform_msrrr_data <- function(object, X, Z = NULL, Y = NULL) {
  if (is.null(object$preprocessing)) {
    stop("object does not contain preprocessing parameters from msrrr_cv().")
  }
  X <- as.matrix(X)
  Z <- .ensure_intercept(Z, nrow(X))
  .apply_msrrr_preprocessing(
    object$preprocessing, X = X, Z = Z, Y = Y
  )
}

predict.msrrr <- function(object, X.new, Z.new = NULL, Y.new = NULL,
                          family = NULL, familygroup = NULL,
                          type = "response", cv.criteria = "pMSE") {
  # 11) Validate prediction options and allow response- or link-scale output.
  type <- match.arg(type, c("response", "link"))
  cv.criteria <- match.arg(cv.criteria, c("pMSE", "deviance"))
  X.new <- as.matrix(X.new)
  q = ncol(object$B)
  n = nrow(X.new)
  if (is.null(family)) family = object$family
  if (is.null(familygroup)) familygroup = object$familygroup
  Y.for.validation <- if (is.null(Y.new)) {
    matrix(NA_real_, nrow = n, ncol = q)
  } else {
    Y.new <- as.matrix(Y.new)
    if (nrow(Y.new) != n || ncol(Y.new) != q) {
      stop("Y.new must have one row per prediction and one column per outcome.")
    }
    Y.new
  }
  spec <- .validate_family_spec(
    Y.for.validation, family, familygroup,
    check_response = !is.null(Y.new)
  )
  familygroup <- spec$familygroup
  fg <- spec$fg

  # Apply the same intercept handling and, for msrrr_cv() objects, the exact
  # whole-development preprocessing learned during the final refit.
  Z.new <- .ensure_intercept(Z.new, n)
  if (!is.null(object$preprocessing)) {
    transformed <- .apply_msrrr_preprocessing(
      object$preprocessing, X = X.new, Z = Z.new, Y = Y.new
    )
    X.new <- transformed$X
    Z.new <- transformed$Z
    if (!is.null(Y.new)) Y.new <- transformed$Y
  }
  if (ncol(X.new) != nrow(object$B)) {
    stop("X.new does not have the number of columns expected by the fitted model.")
  }
  if (ncol(Z.new) != nrow(object$C)) {
    stop(
      "Z.new does not have the number of covariate columns expected by the ",
      "fitted model (including the intercept)."
    )
  }
  eta.new = X.new %*% object$B + Z.new %*% object$C
  mu.new = matrix(NA, n, q)
  
  for (iq in 1:q) {
    fm = fg[[iq]]
    mu.new[, iq] = fm$linkinv(eta.new[, iq])
  }
  
  pred.perf = NULL
  if (!is.null(Y.new)) {
    phi <- if (!is.null(object$fit$phi)) object$fit$phi else object$phi

    if (cv.criteria == "pMSE") {
      pred.perf <- pmse(Y.new, mu.new, phi, family, familygroup)
    } else {
      pred.perf <- 2 * objFun(Y.new, mu.new, phi, family, familygroup)
    }
  }
  
  pred <- if (type == "link") eta.new else mu.new
  return(list(fit = pred, pred.perf = pred.perf))
}


## 5) Reproduce the final msrrr() refit with the same lambda, rank,
## intercept handling, control settings and initial values.
## This is the recommended comparison with object$fit.
refit_msrrr <- function(object, Y, X, Z = NULL, init = object$refit.init,
                        control = object$control) {
  if (is.null(object$nrank) || is.null(object$lam.opt)) {
    stop("object must be returned by msrrr() or msrrr.tuning().")
  }
  if (is.null(object$family) || is.null(object$familygroup)) {
    stop("family and familygroup are missing from object.")
  }
  if (is.null(control)) {
    control <- list(epsilon = 1e-4, maxit = 200, trace = FALSE)
  }

  preprocessing <- NULL
  ensure_intercept <- TRUE
  if (!is.null(object$preprocessing)) {
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    Z <- .ensure_intercept(Z, nrow(Y))
    X_no_scale <- object$preprocessing$X$no_scale_names
    Z_no_scale <- object$preprocessing$Z$no_scale_names
    transformed <- .fit_msrrr_preprocessing(
      Y, X, Z, object$family, object$familygroup,
      X_no_scale = X_no_scale, Z_no_scale = Z_no_scale
    )
    Y <- transformed$Y
    X <- transformed$X
    Z <- transformed$Z
    preprocessing <- transformed$parameters
    ensure_intercept <- FALSE
  }

  fit <- msrrr.fit(
    Y = Y,
    X = X,
    Z = Z,
    family = object$family,
    familygroup = object$familygroup,
    nrank = object$nrank,
    lambda = object$lam.opt,
    init = init,
    control = control,
    ensure_intercept = ensure_intercept
  )

  # Preserve the tuning metadata needed by optional downstream refits.
  fit$nrank <- object$nrank
  fit$lambda <- object$lam.opt
  fit$family <- object$family
  fit$familygroup <- object$familygroup
  fit$control <- control
  fit$preprocessing <- preprocessing
  class(fit) <- unique(c("msrrr_refit", class(fit)))
  fit
}


## 16) Optional post-selection refit without penalisation.
## The selected exposure set and rank are treated as fixed for this refit.
refit_msrrr_unpenalized <- function(
    object, Y, X, Z = NULL,
    tol = .msrrr_selection_tol,
    nrank = object$nrank,
    init = NULL,
    control = object$control) {

  penalized_fit <- if (!is.null(object$B)) {
    object
  } else if (!is.null(object$fit) && !is.null(object$fit$B)) {
    object$fit
  } else {
    stop(
      "object must be returned by refit_msrrr(), msrrr(), or msrrr.tuning()."
    )
  }

  object_family <- if (!is.null(object$family)) {
    object$family
  } else {
    penalized_fit$family
  }
  object_familygroup <- if (!is.null(object$familygroup)) {
    object$familygroup
  } else {
    penalized_fit$familygroup
  }
  if (is.null(object_family) || is.null(object_familygroup)) {
    stop("family and familygroup are missing from object.")
  }
  if (length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol < 0) {
    stop("tol must be one finite non-negative value.")
  }

  Y <- as.matrix(Y)
  X <- as.matrix(X)
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.")
  }
  if (ncol(X) != nrow(penalized_fit$B)) {
    stop(
      paste(
        "X must contain the same exposure columns, in the same order, as",
        "the penalised msrrr fit."
      )
    )
  }
  if (!is.null(Z) && nrow(as.matrix(Z)) != nrow(Y)) {
    stop("Z and Y must have the same number of rows.")
  }

  selected <- sqrt(rowSums(penalized_fit$B^2)) > tol
  if (!any(selected)) {
    reason <- paste(
      "The penalised model selected no exposures at tolerance", tol,
      "; the conditional unpenalised refit is not applicable."
    )
    message("    ! ", reason)
    exposure_names <- colnames(X)
    if (is.null(exposure_names)) {
      exposure_names <- paste0("X", seq_len(ncol(X)))
    }
    out <- list(
      fit = NULL,
      X_selected = X[, FALSE, drop = FALSE],
      selected = selected,
      selected_indices = integer(0),
      selected_exposures = character(0),
      requested_rank = as.integer(nrank)[1L],
      nrank = NA_integer_,
      lambda = 0,
      selection_tol = tol,
      family = object_family,
      familygroup = object_familygroup,
      control = control,
      preprocessing = NULL,
      init = init,
      penalized_lambda = if (!is.null(object$lambda)) object$lambda else object$lam.opt,
      available = FALSE,
      status = "not_applicable_no_selected_exposures",
      reason = reason
    )
    class(out) <- "msrrr_unpenalized_refit"
    return(out)
  }

  X_selected <- X[, selected, drop = FALSE]
  requested_rank <- nrank
  if (length(requested_rank) != 1L || is.na(requested_rank) ||
      !is.finite(requested_rank) || requested_rank < 1 ||
      requested_rank %% 1 != 0) {
    stop("nrank must be one positive integer.")
  }
  requested_rank <- as.integer(requested_rank)
  maximum_rank <- min(ncol(X_selected), ncol(Y))
  refit_rank <- min(requested_rank, maximum_rank)

  if (refit_rank < requested_rank) {
    warning(
      "The selected rank ", requested_rank,
      " exceeds the maximum feasible rank after exposure selection (",
      maximum_rank, "). The unpenalised refit uses rank ", refit_rank, ".",
      call. = FALSE
    )
  }

  if (is.null(control)) {
    control <- list(epsilon = 1e-4, maxit = 200, trace = FALSE)
  }

  preprocessing <- NULL
  ensure_intercept <- TRUE
  if (!is.null(object$preprocessing)) {
    Z <- .ensure_intercept(Z, nrow(Y))
    stored_X_no_scale <- object$preprocessing$X$no_scale_names
    X_no_scale <- intersect(stored_X_no_scale, colnames(X_selected))
    Z_no_scale <- object$preprocessing$Z$no_scale_names
    transformed <- .fit_msrrr_preprocessing(
      Y, X_selected, Z, object_family, object_familygroup,
      X_no_scale = X_no_scale, Z_no_scale = Z_no_scale
    )
    Y <- transformed$Y
    X_selected <- transformed$X
    Z <- transformed$Z
    preprocessing <- transformed$parameters
    ensure_intercept <- FALSE
  }

  fit <- msrrr.fit(
    Y = Y,
    X = X_selected,
    Z = Z,
    family = object_family,
    familygroup = object_familygroup,
    nrank = refit_rank,
    lambda = 0,
    init = init,
    control = control,
    ensure_intercept = ensure_intercept
  )

  exposure_names <- colnames(X)
  if (is.null(exposure_names)) {
    exposure_names <- paste0("X", seq_len(ncol(X)))
  }
  selected_exposures <- exposure_names[selected]
  rownames(fit$B) <- selected_exposures
  if (!is.null(colnames(Y))) colnames(fit$B) <- colnames(Y)

  out <- list(
    fit = fit,
    X_selected = X_selected,
    selected = selected,
    selected_indices = which(selected),
    selected_exposures = selected_exposures,
    requested_rank = requested_rank,
    nrank = refit_rank,
    lambda = 0,
    selection_tol = tol,
    family = object_family,
    familygroup = object_familygroup,
    control = control,
    preprocessing = preprocessing,
    init = init,
    penalized_lambda = if (!is.null(object$lambda)) {
      object$lambda
    } else {
      object$lam.opt
    },
    available = TRUE,
    status = "fitted",
    reason = NA_character_
  )
  class(out) <- "msrrr_unpenalized_refit"
  out
}
