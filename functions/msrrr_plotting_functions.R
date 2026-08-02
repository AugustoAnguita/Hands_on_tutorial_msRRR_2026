## =============================================================================
## msRRR plotting functions
## VERSION: v1 (31 July 2026)
## =============================================================================
##
## Reusable plotting helpers for the msRRR tutorial. Keeping these functions in
## one sourced file ensures that the R script and R Markdown tutorial use the
## same palettes, dimensions, ordering and filtering rules.
##
## Statistical fitting belongs to msrrr_v4.R and bootstrap estimation belongs
## to bootstrapping_functions_v_31_07_26.r. This file contains presentation
## and graphical-export functions only.
##
## plot_final_heatmap() is shared by CV and whole-sample models. It accepts
## coefficient matrices stored either in object$fit$B or object$B, disables
## clustering, fixes cell dimensions and uses the common numerical tolerance.
## =============================================================================

generate_corrplot <- function(
    data_mat, plot_title, file_name,
    label_cex = 0.6, number_cex = 0.4,
    display_in_current_device = FALSE) {

  message(paste("    -> Processing:", plot_title))

  # Calculate Correlation (Pairwise Complete)
  M <- cor(data_mat, use = "pairwise.complete.obs")

  col_palette <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

  draw_corrplot <- function() {
    corrplot::corrplot(
      M,
      method = "color",
      col = col_palette(200),
      type = "upper",
      order = "hclust",
      addCoef.col = "black",
      tl.col = "black",
      tl.cex = label_cex,
      number.cex = number_cex,
      diag = FALSE,
      title = plot_title,
      mar = c(0, 0, 2, 0)
    )
  }

  # Save PDF.
  grDevices::pdf(file = file_name, width = 12, height = 10)
  draw_corrplot()
  dev.off()

  message(paste("       Saved to:", file_name))
  if (display_in_current_device) draw_corrplot()
  invisible(M)
}


generate_diagnostics <- function(model_obj, scenario_name, cv.crit) {

  # 1. Rank vs selected CV metric plot
  rank_data <- data.frame(
    Rank = model_obj$nrankseq,
    CV_metric = unlist(model_obj$tunepath.opt)
  )

  p_rank <- ggplot(rank_data, aes(x = Rank, y = CV_metric)) +
    geom_line(color = "grey") + 
    geom_point(size = 3, color = "#08737f") +
    theme_minimal() + 
    labs(
      title = paste(scenario_name, "Model: Rank vs", cv.crit),
      subtitle = "Optimization of latent factors",
      x = "Rank",
      y = paste("CV", cv.crit)
    )

  # 2. Lambda boxplot for the optimal rank
  opt_rank_name <- paste0("nrank_", model_obj$nrank)

  # Extract fold data, removing the summary column
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
    theme(
      axis.text.x = element_text(angle = 90, size = 6, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      x = "Lambda (log scale)",
      y = paste("CV", cv.crit),
      title = paste("Regularization path:", cv.crit, "across folds")
    )

  # Return as grob, without drawing immediately
  return(gridExtra::arrangeGrob(p_rank, p_lambda, ncol = 1))
}


plot_coeffs <- function(model_obj, outcome_idx, X_matrix, scenario_label) {
  toplot <- data.frame(
    Coefficient = model_obj$fit$B[, outcome_idx],
    Exposure    = colnames(X_matrix)
  )

  toplot$vargroup <- ifelse(
    toplot$Coefficient > selection_tol, "Positive",
    ifelse(
      toplot$Coefficient < -selection_tol, "Negative", "Null"
    )
  )

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


get_pro_breaks <- function(mat) {
  limit <- max(abs(mat), na.rm = TRUE)
  if(limit == 0) limit <- 0.1

  # Build sequences from each extreme to a value very close to zero.
  # This gives the central grey colour an almost invisible interval, ensuring
  # that every non-zero value receives either a green or a red colour.
  eps <- 1e-99
  breaks <- c(
    seq(-limit, -eps, length.out = 51), # 50 negative intervals
    0,                                  # Exact zero
    seq(eps, limit, length.out = 51)    # 50 positive intervals
  )
  return(unique(sort(breaks)))
}


plot_final_heatmap <- function(
    model_obj, X_matrix, scenario_name,
    outcome_names = NULL, vars_info_table = NULL,
    pregnancy_exposure_names = NULL,
    tol = get0("selection_tol", ifnotfound = 1e-8, inherits = TRUE),
    draw = TRUE) {

  # Accept both msrrr_cv()/msrrr() objects ($fit$B) and direct refits ($B).
  if (!is.null(model_obj$fit)) {
    B_mat <- model_obj$fit$B
  } else {
    B_mat <- model_obj$B
  }
  if (is.null(B_mat)) {
    stop("model_obj does not contain a coefficient matrix in $fit$B or $B.")
  }

  rownames(B_mat) <- colnames(X_matrix)
  if (is.null(outcome_names)) {
    outcome_names <- colnames(B_mat)
    if (is.null(outcome_names) && exists("Y_train", inherits = TRUE)) {
      outcome_names <- colnames(get("Y_train", inherits = TRUE))
    }
  }
  if (length(outcome_names) != ncol(B_mat)) {
    stop("outcome_names must contain one name per coefficient-matrix column.")
  }
  colnames(B_mat) <- outcome_names

  # Retain exposures with at least one coefficient above the shared tolerance.
  active_rows <- rowSums(abs(B_mat) > tol) > 0
  coef_subset <- B_mat[active_rows, , drop = FALSE]

  if (nrow(coef_subset) == 0) {
    message("    ! No coefficients available for the heatmap: ", scenario_name)
    empty_grob <- grid::textGrob(
      paste0(
        scenario_name,
        "\nNo exposure was selected at tolerance ",
        format(tol, scientific = TRUE), "."
      ),
      gp = grid::gpar(fontsize = 13)
    )
    if (draw) grid::grid.draw(empty_grob)
    return(invisible(list(gtable = empty_grob, matrix = coef_subset)))
  }

  # Build exposure-family annotations.
  if (is.null(vars_info_table)) {
    vars_info_table <- get0("vars_info", ifnotfound = NULL, inherits = TRUE)
  }
  if (is.null(vars_info_table)) {
    stop("vars_info_table must be supplied (or vars_info must exist).")
  }
  clean_names <- gsub("_None|_Ter_2|_Ter_3", "", rownames(coef_subset))
  row_info <- vars_info_table %>% 
    dplyr::filter(variable_name %in% clean_names) %>%
    dplyr::distinct(variable_name, .keep_all = TRUE)

  df_row_ann <- data.frame(
    Family = row_info$family[match(clean_names, row_info$variable_name)],
    row.names = rownames(coef_subset)
  )

  # Add the life-period annotation to Combined CV and whole-sample heatmaps.
  if (grepl("^Combined", scenario_name)) {
    if (is.null(pregnancy_exposure_names)) {
      pregnancy_matrix_name <- if (grepl("Whole Sample", scenario_name)) {
        "X_preg_whole"
      } else {
        "X_preg_train"
      }
      pregnancy_matrix <- get0(
        pregnancy_matrix_name, ifnotfound = NULL, inherits = TRUE
      )
      if (is.null(pregnancy_matrix)) {
        stop(
          "pregnancy_exposure_names must be supplied for a Combined heatmap."
        )
      }
      pregnancy_exposure_names <- colnames(pregnancy_matrix)
    }
    df_row_ann$Period <- ifelse(
      rownames(coef_subset) %in% pregnancy_exposure_names,
      "Pregnancy", "Childhood"
    )
  }

  # Use separate, explicit palettes for exposure family and life period. This
  # prevents (for example) Childhood and PFAS receiving the same annotation
  # colour merely because pheatmap generated both palettes automatically.
  observed_families <- unique(stats::na.omit(df_row_ann$Family))
  annotation_colors <- list(
    Family = stats::setNames(
      grDevices::hcl.colors(
        max(3L, length(observed_families)), palette = "Dark 3"
      )[seq_along(observed_families)],
      observed_families
    )
  )
  if ("Period" %in% names(df_row_ann)) {
    annotation_colors$Period <- c(
      Pregnancy = "#6A3D9A",
      Childhood = "#00A6D6"
    )
  }

  p <- pheatmap::pheatmap(
    coef_subset,
    main = paste("Environmental Signature:", scenario_name),
    color = custom_colors_pro,
    breaks = get_pro_breaks(coef_subset),
    annotation_row = df_row_ann,
    annotation_colors = annotation_colors,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color = "white",
    fontsize_row = 8,
    cellwidth = 32,
    cellheight = 13,
    silent = TRUE
  )

  if (draw) grid::grid.draw(p$gtable)
  invisible(p)
}


plot_penalized_loco_summary <- function(
  loco_summary,
  scenario_name,
  selection_definition = "selected in the penalised LOCO fit",
  method_label = "penalised LOCO",
  tol = selection_tol,
  rows_to_plot = NULL,
  cellwidth = 32,
  cellheight = 13,
  draw = TRUE) {

  n_selected <- loco_summary$n_selected
  n_positive <- loco_summary$n_positive
  n_negative <- loco_summary$n_negative
  if (is.null(rows_to_plot)) {
    rows_to_plot <- rowSums(n_selected > 0L) > 0L
  } else {
    if (is.character(rows_to_plot)) {
      rows_to_plot <- rownames(n_selected) %in% rows_to_plot
    }
    if (!is.logical(rows_to_plot) || length(rows_to_plot) != nrow(n_selected)) {
      stop("rows_to_plot must identify rows of the LOCO summary matrices.")
    }
    rows_to_plot[is.na(rows_to_plot)] <- FALSE
  }

  if (!any(rows_to_plot)) {
    empty_grob <- grid::textGrob(
      paste0(
        scenario_name,
        ": no exposure-outcome association was ", selection_definition,
        " in any LOCO analysis."
      )
    )
    if (draw) grid::grid.draw(empty_grob)
    return(invisible(empty_grob))
  }

  n_selected <- n_selected[rows_to_plot, , drop = FALSE]
  n_positive <- n_positive[rows_to_plot, , drop = FALSE]
  n_negative <- n_negative[rows_to_plot, , drop = FALSE]

  score <- matrix(
    0,
    nrow = nrow(n_selected),
    ncol = ncol(n_selected),
    dimnames = dimnames(n_selected)
  )
  positive_dominant <- n_positive > n_negative
  negative_dominant <- n_negative > n_positive
  sign_tie <- n_positive == n_negative & n_selected > 0L
  score[positive_dominant] <- n_positive[positive_dominant]
  score[negative_dominant] <- -n_negative[negative_dominant]
  score[sign_tie] <- 0.5

  negative_colors <- RColorBrewer::brewer.pal(9, "Greens")[4:9]
  positive_colors <- RColorBrewer::brewer.pal(9, "Reds")[4:9]
  # Zero/not selected is grey; a selected coefficient with tied direction is
  # yellow. These special values occupy the intervals around 0 and 0.5.
  heatmap_colors <- c(
    rev(negative_colors), "#D9D9D9", "#FFF7BC", positive_colors
  )
  heatmap_breaks <- c(
    seq(-6.5, -0.5, by = 1), 0.25, 0.75, seq(1.5, 6.5, by = 1)
  )

  title_text <- paste0(
    scenario_name, " — ", method_label, "\n",
    selection_definition, "; numbers = count out of ",
    length(cohort_levels)
  )
  # Rebuild with ASCII punctuation so titles render identically on Windows,
  # Linux and in knitted HTML documents.
  title_text <- paste0(
    scenario_name, " - ", method_label, "\n",
    selection_definition, "; numbers = count out of ",
    length(cohort_levels)
  )
  p <- pheatmap::pheatmap(
    score,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    display_numbers = n_selected,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "white",
    fontsize_row = 7,
    fontsize_col = 8,
    fontsize_number = 7,
    cellwidth = cellwidth,
    cellheight = cellheight,
    main = title_text,
    legend = FALSE,
    silent = TRUE
  )

  # Replace the ambiguous continuous -6,...,+6 scale with an explicit legend.
  n_loco <- max(length(cohort_levels), max(n_selected, na.rm = TRUE))
  positive_cols <- positive_colors[seq_len(n_loco)]
  negative_cols <- rev(negative_colors)[seq_len(n_loco)]
  legend_labels <- c(
    "Positive sign", as.character(n_loco:1L),
    "Tie in direction", "0 (never selected)", "",
    "Negative sign", as.character(n_loco:1L)
  )
  legend_fills <- c(
    NA_character_, rev(positive_cols), "#FFF7BC", "#D9D9D9",
    NA_character_, NA_character_, negative_cols
  )
  headings <- legend_labels %in% c("Positive sign", "Negative sign")
  legend_rows <- lapply(seq_along(legend_labels), function(i) {
    key <- if (is.na(legend_fills[i])) {
      grid::nullGrob()
    } else {
      grid::rectGrob(
        width = grid::unit(0.20, "in"), height = grid::unit(0.20, "in"),
        gp = grid::gpar(fill = legend_fills[i], col = "grey25")
      )
    }
    label <- grid::textGrob(
      legend_labels[i], x = 0, just = "left",
      gp = grid::gpar(
        fontsize = 9,
        fontface = if (headings[i]) "bold" else "plain"
      )
    )
    gridExtra::arrangeGrob(
      key, label, ncol = 2,
      widths = grid::unit.c(grid::unit(0.28, "in"), grid::unit(1, "null"))
    )
  })
  # Do not let arrangeGrob distribute the legend entries over the full heatmap
  # height. Give every entry a compact physical height, anchor the complete
  # legend at the top and leave any remaining panel height blank underneath.
  legend_row_heights <- ifelse(
    legend_labels == "", 0.10,
    ifelse(headings, 0.24, 0.22)
  )
  legend_content <- gridExtra::arrangeGrob(
    grobs = legend_rows, ncol = 1,
    heights = grid::unit(legend_row_heights, "in")
  )
  legend_grob <- gridExtra::arrangeGrob(
    legend_content, grid::nullGrob(), ncol = 1,
    heights = grid::unit.c(
      grid::unit(sum(legend_row_heights), "in"),
      grid::unit(1, "null")
    )
  )
  final_grob <- gridExtra::arrangeGrob(
    p$gtable, legend_grob, ncol = 2,
    widths = grid::unit.c(grid::unit(1, "null"), grid::unit(1.55, "in"))
  )
  if (draw) grid::grid.draw(final_grob)
  invisible(final_grob)
}


# Export pre-built grobs on isolated devices. This prevents LOCO plots from
# being drawn into a still-open whole-sample PDF when SVG files are generated.
save_grobs_pdf <- function(grobs, file_name, width, height) {
  grobs <- Filter(Negate(is.null), grobs)
  if (!length(grobs)) stop("No plots were supplied for PDF export.")
  grDevices::pdf(file_name, width = width, height = height, onefile = TRUE)
  device_id <- grDevices::dev.cur()
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(device_id)
    }
  }, add = TRUE)
  for (i in seq_along(grobs)) {
    if (i > 1L) grid::grid.newpage()
    grid::grid.draw(grobs[[i]])
  }
  grDevices::dev.off(device_id)
  invisible(file_name)
}


save_grob_svg <- function(grob, file_name, width, height) {
  if (is.null(grob)) stop("No plot was supplied for SVG export.")
  grDevices::svg(file_name, width = width, height = height)
  device_id <- grDevices::dev.cur()
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(device_id)
    }
  }, add = TRUE)
  grid::grid.draw(grob)
  grDevices::dev.off(device_id)
  invisible(file_name)
}


export_loco_summary_svg <- function(
    file_name, loco_summary, scenario_name, selection_definition,
    method_label) {
  n_rows <- sum(rowSums(loco_summary$n_selected > 0L) > 0L)
  svg_height <- min(18, max(7, 4.5 + 0.20 * n_rows))
  loco_grob <- plot_penalized_loco_summary(
    loco_summary = loco_summary,
    scenario_name = scenario_name,
    selection_definition = selection_definition,
    method_label = method_label,
    draw = FALSE
  )
  save_grob_svg(loco_grob, file_name, width = 11.5, height = svg_height)
}


plot_stable_signature <- function(fit_whole, boot_res, X_mat,
                                  scenario_name,
                                  sel_prob_threshold = 0.90,
                                  sign_consistency_threshold = 0.80,
                                  require_sign_consistency = FALSE) {

  filter_label <- if(require_sign_consistency) {
    paste0("sel_prob >= ", sel_prob_threshold,
           " & sign_consistency >= ", sign_consistency_threshold)
  } else {
    paste0("sel_prob >= ", sel_prob_threshold)
  }

  message(paste0(">>> Generating Stable Heatmap for: ", scenario_name))
  message(paste0("    -> Filter: ", filter_label))

  # A. EXTRACT ORIGINAL BETAS
  # ----------------------------------------------------------------------------
  # Get the B matrix from the whole sample fit
  B_orig <- fit_whole$B
  rownames(B_orig) <- colnames(X_mat)
  colnames(B_orig) <- colnames(Y_whole)

  # B. PREPARE THE SELECTION-PROBABILITY MASK
  # ----------------------------------------------------------------------------
  sel_prob_mask <- boot_res %>%
    dplyr::select(exposure, outcome, sel_prob) %>%
    tidyr::pivot_wider(names_from = outcome, values_from = sel_prob) %>%
    tibble::column_to_rownames("exposure")

  # Align rows and columns with B matrix
  sel_prob_mask <- sel_prob_mask[rownames(B_orig), colnames(B_orig)]
  stability_mask <- as.matrix(sel_prob_mask) >= sel_prob_threshold

  # C. ADD DIRECTIONAL CONSISTENCY MASK IF REQUESTED
  # ----------------------------------------------------------------------------
  if(require_sign_consistency) {
    sign_consistency_mask <- boot_res %>%
      dplyr::select(exposure, outcome, sign_consistency) %>%
      tidyr::pivot_wider(names_from = outcome, values_from = sign_consistency) %>%
      tibble::column_to_rownames("exposure")

    # Align rows and columns with B matrix
    sign_consistency_mask <- sign_consistency_mask[rownames(B_orig), colnames(B_orig)]
    stability_mask <- stability_mask &
      (as.matrix(sign_consistency_mask) >= sign_consistency_threshold)
  }

  stability_mask[is.na(stability_mask)] <- FALSE

  # D. APPLY THE FILTER (MASKING)
  # ----------------------------------------------------------------------------
  B_stable <- B_orig
  B_stable[!stability_mask] <- 0

  # E. POST-FILTER CLEANUP
  # ----------------------------------------------------------------------------
  # Remove rows that are entirely zero after masking (cleaner visualization)
  active_rows <- rowSums(B_stable != 0) > 0
  if(sum(active_rows) == 0) {
    message("    ! WARNING: No associations survived the stability filter.")
    return(NULL)
  }
  coef_subset <- B_stable[active_rows, , drop = FALSE]

  # F. ROW ANNOTATIONS (Families)
  # ----------------------------------------------------------------------------
  clean_names <- gsub("_None|_Ter_2|_Ter_3", "", rownames(coef_subset))
  row_info <- vars_info %>% 
    dplyr::filter(variable_name %in% clean_names) %>%
    dplyr::distinct(variable_name, .keep_all = TRUE)

  df_row_ann <- data.frame(
    Family = row_info$family[match(clean_names, row_info$variable_name)],
    row.names = rownames(coef_subset)
  )

  # G. RENDER PHEATMAP
  # ----------------------------------------------------------------------------
  p <- pheatmap(
    coef_subset,
    main              = paste0(scenario_name, " (Stable Signature: ", filter_label, ")"),
    color             = custom_colors_pro,          # Defined in Block 8
    breaks            = get_pro_breaks(coef_subset),# Defined in Block 8
    annotation_row    = df_row_ann,
    cluster_cols      = FALSE,
    cluster_rows      = TRUE, # Set to FALSE to keep the same predictor order across plots
    border_color      = "white",
    fontsize_row      = 8,
    cellwidth         = 30,
    silent            = TRUE
  )

  print(p)
}


plot_original_vs_stable <- function(
    fit_whole, boot_res, X_mat, scenario_name,
    family_filter = NULL,
    cellwidth = 24, cellheight = 10,
    draw = TRUE) {

  B_orig <- as.matrix(fit_whole$B)
  rownames(B_orig) <- colnames(X_mat)
  colnames(B_orig) <- colnames(Y_whole)

  stable_mask <- boot_res %>%
    dplyr::select(exposure, outcome, passes_stability_filter) %>%
    tidyr::pivot_wider(
      names_from = outcome, values_from = passes_stability_filter
    ) %>%
    tibble::column_to_rownames("exposure")
  stable_mask <- as.matrix(
    stable_mask[rownames(B_orig), colnames(B_orig), drop = FALSE]
  )
  stable_mask[is.na(stable_mask)] <- FALSE

  B_stable <- B_orig
  B_stable[!stable_mask] <- 0

  # Use the common selection tolerance in the displayed matrices. Grey is
  # therefore reserved exclusively for coefficients treated as not selected.
  B_orig[abs(B_orig) <= selection_tol] <- 0
  B_stable[abs(B_stable) <= selection_tol] <- 0

  # Use the original active rows in both panels. Associations that fail the
  # bootstrap filter remain as grey/zero cells on the right.
  rows_to_plot <- rowSums(abs(B_orig) > selection_tol, na.rm = TRUE) > 0
  if (!any(rows_to_plot)) {
    message("    ! No selected coefficients in: ", scenario_name)
    return(invisible(NULL))
  }
  B_orig <- B_orig[rows_to_plot, , drop = FALSE]
  B_stable <- B_stable[rows_to_plot, , drop = FALSE]

  clean_names <- gsub("_None|_Ter_2|_Ter_3", "", rownames(B_orig))
  row_info <- vars_info %>%
    dplyr::filter(variable_name %in% clean_names) %>%
    dplyr::distinct(variable_name, .keep_all = TRUE)
  row_family <- row_info$family[match(clean_names, row_info$variable_name)]

  if (!is.null(family_filter)) {
    keep_family <- !is.na(row_family) & row_family == family_filter
    B_orig <- B_orig[keep_family, , drop = FALSE]
    B_stable <- B_stable[keep_family, , drop = FALSE]
    row_family <- row_family[keep_family]
    if (nrow(B_orig) == 0L) return(invisible(NULL))
  }

  row_annotation <- data.frame(
    Family = row_family,
    row.names = rownames(B_orig)
  )

  common_limit <- max(abs(c(B_orig, B_stable)), na.rm = TRUE)
  if (!is.finite(common_limit) || common_limit <= selection_tol) {
    common_limit <- max(0.1, selection_tol * 10)
  }

  # 50 negative intervals, one grey interval [-tol, +tol], and 50
  # positive intervals. This produces exactly one more break than colours.
  negative_breaks <- seq(
    -common_limit, -selection_tol, length.out = 51L
  )
  positive_breaks <- seq(
    selection_tol, common_limit, length.out = 51L
  )
  common_breaks <- c(
    negative_breaks,
    selection_tol,
    positive_breaks[-1L]
  )

  title_suffix <- if (is.null(family_filter)) {
    scenario_name
  } else {
    paste0(scenario_name, " - ", family_filter)
  }

  common_args <- list(
    color = custom_colors_pro,
    breaks = common_breaks,
    annotation_row = row_annotation,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color = "white",
    fontsize_row = 8,
    cellwidth = cellwidth,
    cellheight = cellheight,
    silent = TRUE
  )

  p_original <- do.call(
    pheatmap,
    c(
      list(mat = B_orig, main = paste0(title_suffix, ": original")),
      common_args
    )
  )
  p_stable <- do.call(
    pheatmap,
    c(
      list(
        mat = B_stable,
        main = paste0(
          title_suffix, ": bootstrap stable\nsel_prob >= ",
          sel_prob_threshold, "; sign_consistency >= ",
          sign_consistency_threshold
        )
      ),
      common_args
    )
  )

  # Draw on the page already opened by the export loop. grid.arrange() opens
  # a new page by default, which previously created blank pages between plots.
  comparison_grob <- gridExtra::arrangeGrob(
    p_original$gtable, p_stable$gtable, ncol = 2
  )
  if (draw) grid::grid.draw(comparison_grob)
  invisible(comparison_grob)
}


make_heatmap_breaks <- function(left_matrix, right_matrix) {
  values <- c(left_matrix, right_matrix)
  finite_values <- values[is.finite(values)]
  common_limit <- if (length(finite_values) == 0L) {
    NA_real_
  } else {
    max(abs(finite_values))
  }
  if (!is.finite(common_limit) || common_limit <= selection_tol) {
    common_limit <- max(0.1, selection_tol * 10)
  }
  negative_breaks <- seq(
    -common_limit, -selection_tol, length.out = 51L
  )
  positive_breaks <- seq(
    selection_tol, common_limit, length.out = 51L
  )
  c(negative_breaks, selection_tol, positive_breaks[-1L])
}


get_heatmap_row_family <- function(exposure_names) {
  clean_names <- gsub("_None|_Ter_2|_Ter_3", "", exposure_names)
  row_info <- vars_info %>%
    dplyr::filter(variable_name %in% clean_names) %>%
    dplyr::distinct(variable_name, .keep_all = TRUE)
  row_info$family[match(clean_names, row_info$variable_name)]
}


draw_heatmap_pair <- function(
    left_matrix,
    right_matrix,
    left_title,
    right_title,
    row_family,
    cellwidth = heatmap_cellwidth,
    cellheight = heatmap_cellheight,
    overall_title = NULL,
    draw = TRUE) {

  stopifnot(
    identical(rownames(left_matrix), rownames(right_matrix)),
    identical(colnames(left_matrix), colnames(right_matrix))
  )
  breaks <- make_heatmap_breaks(left_matrix, right_matrix)
  annotation <- data.frame(
    Family = row_family,
    row.names = rownames(left_matrix)
  )
  common_args <- list(
    color = custom_colors_pro,
    breaks = breaks,
    annotation_row = annotation,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color = "white",
    fontsize_row = 8,
    cellwidth = cellwidth,
    cellheight = cellheight,
    na_col = "grey85",
    silent = TRUE
  )
  left_plot <- do.call(
    pheatmap,
    c(list(mat = left_matrix, main = left_title), common_args)
  )
  right_plot <- do.call(
    pheatmap,
    c(list(mat = right_matrix, main = right_title), common_args)
  )
  pair_args <- list(left_plot$gtable, right_plot$gtable, ncol = 2)
  if (!is.null(overall_title)) {
    pair_args$top <- grid::textGrob(
      overall_title,
      gp = grid::gpar(fontface = "bold", fontsize = 14)
    )
  }
  pair_grob <- do.call(gridExtra::arrangeGrob, pair_args)
  if (draw) grid::grid.draw(pair_grob)
  invisible(pair_grob)
}


plot_unpenalized_pre_post <- function(
    unpenalized_object,
    inference_results,
    scenario_name,
    family_filter = NULL,
    draw = TRUE) {

  B_unpenalized <- as.matrix(unpenalized_object$fit$B)
  rownames(B_unpenalized) <- colnames(unpenalized_object$X_selected)
  colnames(B_unpenalized) <- colnames(Y_whole)

  significance_mask <- inference_results %>%
    dplyr::select(exposure, outcome, significant) %>%
    tidyr::pivot_wider(
      names_from = outcome, values_from = significant
    ) %>%
    tibble::column_to_rownames("exposure")
  significance_mask <- as.matrix(significance_mask[
    rownames(B_unpenalized), colnames(B_unpenalized), drop = FALSE
  ])
  significance_mask[is.na(significance_mask)] <- FALSE

  B_significant_plot <- B_unpenalized
  B_significant_plot[!significance_mask] <- NA_real_
  row_family <- get_heatmap_row_family(rownames(B_unpenalized))

  title_suffix <- scenario_name
  if (!is.null(family_filter)) {
    keep <- !is.na(row_family) & row_family == family_filter
    if (!any(keep)) return(invisible(NULL))
    B_unpenalized <- B_unpenalized[keep, , drop = FALSE]
    B_significant_plot <- B_significant_plot[keep, , drop = FALSE]
    row_family <- row_family[keep]
    title_suffix <- paste0(scenario_name, " - ", family_filter)
  }

  draw_heatmap_pair(
    left_matrix = B_unpenalized,
    right_matrix = B_significant_plot,
    left_title = paste0(title_suffix, ": unpenalised refit"),
    right_title = paste0(
      title_suffix,
      ": post-bootstrap\nOnly p_value < 0.05 associations are coloured"
    ),
    row_family = row_family,
    draw = draw
  )
}
