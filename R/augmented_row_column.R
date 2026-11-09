#' Calculate augmented row-column design capacity
#'
#' Reports the plot counts, error degrees of freedom and feasibility of an
#' augmented row-column layout with superimposed blocks (Piepho and
#' Williams, 2016). Rows are grouped into row groups and columns into column
#' groups; their intersections are the blocks. Groups can have equal sizes
#' (`rows_per_block`, `cols_per_block`) or unequal sizes
#' (`rows_per_rowgroup`, `cols_per_colgroup`), and blocks can hold any number
#' of check plots (`checks_per_block`).
#'
#' @param rows_in_field Number of rows in the whole field.
#' @param cols_in_field Number of columns in the whole field.
#' @param rows_per_block Number of rows in each row group (equal groups).
#' @param cols_per_block Number of columns in each column group (equal
#'   groups).
#' @param n_controls Number of replicated checks/controls (`vc`).
#' @param rows_per_rowgroup Optional vector with the number of rows of each
#'   row group (`kb(i)`), for unequal groups. Replaces `rows_per_block`.
#' @param cols_per_colgroup Optional vector with the number of columns of
#'   each column group (`sb(j)`). Replaces `cols_per_block`.
#' @param checks_per_block Optional number of check plots per block
#'   (`pb(i, j)`): a single number, a matrix with one row per row group and
#'   one column per column group, or a vector by row group. Defaults to
#'   `n_controls` (every block a complete block of checks). Smaller values
#'   give incomplete blocks; larger values give extra check replication.
#'
#' @return A one-row tibble with the field size, number of groups and blocks,
#' check and entry plots, `error_df`, the error degrees of freedom of the
#' row-column ANOVA,
#' `check_plots - n_controls - (rows - 1) - (cols - 1) -
#' (row groups - 1)(column groups - 1)` (with fixed rows, columns and
#' blocks a design needs `error_df > 0`; with `error_df = 0`, as in Fig. 2 of
#' Piepho and Williams (2016), analyze with random blocks or random entries;
#' Clarke and Stefanova (2011) recommend about 50),
#' `check_fraction`, and the smallest slack of the row-group and column-group
#' conditions of Piepho and Williams (2016), which must be non-negative for
#' the row-group and column-group sub-designs to be estimable on their own.
#'
#' @references Piepho, H.-P. and Williams, E. R. (2016). Augmented row-column
#'   designs for a small number of checks. Agronomy Journal 108, 2256-2262.
#'   \doi{10.2134/agronj2016.06.0325}
#'
#' @examples
#' augmented_row_column_capacity(
#'   rows_in_field = 12,
#'   cols_in_field = 12,
#'   rows_per_block = 4,
#'   cols_per_block = 4,
#'   n_controls = 4
#' )
#'
#' # Example 3 of Piepho and Williams (2016): unequal groups and blocks.
#' augmented_row_column_capacity(
#'   n_controls = 4,
#'   rows_per_rowgroup = c(6, 5, 6, 5),
#'   cols_per_colgroup = c(3, 4, 3),
#'   checks_per_block = rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
#' )
#'
#' @export
augmented_row_column_capacity <- function(
    rows_in_field = NULL,
    cols_in_field = NULL,
    rows_per_block = NULL,
    cols_per_block = NULL,
    n_controls,
    rows_per_rowgroup = NULL,
    cols_per_colgroup = NULL,
    checks_per_block = NULL
) {
  spec <- arc_group_spec(
    rows_in_field, cols_in_field, rows_per_block, cols_per_block,
    rows_per_rowgroup, cols_per_colgroup, checks_per_block, n_controls
  )
  arc_spec_capacity(spec)
}

#' Generate an augmented row-column field design
#'
#' Randomizes unreplicated treatments and replicated controls (checks) in an
#' augmented row-column design with superimposed blocks (Piepho and
#' Williams, 2016). Adjacent rows form row groups and adjacent columns form
#' column groups; each intersection is a block that receives a given number
#' of check plots. Row, column and block effects are all estimable from the
#' checks, and the layout must leave positive error degrees of freedom (see
#' [augmented_row_column_capacity()]).
#'
#' Two check-allocation methods are available:
#' \describe{
#'   \item{`"latin"`}{For equal groups where every block is a complete block
#'     of checks (one copy of each). Checks are placed with covering Latin
#'     rectangles, so controls never share a row or column within a block,
#'     the same control never appears twice in a field row or column, and
#'     every field row and column has a check. Candidates are scored for
#'     spatial spread and row/column balance.}
#'   \item{`"optimal"`}{The general strategy of Piepho and Williams (2016):
#'     unequal group sizes and any number of check plots per block
#'     (incomplete blocks or extra replication of checks). Check labels are
#'     first allocated to blocks balancing row groups, column groups and
#'     blocks (Step 2); check plots are then placed within blocks and the
#'     layout improved by a local search that maximizes the D-criterion of
#'     the check-plot model `row + column + rowgroup:colgroup + check`
#'     (Steps 3-5). Check plots next to each other (same column in
#'     consecutive rows, or same row in consecutive columns) are avoided,
#'     inside blocks above all, and spatial spread breaks ties. Entries fill the
#'     remaining plots at random (Step 6).}
#' }
#' `method = "auto"` uses `"latin"` when it applies and `"optimal"`
#' otherwise.
#'
#' @param treatments Character vector of unreplicated treatment/entry names.
#'   Its length must equal the number of entry plots.
#' @param controls Character vector of replicated check/control names.
#' @param rows_in_field Number of rows in the whole field.
#' @param cols_in_field Number of columns in the whole field.
#' @param rows_per_block Number of rows in each row group (equal groups).
#' @param cols_per_block Number of columns in each column group (equal
#'   groups).
#' @param plot_type Plot numbering order. Use `"serpentine"` or `"cartesian"`.
#' @param n_candidates Number of candidate layouts (`"latin"`) or random
#'   restarts of the local search (`"optimal"`, at most 50; the search stops
#'   after 3 restarts without improvement).
#' @param seed Optional random seed.
#' @param output_file Optional path for writing the complete design as a
#'   tab-delimited file.
#' @param trial_name Optional Breedbase trial name. When given, a
#'   `plot_name` column with the plot names used by [save_breedbase()] is
#'   added to the design (see [add_plot_names()]), so the names can be
#'   checked with [breedbase_check()] and changed before saving.
#' @param location Breedbase location of the trial, for the `{location}`
#'   field of the plot names. Used only with `trial_name`.
#' @param plot_name_pattern Pattern for the plot names, as in
#'   [add_plot_names()]. Used only with `trial_name`.
#' @param rows_per_rowgroup,cols_per_colgroup Optional vectors of unequal
#'   group sizes, as in [augmented_row_column_capacity()].
#' @param checks_per_block Optional number of check plots per block, as in
#'   [augmented_row_column_capacity()].
#' @param method Check allocation method: `"auto"`, `"latin"` or
#'   `"optimal"`.
#'
#' @return A list with:
#' \describe{
#'   \item{input}{Validated design inputs.}
#'   \item{summary}{Design size, error df, method, and the D-criterion
#'   `d_criterion` (`det(X'X)^(1/p)` of the check-plot model; higher is
#'   better), and the number of pairs of neighbouring check plots (same
#'   column in consecutive rows, or same row in consecutive columns) inside a
#'   block and across block boundaries.}
#'   \item{design}{Complete plot-level design with row, column, block, entry,
#'   control flag, and plot order.}
#'   \item{breedbase_design}{Breedbase-compatible columns:
#'   `plots`, `block`, `all_entries`, `rep`, and `is_control`.}
#' }
#'
#' @references Piepho, H.-P. and Williams, E. R. (2016). Augmented row-column
#'   designs for a small number of checks. Agronomy Journal 108, 2256-2262.
#'   \doi{10.2134/agronj2016.06.0325}
#'
#' @examples
#' treatments <- paste0("G", seq_len(108))
#' controls <- paste0("Check", seq_len(4))
#'
#' design <- augmented_row_column_design(
#'   treatments = treatments,
#'   controls = controls,
#'   rows_in_field = 12,
#'   cols_in_field = 12,
#'   rows_per_block = 4,
#'   cols_per_block = 4,
#'   n_candidates = 50,
#'   seed = 123
#' )
#'
#' head(design$design)
#' design$summary
#'
#' # Example 3 of Piepho and Williams (2016): 173 entries, 4 checks, unequal
#' # row and column groups, and 3 to 5 check plots per block.
#' pb <- rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
#' ex3 <- augmented_row_column_design(
#'   treatments = paste0("E", 1:173),
#'   controls = c("A", "B", "C", "D"),
#'   rows_per_rowgroup = c(6, 5, 6, 5),
#'   cols_per_colgroup = c(3, 4, 3),
#'   checks_per_block = pb,
#'   n_candidates = 3,
#'   seed = 1
#' )
#' ex3$summary[, c("method", "check_plots", "entry_plots", "error_df", "d_criterion")]
#'
#' @export
augmented_row_column_design <- function(
    treatments,
    controls,
    rows_in_field = NULL,
    cols_in_field = NULL,
    rows_per_block = NULL,
    cols_per_block = NULL,
    plot_type = c("serpentine", "cartesian"),
    n_candidates = 1000,
    seed = NULL,
    output_file = NULL,
    rows_per_rowgroup = NULL,
    cols_per_colgroup = NULL,
    checks_per_block = NULL,
    method = c("auto", "latin", "optimal"),
    trial_name = NULL,
    location = NULL,
    plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}"
) {
  plot_type <- match.arg(plot_type)
  method <- match.arg(method)
  n_candidates <- arc_positive_integer(n_candidates, "n_candidates")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  inputs <- arc_validate_inputs(
    treatments = treatments,
    controls = controls,
    rows_in_field = rows_in_field,
    cols_in_field = cols_in_field,
    rows_per_block = rows_per_block,
    cols_per_block = cols_per_block,
    rows_per_rowgroup = rows_per_rowgroup,
    cols_per_colgroup = cols_per_colgroup,
    checks_per_block = checks_per_block,
    n_candidates = n_candidates,
    plot_type = plot_type,
    method = method
  )
  spec <- inputs$spec
  capacity <- arc_spec_capacity(spec)

  if (inputs$method == "optimal") {
    searched <- arc_flex_search(spec, inputs$n_candidates)
    if (searched$score$rank < searched$score$p) {
      stop(
        "Unable to generate an augmented row-column design in which all row, ",
        "column, and block effects are estimable from the checks. ",
        "Try increasing `n_candidates`, adding check plots, or changing the groups.",
        call. = FALSE
      )
    }
    design_full <- arc_flex_build_design(
      searched$state, spec, inputs$treatments, inputs$controls, inputs$plot_type
    )
    n_valid_candidates <- searched$n_evaluated
    best_score <- -searched$score$value
  } else {
    dims <- arc_design_dims(
      rows_in_field = inputs$rows_in_field,
      cols_in_field = inputs$cols_in_field,
      rows_per_block = inputs$rows_per_block,
      cols_per_block = inputs$cols_per_block,
      n_controls = length(inputs$controls)
    )

    # Candidates only differ in check positions, so each candidate is a cheap
    # integer vector of check coordinates. Entries and plot numbering are
    # assigned once, for the selected candidate.
    best_slots <- NULL
    best_score <- Inf
    n_valid_candidates <- 0L
    for (i in seq_len(inputs$n_candidates)) {
      slots <- arc_make_check_slots(dims)
      candidate_score <- arc_score_check_slots(slots, dims)
      # The estimability check (a small QR) only runs for candidates that
      # would replace the current best.
      if (candidate_score < best_score && arc_checks_estimable(slots, dims)) {
        best_score <- candidate_score
        best_slots <- slots
      }
      n_valid_candidates <- n_valid_candidates + 1L
    }

    if (is.null(best_slots)) {
      stop(
        "Unable to generate an augmented row-column design in which all row, ",
        "column, and block effects are estimable from the checks. ",
        "Try increasing `n_candidates`, adding checks, or changing the block dimensions.",
        call. = FALSE
      )
    }

    best_slots <- arc_polish_latin(best_slots, dims)
    design_full <- arc_build_design(
      slots = best_slots,
      dims = dims,
      treatments = inputs$treatments,
      controls = inputs$controls,
      plot_type = inputs$plot_type
    )
  }

  design <- design_full |>
    dplyr::transmute(
      plots = .data$plots,
      row = .data$row,
      col = .data$col,
      block = .data$block,
      rowgroup = .data$rowgroup,
      colgroup = .data$colgroup,
      all_entries = .data$trt,
      type = .data$type,
      rep = .data$rep,
      is_control = .data$is_control
    )

  breedbase_design <- design |>
    dplyr::select(dplyr::all_of(c("plots", "block", "all_entries", "rep", "is_control")))

  summary <- capacity |>
    dplyr::mutate(
      method = inputs$method,
      plot_type = inputs$plot_type,
      n_candidates = inputs$n_candidates,
      n_valid_candidates = n_valid_candidates,
      design_score = best_score,
      d_criterion = arc_design_d_criterion(design, spec),
      adjacent_checks_in_blocks = arc_design_adjacent(design, spec)[["inside"]],
      adjacent_checks_across_blocks = arc_design_adjacent(design, spec)[["across"]]
    )

  result <- list(
    input = list(
      treatments = inputs$treatments,
      controls = inputs$controls,
      rows_in_field = inputs$rows_in_field,
      cols_in_field = inputs$cols_in_field,
      rows_per_block = inputs$rows_per_block,
      cols_per_block = inputs$cols_per_block,
      rows_per_rowgroup = spec$kb,
      cols_per_colgroup = spec$sb,
      checks_per_block = spec$pb,
      method = inputs$method,
      plot_type = inputs$plot_type,
      n_candidates = inputs$n_candidates,
      seed = seed
    ),
    summary = summary,
    design = design,
    breedbase_design = breedbase_design
  )

  class(result) <- c("augmented_row_column_design", class(result))
  if (!is.null(trial_name)) result <- add_plot_names(result, trial_name, location, plot_name_pattern)

  if (!is.null(output_file)) {
    utils::write.table(
      result$design,
      file = output_file,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
  }

  result
}

#' Plot an augmented row-column design grid
#'
#' Draws a field grid for an augmented row-column design, with plot cells,
#' row and column labels, and thick outlines around row-column blocks.
#'
#' @param design An object returned by `augmented_row_column_design()` or a
#'   design data frame with row, column, and block columns.
#' @param label Cell label to show. Use `"all_entries"`, `"type"`, `"block"`,
#'   or `"none"`.
#' @param fill Cell coloring variable. Use `"type"`, `"block"`, `"rowgroup"`,
#'   `"colgroup"`, or `"is_control"`.
#' @param show_plot_labels Logical. If `TRUE`, labels are printed inside cells
#'   unless `label = "none"`.
#' @param show_block_labels Logical. If `TRUE`, block numbers are printed.
#'   Defaults to `FALSE` so block labels do not overlap accession labels.
#' @param block_label_position Position for block labels when
#'   `show_block_labels = TRUE`. Use `"background"` for a faint watermark
#'   behind accession labels or `"center"` for a regular center label.
#' @param label_cex Optional text size for plot-cell labels.
#' @param block_label_cex Optional text size for block labels.
#' @param block_label_col Color for block labels.
#' @param block_lwd Line width for block boundaries.
#' @param main Optional plot title.
#' @param legend Logical. If `TRUE`, draw a legend for cell colors.
#' @param ... Additional arguments passed to `graphics::plot.window()`.
#'
#' @return Invisibly returns the plotted design data frame with internal
#'   plotting coordinates.
#'
#' @examples
#' design <- augmented_row_column_design(
#'   treatments = paste0("G", seq_len(108)),
#'   controls = paste0("Check", seq_len(4)),
#'   rows_in_field = 12,
#'   cols_in_field = 12,
#'   rows_per_block = 4,
#'   cols_per_block = 4,
#'   n_candidates = 20,
#'   seed = 123
#' )
#'
#' plot_augmented_row_column_design(design, label = "all_entries")
#'
#' @export
plot_augmented_row_column_design <- function(
    design,
    label = c("all_entries", "type", "block", "none"),
    fill = c("type", "block", "rowgroup", "colgroup", "is_control"),
    show_plot_labels = TRUE,
    show_block_labels = FALSE,
    block_label_position = c("background", "center"),
    label_cex = NULL,
    block_label_cex = NULL,
    block_label_col = NULL,
    block_lwd = 2,
    main = NULL,
    legend = TRUE,
    ...
) {
  label <- match.arg(label)
  fill <- match.arg(fill)
  block_label_position <- match.arg(block_label_position)

  if (is.list(design) && !is.data.frame(design) && !is.null(design$design)) {
    design <- design$design
  }

  if (!is.data.frame(design)) {
    stop("`design` must be a design data frame or an augmented row-column design object.", call. = FALSE)
  }

  required_cols <- c("row", "col", "block")
  missing_cols <- setdiff(required_cols, names(design))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (label != "none" && !label %in% names(design)) {
    stop("Column `", label, "` was not found in `design`.", call. = FALSE)
  }

  if (!fill %in% names(design)) {
    fill <- "block"
  }

  design <- as.data.frame(design, stringsAsFactors = FALSE)
  design$row <- suppressWarnings(as.numeric(design$row))
  design$col <- suppressWarnings(as.numeric(design$col))

  if (any(is.na(design$row)) || any(is.na(design$col))) {
    stop("`row` and `col` columns must be numeric or coercible to numeric.", call. = FALSE)
  }

  rows <- sort(unique(design$row))
  cols <- sort(unique(design$col))
  design$.plot_row <- match(design$row, rows)
  design$.plot_col <- match(design$col, cols)
  design$.fill_value <- as.character(design[[fill]])
  fill_levels <- sort(unique(design$.fill_value))
  fill_colors <- grDevices::hcl.colors(length(fill_levels), palette = "Set 3")
  names(fill_colors) <- fill_levels

  if (is.null(label_cex)) {
    label_cex <- max(0.25, min(0.8, 7 / max(length(rows), length(cols))))
  }

  if (is.null(main)) {
    main <- "Augmented row-column design"
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 4, if (isTRUE(legend)) 7 else 2), xpd = NA)

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0.5, length(cols) + 0.5),
    ylim = c(length(rows) + 0.5, 0.5),
    xaxs = "i",
    yaxs = "i",
    asp = 1,
    ...
  )

  graphics::rect(
    xleft = design$.plot_col - 0.5,
    ybottom = design$.plot_row - 0.5,
    xright = design$.plot_col + 0.5,
    ytop = design$.plot_row + 0.5,
    col = fill_colors[design$.fill_value],
    border = "grey70",
    lwd = 0.7
  )

  graphics::abline(v = seq(0.5, length(cols) + 0.5, by = 1), col = "grey80", lwd = 0.5)
  graphics::abline(h = seq(0.5, length(rows) + 0.5, by = 1), col = "grey80", lwd = 0.5)

  block_split <- split(design, design$block)
  for (block_design in block_split) {
    min_col <- min(block_design$.plot_col)
    max_col <- max(block_design$.plot_col)
    min_row <- min(block_design$.plot_row)
    max_row <- max(block_design$.plot_row)

    graphics::rect(
      xleft = min_col - 0.5,
      ybottom = min_row - 0.5,
      xright = max_col + 0.5,
      ytop = max_row + 0.5,
      border = "black",
      lwd = block_lwd
    )

    if (isTRUE(show_block_labels)) {
      block_cex <- block_label_cex
      block_col <- block_label_col
      if (is.null(block_cex)) {
        block_cex <- if (block_label_position == "background") {
          max(1.2, min(4, 0.55 * min(max_col - min_col + 1, max_row - min_row + 1)))
        } else {
          0.9
        }
      }
      if (is.null(block_col)) {
        block_col <- if (block_label_position == "background") {
          grDevices::adjustcolor("black", alpha.f = 0.18)
        } else {
          "black"
        }
      }
      graphics::text(
        x = mean(c(min_col, max_col)),
        y = mean(c(min_row, max_row)),
        labels = paste0("B", unique(block_design$block)[1]),
        cex = block_cex,
        font = 2,
        col = block_col
      )
    }
  }

  if (isTRUE(show_plot_labels) && label != "none") {
    graphics::text(
      x = design$.plot_col,
      y = design$.plot_row,
      labels = as.character(design[[label]]),
      cex = label_cex,
      col = "black"
    )
  }

  graphics::axis(3, at = seq_along(cols), labels = cols, las = 2, tick = FALSE)
  graphics::axis(2, at = seq_along(rows), labels = rows, las = 1, tick = FALSE)
  graphics::mtext("Column", side = 3, line = 2.7)
  graphics::mtext("Row", side = 2, line = 3)
  graphics::title(main = main)
  graphics::box()

  if (isTRUE(legend)) {
    graphics::legend(
      x = length(cols) + 1,
      y = 0.5,
      legend = fill_levels,
      fill = fill_colors,
      border = "grey60",
      title = fill,
      bty = "n",
      xjust = 0,
      yjust = 0,
      cex = 0.8
    )
  }

  invisible(tibble::as_tibble(design))
}

#' Analyze an augmented row-column design
#'
#' Fits the fixed-effects model:
#' `y = mean + rowgroup + colgroup + rowgroup:colgroup + rowgroup:row +
#' colgroup:col + genotype + error`.
#' Sums of squares are sequential, in the order above. Degrees of freedom are
#' the ranks actually estimable from the observed (non-missing) plots. For a
#' complete design they equal the augmented row-column formulas, reported in
#' `anova_df`, which use field rows (`k`), field columns (`s`), row groups
#' (`gk`), column groups (`gs`), check cultivars (`vc`), and unreplicated
#' entries (`ve`):
#' row groups `gk - 1`, rows within row groups `k - gk`, column groups
#' `gs - 1`, columns within column groups `s - gs`, blocks
#' `(gk - 1)(gs - 1)`, genotypes `vc + ve - 1`, and error as the remainder of
#' `k * s - 1`.
#'
#' Unreplicated entries each fit their own plot exactly, so the model is fitted
#' without a dense genotype design matrix: the residual comes from the
#' replicated (check) plots alone. `model` is the corresponding `lm` fit on the
#' check plots, with aliased columns removed. This keeps the analysis fast for
#' large trials.
#'
#' @param data A complete augmented row-column design data frame, usually
#'   `design$design` from `augmented_row_column_design()`, with a response
#'   column added.
#' @param response Name of the numeric response column.
#' @param genotype_col Name of the genotype/treatment column.
#' @param rowgroup_col Name of the row-group column.
#' @param colgroup_col Name of the column-group column.
#' @param row_col Name of the field row column.
#' @param col_col Name of the field column column.
#' @param type_col Optional entry type column. Used only for genotype summaries.
#' @param control_col Optional control flag column. Identifies the checks
#'   for `entries = "random"`; otherwise used only for genotype summaries.
#' @param random_blocks Logical. If `TRUE`, rows within row groups, columns
#'   within column groups and blocks are random effects, which recovers
#'   inter-block information (Wolfinger et al., 1997; Möhring et al., 2015).
#'   Row and column groups stay fixed. Fitted by REML with 'lme4'.
#' @param entries `"fixed"` (default) or `"random"`. With `"random"`, checks
#'   are fixed and the unreplicated entries are random effects from a common
#'   breeding population (Piepho, Williams and Fleck, 2006); entry values are
#'   then BLUPs. Fitted by REML with 'lme4'.
#' @param spatial Logical. If `TRUE`, the plot errors get a spatial add-on: a
#'   two-dimensional P-spline surface over rows and columns fitted with the
#'   'SpATS' package, with row and column groups fixed and rows, columns and
#'   blocks random. Genotypes are fixed.
#' @param spatial_nseg Optional number of P-spline segments for columns and
#'   rows, `c(col, row)`. Defaults to half the number of columns and rows
#'   (at least 3, at most 40).
#'
#' @return An object of class `augmented_row_column_analysis`: a list with
#' \describe{
#'   \item{input, formula, model}{The column mapping, model formula, and an
#'     `lm` fit on the check plots.}
#'   \item{anova, anova_df, anova_formatted}{Sequential ANOVA (estimable df),
#'     design df formulas, and a formatted ANOVA table.}
#'   \item{plot_results}{The input data with fitted values and residuals.}
#'   \item{genotype_summary}{One row per genotype with raw, fitted, and
#'     adjusted means and the standard error of the adjusted mean.}
#'   \item{adjusted_means}{Genotypes ranked by adjusted mean (best first).}
#'   \item{diagnostics}{Error df, residual standard deviation, R-squared,
#'     `mean_sed` (mean standard error of a difference between two adjusted
#'     means), `lsd_5pct` (the matching least significant difference), and
#'     flags for rank deficiency or a saturated fit.}
#'   \item{variance_components, mixed_model}{For mixed or spatial fits: the
#'     REML variance components and the `lmerMod` or `SpATS` fit; otherwise
#'     `NULL`.}
#'   \item{message}{A note when F-tests or adjusted means are not estimable.}
#' }
#'
#' @section Mixed and spatial models:
#' The ANOVA table always comes from the fixed-effects model of Piepho and
#' Williams (2016, Eq. 1). With `random_blocks`, `entries = "random"` or
#' `spatial`, the adjusted means, their standard errors, `mean_sed` and
#' `lsd_5pct` come from the corresponding REML fit instead. Standard errors
#' are exact prediction standard errors from the mixed-model equations for
#' the estimated variance components; `lsd_5pct` uses the fixed-model error
#' df as an approximation.
#'
#' @section Adjusted means:
#' With the default fixed-effects model, genotype effects are fixed, so the
#' adjusted means are best linear unbiased estimates (BLUEs). Each is the genotype effect plus the average row,
#' column, and block effect over all observed plots, i.e. the expected value
#' of the genotype in an average plot of the field. For an unreplicated entry
#' this is its observed value corrected for the row, column, and block
#' effects of its plot, estimated from the checks; its standard error
#' combines the plot error and the uncertainty of that correction. Use
#' `adjusted_mean` (not `raw_mean`) to rank and select entries.
#'
#' @examples
#' design <- augmented_row_column_design(
#'   treatments = paste0("G", seq_len(108)),
#'   controls = paste0("Check", seq_len(4)),
#'   rows_in_field = 12,
#'   cols_in_field = 12,
#'   rows_per_block = 4,
#'   cols_per_block = 4,
#'   n_candidates = 3,
#'   seed = 123
#' )
#'
#' design_data <- design$design
#' design_data$yield <- rnorm(nrow(design_data), mean = 100, sd = 10)
#' analysis <- analyze_augmented_row_column_design(design_data, response = "yield")
#' analysis
#' head(analysis$adjusted_means)
#'
#' @export
analyze_augmented_row_column_design <- function(
    data,
    response,
    genotype_col = "all_entries",
    rowgroup_col = "rowgroup",
    colgroup_col = "colgroup",
    row_col = "row",
    col_col = "col",
    type_col = "type",
    control_col = "is_control",
    random_blocks = FALSE,
    entries = c("fixed", "random"),
    spatial = FALSE,
    spatial_nseg = NULL
) {
  entries <- match.arg(entries)
  if (is.list(data) && !is.data.frame(data) && !is.null(data$design)) {
    data <- data$design
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (isTRUE(spatial) && entries == "random") {
    stop("`spatial = TRUE` is only available with `entries = \"fixed\"`.", call. = FALSE)
  }

  required_cols <- c(response, genotype_col, rowgroup_col, colgroup_col, row_col, col_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  design_structure_data <- tibble::as_tibble(data) |>
    dplyr::transmute(
      genotype = as.factor(.data[[genotype_col]]),
      rowgroup = as.factor(.data[[rowgroup_col]]),
      colgroup = as.factor(.data[[colgroup_col]]),
      row = as.factor(.data[[row_col]]),
      col = as.factor(.data[[col_col]]),
      row_nested = interaction(.data[[rowgroup_col]], .data[[row_col]], drop = TRUE),
      col_nested = interaction(.data[[colgroup_col]], .data[[col_col]], drop = TRUE),
      type = if (type_col %in% names(data)) {
        as.character(.data[[type_col]])
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(data)) {
        as.integer(.data[[control_col]])
      } else {
        NA_integer_
      }
    ) |>
    dplyr::filter(
      !is.na(.data$genotype),
      !is.na(.data$rowgroup),
      !is.na(.data$colgroup),
      !is.na(.data$row),
      !is.na(.data$col)
    )

  analysis_data <- tibble::as_tibble(data) |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    dplyr::transmute(
      .row_id = .data$.row_id,
      response = suppressWarnings(as.numeric(.data[[response]])),
      genotype = as.factor(.data[[genotype_col]]),
      rowgroup = as.factor(.data[[rowgroup_col]]),
      colgroup = as.factor(.data[[colgroup_col]]),
      row = as.factor(.data[[row_col]]),
      col = as.factor(.data[[col_col]]),
      row_nested = interaction(.data[[rowgroup_col]], .data[[row_col]], drop = TRUE),
      col_nested = interaction(.data[[colgroup_col]], .data[[col_col]], drop = TRUE),
      type = if (type_col %in% names(data)) {
        as.character(.data[[type_col]])
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(data)) {
        as.integer(.data[[control_col]])
      } else {
        NA_integer_
      }
    ) |>
    dplyr::filter(
      !is.na(.data$response),
      !is.na(.data$genotype),
      !is.na(.data$rowgroup),
      !is.na(.data$colgroup),
      !is.na(.data$row),
      !is.na(.data$col)
    )

  if (nrow(analysis_data) == 0) {
    stop("No complete observations are available for analysis.", call. = FALSE)
  }

  if (dplyr::n_distinct(analysis_data$genotype) < 2) {
    stop("At least two genotypes are required for analysis.", call. = FALSE)
  }

  model_formula <- stats::as.formula(
    "response ~ rowgroup + row_nested + colgroup + col_nested + rowgroup:colgroup + genotype"
  )

  design_df <- arc_augmented_design_df(design_structure_data)
  fit_parts <- arc_fit_augmented_row_column(analysis_data)
  fit <- fit_parts$model
  anova_table <- fit_parts$anova
  anova_formatted <- arc_format_augmented_anova(
    anova_table = anova_table,
    response = analysis_data$response
  )

  plot_results <- tibble::as_tibble(data) |>
    dplyr::mutate(
      .analysis_response = suppressWarnings(as.numeric(.data[[response]]))
    )

  plot_results$fitted <- NA_real_
  plot_results$residual <- NA_real_
  plot_results$fitted[analysis_data$.row_id] <- fit_parts$fitted
  plot_results$residual[analysis_data$.row_id] <- analysis_data$response - fit_parts$fitted

  genotype_summary <- plot_results |>
    dplyr::group_by(.data[[genotype_col]]) |>
    dplyr::summarise(
      genotype = dplyr::first(as.character(.data[[genotype_col]])),
      type = if (type_col %in% names(plot_results)) {
        dplyr::first(as.character(.data[[type_col]]))
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(plot_results)) {
        dplyr::first(as.integer(.data[[control_col]]))
      } else {
        NA_integer_
      },
      n_observations = sum(!is.na(.data$.analysis_response)),
      raw_mean = mean(.data$.analysis_response, na.rm = TRUE),
      fitted_mean = mean(.data$fitted, na.rm = TRUE),
      residual_mean = mean(.data$residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::select(-dplyr::all_of(genotype_col))

  genotype_summary$raw_mean[is.nan(genotype_summary$raw_mean)] <- NA_real_
  genotype_summary$fitted_mean[is.nan(genotype_summary$fitted_mean)] <- NA_real_
  genotype_summary$residual_mean[is.nan(genotype_summary$residual_mean)] <- NA_real_

  # Mixed-model and spatial alternatives (Piepho and Williams, 2016,
  # Analysis section) replace the fixed-effects adjusted means; the ANOVA
  # always comes from the fixed-effects model.
  model_type <- "fixed"
  adjusted_source <- fit_parts$adjusted_means
  adjusted_vcov <- fit_parts$vcov
  precision <- fit_parts$precision
  variance_components <- NULL
  mixed_model <- NULL
  spatial_effective_dim <- NA_real_
  if (isTRUE(spatial) || isTRUE(random_blocks) || entries == "random") {
    analysis_data$is_check <- if (control_col %in% names(data) && !all(is.na(analysis_data$is_control))) {
      analysis_data$is_control == 1L
    } else {
      analysis_data$genotype %in% names(which(table(analysis_data$genotype) > 1))
    }
    alt <- if (isTRUE(spatial)) {
      model_type <- "spatial"
      arc_rc_spatial(analysis_data, spatial_nseg, fit_parts$error_df)
    } else {
      model_type <- "mixed"
      arc_rc_mixed(analysis_data, isTRUE(random_blocks), entries == "random", fit_parts$error_df)
    }
    adjusted_source <- alt$adjusted_means
    adjusted_vcov <- alt$vcov
    precision <- alt$precision
    # Plot fitted values and residuals of the model that gives the adjusted
    # means (the ANOVA stays the fixed-effects one).
    if (length(alt$fitted) == nrow(analysis_data)) {
      plot_results$fitted[analysis_data$.row_id] <- alt$fitted
      plot_results$residual[analysis_data$.row_id] <- analysis_data$response - alt$fitted
      by_genotype <- split(seq_len(nrow(plot_results)), as.character(plot_results[[genotype_col]]))
      gs_rows <- match(genotype_summary$genotype, names(by_genotype))
      genotype_summary$fitted_mean <- vapply(gs_rows, function(i) {
        v <- plot_results$fitted[by_genotype[[i]]]
        if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
      }, numeric(1))
      genotype_summary$residual_mean <- vapply(gs_rows, function(i) {
        v <- plot_results$residual[by_genotype[[i]]]
        if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
      }, numeric(1))
    }
    variance_components <- alt$variance_components
    mixed_model <- alt$model
    if (!is.null(alt$spatial_effective_dim)) spatial_effective_dim <- alt$spatial_effective_dim
  }

  adjusted_parts <- arc_attach_adjusted_means(genotype_summary, adjusted_source)
  genotype_summary <- adjusted_parts$genotype_summary

  error_df <- fit_parts$error_df
  rss <- fit_parts$rss
  tss <- sum((analysis_data$response - mean(analysis_data$response))^2)
  n_obs <- nrow(analysis_data)
  sigma <- if (error_df > 0) sqrt(rss / error_df) else NA_real_
  r_squared <- if (tss > 0) 1 - rss / tss else NA_real_
  adjusted_r_squared <- if (error_df > 0 && tss > 0) {
    1 - (rss / error_df) / (tss / (n_obs - 1))
  } else {
    NA_real_
  }

  diagnostics <- tibble::tibble(
    n_observations = n_obs,
    n_genotypes = dplyr::n_distinct(analysis_data$genotype),
    model_rank = fit_parts$rank,
    residual_df = error_df,
    anova_error_df = error_df,
    design_error_df = design_df$display_df[design_df$term == "Residuals"],
    aliased_coefficients_dropped = fit_parts$aliased,
    sigma = sigma,
    r_squared = r_squared,
    adjusted_r_squared = adjusted_r_squared,
    rank_deficient = fit_parts$rank < fit_parts$expected_rank,
    saturated_or_near_perfect_fit = error_df == 0 ||
      isTRUE(!is.na(sigma) && sigma < sqrt(.Machine$double.eps)),
    mean_sed = precision$mean_sed,
    lsd_5pct = precision$lsd_5pct,
    anova_warning = NA_character_,
    model_type = model_type,
    spatial_effective_dim = spatial_effective_dim,
    singular_fit = if (inherits(mixed_model, "lmerMod")) lme4::isSingular(mixed_model) else NA
  )

  analysis_message <- arc_augmented_analysis_message(diagnostics)
  analysis_message <- arc_append_estimability_message(analysis_message, adjusted_source)
  if (entries == "random") {
    note <- "Entry values are BLUPs (entries random); checks are BLUEs."
    analysis_message <- if (is.null(analysis_message) || is.na(analysis_message) || !nzchar(analysis_message)) note else paste(analysis_message, note)
  }

  result <- list(
    input = list(
      response = response,
      genotype_col = genotype_col,
      rowgroup_col = rowgroup_col,
      colgroup_col = colgroup_col,
      row_col = row_col,
      col_col = col_col,
      type_col = type_col,
      control_col = control_col,
      random_blocks = isTRUE(random_blocks),
      entries = entries,
      spatial = isTRUE(spatial)
    ),
    formula = model_formula,
    model = fit,
    anova = anova_table,
    anova_df = design_df,
    anova_formatted = anova_formatted,
    plot_results = plot_results |> dplyr::select(-dplyr::all_of(".analysis_response")),
    genotype_summary = genotype_summary,
    adjusted_means = adjusted_parts$adjusted_means,
    adjusted_means_vcov = adjusted_vcov,
    diagnostics = diagnostics,
    variance_components = variance_components,
    mixed_model = mixed_model,
    message = analysis_message
  )

  class(result) <- c("augmented_row_column_analysis", class(result))
  result
}

#' Print augmented row-column analysis
#'
#' @param x An object returned by `analyze_augmented_row_column_design()`.
#' @param ... Additional arguments passed to `print()`.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.augmented_row_column_analysis <- function(x, ...) {
  cat("Augmented row-column design analysis\n\n")
  print(x$anova_formatted, ...)
  arc_print_adjusted_means(x)
  if (!is.null(x$message) && !is.na(x$message) && nzchar(x$message)) {
    cat("\n", x$message, "\n", sep = "")
  }
  invisible(x)
}

arc_print_adjusted_means <- function(x, n = 10) {
  means <- x$adjusted_means
  if (is.null(means) || all(is.na(means$adjusted_mean))) {
    return(invisible(NULL))
  }
  cat("\nTop adjusted means (", sum(!is.na(means$adjusted_mean)), " genotypes; ",
      "see $adjusted_means):\n", sep = "")
  top <- utils::head(means[!is.na(means$adjusted_mean), ], n)
  print(as.data.frame(top[, intersect(c("rank", "genotype", "n_observations", "adjusted_mean", "std_error"),
                                      names(top))]), row.names = FALSE, digits = 4)
  if (!is.na(x$diagnostics$mean_sed)) {
    cat(sprintf("Mean SED = %.4g; LSD (5%%) = %.4g\n", x$diagnostics$mean_sed, x$diagnostics$lsd_5pct))
  }
  invisible(NULL)
}

arc_positive_integer <- function(x, label) {
  if (length(x) != 1 || is.na(x)) {
    stop("`", label, "` must be a single positive integer.", call. = FALSE)
  }

  x <- as.integer(x)

  if (is.na(x) || x < 1) {
    stop("`", label, "` must be a single positive integer.", call. = FALSE)
  }

  x
}

arc_validate_inputs <- function(
    treatments,
    controls,
    rows_in_field,
    cols_in_field,
    rows_per_block,
    cols_per_block,
    n_candidates,
    plot_type,
    rows_per_rowgroup = NULL,
    cols_per_colgroup = NULL,
    checks_per_block = NULL,
    method = "auto"
) {
  treatments <- as.character(treatments)
  controls <- as.character(controls)

  treatments <- treatments[!is.na(treatments) & treatments != ""]
  controls <- controls[!is.na(controls) & controls != ""]
  n_candidates <- arc_positive_integer(n_candidates, "n_candidates")

  if (length(controls) == 0) {
    stop("At least one control/check is required.", call. = FALSE)
  }
  if (length(treatments) == 0) {
    stop("At least one treatment/entry is required.", call. = FALSE)
  }
  if (anyDuplicated(controls)) {
    stop("`controls` contains duplicated names.", call. = FALSE)
  }
  if (anyDuplicated(treatments)) {
    stop("`treatments` contains duplicated names.", call. = FALSE)
  }
  overlap <- intersect(treatments, controls)
  if (length(overlap) > 0) {
    stop(
      "`treatments` and `controls` must not contain the same name: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }

  spec <- arc_group_spec(
    rows_in_field, cols_in_field, rows_per_block, cols_per_block,
    rows_per_rowgroup, cols_per_colgroup, checks_per_block, length(controls)
  )
  capacity <- arc_spec_capacity(spec)
  n_controls <- length(controls)

  if (length(treatments) != capacity$entry_plots) {
    stop(
      "The number of treatments does not match the available entry plots. ",
      "Expected ",
      capacity$entry_plots,
      " treatments, but received ",
      length(treatments),
      ".",
      call. = FALSE
    )
  }

  # The Latin construction needs equal groups, complete blocks of checks and
  # room for every check in distinct rows and columns of each block.
  latin_ok <- spec$regular &&
    spec$gs <= spec$kb[1] && spec$gk <= spec$sb[1] &&
    n_controls <= spec$kb[1] && n_controls <= spec$sb[1]
  if (method == "latin" && !spec$regular) {
    stop(
      "Method \"latin\" needs equal row and column groups and one copy of every check ",
      "per block. Use method \"optimal\".",
      call. = FALSE
    )
  }
  if (method == "latin") {
    if (n_controls > spec$kb[1] || n_controls > spec$sb[1]) {
      stop(
        "`n_controls` must fit in each block without duplicated check rows or columns.",
        call. = FALSE
      )
    }
    if (spec$gs > spec$kb[1]) {
      stop(
        "`rows_per_block` must be at least the number of column groups ",
        "to avoid repeating the same control in a field row.",
        call. = FALSE
      )
    }
    if (spec$gk > spec$sb[1]) {
      stop(
        "`cols_per_block` must be at least the number of row groups ",
        "to avoid repeating the same control in a field column.",
        call. = FALSE
      )
    }
  }
  if (method == "auto") {
    method <- if (latin_ok) "latin" else "optimal"
  }

  row_checks <- rowSums(spec$pb)
  if (any(row_checks < spec$kb)) {
    i <- which(row_checks < spec$kb)[1]
    stop(
      "Not every field row can contain a check: row group ", i, " has only ",
      row_checks[i], " check plots for ", spec$kb[i],
      " rows. Use more checks, more column groups, or fewer rows per row group.",
      call. = FALSE
    )
  }
  col_checks <- colSums(spec$pb)
  if (any(col_checks < spec$sb)) {
    j <- which(col_checks < spec$sb)[1]
    stop(
      "Not every field column can contain a check: column group ", j, " has only ",
      col_checks[j], " check plots for ", spec$sb[j],
      " columns. Use more checks, more row groups, or fewer columns per column group.",
      call. = FALSE
    )
  }
  if (capacity$error_df < 0) {
    stop(
      "This layout has fewer check plots than row, column and block effects to estimate ",
      "(error df = ", capacity$error_df, "). ",
      "Use more checks or larger blocks (fewer blocks per field).",
      call. = FALSE
    )
  }
  # error_df == 0 (e.g. Piepho and Williams, 2016, Fig. 2) is accepted
  # silently: the effects are estimable and the layout is analyzed with
  # random blocks or random entries. summary$error_df and check_design()
  # report it.

  list(
    treatments = treatments,
    controls = controls,
    rows_in_field = spec$k,
    cols_in_field = spec$s,
    rows_per_block = if (length(unique(spec$kb)) == 1) spec$kb[1] else NA_integer_,
    cols_per_block = if (length(unique(spec$sb)) == 1) spec$sb[1] else NA_integer_,
    n_candidates = n_candidates,
    plot_type = plot_type,
    method = method,
    spec = spec
  )
}

arc_make_field_template <- function(rows_in_field, cols_in_field, rows_per_block, cols_per_block) {
  block_cols <- cols_in_field / cols_per_block

  field <- expand.grid(
    row = seq_len(rows_in_field),
    col = seq_len(cols_in_field),
    KEEP.OUT.ATTRS = FALSE
  )

  field$rowgroup <- ((field$row - 1L) %/% rows_per_block) + 1L
  field$colgroup <- ((field$col - 1L) %/% cols_per_block) + 1L
  field$block <- ((field$rowgroup - 1L) * block_cols) + field$colgroup
  field$trt <- NA_character_
  field$type <- "entry"
  field
}

arc_lm_anova_table <- function(fit) {
  if (!is.null(fit$x) && !is.null(fit$y)) {
    out <- arc_sequential_anova_from_fit(fit)
    warning_text <- NULL
  } else {
    warning_text <- NULL
    tab <- withCallingHandlers(
      stats::anova(fit),
      warning = function(w) {
        warning_text <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )

    out <- tibble::tibble(
      term = rownames(tab),
      df = as.numeric(tab[["Df"]]),
      sum_sq = as.numeric(tab[["Sum Sq"]]),
      mean_sq = as.numeric(tab[["Mean Sq"]]),
      f_value = as.numeric(tab[["F value"]]),
      p_value = as.numeric(tab[["Pr(>F)"]])
    )
    rownames(out) <- NULL
  }

  out$term[out$term == "row_nested"] <- "rowgroup:row"
  out$term[out$term == "col_nested"] <- "colgroup:col"
  out <- arc_clean_anova_numbers(out)

  list(
    table = out,
    warning = warning_text
  )
}

# Sequential (type I) ANOVA from a single QR decomposition: the effects
# Q'y of the fitted QR, grouped by the term of each pivoted column, give the
# sequential sums of squares (the same computation as `anova.lm()`).
arc_sequential_anova_from_fit <- function(fit) {
  term_labels <- attr(fit$terms, "term.labels")
  rank <- fit$rank
  kept <- seq_len(rank)
  term_of_effect <- fit$assign[fit$qr$pivot[kept]]
  effects <- fit$effects[kept]

  term_df <- vapply(seq_along(term_labels), function(i) sum(term_of_effect == i), numeric(1))
  term_ss <- vapply(seq_along(term_labels), function(i) sum(effects[term_of_effect == i]^2), numeric(1))
  residual_df <- length(fit$y) - rank
  residual_ss <- sum(fit$residuals^2)
  residual_mean_sq <- if (residual_df > 0) residual_ss / residual_df else NA_real_

  out <- tibble::tibble(
    term = term_labels,
    df = term_df,
    sum_sq = term_ss,
    mean_sq = ifelse(term_df > 0, term_ss / term_df, NA_real_),
    f_value = NA_real_,
    p_value = NA_real_
  )

  if (!is.na(residual_mean_sq) && residual_mean_sq > 0) {
    testable <- out$df > 0
    out$f_value[testable] <- out$mean_sq[testable] / residual_mean_sq
    out$p_value[testable] <- stats::pf(out$f_value[testable], df1 = out$df[testable],
                                       df2 = residual_df, lower.tail = FALSE)
  }

  dplyr::bind_rows(
    out,
    tibble::tibble(
      term = "Residuals",
      df = residual_df,
      sum_sq = residual_ss,
      mean_sq = residual_mean_sq,
      f_value = NA_real_,
      p_value = NA_real_
    )
  )
}

arc_lm_drop_aliased_terms <- function(formula, data) {
  model_frame <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  response <- stats::model.response(model_frame)
  terms_obj <- stats::terms(formula, data = data, keep.order = TRUE)
  model_matrix <- stats::model.matrix(terms_obj, model_frame)
  full_model_matrix <- model_matrix
  full_assign <- attr(full_model_matrix, "assign")

  qr_obj <- qr(full_model_matrix)
  selected <- sort(qr_obj$pivot[seq_len(qr_obj$rank)])

  model_matrix <- full_model_matrix[, selected, drop = FALSE]
  dropped <- setdiff(seq_len(ncol(full_model_matrix)), selected)

  fit <- stats::lm.fit(x = model_matrix, y = response)
  fit$terms <- terms_obj
  fit$model <- model_frame
  fit$x <- model_matrix
  fit$y <- response
  fit$call <- match.call()
  fit$assign <- full_assign[selected]
  fit$contrasts <- attr(full_model_matrix, "contrasts")
  fit$xlevels <- stats::.getXlevels(terms_obj, model_frame)
  fit$formula <- formula
  fit$na.action <- attr(model_frame, "na.action")
  class(fit) <- "lm"
  attr(fit, "aliased_coefficients_dropped") <- length(dropped)
  attr(fit, "dropped_coefficients") <- colnames(full_model_matrix)[dropped]

  fit
}

arc_format_augmented_anova <- function(anova_table, response) {
  source_map <- data.frame(
    term = c(
      "rowgroup",
      "rowgroup:row",
      "colgroup",
      "colgroup:col",
      "rowgroup:colgroup",
      "genotype",
      "Residuals"
    ),
    source = c(
      "Row groups",
      "Rows, nested within row groups",
      "Column groups",
      "Columns, nested within column groups",
      "Row groups x column groups (blocks)",
      "Genotypes",
      "Error"
    ),
    order = seq_len(7),
    stringsAsFactors = FALSE
  )

  out <- source_map |>
    dplyr::left_join(anova_table, by = "term") |>
    dplyr::arrange(.data$order) |>
    dplyr::transmute(
      Source = .data$source,
      Df = .data$df,
      `Sum Sq` = .data$sum_sq,
      `Mean Sq` = .data$mean_sq,
      `F value` = .data$f_value,
      `Pr(>F)` = .data$p_value,
      Significance = arc_significance_stars(.data$p_value)
    )

  out <- arc_clean_anova_numbers(out)

  corrected_total <- tibble::tibble(
    Source = "Corrected total",
    Df = length(response) - 1,
    `Sum Sq` = sum((response - mean(response))^2),
    `Mean Sq` = NA_real_,
    `F value` = NA_real_,
    `Pr(>F)` = NA_real_,
    Significance = ""
  )

  dplyr::bind_rows(out, corrected_total)
}

# Fit the augmented row-column model
#   y = rowgroup + rowgroup:row + colgroup + colgroup:col + rowgroup:colgroup
#       + genotype + error
# without building a dense genotype design matrix.
#
# An unreplicated genotype has one parameter for one plot, so it fits that plot
# exactly and contributes nothing else to the fit. The full-model residual sum
# of squares and rank therefore come from a small regression on the plots of
# replicated genotypes (the checks) only. The sequential sums of squares of the
# design terms come from one pivoted QR of the design columns on all plots
# (the same decomposition `anova.lm()` uses), and the genotype sum of squares
# is the remaining reduction in residual sum of squares.
arc_fit_augmented_row_column <- function(data) {
  data <- droplevels(data)
  y <- data$response
  n <- length(y)
  data$block_nested <- interaction(data$rowgroup, data$colgroup, drop = TRUE)

  term_labels <- c("rowgroup", "rowgroup:row", "colgroup", "colgroup:col", "rowgroup:colgroup")
  x_design <- arc_factor_model_matrix(
    data[c("rowgroup", "row_nested", "colgroup", "col_nested", "block_nested")]
  )
  design_assign <- attr(x_design, "assign")

  qr_design <- qr(x_design)
  rank_design <- qr_design$rank
  effects <- qr.qty(qr_design, y)
  kept_terms <- design_assign[qr_design$pivot[seq_len(rank_design)]]
  kept_effects <- effects[seq_len(rank_design)]
  design_ss <- vapply(seq_along(term_labels), function(i) sum(kept_effects[kept_terms == i]^2), numeric(1))
  design_df <- vapply(seq_along(term_labels), function(i) sum(kept_terms == i), numeric(1))
  rss_design <- sum(qr.resid(qr_design, y)^2)

  genotype_counts <- table(data$genotype)
  replicated <- data$genotype %in% names(genotype_counts)[genotype_counts > 1]
  n_unreplicated <- sum(!replicated)
  fitted <- y
  rank_replicated <- 0L
  rss <- 0
  model <- NULL

  if (any(replicated)) {
    replicated_genotype <- droplevels(data$genotype[replicated])
    genotype_x <- arc_dummy_matrix(replicated_genotype, prefix = "genotype")
    x_replicated <- cbind(x_design[replicated, , drop = FALSE], genotype_x)
    qr_replicated <- qr(x_replicated)
    rank_replicated <- qr_replicated$rank
    resid_replicated <- qr.resid(qr_replicated, y[replicated])
    rss <- sum(resid_replicated^2)
    fitted[replicated] <- y[replicated] - resid_replicated

    keep <- sort(qr_replicated$pivot[seq_len(rank_replicated)])
    x_keep <- x_replicated[, keep, drop = FALSE]
    y_replicated <- y[replicated]
    model <- stats::lm(y_replicated ~ 0 + x_keep)
    names(model$coefficients) <- colnames(x_keep)
  }

  rank_full <- rank_replicated + n_unreplicated
  error_df <- n - rank_full

  # Adjusted genotype means: the genotype effect plus the average design
  # effect over all observed plots (the prediction for an "average plot").
  # For a replicated genotype this is the estimable function
  # (mean design row, genotype indicator) of the check-plot fit. For an
  # unreplicated entry observed at plot p it is y_p minus the estimated
  # design effect of plot p relative to the average plot, so its variance is
  # sigma^2 plus the variance of that estimated design adjustment.
  adjusted <- tibble::tibble(
    genotype = character(0),
    adjusted_mean = numeric(0),
    std_error = numeric(0),
    estimable = logical(0)
  )
  precision <- list(mean_sed = NA_real_, lsd_5pct = NA_real_)
  vcov <- NULL

  if (!is.null(model)) {
    sigma2 <- if (error_df > 0) rss / error_df else NA_real_
    x_mean <- colMeans(x_design)
    n_design_cols <- ncol(x_design)
    n_check <- ncol(genotype_x)
    check_l <- cbind(
      matrix(x_mean, nrow = n_check, ncol = n_design_cols, byrow = TRUE),
      diag(1, n_check)
    )
    entry_l <- cbind(
      -sweep(x_design[!replicated, , drop = FALSE], 2, x_mean),
      matrix(0, nrow = n_unreplicated, ncol = n_check)
    )
    est <- arc_linear_estimates(
      qr_replicated,
      qr.coef(qr_replicated, y[replicated]),
      rbind(check_l, entry_l)
    )
    extra_var <- c(rep(0, n_check), rep(1, n_unreplicated))
    estimate <- est$estimate + c(rep(0, n_check), y[!replicated])
    std_error <- sqrt(sigma2 * (colSums(est$W^2) + extra_var))
    estimate[!est$estimable] <- NA_real_
    std_error[!est$estimable] <- NA_real_

    adjusted <- tibble::tibble(
      genotype = c(levels(replicated_genotype), as.character(data$genotype[!replicated])),
      adjusted_mean = estimate,
      std_error = std_error,
      estimable = est$estimable
    )
    precision <- arc_mean_precision(
      est$W[, est$estimable, drop = FALSE],
      extra_var[est$estimable],
      sigma2,
      error_df
    )
    # Covariance of the adjusted means (for two-stage weighting):
    # sigma^2 (W'W + diag(extra)), estimable genotypes only.
    ok <- which(est$estimable)
    W <- est$W[, ok, drop = FALSE]
    vcov <- sigma2 * (crossprod(W) + diag(extra_var[ok], length(ok)))
    dimnames(vcov) <- list(adjusted$genotype[ok], adjusted$genotype[ok])
  }
  genotype_df <- rank_full - rank_design
  genotype_ss <- max(rss_design - rss, 0)

  table_out <- tibble::tibble(
    term = c(term_labels, "genotype", "Residuals"),
    df = c(design_df, genotype_df, error_df),
    sum_sq = c(design_ss, genotype_ss, rss)
  )
  table_out$mean_sq <- ifelse(table_out$df > 0, table_out$sum_sq / table_out$df, NA_real_)
  error_mean_sq <- table_out$mean_sq[table_out$term == "Residuals"]
  testable <- table_out$term != "Residuals" & table_out$df > 0 &
    !is.na(error_mean_sq) & error_mean_sq > 0
  table_out$f_value <- ifelse(testable, table_out$mean_sq / error_mean_sq, NA_real_)
  table_out$p_value <- ifelse(
    testable,
    stats::pf(table_out$f_value, df1 = table_out$df, df2 = error_df, lower.tail = FALSE),
    NA_real_
  )

  n_row <- nlevels(data$row_nested)
  n_col <- nlevels(data$col_nested)
  n_rowgroup <- nlevels(data$rowgroup)
  n_colgroup <- nlevels(data$colgroup)
  n_genotype <- nlevels(data$genotype)
  expected_rank <- n_row + n_col - 1 + (n_rowgroup - 1) * (n_colgroup - 1) + n_genotype - 1

  list(
    anova = arc_clean_anova_numbers(table_out),
    model = model,
    fitted = fitted,
    rss = rss,
    rank = rank_full,
    expected_rank = expected_rank,
    error_df = error_df,
    adjusted_means = adjusted,
    vcov = vcov,
    precision = precision,
    # Model columns (intercept, design contrasts, one dummy per genotype) that
    # are aliased and therefore dropped from the fit.
    aliased = as.integer(ncol(x_design) + n_genotype - rank_full)
  )
}

# Estimates of the linear functions `L %*% beta` from the pivoted QR of a
# (possibly rank-deficient) model matrix X. Aliased coefficients are set to
# zero, which is valid for estimable functions. A row l of L is estimable when
# it lies in the row space of X, i.e. in the row space of R = [R11 R12]:
# l = a'R with a = R11^-T l_K must also reproduce the aliased part of l.
# W = R11^-T L_K' gives Cov(L beta_hat) = sigma^2 W'W.
arc_linear_estimates <- function(qr_obj, coef, L, tol = 1e-6) {
  rank <- qr_obj$rank
  p <- ncol(qr_obj$qr)
  kept <- seq_len(rank)
  L_pivoted <- L[, qr_obj$pivot, drop = FALSE]
  R11 <- qr_obj$qr[kept, kept, drop = FALSE]
  W <- backsolve(R11, t(L_pivoted[, kept, drop = FALSE]), transpose = TRUE)
  W <- matrix(W, nrow = rank)

  estimable <- rep(TRUE, nrow(L))
  if (rank < p) {
    aliased <- (rank + 1):p
    implied <- crossprod(W, qr_obj$qr[kept, aliased, drop = FALSE])
    discrepancy <- rowSums(abs(L_pivoted[, aliased, drop = FALSE] - implied))
    estimable <- discrepancy <= tol * pmax(1, rowSums(abs(L)))
  }

  coef[is.na(coef)] <- 0
  list(estimate = drop(L %*% coef), W = W, estimable = estimable)
}

# Precision of a set of adjusted means with covariance
# sigma^2 * (W'W + diag(extra_var)): the root mean variance of all pairwise
# differences (mean standard error of a difference) and the matching 5% LSD.
arc_mean_precision <- function(W, extra_var, sigma2, error_df) {
  m <- ncol(W)
  if (m < 2 || is.na(sigma2)) {
    return(list(mean_sed = NA_real_, lsd_5pct = NA_real_))
  }
  trace_cov <- sum(W^2) + sum(extra_var)
  sum_cov <- sum(rowSums(W)^2) + sum(extra_var)
  mean_var_diff <- sigma2 * 2 * (m * trace_cov - sum_cov) / (m * (m - 1))
  mean_sed <- sqrt(max(mean_var_diff, 0))
  list(
    mean_sed = mean_sed,
    lsd_5pct = if (error_df > 0) stats::qt(0.975, error_df) * mean_sed else NA_real_
  )
}

# Adds adjusted means to a genotype summary and returns the ranked table.
arc_attach_adjusted_means <- function(genotype_summary, adjusted) {
  matched <- match(genotype_summary$genotype, adjusted$genotype)
  genotype_summary$adjusted_mean <- adjusted$adjusted_mean[matched]
  genotype_summary$std_error <- adjusted$std_error[matched]

  ranked <- genotype_summary[order(-genotype_summary$adjusted_mean, na.last = TRUE), ]
  ranked <- ranked[, setdiff(names(ranked), c("fitted_mean", "residual_mean"))]
  ranked$rank <- ifelse(
    is.na(ranked$adjusted_mean),
    NA_integer_,
    rank(-ranked$adjusted_mean, ties.method = "min", na.last = "keep")
  )
  ranked <- ranked[, c("rank", setdiff(names(ranked), "rank"))]

  list(genotype_summary = genotype_summary, adjusted_means = ranked)
}

# Treatment-contrast model matrix for an intercept plus a list of factors, in
# order, with an `assign` attribute mapping columns to factors (0 = intercept).
arc_factor_model_matrix <- function(factors) {
  n <- length(factors[[1]])
  blocks <- lapply(factors, function(f) arc_dummy_matrix(f, drop_first = TRUE))
  x <- do.call(cbind, c(list(`(Intercept)` = rep(1, n)), blocks))
  colnames(x)[1] <- "(Intercept)"
  attr(x, "assign") <- c(0L, rep(seq_along(blocks), vapply(blocks, ncol, integer(1))))
  x
}

arc_dummy_matrix <- function(f, drop_first = FALSE, prefix = NULL) {
  f <- as.factor(f)
  lv <- levels(f)
  x <- matrix(0, nrow = length(f), ncol = length(lv))
  x[cbind(seq_along(f), as.integer(f))] <- 1
  colnames(x) <- paste0(prefix %||% "", lv)
  if (drop_first) x[, -1, drop = FALSE] else x
}

arc_augmented_design_df <- function(data) {
  k <- dplyr::n_distinct(data$row)
  s <- dplyr::n_distinct(data$col)
  gk <- dplyr::n_distinct(data$rowgroup)
  gs <- dplyr::n_distinct(data$colgroup)

  if ("is_control" %in% names(data) && any(!is.na(data$is_control))) {
    vc <- dplyr::n_distinct(data$genotype[data$is_control == 1])
    ve <- dplyr::n_distinct(data$genotype[data$is_control != 1 | is.na(data$is_control)])
  } else if ("type" %in% names(data) && any(!is.na(data$type))) {
    vc <- dplyr::n_distinct(data$genotype[data$type == "check"])
    ve <- dplyr::n_distinct(data$genotype[data$type != "check" | is.na(data$type)])
  } else {
    vc <- 0
    ve <- dplyr::n_distinct(data$genotype)
  }

  genotype_df <- vc + ve - 1
  block_df <- (gk - 1) * (gs - 1)
  corrected_total_df <- k * s - 1
  source_df <- c(
    gk - 1,
    k - gk,
    gs - 1,
    s - gs,
    block_df,
    genotype_df
  )
  error_df <- corrected_total_df - sum(source_df)

  tibble::tibble(
    term = c(
      "rowgroup",
      "rowgroup:row",
      "colgroup",
      "colgroup:col",
      "rowgroup:colgroup",
      "genotype",
      "Residuals"
    ),
    display_df = c(
      source_df,
      error_df
    ),
    k = k,
    s = s,
    gk = gk,
    gs = gs,
    vc = vc,
    ve = ve,
    corrected_total_df = corrected_total_df,
    error_df_formula = paste0(
      "(",
      k,
      " * ",
      s,
      " - 1) - [(",
      gk,
      " - 1) + (",
      k,
      " - ",
      gk,
      ") + (",
      gs,
      " - 1) + (",
      s,
      " - ",
      gs,
      ") + (",
      gk,
      " - 1)(",
      gs,
      " - 1) + (",
      vc,
      " + ",
      ve,
      " - 1)]"
    ),
    previous_error_df_formula = paste0(
      "(",
      k,
      " - 1)(",
      s,
      " - 1) - (",
      gk,
      " - 1)(",
      gs,
      " - 1) - (",
      vc,
      " + ",
      ve,
      " - 1)"
    )
  )
}

arc_significance_stars <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ "",
    p_value <= 0.01 ~ "**",
    p_value <= 0.05 ~ "*",
    TRUE ~ ""
  )
}

arc_clean_anova_numbers <- function(x) {
  numeric_cols <- vapply(x, is.numeric, logical(1))
  x[numeric_cols] <- lapply(x[numeric_cols], function(col) {
    col[is.nan(col)] <- NA_real_
    col
  })
  x
}

arc_augmented_analysis_message <- function(diagnostics) {
  if (isTRUE(diagnostics$anova_error_df <= 0)) {
    return(
      "F-tests are not estimable because the formatted ANOVA has no positive error degrees of freedom. The design/model is saturated for these data."
    )
  }

  if (isTRUE(diagnostics$saturated_or_near_perfect_fit)) {
    return(
      "F-tests should be interpreted cautiously because the model has an essentially perfect fit."
    )
  }

  NA_character_
}

arc_append_estimability_message <- function(message, adjusted) {
  n_missing <- sum(!adjusted$estimable)
  if (n_missing == 0) {
    return(message)
  }
  note <- paste0(
    "Adjusted means are not estimable for ", n_missing, " of ", nrow(adjusted),
    " genotypes because some design effects cannot be estimated from the observed ",
    "plots (for example, a field row or column, or a block, with no observed check). ",
    "These genotypes have NA adjusted means."
  )
  if (is.na(message) || !nzchar(message)) note else paste(message, note)
}

arc_make_plot_order <- function(df, plot_type = "serpentine") {
  if (plot_type == "serpentine") {
    df$plot_order_col <- ifelse(df$row %% 2 == 1, df$col, -df$col)
    df <- df[order(df$row, df$plot_order_col), ]
    df$plot_order_col <- NULL
  } else {
    df <- df[order(df$row, df$col), ]
  }

  rownames(df) <- NULL
  df$plots <- seq_len(nrow(df))
  df
}

arc_validate_check_layout <- function(df) {
  check_df <- df[df$type == "check", c("block", "row", "col", "trt")]

  for (b in sort(unique(check_df$block))) {
    x <- check_df[check_df$block == b, ]

    if (anyDuplicated(x$row)) {
      stop("Invalid layout: block ", b, " has more than one check in the same row.", call. = FALSE)
    }

    if (anyDuplicated(x$col)) {
      stop("Invalid layout: block ", b, " has more than one check in the same column.", call. = FALSE)
    }

    if (anyDuplicated(x$trt)) {
      stop("Invalid layout: block ", b, " has duplicated check labels.", call. = FALSE)
    }
  }

  if (any(table(check_df$trt, check_df$row) > 1)) {
    stop("Invalid layout: the same check appears more than once in the same field row.", call. = FALSE)
  }

  if (any(table(check_df$trt, check_df$col) > 1)) {
    stop("Invalid layout: the same check appears more than once in the same field column.", call. = FALSE)
  }

  invisible(TRUE)
}

arc_design_dims <- function(rows_in_field, cols_in_field, rows_per_block, cols_per_block, n_controls) {
  n_block_rows <- rows_in_field %/% rows_per_block
  n_block_cols <- cols_in_field %/% cols_per_block
  n_slots <- n_controls * n_block_cols * n_block_rows

  # Fixed slot structure: one slot per (control, column group, row group).
  slot_control <- rep.int(seq_len(n_controls), n_block_cols * n_block_rows)
  slot_colgroup <- rep.int(rep(seq_len(n_block_cols), each = n_controls), n_block_rows)
  slot_rowgroup <- rep(seq_len(n_block_rows), each = n_controls * n_block_cols)

  list(
    rows_in_field = rows_in_field,
    cols_in_field = cols_in_field,
    rows_per_block = rows_per_block,
    cols_per_block = cols_per_block,
    n_controls = n_controls,
    n_block_rows = n_block_rows,
    n_block_cols = n_block_cols,
    n_slots = n_slots,
    slot_control = slot_control,
    slot_rowgroup = slot_rowgroup,
    slot_colgroup = slot_colgroup,
    slot_block = (slot_rowgroup - 1L) * n_block_cols + slot_colgroup,
    row_index = cbind(slot_control, slot_colgroup, slot_rowgroup),
    col_index = cbind(slot_control, slot_rowgroup, slot_colgroup),
    # Rank of [row, column, block, check] dummies on the check plots when every
    # row, column, and block x block interaction effect is estimable.
    expected_check_rank = rows_in_field + cols_in_field - 1L +
      (n_block_rows - 1L) * (n_block_cols - 1L) + n_controls - 1L
  )
}

# Random check coordinates for one candidate layout.
#
# For row group r, `row_layouts[c, g, r]` is the local row of control c in the
# block at column group g. Each layout is a Latin rectangle (distinct entries in
# every row and column), so within a block the checks use different rows, and a
# control never reuses a field row across the blocks of a row group. Columns
# work the same way with `col_layouts[c, r, g]`. The rectangles also use every
# local row/column, so every field row and column receives a check.
arc_make_check_slots <- function(dims) {
  nc <- dims$n_controls
  nbr <- dims$n_block_rows
  nbc <- dims$n_block_cols

  row_layouts <- unlist(lapply(
    seq_len(nbr),
    function(i) arc_covering_latin_rectangle(nc, nbc, dims$rows_per_block)
  ))
  col_layouts <- unlist(lapply(
    seq_len(nbc),
    function(i) arc_covering_latin_rectangle(nc, nbr, dims$cols_per_block)
  ))
  dim(row_layouts) <- c(nc, nbc, nbr)
  dim(col_layouts) <- c(nc, nbr, nbc)

  list(
    row = (dims$slot_rowgroup - 1L) * dims$rows_per_block + row_layouts[dims$row_index],
    col = (dims$slot_colgroup - 1L) * dims$cols_per_block + col_layouts[dims$col_index]
  )
}

# Random n_rows x n_cols Latin rectangle on symbols 1..n_symbols that uses
# every symbol (requires n_rows, n_cols <= n_symbols <= n_rows * n_cols).
#
# Built as a sub-array of the cyclic Latin square (i + j) mod n with row offsets
# {0, ..., a - 1} and column offsets {0, a, 2a, ...} (plus random extra
# columns), whose sums cover every residue mod n. Random symbol, row, and
# column permutations then randomize it while keeping both properties.
# Latin rectangle of symbols (rows within a block) for checks (rows) by
# blocks (columns): symbols are distinct within each column and within each
# row. Each of the longer dimension's lines takes a cyclic run of symbols
# starting at evenly spaced offsets, so every symbol is used either
# floor(n / n_symbols) or ceiling(n / n_symbols) times: rows (or columns)
# within a group receive equal numbers of check plots. Random rotation,
# symbol relabelling and row/column order keep the layout randomized.
arc_covering_latin_rectangle <- function(n_rows, n_cols, n_symbols) {
  transpose <- n_rows > n_cols
  a <- if (transpose) n_cols else n_rows
  b <- if (transpose) n_rows else n_cols

  starts <- (sample.int(n_symbols, 1L) - 1L + floor((seq_len(b) - 1L) * n_symbols / b)) %% n_symbols
  offsets <- outer(seq_len(a) - 1L, starts, "+") %% n_symbols

  if (transpose) {
    offsets <- t(offsets)
  }

  symbols <- sample.int(n_symbols)
  out <- matrix(symbols[offsets + 1L], nrow = n_rows, ncol = n_cols)
  out[sample.int(n_rows), sample.int(n_cols), drop = FALSE]
}

# Lower is better: penalize adjacent checks and unequal check counts across
# field rows and columns. Row/column/block constraints hold by construction.
# Moves checks apart by swapping two rows within a row group or two columns
# within a column group. Such swaps keep every property of the Latin
# construction (one of each check per block in distinct rows and columns, no
# check twice in a field row or column, estimability) and only change which
# plots are neighbours; a swap is kept when it does not worsen the score.
arc_polish_latin <- function(slots, dims, n_iter = 2000L) {
  score <- arc_score_check_slots(slots, dims)
  rb <- dims$rows_per_block
  cb <- dims$cols_per_block
  for (it in seq_len(n_iter)) {
    trial <- slots
    if (stats::runif(1) < 0.5 && rb > 1) {
      g <- sample.int(dims$n_block_rows, 1L) - 1L
      ab <- g * rb + sample.int(rb, 2L)
      hit_a <- slots$row == ab[1]
      hit_b <- slots$row == ab[2]
      trial$row[hit_a] <- ab[2]
      trial$row[hit_b] <- ab[1]
    } else if (cb > 1) {
      g <- sample.int(dims$n_block_cols, 1L) - 1L
      ab <- g * cb + sample.int(cb, 2L)
      hit_a <- slots$col == ab[1]
      hit_b <- slots$col == ab[2]
      trial$col[hit_a] <- ab[2]
      trial$col[hit_b] <- ab[1]
    } else {
      next
    }
    trial_score <- arc_score_check_slots(trial, dims)
    if (trial_score <= score) {
      slots <- trial
      score <- trial_score
    }
  }
  slots
}

arc_score_check_slots <- function(slots, dims) {
  nr <- dims$rows_in_field
  nc <- dims$cols_in_field
  occupied <- matrix(FALSE, nrow = nr, ncol = nc)
  occupied[cbind(slots$row, slots$col)] <- TRUE

  pairs_at <- function(dr, dc) {
    if (dr >= nr || abs(dc) >= nc) {
      return(0)
    }
    a <- occupied[seq_len(nr - dr), max(1, 1 - dc):min(nc, nc - dc), drop = FALSE]
    b <- occupied[(1 + dr):nr, max(1, 1 + dc):min(nc, nc + dc), drop = FALSE]
    sum(a & b)
  }

  distance_1 <- pairs_at(1, 0) + pairs_at(0, 1)
  distance_2 <- pairs_at(2, 0) + pairs_at(0, 2) + pairs_at(1, 1) + pairs_at(1, -1)

  row_counts <- tabulate(slots$row, nbins = nr)
  col_counts <- tabulate(slots$col, nbins = nc)

  500 * distance_1 +
    100 * distance_2 +
    50 * sum((row_counts - mean(row_counts))^2) +
    50 * sum((col_counts - mean(col_counts))^2)
}

# TRUE when row, column, and block effects are all estimable from the checks,
# i.e. the check-only model has the full rank implied by the design.
arc_checks_estimable <- function(slots, dims) {
  n <- dims$n_slots
  idx <- seq_len(n)
  n_blocks <- dims$n_block_rows * dims$n_block_cols
  offsets <- cumsum(c(0L, dims$rows_in_field, dims$cols_in_field, n_blocks))
  x <- matrix(0, nrow = n, ncol = offsets[4] + dims$n_controls)
  x[cbind(idx, slots$row)] <- 1
  x[cbind(idx, offsets[2] + slots$col)] <- 1
  x[cbind(idx, offsets[3] + dims$slot_block)] <- 1
  x[cbind(idx, offsets[4] + dims$slot_control)] <- 1
  qr(x)$rank == dims$expected_check_rank
}

arc_build_design <- function(slots, dims, treatments, controls, plot_type) {
  field <- arc_make_field_template(
    rows_in_field = dims$rows_in_field,
    cols_in_field = dims$cols_in_field,
    rows_per_block = dims$rows_per_block,
    cols_per_block = dims$cols_per_block
  )

  check_idx <- (slots$col - 1L) * dims$rows_in_field + slots$row
  control_labels <- sample(controls)
  field$trt[check_idx] <- control_labels[dims$slot_control]
  field$type[check_idx] <- "check"

  entry_idx <- which(field$type == "entry")
  field$trt[entry_idx] <- treatments[sample.int(length(treatments))]

  arc_validate_check_layout(field)

  field <- arc_make_plot_order(field, plot_type = plot_type)
  field$rep <- ifelse(field$type == "check", field$block, 1L)
  field$is_control <- ifelse(field$type == "check", 1L, 0L)

  field
}
