#' Suggest trial layouts from the number of entries and available field space
#'
#' Lists parameter options for the package's three design families, starting
#' with the information the breeder knows. Rows, columns and check counts
#' can be left unspecified. No randomization or optimization is performed.
#'
#' @param n_entries Number of test entries, excluding checks.
#' @param design Design families to consider: `"augmented_row_column"` (augmented
#'   row-column), `"prep"` (partially replicated), and `"resolvable"`
#'   (fully replicated row-column). Defaults to all three.
#' @param n_checks Distinct check genotypes to consider, a non-negative
#'   integer or vector. `NULL` tries 2 through 6 for augmented row-column designs and
#'   zero for the other families. p-rep uses no separate checks; it is
#'   excluded when only positive check counts are requested. Resolvable
#'   designs replicate checks just like the other genotypes.
#' @param rows_in_field,cols_in_field Known dimensions of the whole field at
#'   one location (plot rows and plot columns, including all checks and
#'   replicates), not of a block. Block dimensions are suggested in the
#'   result as `rows_per_block` and `cols_per_block`. Each can be a positive
#'   integer, a vector of dimensions to consider, or `NULL` to search
#'   automatically. Supplied dimensions are exact candidate sizes, not upper
#'   bounds. When only one is supplied, the square field (the same number of
#'   rows and columns) is also considered. For p-rep, the last row may be
#'   partially occupied, but wholly empty final rows are not proposed.
#' @param n_locations Number of locations. Augmented row-column and resolvable layouts
#'   describe one field to randomize independently at each location. A
#'   p-rep option describes the combined multi-location design.
#' @param n_reps Full replication counts to try for resolvable designs
#'   (at least 2). Other design families use their own replication scheme.
#' @param block_size Maximum incomplete block sizes to try for p-rep.
#' @param prep_prop Proportions of entries to duplicate at each p-rep
#'   location. Counts use `round(n_entries * prep_prop)`, as in the generator.
#'   An entry may be duplicated at only one location, so the sum of these
#'   counts over locations must not exceed `n_entries`. `NULL` tries 15%,
#'   20% and 25%, and, for several locations, a balanced share of entries
#'   per location rounded down to respect this budget.
#' @param max_plots Largest field area (rows times columns) per location.
#'   Defaults to four times the number of genotypes (increased when more
#'   full replicates are requested), or the largest supplied field area
#'   when both dimensions are given.
#' @param max_dimension Largest automatically searched row or column
#'   dimension. Defaults to twice the square root of `max_plots`, with a
#'   minimum of 20. Explicit dimensions are not restricted by this bound.
#'   When neither dimension is supplied, fields with an aspect ratio above
#'   3 are excluded. Supply a dimension explicitly to consider longer fields.
#' @param min_error_df Minimum nominal error degrees of freedom. These are
#'   count-based values for the design model, assuming its effects are
#'   estimable, not degrees of freedom of a later spatial model. `NULL`
#'   (default) requires at least 1, except for the layout of Fig. 2 of
#'   Piepho and Williams (2016): blocks of `n_checks` x `n_checks` plots with
#'   each check exactly once, such as 36 entries and three checks in a 9 x 6
#'   field. That layout is saturated with fixed rows, columns and blocks
#'   (error df = 0) and is analyzed with
#'   `analyze_augmented_row_column_design(random_blocks = TRUE)` or
#'   `entries = "random"`. A supplied value applies to every layout; 0 also
#'   admits other saturated layouts.
#' @param n_options Maximum options returned per design family. The
#'   shortlist includes different check, full-replication or p-rep counts
#'   before adding further layouts with the same counts.
#' @param n_candidates Search candidates to pass to the augmented row-column
#'   and p-rep generators. Defaults to 200. This is a user-selected search
#'   effort, not a dimension inferred from entry counts.
#' @param rows,cols Deprecated names of `rows_in_field` and `cols_in_field`.
#'
#' @details
#' This is a parameter planner. Candidates pass plot-count, block-size and
#' nominal error-df checks. Augmented row-column candidates also pass the check
#' coverage and row/column group conditions used by the generator. p-rep
#' candidates have enough repeated plots to connect their blocks and a
#' pair-count budget compatible with (0,1) concurrence. These necessary
#' conditions do not certify a randomized design's connectivity or
#' efficiency: generate the chosen option, inspect its summary and run
#' [check_design()]. The search is bounded and not exhaustive.
#'
#' Augmented row-column options with equal row and column groups have
#' blocks of at least `n_checks` rows and at least `n_checks` columns,
#' leaving room for one copy of each check in a distinct row and column. For
#' example, four checks require blocks of at least 4 rows by 4 columns. Every
#' such block contains all checks at least once; extra positions are assigned
#' to extra check replication, so the number of test-entry positions is
#' exactly `n_entries`. Flexible allocations use `method = "optimal"`.
#'
#' Options with unequal row or column groups (Piepho and Williams, 2016,
#' Example 3 and Fig. 4) split the field into near-equal groups, such as
#' rows 6, 5, 6, 5 and columns 3, 4, 3 for a 22 x 10 field. Their blocks hold
#' `n_checks - 1` to `n_checks + 1` check plots, in proportion to block size,
#' so some blocks may be incomplete for the checks. These options report
#' `rows_per_rowgroup`, `cols_per_colgroup` and `checks_per_block`, with
#' `rows_per_block` and `cols_per_block` set to `NA`.
#'
#' Resolvable options may contain filler positions, which are reported in
#' `unused_plots`. Include both entries and checks in the generator's
#' `treatments`. Its native generator treats them all as treatments; if
#' checks need to be marked in Breedbase, set `is_control` in both the
#' generated `design` and `breedbase_design` tables before saving.
#'
#' Square blocks are preferred but not required: for augmented row-column
#' designs, when a layout with blocks of `n_checks` rows by `n_checks`
#' columns containing each check exactly once exists (for example 3 x 3
#' blocks with checks A, B and C, as in Fig. 2 of Piepho and Williams, 2016),
#' the smallest one is the first option and `square_check_blocks` flags it.
#' Blocks of that size with extra check plots are not flagged. The other
#' options follow in the usual order: within each
#' resource-count group, smaller field areas and fewer unused positions come
#' first. Within a field, augmented row-column layouts whose blocks hold about
#' one check plot per check (complete blocks for the checks, as in Piepho and
#' Williams, 2016) come first, then Latin constructions, squarer blocks and
#' less elongated fields.
#' This ordering does not estimate statistical precision or choose the
#' biologically best design.
#'
#' @return A tibble with option id, design family, entry/check counts,
#'   field dimensions, a `square_check_blocks` flag, location and replication counts, active plots per
#'   location and in total, field area, unused positions, check plots,
#'   repeated entries, blocks, nominal `error_df`, method and notes.
#'   Field dimensions are named `rows_in_field` and `cols_in_field`, matching
#'   [augmented_row_column_design()]. Block and replicate dimensions are
#'   included where applicable. `n_candidates` reports the chosen search
#'   effort. `checks_per_block` is a list-column with each augmented
#'   row-column option's check-allocation matrix, or `NULL` for the default
#'   complete check blocks. These columns support an explicit generator call.
#'   `rows_per_rowgroup` and `cols_per_colgroup` are list-columns with the
#'   group sizes of unequal-group options, or `NULL`.
#'   `generator` names the function to call; `arguments` is a list-column of
#'   its layout arguments. Add your own `treatments`, `controls` (augmented row-column)
#'   and `seed` before calling it with `do.call()`. p-rep arguments include
#'   placeholder location names, which can be replaced with your own names.
#'   Attribute `diagnostics` reports the number of shortlisted options and
#'   the reason for any excluded family. No matches return a typed empty
#'   tibble and a message suggesting which constraints to relax.
#'
#' @seealso [augmented_row_column_layouts()], [augmented_prep_parameters()],
#'   [augmented_row_column_design()], [augmented_prep_design()],
#'   [resolvable_row_column_design()]
#' @examples
#' # Only the number of entries is known:
#' options <- design_options(100)
#' options[, c("option", "design", "n_checks", "rows_in_field", "cols_in_field", "plots_per_location")]
#'
#' # A known field, four checks and exactly 108 test entries:
#' options <- design_options(108, design = "augmented_row_column", n_checks = 4,
#'                           rows_in_field = 12, cols_in_field = 12, n_candidates = 5)
#' selected <- options[1, ]
#' trial <- augmented_row_column_design(
#'   treatments = paste0("LINE_", seq_len(108)),
#'   controls = paste0("CHECK_", seq_len(selected$n_checks)),
#'   rows_in_field = selected$rows_in_field,
#'   cols_in_field = selected$cols_in_field,
#'   rows_per_block = selected$rows_per_block,
#'   cols_per_block = selected$cols_per_block,
#'   checks_per_block = selected$checks_per_block[[1]],
#'   method = selected$method,
#'   n_candidates = selected$n_candidates,
#'   seed = 2026
#' )
#'
#' # Fully replicated options with one known dimension:
#' design_options(42, design = "resolvable", cols_in_field = 7, n_reps = 3)
#'
#' # Only the number of field rows is known. With three checks, the first
#' # option is the layout of Piepho and Williams (2016, Fig. 2): a 9 x 6 field
#' # of six 3 x 3 blocks, each check once per block (error df = 0).
#' design_options(36, design = "augmented_row_column", n_checks = 3, rows_in_field = 9)
#' # Require positive error df for every layout:
#' design_options(36, design = "augmented_row_column", n_checks = 3, rows_in_field = 9,
#'                min_error_df = 1)
#' # Partial replication at three locations:
#' design_options(100, design = "prep", n_locations = 3)
#' @export
design_options <- function(n_entries, design = c("augmented_row_column", "prep", "resolvable"),
                           n_checks = NULL, rows_in_field = NULL, cols_in_field = NULL,
                           n_locations = 1L, n_reps = 2:3,
                           block_size = c(5L, 8L, 10L, 12L, 16L, 20L),
                           prep_prop = NULL,
                           max_plots = NULL, max_dimension = NULL,
                           min_error_df = NULL, n_options = 6L, n_candidates = 200L,
                           rows = NULL, cols = NULL) {
  integers <- function(x, name, minimum = 1L, scalar = FALSE) {
    if (!is.numeric(x) || !length(x) || (scalar && length(x) != 1L) ||
        anyNA(x) || any(!is.finite(x) | x < minimum | x > .Machine$integer.max | x != floor(x))) {
      stop("`", name, "` must contain ", if (scalar) "one ", "integer",
           if (!scalar) "s", " >= ", minimum, ".", call. = FALSE)
    }
    sort(unique(as.integer(x)))
  }
  n_entries <- integers(n_entries, "n_entries", scalar = TRUE)
  if ("augmented" %in% design) {
    stop("Use `design = \"augmented_row_column\"` for augmented row-column designs.", call. = FALSE)
  }
  design <- match.arg(design, c("augmented_row_column", "prep", "resolvable"), several.ok = TRUE)
  if (!is.null(n_checks)) n_checks <- integers(n_checks, "n_checks", minimum = 0L)
  if (!is.null(rows) || !is.null(cols)) {
    if ((!is.null(rows) && !is.null(rows_in_field)) || (!is.null(cols) && !is.null(cols_in_field))) {
      stop("Use `rows_in_field` and `cols_in_field` only; `rows` and `cols` are deprecated names.",
           call. = FALSE)
    }
    warning("`rows` and `cols` are deprecated; use `rows_in_field` and `cols_in_field` ",
            "(dimensions of the whole field, not of a block).", call. = FALSE)
    rows_in_field <- rows_in_field %||% rows
    cols_in_field <- cols_in_field %||% cols
  }
  rows <- if (!is.null(rows_in_field)) integers(rows_in_field, "rows_in_field")
  cols <- if (!is.null(cols_in_field)) integers(cols_in_field, "cols_in_field")
  n_locations <- integers(n_locations, "n_locations", scalar = TRUE)
  n_reps <- integers(n_reps, "n_reps", minimum = 2L)
  block_size <- integers(block_size, "block_size", minimum = 2L)
  # Default: at least 1 error df, except the saturated layout of Piepho and
  # Williams (2016, Fig. 2), with each check exactly once per n_checks x
  # n_checks block, which is analyzed with random blocks.
  allow_fig2 <- is.null(min_error_df)
  if (allow_fig2) min_error_df <- 1L
  min_error_df <- integers(min_error_df, "min_error_df", scalar = TRUE, minimum = 0L)
  n_options <- integers(n_options, "n_options", scalar = TRUE)
  n_candidates <- integers(n_candidates, "n_candidates", scalar = TRUE)
  if (is.null(prep_prop)) {
    prep_prop <- c(0.15, 0.20, 0.25)
    if (n_locations > 1L && n_entries >= n_locations) {
      prep_prop <- unique(c(prep_prop, floor(n_entries / n_locations) / n_entries))
    }
  }
  if (!is.numeric(prep_prop) || !length(prep_prop) || anyNA(prep_prop) ||
      any(!is.finite(prep_prop) | prep_prop <= 0 | prep_prop > 1)) {
    stop("`prep_prop` must contain proportions greater than 0 and no larger than 1.", call. = FALSE)
  }
  max_checks <- if (is.null(n_checks)) 6L else max(n_checks)
  default_plots <- is.null(max_plots)
  if (is.null(max_plots)) {
    multiplier <- if ("resolvable" %in% design) max(4, as.double(n_reps) + 1) else 4
    max_plots <- if (!is.null(rows) && !is.null(cols)) as.double(max(rows)) * max(cols) else
      multiplier * (as.double(n_entries) + max_checks)
  }
  max_plots <- integers(max_plots, "max_plots", scalar = TRUE)
  if (xor(is.null(rows), is.null(cols)) && default_plots) {
    # With one known dimension, always test the square field too.
    known <- as.double(max(c(rows, cols)))
    if (known * known >= n_entries) max_plots <- max(max_plots, as.integer(known * known))
  }
  if (is.null(max_dimension)) max_dimension <- max(20, ceiling(2 * sqrt(max_plots)))
  max_dimension <- integers(max_dimension, "max_dimension", scalar = TRUE)
  fields <- do_option_fields(rows, cols, n_entries, max_plots, max_dimension)
  pieces <- list()
  diagnostics <- list()
  for (kind in unique(design)) {
    checks <- if (is.null(n_checks)) if (kind == "augmented_row_column") 2:6 else 0L else n_checks
    why <- "No layout within the supplied dimensions, field-area budget and error-df limit."
    choices <- switch(kind,
      augmented_row_column = if (any(checks > 0L)) {
        why <- paste0("No augmented row-column layout meets the supplied limits with at least ",
                       "n_checks rows and n_checks columns in every block.")
        do_option_augmented(n_entries, checks[checks > 0L], fields, n_locations, min_error_df, n_options,
                            allow_fig2)
      } else {
        why <- "Augmented row-column designs require at least one distinct check."
        NULL
      },
      prep = if (0L %in% checks && n_entries >= 2L) {
        do_option_prep(n_entries, fields, n_locations, block_size, prep_prop, min_error_df, n_options)
      } else {
        why <- "p-rep needs at least two entries and uses no separate checks; allow n_checks = 0."
        NULL
      },
      resolvable = do_option_resolvable(n_entries, checks, fields, n_locations, n_reps, min_error_df, n_options)
    )
    choices <- dplyr::bind_rows(choices)
    if (nrow(choices)) {
      if (kind %in% c("augmented_row_column", "prep")) {
        choices$n_candidates <- n_candidates
        choices$arguments <- lapply(choices$arguments, function(args) {
          args$n_candidates <- n_candidates
          args
        })
      }
      group <- paste(choices$n_checks, choices$n_reps, choices$n_repeated, sep = "/")
      # Within a field: blocks holding about one plot per check (complete
      # blocks for the checks), Latin constructions, then squarer blocks.
      per_block <- if (kind == "augmented_row_column") abs(choices$check_plots / choices$n_blocks - choices$n_checks) else rep(0, nrow(choices))
      block_rows <- choices$rows_in_field / vapply(choices$arguments, function(a)
        if (!is.null(a$rows_per_rowgroup)) length(a$rows_per_rowgroup) else
          if (!is.null(a$rows_per_block)) a$rows_in_field / a$rows_per_block else 1, numeric(1))
      block_cols <- choices$cols_in_field / vapply(choices$arguments, function(a)
        if (!is.null(a$cols_per_colgroup)) length(a$cols_per_colgroup) else
          if (!is.null(a$cols_per_block)) a$cols_in_field / a$cols_per_block else 1, numeric(1))
      choices <- choices[order(choices$field_plots, choices$unused_plots, per_block,
        choices$method != "latin", abs(log(block_rows / block_cols)),
        abs(log(choices$rows_in_field / choices$cols_in_field)), -choices$error_df), ]
      # Preference, not a requirement: when possible, the best option has
      # square blocks of n_checks rows by n_checks columns.
      group <- paste(choices$n_checks, choices$n_reps, choices$n_repeated, sep = "/")
      best <- choices$square_check_blocks & !duplicated(paste(group, choices$square_check_blocks))
      choices <- choices[order(!best), ]
      group <- paste(choices$n_checks, choices$n_reps, choices$n_repeated, sep = "/")
      within_group <- stats::ave(seq_len(nrow(choices)), group, FUN = seq_along)
      choices <- choices[utils::head(order(within_group, seq_len(nrow(choices))), n_options), ]
      pieces[[kind]] <- choices
    }
    diagnostics[[kind]] <- tibble::tibble(design = kind, options = nrow(choices),
                                          reason = if (nrow(choices)) "" else why)
  }
  out <- dplyr::bind_rows(pieces)
  if (!nrow(out)) {
    out <- do_option_row()[FALSE, ]
    message("No design options found. Try other field dimensions, a larger `max_plots` or `max_dimension`, ",
            "more replication, or a smaller `min_error_df`. See attr(result, 'diagnostics').")
  }
  out$option <- seq_len(nrow(out))
  out <- out[, c("option", setdiff(names(out), "option"))]
  attr(out, "diagnostics") <- dplyr::bind_rows(diagnostics)
  out
}

do_option_fields <- function(rows, cols, n_entries, max_plots, max_dimension) {
  both_unknown <- is.null(rows) && is.null(cols)
  if (is.null(rows)) rows <- seq_len(min(max_dimension, max_plots))
  choices <- lapply(rows, function(r) {
    candidates <- cols
    if (is.null(candidates)) {
      lo <- max(1L, ceiling(n_entries / r))
      hi <- min(max_dimension, floor(max_plots / r))
      if (lo > hi) return(NULL)
      candidates <- seq.int(lo, hi)
    }
    keep <- as.double(r) * candidates >= n_entries & as.double(r) * candidates <= max_plots
    if (both_unknown) keep <- keep & pmax(r / candidates, candidates / r) <= 3
    if (!any(keep)) return(NULL)
    data.frame(rows = r, cols = candidates[keep])
  })
  out <- dplyr::bind_rows(choices)
  if (!nrow(out)) return(data.frame(rows = integer(), cols = integer()))
  out[order(as.double(out$rows) * out$cols, abs(log(out$rows / out$cols))), ]
}

do_option_row <- function(design = "", n_entries = 0L, n_checks = 0L, rows = 0L, cols = 0L,
                           n_locations = 1L, n_reps = NA_integer_, plots = 0L, check_plots = 0L,
                           repeated = NA_integer_, blocks = 0L, error_df = 0L, method = "",
                           generator = "", arguments = list(), note = "") {
  tibble::tibble(design = design, n_entries = n_entries, n_checks = n_checks,
    rows_in_field = rows, cols_in_field = cols,
    square_check_blocks = isTRUE(n_checks > 0L && identical(as.integer(arguments$rows_per_block), as.integer(n_checks)) &&
                                   identical(as.integer(arguments$cols_per_block), as.integer(n_checks)) &&
                                   is.null(arguments$checks_per_block)),
    n_locations = n_locations, n_reps = n_reps,
    rows_per_block = arguments$rows_per_block %||% NA_integer_,
    rows_per_rowgroup = list(arguments$rows_per_rowgroup),
    cols_per_colgroup = list(arguments$cols_per_colgroup),
    cols_per_block = arguments$cols_per_block %||% NA_integer_,
    rows_per_rep = arguments$rows_per_rep %||% NA_integer_,
    cols_per_rep = arguments$cols_per_rep %||% NA_integer_,
    block_size = arguments$block_size %||% NA_integer_,
    n_candidates = NA_integer_, checks_per_block = list(arguments$checks_per_block),
    plots_per_location = plots, total_plots = as.double(plots) * n_locations,
    field_plots = as.double(rows) * cols, unused_plots = as.double(rows) * cols - plots,
    check_plots = check_plots, n_repeated = repeated, n_blocks = blocks,
    error_df = error_df, method = method, generator = generator,
    arguments = list(arguments), note = note)
}

do_option_divisors <- function(n) {
  small <- seq_len(floor(sqrt(n)))
  small <- small[n %% small == 0L]
  sort(unique(c(small, n %/% small)))
}

do_option_augmented <- function(n, checks, fields, locations, min_df, limit, allow_fig2 = FALSE) {
  if (!length(checks)) return(NULL)
  choices <- list()
  seen_checks <- integer()
  seen_square <- integer()
  stop_area <- NA_real_
  for (i in seq_len(nrow(fields))) {
    r <- fields$rows[i]; c <- fields$cols[i]
    p <- as.double(r) * c - n
    if (p < r + c - 2 + min(checks) + if (allow_fig2) 0 else min_df) next
    rb <- do_option_divisors(r); rb <- rb[rb >= 2L]
    cb <- do_option_divisors(c); cb <- cb[cb >= 2L]
    grid <- expand.grid(rb = rb, cb = cb, checks = checks)
    if (!is.na(stop_area)) {
      # Past the usual search: only n_checks x n_checks blocks still missing.
      grid <- grid[grid$rb == grid$checks & grid$cb == grid$checks &
                     !grid$checks %in% seen_square, , drop = FALSE]
    }
    nr <- r / grid$rb; nc <- c / grid$cb; b <- nr * nc
    df <- p - grid$checks - (r - 1) - (c - 1) - (nr - 1) * (nc - 1)
    # Each check exactly once in every n_checks x n_checks block.
    fig2 <- grid$rb == grid$checks & grid$cb == grid$checks & p == b * grid$checks
    grid <- grid[grid$rb >= grid$checks & grid$cb >= grid$checks &
                   (df >= min_df | (allow_fig2 & fig2 & df >= 0)) & floor(p / b) >= grid$checks &
                   ceiling(p / b) < grid$rb * grid$cb, , drop = FALSE]
    for (j in seq_len(nrow(grid))) {
      g <- grid[j, ]; nr <- r %/% g$rb; nc <- c %/% g$cb; b <- nr * nc
      pb <- matrix(p %/% b, nr, nc)
      extra <- matrix(0L, nr, nc)
      for (k in seq_len(p %% b)) {
        slots <- which(extra == 0L, arr.ind = TRUE)
        score <- rowSums(extra)[slots[, 1]] * (b + 1) + colSums(extra)[slots[, 2]]
        pos <- slots[which.min(score), ]
        extra[pos[1], pos[2]] <- 1L
      }
      pb <- pb + extra
      cap <- augmented_row_column_capacity(r, c, g$rb, g$cb, g$checks, checks_per_block = pb)
      if (!cap$every_row_has_check || !cap$every_col_has_check ||
          cap$min_rowgroup_slack < 0 || cap$min_colgroup_slack < 0) next
      regular <- all(pb == g$checks)
      latin <- regular && nc <= g$rb && nr <= g$cb && g$checks <= min(g$rb, g$cb)
      method <- if (latin) "latin" else "optimal"
      args <- list(rows_in_field = r, cols_in_field = c, rows_per_block = g$rb,
                    cols_per_block = g$cb, method = method)
      if (!regular) args$checks_per_block <- pb
      choices[[length(choices) + 1L]] <- do_option_row("augmented_row_column", n, g$checks, r, c, locations,
        1L, r * c, p, blocks = b, error_df = cap$error_df, method = method,
        generator = "augmented_row_column_design", arguments = args,
        note = if (regular && cap$error_df == 0) paste0("Every check once per block. Error df = 0 with fixed ",
          "rows, columns and blocks (Piepho and Williams, 2016, Fig. 2): analyze with random_blocks = TRUE.") else
          if (regular) "Every check once per block." else
          "Extra check replication fills the field exactly; allocation needs optimization.")
      seen_checks <- union(seen_checks, g$checks)
      if (regular && g$rb == g$checks && g$cb == g$checks) seen_square <- union(seen_square, g$checks)
    }
    if (is.na(stop_area)) {
      unequal <- do_option_unequal(n, checks, r, c, p, locations, min_df)
      choices <- c(choices, unequal)
      seen_checks <- union(seen_checks, vapply(unequal, function(x) x$n_checks, integer(1)))
    }
    # Search the smallest fields first; keep a pool larger than the shortlist.
    wanted <- utils::head(checks, limit)
    if (is.na(stop_area) && length(choices) >= 10L * limit && all(wanted %in% seen_checks)) {
      stop_area <- as.double(r) * c
    }
    if (!is.na(stop_area)) {
      # Look a little further for n_checks x n_checks blocks, without
      # promoting a much larger field.
      rest <- fields[-seq_len(i), , drop = FALSE]
      rest <- rest[as.double(rest$rows) * rest$cols <= 1.5 * stop_area, , drop = FALSE]
      # Each check once per k x k block: k check plots in every k^2 plots.
      possible <- vapply(wanted, function(k) any(rest$rows %% k == 0L & rest$cols %% k == 0L &
        as.double(rest$rows) * rest$cols * (k - 1) == as.double(n) * k), logical(1))
      if (all(wanted %in% seen_square | !possible)) break
    }
  }
  choices
}

# Near-equal split of `total` rows (or columns) into `g` groups, larger
# groups spread out, e.g. 22 rows in 4 groups -> 6, 5, 6, 5 and 10 columns
# in 3 groups -> 3, 4, 3, as in Piepho and Williams (2016, Fig. 4).
do_option_split <- function(total, g) {
  sizes <- rep(total %/% g, g)
  extra <- total %% g
  if (extra) {
    step <- g / extra
    sizes[floor((seq_len(extra) - 1) * step + (step - 1) / 2) + 1L] <- sizes[1] + 1L
  }
  as.integer(sizes)
}

# Check plots per block proportional to block area (largest remainder),
# as in Piepho and Williams (2016, Fig. 4).
do_option_allocate <- function(p, rows, cols) {
  area <- outer(rows, cols)
  quota <- p * area / sum(area)
  pb <- floor(quota)
  left <- p - sum(pb)
  if (left > 0) {
    take <- order(-(quota - pb), seq_along(quota))[seq_len(left)]
    pb[take] <- pb[take] + 1
  }
  matrix(as.integer(pb), length(rows), length(cols))
}

# Augmented row-column layouts with unequal row or column groups (Piepho and
# Williams, 2016, Example 3). Blocks hold about n_checks check plots each
# (n_checks - 1 to n_checks + 1), so some may be incomplete for the checks.
do_option_unequal <- function(n, checks, r, c, p, locations, min_df) {
  splits <- function(total, k) {
    out <- list()
    for (g in seq_len(total %/% 2L)) {
      size <- total / g
      if (size < max(2, ceiling(k / 2)) || size > 2 * k) next
      out[[length(out) + 1L]] <- do_option_split(total, g)
    }
    out
  }
  choices <- list()
  for (k in checks) {
    row_splits <- splits(r, k); col_splits <- splits(c, k)
    for (rs in row_splits) for (cs in col_splits) {
      if (length(unique(rs)) == 1L && length(unique(cs)) == 1L) next  # equal groups: searched above
      b <- length(rs) * length(cs)
      if (p / b < max(2, k - 1) || p / b > k + 1) next
      pb <- do_option_allocate(p, rs, cs)
      if (any(pb < 1L) || any(pb >= outer(rs, cs))) next
      cap <- augmented_row_column_capacity(n_controls = k, rows_per_rowgroup = rs,
                                           cols_per_colgroup = cs, checks_per_block = pb)
      if (!cap$every_row_has_check || !cap$every_col_has_check ||
          cap$min_rowgroup_slack < 0 || cap$min_colgroup_slack < 0 || cap$error_df < min_df) next
      args <- list(rows_per_rowgroup = rs, cols_per_colgroup = cs, checks_per_block = pb,
                   method = "optimal")
      choices[[length(choices) + 1L]] <- do_option_row("augmented_row_column", n, k, r, c, locations,
        1L, r * c, p, blocks = b, error_df = cap$error_df, method = "optimal",
        generator = "augmented_row_column_design", arguments = args,
        note = paste0("Unequal groups: rows ", paste(rs, collapse = ","), "; columns ",
                      paste(cs, collapse = ","), ". Check plots per block follow block size."))
    }
  }
  choices
}

do_option_prep <- function(n, fields, locations, sizes, proportions, min_df, limit) {
  choices <- list()
  for (d in unique(as.integer(round(n * proportions)))) {
    if (d < 1L || as.double(d) * locations > n) next
    plots <- n + d
    for (size in sizes[sizes <= n]) {
      blocks <- ceiling(plots / size)
      df <- d - blocks + 1L
      if (blocks < 2L || df < min_df) next
      sizes_here <- prep_balanced_counts(plots, blocks)
      if (sum(sizes_here * (sizes_here - 1) / 2) * locations > as.double(n) * (n - 1) / 2) next
      fits <- fields[fields$rows == ceiling(plots / fields$cols), , drop = FALSE]
      fits <- utils::head(fits, 3L * limit)
      for (i in seq_len(nrow(fits))) {
        r <- fits$rows[i]; c <- fits$cols[i]
        args <- list(locations = paste0("Location_", seq_len(locations)), method = "search",
          block_size = size, n_duplicate_per_location = d, field_cols = c)
        choices[[length(choices) + 1L]] <- do_option_row("prep", n, 0L, r, c, locations,
          plots = plots, repeated = d, blocks = blocks, error_df = df, method = "search",
          generator = "augmented_prep_design", arguments = args,
          note = "Two plots for selected entries; one for the others. Check connectivity and concurrence after search.")
      }
    }
  }
  choices
}

do_option_resolvable <- function(n, checks, fields, locations, reps, min_df, limit) {
  choices <- list()
  for (check in checks) {
    v <- as.double(n) + check
    if (v < 2) next
    for (rep in reps) {
      for (layout in c("vertical", "horizontal")) {
        k <- if (layout == "vertical") fields$rows / rep else fields$rows
        s <- if (layout == "horizontal") fields$cols / rep else fields$cols
        df <- rep * v - rep * (k + s - 1) - (v - 1)
        keep <- k == floor(k) & s == floor(s) & k >= 2 & s >= 2 & k * s >= v & df >= min_df
        for (i in utils::head(which(keep), 3L * limit)) {
          args <- list(n_reps = rep, rows_per_rep = k[i], cols_per_rep = s[i], replicate_layout = layout)
          choices[[length(choices) + 1L]] <- do_option_row("resolvable", n, check,
            fields$rows[i], fields$cols[i], locations, rep, rep * v, check * rep,
            repeated = v, blocks = rep * k[i], error_df = df[i], method = "search",
            generator = "resolvable_row_column_design", arguments = args,
            note = if (check) "Include checks in treatments; mark is_control in the generated field books." else
              "Every entry once per replicate; unused positions are fillers.")
        }
      }
    }
  }
  choices
}
