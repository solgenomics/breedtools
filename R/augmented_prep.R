#' Generate an augmented p-rep design
#'
#' Creates a multi-location augmented partially replicated (p-rep) design
#' (Williams, Piepho and Whitaker, 2011). Instead of replicated checks, a
#' subset of the test entries is replicated at each location: every entry is
#' planted once at every location, and the entries chosen for p-rep at a
#' location get extra plots there. All plots of a location are arranged in
#' incomplete blocks, chosen so that pairs of entries share a block at most
#' once across the whole trial ((0,1)-designs).
#'
#' Three construction methods are available through `method`:
#' \describe{
#'   \item{`"search"` (default)}{A flexible randomized search. Any number of
#'     entries, locations, block sizes, p-rep proportions (`prep_prop`,
#'     `n_duplicate_per_location`), replication levels
#'     (`replications_per_repeated_entry`) and block counts can be used, with
#'     an optional resolvable layout (`resolvable`).}
#'   \item{`"alpha"`}{The partial alpha-array construction of Section 2.1 of
#'     the paper. An alpha-array for `l * m` rows and `2l` replicates is
#'     searched, the `m x 2` rectangles on its diagonal give the two-replicate
#'     design for the p-rep entries of each location, and half of the
#'     remaining elements of each column are deleted (Section 2.2) so that the
#'     design is a (0,1)-design. Every entry is duplicated at exactly one
#'     location, each location is resolvable with two replicates of `s`
#'     blocks, and each block holds `m` p-rep entries; the final block size is
#'     `(l + 1) m / 2` (blocks of the two replicates differ by one plot when
#'     this is not an integer).}
#'   \item{`"alpha2"`}{The augmented alpha(2)-array construction of
#'     Section 2.4, for even `m` and an odd number of blocks per replicate
#'     `s`. Every block has `(l + 1) m / 2` plots.}
#' }
#' For the algebraic methods, `m` is set with `n_prep_per_block` (or chosen as
#' the largest value whose blocks fit in `block_size`), the number of blocks
#' and of duplicated entries follow from the construction, and
#' `replications_per_repeated_entry` must be 2. When the number of entries is
#' not a multiple of the array size, the surplus array positions are dropped
#' as dummy entries, so some blocks are one plot smaller.
#' `augmented_prep_parameters()` lists the feasible parameter sets.
#'
#' The incomplete block structure is independent from the field geometry. By
#' default, `field_cols` lays the ordered plots into rows and columns. For a
#' rectangular block layout, use `rows_per_block` and `blocks_per_row` to place
#' multiple blocks side by side within each location.
#'
#' @param treatments Character vector of entry names.
#' @param locations Character vector of location names.
#' @param method Construction method: `"search"`, `"alpha"` or `"alpha2"`
#'   (see Details).
#' @param block_size Maximum final incomplete block size. Can be a single
#'   integer or a named vector indexed by location (the algebraic methods use
#'   its smallest value). Actual block sizes may differ by one plot when
#'   totals do not divide evenly.
#' @param prep_prop Optional proportion of entries to duplicate per location
#'   (`method = "search"` only). If `NULL`, entries are split across locations
#'   so each entry is duplicated at approximately one location.
#' @param n_duplicate_per_location Optional number of entries duplicated at
#'   each location (`method = "search"` only). Can be a single integer, a
#'   vector matching `locations`, or a named vector indexed by location.
#' @param replications_per_repeated_entry Total number of plots for entries
#'   selected for p-rep replication at a location. The default `2` gives one
#'   base plot plus one extra plot. Can be a single integer or a named vector
#'   indexed by location. Must be 2 for the algebraic methods.
#' @param n_blocks_per_location Optional number of incomplete blocks at each
#'   location (`method = "search"` only). If `NULL`, the number of blocks is
#'   the minimum needed so no block is larger than `block_size`. Can be a
#'   single integer or a named vector indexed by location.
#' @param resolvable Logical (`method = "search"`). If `TRUE`, each location
#'   is laid out as a resolvable design: the blocks are grouped into
#'   `replications_per_repeated_entry` replicates, each replicate holds every
#'   p-rep entry of that location exactly once plus an even share of the
#'   unreplicated entries, and all replicates have the same number of blocks.
#'   If `n_blocks_per_location` is given it must be a multiple of the number
#'   of replicates. Designs from the algebraic methods are always resolvable.
#' @param n_prep_per_block For the algebraic methods, the number `m` of p-rep
#'   entries in every block. If `NULL`, the largest feasible `m` whose blocks
#'   fit in `block_size` is used.
#' @param n_candidates Maximum number of candidate designs (`"search"`) or of
#'   random starting arrays (`"alpha"`, `"alpha2"`) to generate.
#'   `summary$n_candidates_evaluated` reports how many were used.
#' @param n_efficiency_candidates Number of optimal candidates (no entry
#'   repeated within a block, connected locations, and no pair of entries
#'   sharing more than one block) to compare on the average efficiency factor;
#'   the one with the highest efficiency is returned. Use `1` to return the
#'   first optimal candidate. Efficiency is compared for up to 2000 entries.
#' @param seed Optional random seed.
#' @param field_cols Optional number of field columns for row/column layout.
#'   Can be a single integer or a named vector indexed by location.
#' @param rows_per_block Optional number of physical field rows each incomplete
#'   block should occupy. Can be a single integer or a named vector for selected
#'   locations. Must be used with `blocks_per_row` for those locations.
#' @param blocks_per_row Optional number of incomplete blocks to place side by
#'   side in each field row band. Can be a single integer or a named vector for
#'   selected locations. Must be used with `rows_per_block` for those locations.
#' @param serpentine Logical. If `TRUE`, row/column coordinates use serpentine
#'   ordering within each location.
#' @param output_file Optional tab-delimited output file for the complete
#'   design.
#'
#' @param trial_name Optional Breedbase trial name. When given, a
#'   `plot_name` column with the plot names used by [save_breedbase()] is
#'   added to the design (see [add_plot_names()]), so the names can be
#'   checked with [breedbase_check()] and changed before saving. For a
#'   multi-location design the location is added, one trial per location.
#' @param plot_name_pattern Pattern for the plot names, as in
#'   [add_plot_names()]. Used only with `trial_name`; `{location}` is the
#'   design location.
#' @details
#' With `method = "search"`, the p-rep plots of every location are spread
#' evenly over its blocks, and the p-rep entries are placed so that they link
#' all blocks of the location (the location is connected), so entry effects
#' can also be estimated from each location on its own. A warning is given
#' when there are too few p-rep plots to connect the blocks. After the random
#' greedy search, the best candidates are improved by swapping plots between
#' blocks to remove pairs of entries that share more than one block.
#'
#' For all methods, `summary` reports the maximum pairwise concurrence, the
#' number of locations whose blocks are not connected, and the average
#' efficiency factor `efficiency_factor` of the block design (Section 2.3 of
#' the paper; higher is better), which can be used to compare designs.
#' Entry labels, block order within replicates and plot order within blocks
#' are randomized.
#'
#' @references Williams, E., Piepho, H.-P. and Whitaker, D. (2011). Augmented
#'   p-rep designs. Biometrical Journal 53, 19-27. \doi{10.1002/bimj.201000102}
#'
#' @return A list with design inputs, summary diagnostics, complete plot-level
#' design, Breedbase-compatible design columns, block summary, replication
#' summary, repeated-entry groups, and pairwise concurrence matrix. In the
#' complete design, `repeated_at` is the location where the entry receives its
#' extra p-rep plot(s), or `"not_repeated"` when the entry is not selected for
#' extra replication. `replicate` is the replicate (superblock) of the plot
#' for resolvable designs, and 1 otherwise; the Breedbase `rep` column is
#' `replicate` for resolvable designs and the plot's occurrence number
#' otherwise.
#'
#' @seealso [augmented_prep_parameters()], [analyze_augmented_prep_design()],
#'   [plot_augmented_prep_design()]
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(12)),
#'   locations = c("Loc1", "Loc2", "Loc3"),
#'   block_size = c(Loc1 = 5, Loc2 = 5, Loc3 = 6),
#'   n_blocks_per_location = c(Loc1 = 4, Loc2 = 4, Loc3 = 4),
#'   n_duplicate_per_location = c(Loc1 = 4, Loc2 = 2, Loc3 = 4),
#'   replications_per_repeated_entry = c(Loc1 = 2, Loc2 = 3, Loc3 = 2),
#'   n_candidates = 5,
#'   seed = 123,
#'   field_cols = 8
#' )
#'
#' design$summary
#' head(design$design)
#'
#' # Example 1 of Williams, Piepho and Whitaker (2011): 20 entries, two
#' # locations, blocks of 3 built from a partial alpha-array.
#' alpha <- augmented_prep_design(
#'   treatments = paste0("G", 1:20),
#'   locations = c("L1", "L2"),
#'   method = "alpha",
#'   block_size = 3,
#'   seed = 1
#' )
#' alpha$summary[, c("method", "m", "s", "max_pairwise_concurrence", "efficiency_factor")]
#'
#' @export
augmented_prep_design <- function(
    treatments,
    locations,
    block_size = 8,
    prep_prop = NULL,
    n_duplicate_per_location = NULL,
    replications_per_repeated_entry = 2,
    n_blocks_per_location = NULL,
    resolvable = FALSE,
    n_candidates = 1000,
    seed = NULL,
    field_cols = NULL,
    rows_per_block = NULL,
    blocks_per_row = NULL,
    serpentine = TRUE,
    output_file = NULL,
    method = c("search", "alpha", "alpha2"),
    n_prep_per_block = NULL,
    n_efficiency_candidates = 5,
    trial_name = NULL,
    plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}"
) {
  method <- match.arg(method)
  n_efficiency_candidates <- arc_positive_integer(n_efficiency_candidates, "n_efficiency_candidates")

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

  if (method != "search") {
    not_used <- c(
      prep_prop = !is.null(prep_prop),
      n_duplicate_per_location = !is.null(n_duplicate_per_location),
      n_blocks_per_location = !is.null(n_blocks_per_location)
    )
    if (any(not_used)) {
      stop(
        "`", paste(names(not_used)[not_used], collapse = "`, `"), "` cannot be used with method \"",
        method, "\": the construction fixes the p-rep entries and blocks. Use ",
        "`n_prep_per_block` and `block_size`, or method \"search\".",
        call. = FALSE
      )
    }
    if (any(as.integer(replications_per_repeated_entry) != 2L)) {
      stop("Method \"", method, "\" needs `replications_per_repeated_entry = 2`.", call. = FALSE)
    }
  } else if (!is.null(n_prep_per_block)) {
    stop("`n_prep_per_block` is only used by methods \"alpha\" and \"alpha2\".", call. = FALSE)
  }

  inputs <- prep_validate_inputs(
    treatments = treatments,
    locations = locations,
    block_size = block_size,
    prep_prop = prep_prop,
    n_duplicate_per_location = n_duplicate_per_location,
    replications_per_repeated_entry = replications_per_repeated_entry,
    n_blocks_per_location = n_blocks_per_location,
    resolvable = resolvable,
    n_candidates = n_candidates,
    field_cols = field_cols,
    rows_per_block = rows_per_block,
    blocks_per_row = blocks_per_row
  )
  inputs$method <- method

  if (method == "search") {
    setup <- prep_candidate_setup(inputs)
    searched <- prep_search_design(setup, inputs$n_candidates, n_efficiency_candidates, length(inputs$treatments))
    best_candidate <- searched$candidate
    n_evaluated <- searched$n_evaluated
    parameters <- NULL
  } else {
    built <- prep_algebraic_design(method, inputs, n_prep_per_block, n_efficiency_candidates)
    best_candidate <- built$candidate
    setup <- built$setup
    n_evaluated <- built$n_evaluated
    parameters <- built$parameters
    inputs$resolvable <- TRUE
    inputs$n_prep_per_block <- parameters$m
    inputs$n_duplicate_per_location <- vapply(best_candidate$groups, length, integer(1))
  }

  best <- prep_build_candidate_design(
    candidate = best_candidate,
    inputs = inputs,
    setup = setup,
    serpentine = serpentine
  )
  best$summary$n_candidates_evaluated <- n_evaluated
  best$summary$method <- method
  best$summary$m <- if (is.null(parameters)) NA_integer_ else parameters$m
  best$summary$s <- if (is.null(parameters)) NA_integer_ else parameters$s
  best$summary$n_dummies <- if (is.null(parameters)) 0L else parameters$d

  block_summary <- best$design |>
    dplyr::group_by(.data$location, .data$block, .data$block_uid) |>
    dplyr::summarise(n_plots = dplyr::n(), .groups = "drop")

  replication_summary <- best$design |>
    dplyr::group_by(.data$all_entries) |>
    dplyr::summarise(
      total_plots = dplyr::n(),
      locations = paste(unique(.data$location), collapse = ", "),
      repeated_at = paste(unique(.data$repeated_at[.data$repeated_at != "not_repeated"]), collapse = ", "),
      .groups = "drop"
    )
  replication_summary$repeated_at[replication_summary$repeated_at == ""] <- "not_repeated"

  breedbase_design <- best$design |>
    dplyr::select(dplyr::all_of(c("plots", "block", "all_entries", "rep", "is_control")))

  result <- list(
    input = inputs,
    summary = best$summary,
    design = best$design,
    breedbase_design = breedbase_design,
    block_summary = block_summary,
    replication_summary = replication_summary,
    repeated_groups = best$repeated_groups,
    concurrence_matrix = best$concurrence_matrix
  )
  class(result) <- c("augmented_prep_design", class(result))
  if (!is.null(trial_name)) result <- add_plot_names(result, trial_name, plot_name_pattern = plot_name_pattern)

  if (!is.null(output_file)) {
    utils::write.table(result$design, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
  }

  result
}

#' Calculate pairwise concurrence for a p-rep design
#'
#' Counts how many times each pair of entries appears together in the same
#' incomplete block across all locations.
#'
#' @param design A p-rep design data frame from `augmented_prep_design()`.
#' @param treatments Optional treatment names defining matrix order.
#'
#' @return A symmetric integer matrix of pairwise block concurrences.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(12)),
#'   locations = c("Loc1", "Loc2"),
#'   block_size = 4,
#'   n_candidates = 5,
#'   seed = 1
#' )
#' conc <- augmented_prep_concurrence(design$design)
#' table(conc[upper.tri(conc)])
#'
#' @export
augmented_prep_concurrence <- function(design, treatments = NULL) {
  if (is.null(treatments)) {
    treatments <- sort(unique(as.character(design$all_entries)))
  }

  treatments <- as.character(treatments)
  entry_idx <- match(as.character(design$all_entries), treatments)
  if (anyNA(entry_idx)) {
    stop("`design` contains entries not listed in `treatments`.", call. = FALSE)
  }
  blocks <- interaction(design$location, design$block, drop = TRUE)

  # Block x entry presence matrix; crossprod counts shared blocks per pair.
  presence <- matrix(0, nrow = nlevels(blocks), ncol = length(treatments))
  presence[cbind(as.integer(blocks), entry_idx)] <- 1
  mat <- crossprod(presence)
  diag(mat) <- 0
  storage.mode(mat) <- "integer"
  dimnames(mat) <- list(treatments, treatments)
  mat
}


#' Analyze an augmented p-rep design
#'
#' Fits the fixed-effects incomplete-block model:
#' `y = mean + location + block(location) + genotype + error`.
#' For resolvable designs (a `replicate` column with more than one replicate
#' per location) the model is
#' `y = mean + location + replicate(location) + block(replicate, location) +
#' genotype + error`. With `genotype_by_location = "fixed"` (or `TRUE`) a
#' fixed genotype-by-location interaction is added after genotypes; the error
#' is then estimated only from the p-rep entries, as in a per-location
#' analysis with pooled error (Williams, Piepho and Whitaker, 2011,
#' Section 3). With `genotype_by_location = "random"` the interaction is a
#' random effect, as recommended by the paper when locations are a sample of
#' environments, and the model is fitted by REML with 'lme4'.
#' Sums of squares are sequential, in the order above. Degrees of freedom in
#' the formatted ANOVA are those estimable from the observed (non-missing)
#' plots, so F-tests stay valid when plots are missing. For a complete,
#' connected design they equal the design formulas reported in `anova_df`,
#' which use total plots (`n`), locations (`l`), incomplete blocks (`b`), and
#' genotypes (`v`). Terms with a single level in the data (for example the
#' location term in a single-location trial) are left out of the model and
#' shown with 0 df.
#'
#' @param data A complete augmented p-rep design data frame, usually
#'   `design$design` from `augmented_prep_design()`, or an object returned by
#'   `augmented_prep_design()`, with a response column added.
#' @param response Name of the numeric response column.
#' @param genotype_col Name of the genotype/treatment column.
#' @param location_col Name of the location column.
#' @param block_col Name of the incomplete block column. Block numbers may
#'   repeat across locations; blocks are treated as nested within locations.
#' @param entry_type_col Optional entry type column. Used only for genotype
#'   summaries.
#' @param extra_repeat_col Optional extra-repeat indicator column. Used only
#'   for genotype summaries.
#' @param repeated_at_col Optional column identifying where each entry is
#'   repeated. Used only for genotype summaries.
#' @param control_col Optional control flag column. Used only for genotype
#'   summaries.
#' @param replicate_col Optional replicate (superblock) column. When present
#'   and a location has more than one replicate, replicates nested within
#'   locations are fitted, and blocks are nested within replicates. Set to
#'   `NULL` to ignore it.
#' @param genotype_by_location Genotype-by-location interaction in a
#'   multi-location analysis: `"none"` (default, or `FALSE`), `"fixed"` (or
#'   `TRUE`), or `"random"` (needs the 'lme4' package). With `"fixed"`,
#'   adjusted means are averaged over locations. With `"random"`, adjusted
#'   means are the genotype BLUEs from the mixed model, `variance_components`
#'   gives the REML estimates, and the ANOVA table is the fixed-interaction
#'   table, whose genotype-by-location and error mean squares show the size of
#'   the interaction.
#'
#' @return An object of class `augmented_prep_analysis`: a list with
#' \describe{
#'   \item{input, formula, model}{The column mapping, model formula, and the
#'     `lm` fit.}
#'   \item{anova, anova_df, anova_formatted}{Sequential ANOVA (estimable df),
#'     design df formulas, and a formatted ANOVA table.}
#'   \item{plot_results}{The input data with fitted values and residuals.}
#'   \item{genotype_summary}{One row per genotype with raw, fitted, and
#'     adjusted means and the standard error of the adjusted mean.}
#'   \item{adjusted_means}{Genotypes ranked by adjusted mean (best first).}
#'   \item{adjusted_means_vcov}{Covariance matrix of the estimable adjusted
#'     means, used by [analyze_augmented_prep_two_stage()] for weighting.}
#'   \item{diagnostics}{Error df, residual standard deviation, R-squared,
#'     `mean_sed` (mean standard error of a difference between two adjusted
#'     means), `lsd_5pct` (the matching least significant difference), and
#'     flags for rank deficiency or a saturated fit.}
#'   \item{variance_components, mixed_model}{For
#'     `genotype_by_location = "random"`: REML variance components and the
#'     `lmerMod` fit; otherwise `NULL`.}
#'   \item{message}{A note when F-tests or adjusted means are not estimable.}
#' }
#'
#' @section Adjusted means:
#' Genotype effects are fixed, so the adjusted means are best linear unbiased
#' estimates (BLUEs): the genotype effect plus the average location and block
#' effect over all observed plots, i.e. the expected value of the genotype in
#' an average plot of the trial. They remove the advantage or penalty of the
#' blocks and locations where each entry happened to be planted, so use
#' `adjusted_mean` (not `raw_mean`) to rank and select entries. Differences
#' between adjusted means are the same as those from `emmeans` on the fitted
#' model. With `genotype_by_location = "fixed"` the genotype-by-location
#' effects are averaged over locations, weighting each location by its share
#' of the observed plots. With `"random"`, standard errors come from the
#' REML covariance of the fixed effects, which includes the
#' genotype-by-location variance; `lsd_5pct` uses the error df of the
#' fixed-interaction model as an approximation.
#'
#' @section Per-location analysis:
#' Williams, Piepho and Whitaker (2011) suggest analysing each location with
#' an incomplete-block model. Pass the rows of one location, e.g.
#' `analyze_augmented_prep_design(data[data$location == "A", ], "yield")`.
#' `augmented_prep_design()` keeps the blocks of every location connected,
#' so this is estimable whenever the location has enough p-rep plots.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(8)),
#'   locations = c("Loc1", "Loc2"),
#'   block_size = 4,
#'   n_duplicate_per_location = c(Loc1 = 2, Loc2 = 2),
#'   n_candidates = 3,
#'   seed = 123
#' )
#'
#' design_data <- design$design
#' design_data$yield <- rnorm(nrow(design_data), mean = 100, sd = 10)
#' analysis <- analyze_augmented_prep_design(design_data, response = "yield")
#' analysis
#' head(analysis$adjusted_means)
#'
#' @export
analyze_augmented_prep_design <- function(
    data,
    response,
    genotype_col = "all_entries",
    location_col = "location",
    block_col = "block",
    entry_type_col = "entry_type",
    extra_repeat_col = "is_extra_repeat",
    repeated_at_col = "repeated_at",
    control_col = "is_control",
    replicate_col = "replicate",
    genotype_by_location = FALSE
) {
  gxl_mode <- prep_gxl_mode(genotype_by_location)
  genotype_by_location <- gxl_mode != "none"
  if (is.list(data) && !is.data.frame(data) && !is.null(data$design)) {
    data <- data$design
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or an augmented p-rep design object.", call. = FALSE)
  }

  required_cols <- c(response, genotype_col, location_col, block_col)
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
      location = as.factor(.data[[location_col]]),
      block = as.factor(.data[[block_col]]),
      replicate = if (!is.null(replicate_col) && replicate_col %in% names(data)) {
        as.factor(.data[[replicate_col]])
      } else {
        factor(1L)
      },
      replicate_nested = interaction(.data[[location_col]], .data$replicate, drop = TRUE),
      block_nested = interaction(.data[[location_col]], .data$replicate, .data[[block_col]], drop = TRUE),
      entry_type = if (entry_type_col %in% names(data)) {
        as.character(.data[[entry_type_col]])
      } else {
        NA_character_
      },
      is_extra_repeat = if (extra_repeat_col %in% names(data)) {
        as.integer(.data[[extra_repeat_col]])
      } else {
        NA_integer_
      },
      repeated_at = if (repeated_at_col %in% names(data)) {
        as.character(.data[[repeated_at_col]])
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
      !is.na(.data$location),
      !is.na(.data$block),
      !is.na(.data$block_nested)
    )
  has_replicates <- nlevels(droplevels(design_structure_data$replicate_nested)) >
    dplyr::n_distinct(design_structure_data$location)

  analysis_data <- tibble::as_tibble(data) |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    dplyr::transmute(
      .row_id = .data$.row_id,
      response = suppressWarnings(as.numeric(.data[[response]])),
      genotype = as.factor(.data[[genotype_col]]),
      location = as.factor(.data[[location_col]]),
      block = as.factor(.data[[block_col]]),
      replicate = if (!is.null(replicate_col) && replicate_col %in% names(data)) {
        as.factor(.data[[replicate_col]])
      } else {
        factor(1L)
      },
      replicate_nested = interaction(.data[[location_col]], .data$replicate, drop = TRUE),
      block_nested = interaction(.data[[location_col]], .data$replicate, .data[[block_col]], drop = TRUE),
      entry_type = if (entry_type_col %in% names(data)) {
        as.character(.data[[entry_type_col]])
      } else {
        NA_character_
      },
      is_extra_repeat = if (extra_repeat_col %in% names(data)) {
        as.integer(.data[[extra_repeat_col]])
      } else {
        NA_integer_
      },
      repeated_at = if (repeated_at_col %in% names(data)) {
        as.character(.data[[repeated_at_col]])
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
      !is.na(.data$location),
      !is.na(.data$block),
      !is.na(.data$block_nested)
    )

  if (nrow(analysis_data) == 0) {
    stop("No complete observations are available for analysis.", call. = FALSE)
  }

  if (dplyr::n_distinct(analysis_data$genotype) < 2) {
    stop("At least two genotypes are required for analysis.", call. = FALSE)
  }

  analysis_data <- droplevels(analysis_data)

  # Only fit terms that vary in the observed data: a single-location trial has
  # no location term, and a location with one block has no block term.
  multi_location <- nlevels(analysis_data$location) > 1
  model_terms <- c(
    if (multi_location) "location",
    if (nlevels(analysis_data$replicate_nested) > nlevels(analysis_data$location)) "replicate_nested",
    if (nlevels(analysis_data$block_nested) > nlevels(analysis_data$replicate_nested)) "block_nested",
    "genotype",
    if (genotype_by_location && multi_location) "location:genotype"
  )
  fixed_terms <- setdiff(model_terms, "location:genotype")
  model_formula <- stats::as.formula(paste("response ~", paste(model_terms, collapse = " + ")))
  design_df <- prep_augmented_prep_design_df(
    design_structure_data,
    genotype_by_location = genotype_by_location
  )
  fit <- arc_lm_drop_aliased_terms(model_formula, analysis_data)
  anova_result <- arc_lm_anova_table(fit)
  anova_table <- anova_result$table
  anova_formatted <- prep_format_augmented_prep_anova(
    anova_table = anova_table,
    response = analysis_data$response,
    has_replicates = has_replicates,
    genotype_by_location = genotype_by_location
  )

  plot_results <- tibble::as_tibble(data) |>
    dplyr::mutate(
      .analysis_response = suppressWarnings(as.numeric(.data[[response]]))
    )

  plot_results$fitted <- NA_real_
  plot_results$residual <- NA_real_
  plot_results$fitted[analysis_data$.row_id] <- stats::fitted(fit)
  plot_results$residual[analysis_data$.row_id] <- stats::residuals(fit)

  genotype_summary <- plot_results |>
    dplyr::group_by(.data[[genotype_col]]) |>
    dplyr::summarise(
      genotype = dplyr::first(as.character(.data[[genotype_col]])),
      entry_type = if (entry_type_col %in% names(plot_results)) {
        paste(sort(unique(as.character(.data[[entry_type_col]]))), collapse = ", ")
      } else {
        NA_character_
      },
      repeated_at = if (repeated_at_col %in% names(plot_results)) {
        paste(sort(unique(as.character(.data[[repeated_at_col]]))), collapse = ", ")
      } else {
        NA_character_
      },
      is_control = if (control_col %in% names(plot_results)) {
        dplyr::first(as.integer(.data[[control_col]]))
      } else {
        NA_integer_
      },
      n_observations = sum(!is.na(.data$.analysis_response)),
      n_locations = if (location_col %in% names(plot_results)) {
        dplyr::n_distinct(.data[[location_col]][!is.na(.data$.analysis_response)])
      } else {
        NA_integer_
      },
      n_extra_repeats = if (extra_repeat_col %in% names(plot_results)) {
        sum(as.integer(.data[[extra_repeat_col]]) == 1L & !is.na(.data$.analysis_response), na.rm = TRUE)
      } else {
        NA_integer_
      },
      raw_mean = mean(.data$.analysis_response, na.rm = TRUE),
      fitted_mean = mean(.data$fitted, na.rm = TRUE),
      residual_mean = mean(.data$residual, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::select(-dplyr::all_of(genotype_col))

  genotype_summary$raw_mean[is.nan(genotype_summary$raw_mean)] <- NA_real_
  genotype_summary$fitted_mean[is.nan(genotype_summary$fitted_mean)] <- NA_real_
  genotype_summary$residual_mean[is.nan(genotype_summary$residual_mean)] <- NA_real_

  variance_components <- NULL
  mixed_model <- NULL
  if (gxl_mode == "random" && multi_location) {
    mixed <- prep_fit_random_gxl(analysis_data, fixed_terms, error_df = stats::df.residual(fit))
    adjusted_fit <- mixed$adjusted
    variance_components <- mixed$variance_components
    mixed_model <- mixed$model
  } else {
    adjusted_fit <- prep_adjusted_means(fit, analysis_data)
  }
  adjusted_parts <- arc_attach_adjusted_means(genotype_summary, adjusted_fit$adjusted_means)
  genotype_summary <- adjusted_parts$genotype_summary
  precision <- adjusted_fit$precision

  anova_warning <- anova_result$warning
  if (is.null(anova_warning)) {
    anova_warning <- NA_character_
  }

  aliased_coefficients_dropped <- attr(fit, "aliased_coefficients_dropped")
  if (is.null(aliased_coefficients_dropped)) {
    aliased_coefficients_dropped <- 0L
  }

  diagnostics <- tibble::tibble(
    n_observations = stats::nobs(fit),
    n_locations = dplyr::n_distinct(analysis_data$location),
    n_blocks = dplyr::n_distinct(analysis_data$block_nested),
    n_genotypes = dplyr::n_distinct(analysis_data$genotype),
    model_rank = fit$rank,
    residual_df = stats::df.residual(fit),
    anova_error_df = stats::df.residual(fit),
    design_error_df = design_df$display_df[design_df$term == "Residuals"],
    aliased_coefficients_dropped = aliased_coefficients_dropped,
    sigma = summary(fit)$sigma,
    r_squared = summary(fit)$r.squared,
    adjusted_r_squared = summary(fit)$adj.r.squared,
    rank_deficient = fit$rank < length(stats::coef(fit)),
    saturated_or_near_perfect_fit = isTRUE(stats::df.residual(fit) == 0) ||
      isTRUE(!is.na(summary(fit)$sigma) && summary(fit)$sigma < sqrt(.Machine$double.eps)),
    mean_sed = precision$mean_sed,
    lsd_5pct = precision$lsd_5pct,
    anova_warning = anova_warning
  )

  analysis_message <- arc_augmented_analysis_message(diagnostics)
  analysis_message <- arc_append_estimability_message(analysis_message, adjusted_fit$adjusted_means)

  result <- list(
    input = list(
      response = response,
      genotype_col = genotype_col,
      location_col = location_col,
      block_col = block_col,
      entry_type_col = entry_type_col,
      extra_repeat_col = extra_repeat_col,
      repeated_at_col = repeated_at_col,
      control_col = control_col,
      replicate_col = replicate_col,
      genotype_by_location = gxl_mode
    ),
    formula = model_formula,
    model = fit,
    anova = anova_table,
    anova_df = design_df,
    anova_formatted = anova_formatted,
    plot_results = plot_results |> dplyr::select(-dplyr::all_of(".analysis_response")),
    genotype_summary = genotype_summary,
    adjusted_means = adjusted_parts$adjusted_means,
    adjusted_means_vcov = adjusted_fit$vcov,
    diagnostics = diagnostics,
    variance_components = variance_components,
    mixed_model = mixed_model,
    message = analysis_message
  )

  class(result) <- c("augmented_prep_analysis", class(result))
  result
}

#' Print augmented p-rep analysis
#'
#' @param x An object returned by `analyze_augmented_prep_design()`.
#' @param ... Additional arguments passed to `print()`.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.augmented_prep_analysis <- function(x, ...) {
  cat("Augmented p-rep design analysis\n\n")
  print(x$anova_formatted, ...)
  if (!is.null(x$variance_components)) {
    cat("\nREML variance components (random genotype x location):\n")
    print(as.data.frame(x$variance_components), row.names = FALSE, ...)
  }
  arc_print_adjusted_means(x)
  if (!is.null(x$message) && !is.na(x$message) && nzchar(x$message)) {
    cat("\n", x$message, "\n", sep = "")
  }
  invisible(x)
}


#' Plot an augmented p-rep design grid
#'
#' Draws one field map per location for an augmented p-rep design.
#'
#' A p-rep design has no check plots: at each location some entries are
#' planted more than once (repeated entries), and these replace the checks
#' for estimating the error and the block effects. By default the map
#' colours the repeated entries (orange, names in bold) and the entries with
#' a single plot at that location (light grey). Thick black lines outline the
#' incomplete blocks; in resolvable designs a thick blue line separates the
#' replicates, which are labelled at the right. A caption under each map
#' explains this.
#'
#' @param design An object returned by `augmented_prep_design()` or its
#'   `design` data frame.
#' @param locations Optional character vector of locations to plot. Defaults
#'   to all locations in design order.
#' @param label Cell label to show. Use `"all_entries"`, `"entry_type"`,
#'   `"repeated_at"`, `"block"`, `"replicate"`, or `"none"`.
#' @param fill Cell coloring variable. Use `"entry_type"`,
#'   `"is_extra_repeat"`, `"repeated_at"`, `"block"`, or `"replicate"`.
#' @param show_plot_labels Logical. If `TRUE`, labels are printed inside
#'   cells unless `label = "none"`.
#' @param show_block_labels Logical. If `TRUE`, block numbers are printed.
#'   Defaults to `FALSE` so block labels do not overlap accession labels.
#' @param block_label_position Position for block labels when
#'   `show_block_labels = TRUE`. Use `"background"` for a faint watermark
#'   behind accession labels or `"center"` for a regular center label.
#' @param label_cex Optional text size for plot-cell labels.
#' @param block_label_cex Optional text size for block labels.
#' @param block_label_col Color for block labels.
#' @param block_lwd Line width for block boundaries.
#' @param main Optional plot title. If multiple locations are plotted, a
#'   single value is combined with the location name; a vector matching the
#'   number of locations is used panel by panel.
#' @param legend Logical. If `TRUE`, draw a legend for cell colors.
#' @param caption Logical. If `TRUE`, print a short explanation of the
#'   colours and lines under each map.
#' @param ... Additional arguments passed to `graphics::plot.window()`.
#'
#' @return Invisibly returns the plotted design data frame with internal
#'   plotting coordinates.
#'
#' @examples
#' design <- augmented_prep_design(
#'   treatments = paste0("G", seq_len(12)),
#'   locations = c("Loc1", "Loc2", "Loc3"),
#'   block_size = 5,
#'   n_candidates = 5,
#'   seed = 123,
#'   field_cols = 5
#' )
#'
#' plot_augmented_prep_design(design)
#'
#' @export
plot_augmented_prep_design <- function(
    design,
    locations = NULL,
    label = c("all_entries", "entry_type", "repeated_at", "block", "replicate", "none"),
    fill = c("entry_type", "is_extra_repeat", "repeated_at", "block", "replicate"),
    show_plot_labels = TRUE,
    show_block_labels = FALSE,
    block_label_position = c("background", "center"),
    label_cex = NULL,
    block_label_cex = NULL,
    block_label_col = NULL,
    block_lwd = 2,
    main = NULL,
    legend = TRUE,
    caption = TRUE,
    ...
) {
  label <- match.arg(label)
  fill <- match.arg(fill)
  block_label_position <- match.arg(block_label_position)

  input_locations <- NULL
  if (is.list(design) && !is.data.frame(design) && !is.null(design$design)) {
    input_locations <- design$input$locations
    design <- design$design
  }

  if (!is.data.frame(design)) {
    stop("`design` must be a design data frame or an augmented p-rep design object.", call. = FALSE)
  }

  required_cols <- c("location", "row", "col", "block")
  missing_cols <- setdiff(required_cols, names(design))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (label != "none" && !label %in% names(design)) {
    stop("Column `", label, "` was not found in `design`.", call. = FALSE)
  }
  if (!fill %in% names(design)) {
    stop("Column `", fill, "` was not found in `design`.", call. = FALSE)
  }

  design <- as.data.frame(design, stringsAsFactors = FALSE)
  design$location <- as.character(design$location)
  design$row <- suppressWarnings(as.numeric(design$row))
  design$col <- suppressWarnings(as.numeric(design$col))

  if (any(is.na(design$row)) || any(is.na(design$col))) {
    stop(
      "Plotting augmented p-rep designs requires row and column coordinates. ",
      "Call `augmented_prep_design(..., field_cols = ...)` before plotting.",
      call. = FALSE
    )
  }

  if (is.null(locations)) {
    # Locations in the order they were given to augmented_prep_design().
    locations <- unique(design$location)
    if (!is.null(input_locations)) locations <- intersect(input_locations, locations)
  } else {
    locations <- as.character(locations)
  }

  missing_locations <- setdiff(locations, unique(design$location))
  if (length(missing_locations) > 0) {
    stop("Location(s) not found in `design`: ", paste(missing_locations, collapse = ", "), call. = FALSE)
  }

  design <- design[design$location %in% locations, , drop = FALSE]
  design$.panel_index <- match(design$location, locations)
  design <- design[order(design$.panel_index, design$row, design$col), , drop = FALSE]
  design$.fill_value <- prep_plot_fill_values(design, fill)

  fill_colors <- prep_plot_fill_colors(design$.fill_value, fill)
  plot_args <- list(...)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n_locations <- length(locations)
  panel_cols <- ceiling(sqrt(n_locations))
  panel_rows <- ceiling(n_locations / panel_cols)
  graphics::par(
    mfrow = c(panel_rows, panel_cols),
    mar = c(if (isTRUE(caption)) 5 else 2, 4.5, 5.5, if (isTRUE(legend)) 10 else 4),
    xpd = NA
  )

  plotted <- vector("list", n_locations)
  for (i in seq_along(locations)) {
    location <- locations[i]
    location_design <- design[design$location == location, , drop = FALSE]
    panel_main <- if (is.null(main)) {
      paste("Location", location)
    } else if (length(main) == n_locations) {
      main[i]
    } else {
      paste(main[1], location, sep = " - ")
    }

    plotted[[i]] <- prep_plot_prep_location(
      design = location_design,
      fill_colors = fill_colors,
      label = label,
      show_plot_labels = show_plot_labels,
      show_block_labels = show_block_labels,
      block_label_position = block_label_position,
      label_cex = label_cex,
      block_label_cex = block_label_cex,
      block_label_col = block_label_col,
      block_lwd = block_lwd,
      main = panel_main,
      legend = legend,
      fill_title = prep_plot_fill_title(fill),
      caption = caption,
      plot_args = plot_args
    )
  }

  rownames(design) <- NULL
  invisible(tibble::as_tibble(do.call(rbind, plotted)))
}

prep_plot_fill_values <- function(design, fill) {
  value <- as.character(design[[fill]])
  if (fill == "entry_type") {
    value <- ifelse(
      value == "p_rep_entry_at_this_location",
      "Repeated here",
      "Single plot here"
    )
  }
  if (fill == "is_extra_repeat") {
    value <- ifelse(as.integer(design[[fill]]) == 1L, "Extra plot of a repeated entry", "First plot of the entry")
  }
  if (fill == "replicate") value <- paste("Replicate", value)
  if (fill == "block") value <- paste("Block", value)
  if (fill == "repeated_at") value <- ifelse(value == "not_repeated", "Not repeated", paste("Repeated at", value))
  value
}

prep_plot_fill_title <- function(fill) {
  switch(fill,
    entry_type = "Plot type",
    is_extra_repeat = "Plot",
    repeated_at = "Where the entry is repeated",
    block = "Block",
    replicate = "Replicate"
  )
}

# Strong, distinct colours: repeated entries stand out against light grey.
prep_plot_fill_colors <- function(values, fill) {
  if (fill == "entry_type") {
    out <- c("Repeated here" = "#E69F00", "Single plot here" = "#EEEEEE")
    return(out[intersect(names(out), unique(values))])
  }
  if (fill == "is_extra_repeat") {
    out <- c("Extra plot of a repeated entry" = "#E69F00", "First plot of the entry" = "#EEEEEE")
    return(out[intersect(names(out), unique(values))])
  }
  lv <- unique(values)
  lv <- lv[order(suppressWarnings(as.numeric(sub("^\\D*", "", lv))), lv)]
  cols <- if (length(lv) <= 8) {
    c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#CC79A7", "#D55E00", "#0072B2", "#999999")[seq_along(lv)]
  } else {
    grDevices::hcl.colors(length(lv), palette = "Set 2")
  }
  stats::setNames(cols, lv)
}

prep_plot_prep_location <- function(design, fill_colors, label, show_plot_labels,
                                    show_block_labels, block_label_position,
                                    label_cex, block_label_cex, block_label_col,
                                    block_lwd, main, legend, fill_title,
                                    caption = TRUE, plot_args) {
  rows <- sort(unique(design$row))
  cols <- sort(unique(design$col))
  design$.plot_row <- match(design$row, rows)
  design$.plot_col <- match(design$col, cols)

  if (is.null(label_cex)) {
    label_cex <- max(0.25, min(0.8, 7 / max(length(rows), length(cols))))
  }

  graphics::plot.new()
  do.call(
    graphics::plot.window,
    c(
      list(
        xlim = c(0.5, length(cols) + 0.5),
        ylim = c(length(rows) + 0.5, 0.5),
        xaxs = "i",
        yaxs = "i",
        asp = 1
      ),
      plot_args
    )
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

  prep_plot_prep_block_boundaries(design, block_lwd = block_lwd)
  has_reps <- "replicate" %in% names(design) && length(unique(design$replicate)) > 1
  if (has_reps) prep_plot_replicate_boundaries(design, lwd = block_lwd * 2.5)

  if (isTRUE(show_block_labels)) {
    block_split <- split(design, design$block)
    for (block_design in block_split) {
      block_cex <- block_label_cex
      block_col <- block_label_col
      if (is.null(block_cex)) {
        block_cex <- if (block_label_position == "background") {
          max(1.1, min(3.5, 0.6 * sqrt(nrow(block_design))))
        } else {
          0.85
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
        x = mean(block_design$.plot_col),
        y = mean(block_design$.plot_row),
        labels = paste0("B", unique(block_design$block)[1]),
        cex = block_cex,
        font = 2,
        col = block_col
      )
    }
  }

  repeated <- if ("entry_type" %in% names(design)) {
    design$entry_type == "p_rep_entry_at_this_location"
  } else {
    rep(FALSE, nrow(design))
  }
  if (isTRUE(show_plot_labels) && label != "none") {
    graphics::text(
      x = design$.plot_col,
      y = design$.plot_row,
      labels = as.character(design[[label]]),
      cex = label_cex,
      font = ifelse(repeated, 2, 1),
      col = "black"
    )
  }

  # Frame around the field only.
  graphics::rect(0.5, length(rows) + 0.5, length(cols) + 0.5, 0.5, border = "black", lwd = 1)
  graphics::axis(3, at = seq_along(cols), labels = cols, las = 1, tick = FALSE, line = -0.6, cex.axis = 0.8)
  graphics::axis(2, at = seq_along(rows), labels = rows, las = 1, tick = FALSE, line = -0.6, cex.axis = 0.8)
  graphics::mtext("Column", side = 3, line = 1.6, cex = 0.8)
  graphics::mtext("Row", side = 2, line = 2.2, cex = 0.8)
  graphics::title(main = main, line = 3.2)

  # Replicate labels at the right of the field.
  if (has_reps) {
    for (rp in sort(unique(design$replicate))) {
      cells <- design[design$replicate == rp, ]
      graphics::text(length(cols) + 0.65, mean(range(cells$.plot_row)), paste("Rep", rp),
                     col = "#1F4E79", font = 2, cex = 0.85, adj = 0)
    }
  }

  if (isTRUE(legend)) {
    graphics::legend(
      x = length(cols) + 0.6,
      y = 0.5,
      legend = names(fill_colors),
      fill = fill_colors,
      border = "grey40",
      title = fill_title,
      title.adj = 0,
      bty = "n",
      xjust = 0,
      yjust = 0,
      cex = 0.8
    )
  }

  if (isTRUE(caption)) {
    colour_note <- if (identical(fill_title, "Plot type") || identical(fill_title, "Plot")) "orange, names in bold" else "names in bold"
    lines <- c(
      paste0("No check plots: the repeated entries (", colour_note,
             ") are planted more than once here and replace the checks."),
      paste0("Thick black lines: incomplete blocks", if (has_reps) "; thick blue line: replicates." else ".")
    )
    graphics::mtext(lines, side = 1, line = c(0.8, 1.8), cex = 0.7, adj = 0, at = 0.5)
  }

  design
}

# Thick line between plots of different replicates.
prep_plot_replicate_boundaries <- function(design, lwd) {
  rep_lookup <- stats::setNames(as.character(design$replicate),
                                paste(design$.plot_row, design$.plot_col, sep = ":"))
  for (i in seq_len(nrow(design))) {
    r <- design$.plot_row[i]
    cl <- design$.plot_col[i]
    here <- as.character(design$replicate[i])
    right <- rep_lookup[paste(r, cl + 1L, sep = ":")]
    below <- rep_lookup[paste(r + 1L, cl, sep = ":")]
    if (!is.na(right) && right != here) {
      graphics::segments(cl + 0.5, r - 0.5, cl + 0.5, r + 0.5, lwd = lwd, col = "#1F4E79")
    }
    if (!is.na(below) && below != here) {
      graphics::segments(cl - 0.5, r + 0.5, cl + 0.5, r + 0.5, lwd = lwd, col = "#1F4E79")
    }
  }
}

prep_plot_prep_block_boundaries <- function(design, block_lwd) {
  block_lookup <- stats::setNames(
    as.character(design$block),
    paste(design$.plot_row, design$.plot_col, sep = ":")
  )
  same_block <- function(row, col, block) {
    key <- paste(row, col, sep = ":")
    !is.na(block_lookup[key]) && identical(unname(block_lookup[key]), as.character(block))
  }

  for (i in seq_len(nrow(design))) {
    row <- design$.plot_row[i]
    col <- design$.plot_col[i]
    block <- design$block[i]
    left <- col - 0.5
    right <- col + 0.5
    top <- row - 0.5
    bottom <- row + 0.5

    if (!same_block(row, col - 1L, block)) {
      graphics::segments(left, top, left, bottom, lwd = block_lwd, col = "black")
    }
    if (!same_block(row, col + 1L, block)) {
      graphics::segments(right, top, right, bottom, lwd = block_lwd, col = "black")
    }
    if (!same_block(row - 1L, col, block)) {
      graphics::segments(left, top, right, top, lwd = block_lwd, col = "black")
    }
    if (!same_block(row + 1L, col, block)) {
      graphics::segments(left, bottom, right, bottom, lwd = block_lwd, col = "black")
    }
  }
}


# Adjusted genotype means for the p-rep model: the genotype effect plus the
# average location/block effect over all observed plots (the prediction for
# an "average plot" of the trial). Estimability is checked against the full
# model matrix, so genotypes whose effect is aliased get NA.
prep_adjusted_L <- function(fit, analysis_data) {
  x_full <- stats::model.matrix(fit$terms, fit$model)
  y <- stats::model.response(fit$model)
  assign <- attr(x_full, "assign")
  genotype_term <- match("genotype", attr(fit$terms, "term.labels"))
  genotype_cols <- which(assign == genotype_term)
  genotype_levels <- levels(analysis_data$genotype)

  L <- matrix(colMeans(x_full), nrow = length(genotype_levels), ncol = ncol(x_full), byrow = TRUE)
  L[, genotype_cols] <- 0
  col_level <- match(colnames(x_full)[genotype_cols], paste0("genotype", genotype_levels))
  L[cbind(col_level, genotype_cols)] <- 1

  # Genotype-by-location columns: average genotype g over locations, each
  # location weighted by its share of the observed plots (the same weights
  # the location main effect gets above).
  gxl_term <- match("location:genotype", attr(fit$terms, "term.labels"))
  if (!is.na(gxl_term)) {
    gxl_cols <- which(assign == gxl_term)
    location_levels <- levels(analysis_data$location)
    weight <- as.numeric(table(fit$model$location)[location_levels]) / nrow(x_full)
    names_grid <- outer(
      paste0("genotype", genotype_levels),
      paste0("location", location_levels),
      function(g, l) paste(l, g, sep = ":")
    )
    L[, gxl_cols] <- 0
    hit <- match(colnames(x_full)[gxl_cols], names_grid)
    g_idx <- (hit - 1L) %% length(genotype_levels) + 1L
    l_idx <- (hit - 1L) %/% length(genotype_levels) + 1L
    ok <- !is.na(hit)
    L[cbind(g_idx[ok], gxl_cols[ok])] <- weight[l_idx[ok]]
  }

  list(L = L, x_full = x_full, y = y, genotype_levels = genotype_levels)
}

prep_adjusted_means <- function(fit, analysis_data) {
  parts <- prep_adjusted_L(fit, analysis_data)
  L <- parts$L
  x_full <- parts$x_full
  y <- parts$y
  genotype_levels <- parts$genotype_levels
  qr_full <- qr(x_full)
  est <- arc_linear_estimates(qr_full, qr.coef(qr_full, y), L)
  error_df <- stats::df.residual(fit)
  sigma2 <- if (error_df > 0) sum(stats::residuals(fit)^2) / error_df else NA_real_
  std_error <- sqrt(sigma2 * colSums(est$W^2))
  estimate <- est$estimate
  estimate[!est$estimable] <- NA_real_
  std_error[!est$estimable] <- NA_real_

  list(
    adjusted_means = tibble::tibble(
      genotype = genotype_levels,
      adjusted_mean = estimate,
      std_error = std_error,
      estimable = est$estimable
    ),
    vcov = prep_means_vcov(est$W, est$estimable, sigma2, genotype_levels),
    precision = arc_mean_precision(
      est$W[, est$estimable, drop = FALSE],
      rep(0, sum(est$estimable)),
      sigma2,
      error_df
    )
  )
}

# Covariance matrix of the estimable adjusted means: sigma2 * W'W.
prep_means_vcov <- function(W, estimable, sigma2, genotype_levels) {
  ok <- which(estimable)
  out <- sigma2 * crossprod(W[, ok, drop = FALSE])
  dimnames(out) <- list(genotype_levels[ok], genotype_levels[ok])
  out
}

prep_gxl_mode <- function(x) {
  if (is.logical(x) && length(x) == 1 && !is.na(x)) {
    return(if (x) "fixed" else "none")
  }
  if (is.character(x) && length(x) == 1 && x %in% c("none", "fixed", "random")) {
    return(x)
  }
  stop("`genotype_by_location` must be \"none\", \"fixed\", \"random\", TRUE or FALSE.", call. = FALSE)
}

# Mixed model with random genotype-by-location effects, fitted by REML with
# lme4 on the full-rank fixed-effects matrix of the model without the
# interaction. Genotype adjusted means are estimable functions of the fixed
# effects, so they are computed from the kept columns.
prep_fit_random_gxl <- function(analysis_data, fixed_terms, error_df) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("`genotype_by_location = \"random\"` needs the 'lme4' package.", call. = FALSE)
  }
  formula0 <- stats::as.formula(paste("response ~", paste(fixed_terms, collapse = " + ")))
  fit0 <- arc_lm_drop_aliased_terms(formula0, analysis_data)
  parts <- prep_adjusted_L(fit0, analysis_data)
  qr_full <- qr(parts$x_full)
  estimable <- arc_linear_estimates(qr_full, qr.coef(qr_full, parts$y), parts$L)$estimable

  x_kept <- fit0$x
  kept <- match(colnames(x_kept), colnames(parts$x_full))
  mixed_data <- data.frame(
    response = fit0$y,
    gxl = interaction(analysis_data$genotype, analysis_data$location, drop = TRUE)
  )
  mixed_data$X <- x_kept
  model <- lme4::lmer(
    response ~ 0 + X + (1 | gxl),
    data = mixed_data,
    REML = TRUE,
    control = lme4::lmerControl(check.rankX = "ignore")
  )
  beta <- lme4::fixef(model)
  vcov_beta <- as.matrix(stats::vcov(model))
  L_kept <- parts$L[, kept, drop = FALSE]
  estimate <- drop(L_kept %*% beta)
  cov_means <- L_kept %*% vcov_beta %*% t(L_kept)
  std_error <- sqrt(pmax(diag(cov_means), 0))
  estimate[!estimable] <- NA_real_
  std_error[!estimable] <- NA_real_
  means_vcov <- cov_means[estimable, estimable, drop = FALSE]
  dimnames(means_vcov) <- list(parts$genotype_levels[estimable], parts$genotype_levels[estimable])

  ok <- which(estimable)
  m <- length(ok)
  precision <- list(mean_sed = NA_real_, lsd_5pct = NA_real_)
  if (m >= 2) {
    sub <- cov_means[ok, ok, drop = FALSE]
    mean_var_diff <- 2 * (m * sum(diag(sub)) - sum(sub)) / (m * (m - 1))
    mean_sed <- sqrt(max(mean_var_diff, 0))
    precision <- list(
      mean_sed = mean_sed,
      lsd_5pct = if (error_df > 0) stats::qt(0.975, error_df) * mean_sed else NA_real_
    )
  }

  vc <- as.data.frame(lme4::VarCorr(model))
  list(
    model = model,
    adjusted = list(
      adjusted_means = tibble::tibble(
        genotype = parts$genotype_levels,
        adjusted_mean = estimate,
        std_error = std_error,
        estimable = estimable
      ),
      vcov = means_vcov,
      precision = precision
    ),
    variance_components = tibble::tibble(
      component = c("genotype x location", "residual"),
      variance = vc$vcov[match(c("gxl", "Residual"), vc$grp)]
    )
  )
}

prep_format_augmented_prep_anova <- function(anova_table, response,
                                             has_replicates = FALSE,
                                             genotype_by_location = FALSE) {
  source_map <- data.frame(
    term = c(
      "location",
      "replicate_nested",
      "block_nested",
      "genotype",
      "location:genotype",
      "Residuals"
    ),
    source = c(
      "Locations",
      "Replicates, nested within locations",
      if (has_replicates) "Blocks, nested within replicates" else "Blocks, nested within locations",
      "Genotypes",
      "Genotype x location",
      "Error"
    ),
    order = seq_len(6),
    stringsAsFactors = FALSE
  )
  source_map <- source_map[
    (has_replicates | source_map$term != "replicate_nested") &
      (genotype_by_location | source_map$term != "location:genotype"),
    ,
    drop = FALSE
  ]

  # Degrees of freedom are those estimable from the observed plots. A term that
  # was not fitted (e.g. locations in a single-location trial) has 0 df.
  out <- source_map |>
    dplyr::left_join(anova_table, by = "term") |>
    dplyr::arrange(.data$order) |>
    dplyr::mutate(
      df = dplyr::coalesce(.data$df, 0),
      sum_sq = dplyr::coalesce(.data$sum_sq, 0),
      mean_sq = dplyr::if_else(.data$df > 0, .data$sum_sq / .data$df, NA_real_)
    )

  error_df <- out$df[out$term == "Residuals"]
  error_mean_sq <- out$mean_sq[out$term == "Residuals"]

  out <- out |>
    dplyr::mutate(
      f_value = dplyr::if_else(
        .data$term != "Residuals" &
          !is.na(error_mean_sq) &
          error_mean_sq > 0 &
          .data$df > 0,
        .data$mean_sq / error_mean_sq,
        NA_real_
      ),
      p_value = dplyr::if_else(
        !is.na(.data$f_value) &
          !is.na(error_df) &
          error_df > 0,
        stats::pf(.data$f_value, df1 = .data$df, df2 = error_df, lower.tail = FALSE),
        NA_real_
      )
    ) |>
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

prep_augmented_prep_design_df <- function(data, genotype_by_location = FALSE) {
  n <- nrow(data)
  l <- dplyr::n_distinct(data$location)
  r <- dplyr::n_distinct(data$replicate_nested)
  b <- dplyr::n_distinct(data$block_nested)
  v <- dplyr::n_distinct(data$genotype)
  gxl <- genotype_by_location && l > 1
  # Genotype-by-location cells actually present in the design.
  cells <- dplyr::n_distinct(paste(data$genotype, data$location, sep = "\r"))

  terms <- c("location", if (r > l) "replicate_nested", "block_nested", "genotype", if (gxl) "location:genotype")
  source_df <- c(
    l - 1,
    if (r > l) r - l,
    b - r,
    v - 1,
    if (gxl) cells - v - l + 1
  )
  corrected_total_df <- n - 1
  error_df <- corrected_total_df - sum(source_df)

  formula_parts <- c(
    paste0("(", l, " - 1)"),
    if (r > l) paste0("(", r, " - ", l, ")"),
    paste0("(", b, " - ", r, ")"),
    paste0("(", v, " - 1)"),
    if (gxl) paste0("(", cells, " - ", v, " - ", l, " + 1)")
  )

  tibble::tibble(
    term = c(terms, "Residuals"),
    display_df = c(source_df, error_df),
    n = n,
    l = l,
    b = b,
    v = v,
    corrected_total_df = corrected_total_df,
    error_df_formula = paste0("(", n, " - 1) - [", paste(formula_parts, collapse = " + "), "]")
  )
}


prep_validate_inputs <- function(treatments, locations, block_size, prep_prop,
                                 n_duplicate_per_location,
                                 replications_per_repeated_entry,
                                 n_blocks_per_location, n_candidates,
                                 field_cols, rows_per_block,
                                 blocks_per_row, resolvable = FALSE) {
  treatments <- as.character(treatments)
  locations <- as.character(locations)
  treatments <- treatments[!is.na(treatments) & treatments != ""]
  locations <- locations[!is.na(locations) & locations != ""]
  n_candidates <- arc_positive_integer(n_candidates, "n_candidates")
  if (!is.logical(resolvable) || length(resolvable) != 1 || is.na(resolvable)) {
    stop("`resolvable` must be TRUE or FALSE.", call. = FALSE)
  }

  if (length(treatments) < 2) {
    stop("At least two treatments are required.", call. = FALSE)
  }
  if (length(locations) < 1) {
    stop("At least one location is required.", call. = FALSE)
  }
  if (anyDuplicated(treatments)) {
    stop("`treatments` must be unique.", call. = FALSE)
  }
  if (anyDuplicated(locations)) {
    stop("`locations` must be unique.", call. = FALSE)
  }
  if (!is.null(prep_prop) && (length(prep_prop) != 1 || is.na(prep_prop) ||
    prep_prop <= 0 || prep_prop > 1)) {
    stop("`prep_prop` must be a single number greater than 0 and no larger than 1.", call. = FALSE)
  }

  n_duplicate_per_location <- prep_normalize_duplicate_counts(
    treatments, locations, prep_prop, n_duplicate_per_location
  )
  block_size <- prep_normalize_location_counts(block_size, locations, "block_size")
  replications_per_repeated_entry <- prep_normalize_location_counts(
    replications_per_repeated_entry,
    locations,
    "replications_per_repeated_entry"
  )
  if (any(replications_per_repeated_entry < 2L)) {
    stop("`replications_per_repeated_entry` must be at least 2.", call. = FALSE)
  }
  n_blocks_per_location <- prep_normalize_optional_location_counts(
    n_blocks_per_location,
    locations,
    "n_blocks_per_location"
  )
  field_cols <- prep_normalize_field_cols(field_cols, locations)
  rows_per_block <- prep_normalize_optional_partial_location_counts(
    rows_per_block,
    locations,
    "rows_per_block"
  )
  blocks_per_row <- prep_normalize_optional_partial_location_counts(
    blocks_per_row,
    locations,
    "blocks_per_row"
  )
  prep_validate_block_layout_counts(rows_per_block, blocks_per_row, locations)

  list(
    treatments = treatments,
    locations = locations,
    block_size = block_size,
    prep_prop = prep_prop,
    n_duplicate_per_location = n_duplicate_per_location,
    replications_per_repeated_entry = replications_per_repeated_entry,
    n_blocks_per_location = n_blocks_per_location,
    resolvable = resolvable,
    n_candidates = n_candidates,
    field_cols = field_cols,
    rows_per_block = rows_per_block,
    blocks_per_row = blocks_per_row
  )
}

prep_normalize_duplicate_counts <- function(treatments, locations, prep_prop, n_duplicate_per_location) {
  n_treatments <- length(treatments)
  n_locations <- length(locations)

  if (is.null(n_duplicate_per_location)) {
    if (is.null(prep_prop)) {
      out <- prep_balanced_counts(n_treatments, n_locations)
    } else {
      out <- rep(round(n_treatments * prep_prop), n_locations)
    }
  } else if (length(n_duplicate_per_location) == 1) {
    out <- rep(as.integer(n_duplicate_per_location), n_locations)
  } else if (!is.null(names(n_duplicate_per_location))) {
    out <- as.integer(n_duplicate_per_location[locations])
  } else {
    out <- as.integer(n_duplicate_per_location)
  }

  if (length(out) != n_locations || any(is.na(out)) || any(out < 0)) {
    stop("`n_duplicate_per_location` must be non-negative and match locations.", call. = FALSE)
  }
  if (sum(out) > n_treatments) {
    stop("Total duplicated entries cannot exceed number of treatments.", call. = FALSE)
  }

  out
}


prep_normalize_location_counts <- function(x, locations, label) {
  if (length(x) == 1) {
    out <- rep(arc_positive_integer(x, label), length(locations))
    names(out) <- locations
    return(out)
  }

  if (!is.null(names(x))) {
    out <- as.integer(x[locations])
  } else {
    out <- as.integer(x)
  }

  if (length(out) != length(locations) || any(is.na(out)) || any(out < 1)) {
    stop("`", label, "` must contain positive integers for every location.", call. = FALSE)
  }
  names(out) <- locations
  out
}

prep_normalize_optional_location_counts <- function(x, locations, label) {
  if (is.null(x)) {
    return(NULL)
  }
  prep_normalize_location_counts(x, locations, label)
}

prep_normalize_optional_partial_location_counts <- function(x, locations, label) {
  out <- stats::setNames(rep(NA_integer_, length(locations)), locations)
  if (is.null(x)) {
    return(out)
  }

  if (length(x) == 1 && is.null(names(x))) {
    out[] <- arc_positive_integer(x, label)
    return(out)
  }

  if (!is.null(names(x))) {
    unknown <- setdiff(names(x), locations)
    if (length(unknown) > 0) {
      stop("`", label, "` contains unknown location(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    values <- as.integer(x)
    if (any(is.na(values)) || any(values < 1)) {
      stop("`", label, "` must contain positive integers.", call. = FALSE)
    }
    out[names(x)] <- values
    return(out)
  }

  values <- as.integer(x)
  if (length(values) != length(locations) || any(is.na(values)) || any(values < 1)) {
    stop("`", label, "` must be a positive integer, a named vector by location, or match `locations`.", call. = FALSE)
  }
  out[] <- values
  out
}

prep_validate_block_layout_counts <- function(rows_per_block, blocks_per_row, locations) {
  rows_set <- !is.na(rows_per_block)
  blocks_set <- !is.na(blocks_per_row)
  incomplete <- locations[xor(rows_set, blocks_set)]
  if (length(incomplete) > 0) {
    stop(
      "`rows_per_block` and `blocks_per_row` must both be supplied for location(s): ",
      paste(incomplete, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


prep_normalize_field_cols <- function(field_cols, locations) {
  if (is.null(field_cols)) {
    return(NULL)
  }
  if (length(field_cols) == 1) {
    out <- rep(arc_positive_integer(field_cols, "field_cols"), length(locations))
    names(out) <- locations
    return(out)
  }
  if (is.null(names(field_cols))) {
    stop("If `field_cols` has length greater than 1, it must be named by location.", call. = FALSE)
  }
  out <- as.integer(field_cols[locations])
  if (any(is.na(out)) || any(out < 1)) {
    stop("`field_cols` must contain positive integers for every location.", call. = FALSE)
  }
  names(out) <- locations
  out
}

# Fixed, candidate-independent quantities for the candidate search.
#
# For every location this fixes the block sizes, the replicate (superblock)
# of each block when the design is resolvable, and how many p-rep plots each
# block receives. Spreading the p-rep plots evenly over blocks follows
# Williams, Piepho and Whitaker (2011), where every block holds the same
# number (m) of p-rep entries and is filled up with unreplicated entries.
prep_candidate_setup <- function(inputs) {
  n_treatments <- length(inputs$treatments)
  locations <- inputs$locations
  reps <- as.integer(inputs$replications_per_repeated_entry[locations])
  n_duplicate <- as.integer(inputs$n_duplicate_per_location)
  resolvable <- isTRUE(inputs$resolvable)

  layout <- lapply(seq_along(locations), function(l) {
    loc <- locations[l]
    n_blocks_here <- if (is.null(inputs$n_blocks_per_location)) NULL else inputs$n_blocks_per_location[[loc]]
    prep_location_layout(
      loc = loc,
      n_treatments = n_treatments,
      n_dup = n_duplicate[l],
      reps = reps[l],
      max_block_size = inputs$block_size[[loc]],
      n_blocks = n_blocks_here,
      resolvable = resolvable && n_duplicate[l] > 0
    )
  })
  block_sizes <- lapply(layout, `[[`, "sizes")

  # Each entry has one plot per location plus (reps - 1) extra plots at the
  # location where it is repeated, so the multiset of replication counts (and
  # its variance term in the score) is the same for every candidate.
  replication_counts <- c(
    rep(length(locations), n_treatments - sum(n_duplicate)),
    rep(length(locations) + reps - 1L, n_duplicate)
  )
  replication_variance <- sum((replication_counts - mean(replication_counts))^2)
  any_pairs <- any(unlist(block_sizes) > 1L)

  # Blocks at a location are linked only through entries with more than one
  # plot there. An entry with r plots joins at most r blocks, so at least
  # (blocks - 1) - n_dup * (reps - 1) components cannot be avoided.
  min_extra_components <- vapply(seq_along(locations), function(l) {
    max(0L, (length(block_sizes[[l]]) - 1L) - n_duplicate[l] * (reps[l] - 1L))
  }, integer(1))
  unavoidable <- locations[min_extra_components > 0]
  if (length(unavoidable) > 0) {
    warning(
      "The blocks at location(s) ", paste(unavoidable, collapse = ", "),
      " cannot all be connected by the p-rep entries (too few repeated plots ",
      "for the number of blocks). Entry effects are then not estimable from ",
      "that location alone; a combined analysis across locations still is. ",
      "Increase `n_duplicate_per_location` / `prep_prop` or use fewer, ",
      "larger blocks.",
      call. = FALSE
    )
  }

  list(
    n_treatments = n_treatments,
    n_locations = length(locations),
    reps = reps,
    n_duplicate = n_duplicate,
    block_sizes = block_sizes,
    block_rep = lapply(layout, `[[`, "block_rep"),
    prep_quota = lapply(layout, `[[`, "prep_quota"),
    resolvable = vapply(layout, `[[`, logical(1), "resolvable"),
    min_extra_components = min_extra_components,
    replication_variance = replication_variance,
    score_lower_bound = replication_variance +
      (if (any_pairs) 1000 else 0) +
      sum(min_extra_components) * 1e10
  )
}

# Block sizes, replicate membership and p-rep quota for one location.
prep_location_layout <- function(loc, n_treatments, n_dup, reps, max_block_size,
                                 n_blocks, resolvable) {
  n_plots <- n_treatments + n_dup * (reps - 1L)

  if (!resolvable) {
    sizes <- prep_balanced_block_sizes(n_plots, max_block_size, n_blocks)
    if (n_dup > 0 && reps > length(sizes)) {
      stop(
        "Location `", loc, "` has ", reps,
        " total reps for selected entries but only ", length(sizes),
        " blocks. Increase `n_blocks_per_location` or reduce `replications_per_repeated_entry`.",
        call. = FALSE
      )
    }
    return(list(
      sizes = sizes,
      block_rep = rep(1L, length(sizes)),
      prep_quota = prep_balanced_counts(n_dup * reps, length(sizes)),
      resolvable = FALSE
    ))
  }

  # Resolvable: `reps` replicates, each holding every p-rep entry once plus an
  # even share of the unreplicated entries, with the same number of blocks.
  rep_totals <- n_dup + prep_balanced_counts(n_treatments - n_dup, reps)
  if (is.null(n_blocks)) {
    blocks_per_rep <- max(ceiling(rep_totals / max_block_size))
  } else {
    n_blocks <- arc_positive_integer(n_blocks, "n_blocks_per_location")
    if (n_blocks %% reps != 0) {
      stop(
        "With `resolvable = TRUE`, `n_blocks_per_location` for location `", loc,
        "` must be a multiple of `replications_per_repeated_entry` (", reps, ").",
        call. = FALSE
      )
    }
    blocks_per_rep <- n_blocks %/% reps
  }
  if (blocks_per_rep > min(rep_totals)) {
    stop(
      "Location `", loc, "` has too many blocks per replicate for its plots.",
      call. = FALSE
    )
  }
  sizes <- unlist(lapply(rep_totals, prep_balanced_counts, n_groups = blocks_per_rep))
  if (max(sizes) > max_block_size) {
    stop(
      "`n_blocks_per_location` is too small for `block_size` at location `", loc,
      "`; at least one block would contain ", max(sizes), " plots.",
      call. = FALSE
    )
  }

  list(
    sizes = sizes,
    block_rep = rep(seq_len(reps), each = blocks_per_rep),
    prep_quota = rep(prep_balanced_counts(n_dup, blocks_per_rep), reps),
    resolvable = TRUE
  )
}

# Number of connected components among the blocks of one location, where two
# blocks are linked when they share an entry. `entry` and `block` describe the
# plots of the location.
prep_block_components <- function(entry, block, n_blocks) {
  if (n_blocks <= 1L) {
    return(1L)
  }
  labels <- seq_len(n_blocks)
  multi <- entry[duplicated(entry)]
  for (e in unique(multi)) {
    labs <- unique(labels[block[entry == e]])
    if (length(labs) > 1L) {
      labels[labels %in% labs] <- min(labs)
    }
  }
  length(unique(labels))
}

# One candidate: random repeated-entry groups and a greedy block allocation.
# At each location the p-rep entries are placed first, one entry at a time, so
# that their plots link different blocks (keeping the location connected), and
# then the unreplicated entries fill the remaining space. Placements avoid
# repeating an entry in a block and avoid pairs of entries meeting in more than
# one block (across locations). Works on integer entry ids and keeps the score
# statistics incrementally.
prep_make_candidate <- function(setup) {
  n <- setup$n_treatments
  n_dup <- setup$n_duplicate
  duplicated_ids <- sample.int(n, sum(n_dup))
  groups <- split(duplicated_ids, factor(rep.int(seq_along(n_dup), n_dup), levels = seq_along(n_dup)))

  concurrence <- matrix(0L, n, n)
  same_entry_same_block <- 0
  allocations <- vector("list", setup$n_locations)

  for (l in seq_len(setup$n_locations)) {
    g <- groups[[l]]
    if (length(g) > 1L) g <- g[sample.int(length(g))]
    reps <- setup$reps[l]
    singles <- setdiff(seq_len(n), g)
    if (length(singles) > 1L) singles <- singles[sample.int(length(singles))]

    occ_entry <- c(rep(g, each = reps), singles)
    occ_number <- c(rep.int(seq_len(reps), length(g)), rep.int(1L, length(singles)))
    occ_is_prep <- c(rep(TRUE, length(g) * reps), rep(FALSE, length(singles)))

    block_sizes <- setup$block_sizes[[l]]
    block_rep <- setup$block_rep[[l]]
    resolvable <- setup$resolvable[l]
    n_blocks <- length(block_sizes)
    prep_left <- setup$prep_quota[[l]]
    single_left <- block_sizes - prep_left
    presence <- matrix(0, nrow = n_blocks, ncol = n)
    members <- vector("list", n_blocks)
    component <- seq_len(n_blocks)
    occ_block <- integer(length(occ_entry))

    for (i in seq_along(occ_entry)) {
      e <- occ_entry[i]
      is_prep <- occ_is_prep[i]
      open <- if (is_prep) which(prep_left > 0L) else which(single_left > 0L)
      if (is_prep && resolvable) {
        in_rep <- open[block_rep[open] == occ_number[i]]
        if (length(in_rep) > 0) open <- in_rep
      }
      if (length(open) == 0) {
        open <- which((prep_left + single_left) > 0L)
      }
      free <- presence[open, e] == 0
      if (any(free)) {
        open <- open[free]
      }

      filled <- block_sizes[open] - prep_left[open] - single_left[open]
      penalty <- filled / block_sizes[open] + stats::runif(length(open))
      ce <- concurrence[, e]
      met <- which(ce > 0L)
      if (length(met) > 0) {
        # sum over block members m of: 1e4 * [conc >= 1] + 100 * conc^2
        penalty <- penalty + as.vector(presence[open, met, drop = FALSE] %*% (1e4 + 100 * ce[met]^2))
      }
      if (is_prep && occ_number[i] > 1L) {
        # Prefer a block not yet linked to the blocks already holding `e`.
        own <- unique(component[which(presence[, e] > 0)])
        penalty <- penalty + 1e5 * (component[open] %in% own)
      }
      b <- open[which.min(penalty)]

      if (presence[b, e] == 0) {
        mem <- members[[b]]
        if (length(mem) > 0) {
          new_conc <- concurrence[mem, e] + 1L
          concurrence[mem, e] <- new_conc
          concurrence[e, mem] <- new_conc
        }
        if (is_prep && occ_number[i] > 1L) {
          own <- unique(component[which(presence[, e] > 0)])
          merged <- c(own, component[b])
          component[component %in% merged] <- min(merged)
        }
        presence[b, e] <- 1
        members[[b]] <- c(mem, e)
      } else {
        same_entry_same_block <- same_entry_same_block + 1
      }

      if (is_prep && prep_left[b] > 0L) {
        prep_left[b] <- prep_left[b] - 1L
      } else if (!is_prep && single_left[b] > 0L) {
        single_left[b] <- single_left[b] - 1L
      } else if (prep_left[b] > 0L) {
        prep_left[b] <- prep_left[b] - 1L
      } else {
        single_left[b] <- single_left[b] - 1L
      }
      occ_block[i] <- b
    }

    allocations[[l]] <- list(
      entry = occ_entry,
      occurrence = occ_number,
      block = occ_block,
      is_prep = occ_is_prep,
      extra_components = length(unique(component)) - 1L
    )
  }

  candidate <- list(
    groups = groups,
    allocations = allocations,
    concurrence = concurrence,
    same_entry_same_block = same_entry_same_block
  )
  candidate$score <- prep_candidate_score(candidate, setup)
  candidate
}

prep_candidate_score <- function(candidate, setup) {
  upper <- candidate$concurrence[upper.tri(candidate$concurrence)]
  extra_components <- sum(vapply(candidate$allocations, `[[`, integer(1), "extra_components"))
  # Excess concurrence is penalized quadratically, so one pair meeting in
  # every location costs more than several pairs meeting twice.
  candidate$same_entry_same_block * 1e12 +
    extra_components * 1e10 +
    sum(pmax(upper - 1L, 0L)^2) * 1e7 +
    (if (length(upper) > 0) max(upper) else 0) * 1000 +
    setup$replication_variance
}

# Local improvement of a candidate by swapping two plots between blocks of the
# same location. Unreplicated plots are swapped with unreplicated plots and
# p-rep plots with p-rep plots (within the same replicate when resolvable), so
# block sizes, the p-rep quota of every block and resolvability are kept.
# Swaps that lower the number of pairs concurring more than once are kept, and
# neutral swaps are taken at random to move across plateaus; p-rep swaps that
# would disconnect the blocks of a location are rejected.
prep_improve_candidate <- function(candidate, setup, max_tries = NULL) {
  if (candidate$same_entry_same_block > 0) {
    candidate <- prep_repair_duplicates(candidate, setup)
    if (candidate$same_entry_same_block > 0) return(candidate)
  }
  conc <- candidate$concurrence
  allocs <- candidate$allocations
  n <- nrow(conc)
  n_plots <- sum(vapply(allocs, function(a) length(a$entry), integer(1)))
  if (is.null(max_tries)) {
    max_tries <- min(15000L, 8L * n_plots)
  }
  pair_cost <- function(x) 1e7 * pmax(x - 1L, 0L)^2

  # Per location: block x entry presence and the entries of every block.
  presence <- lapply(seq_along(allocs), function(l) {
    m <- matrix(FALSE, length(setup$block_sizes[[l]]), n)
    m[cbind(allocs[[l]]$block, allocs[[l]]$entry)] <- TRUE
    m
  })
  block_entries <- lapply(seq_along(allocs), function(l) {
    split(allocs[[l]]$entry, factor(allocs[[l]]$block, levels = seq_along(setup$block_sizes[[l]])))
  })

  best_score <- candidate$score
  best_state <- list(conc = conc, allocs = allocs)
  tries <- 0L
  stale <- 0L
  while (tries < max_tries && stale < 3000L) {
    bad <- which(conc >= 2L, arr.ind = TRUE)
    bad <- bad[bad[, 1] < bad[, 2], , drop = FALSE]
    if (nrow(bad) == 0) break
    tries <- tries + 1L
    stale <- stale + 1L

    pick <- bad[sample.int(nrow(bad), 1L), ]
    if (stats::runif(1) < 0.5) pick <- pick[2:1]
    a <- pick[1]
    x <- pick[2]

    # A plot of `a` in a block that also holds `x`.
    loc_ok <- which(vapply(presence, function(m) any(m[, a] & m[, x]), logical(1)))
    l <- loc_ok[sample.int(length(loc_ok), 1L)]
    al <- allocs[[l]]
    pres <- presence[[l]]
    shared <- which(pres[, a] & pres[, x])
    b1 <- shared[sample.int(length(shared), 1L)]
    p <- which(al$entry == a & al$block == b1)[1]
    block_rep <- setup$block_rep[[l]]

    same_type <- al$is_prep == al$is_prep[p] & al$block != b1
    if (al$is_prep[p] && setup$resolvable[l]) {
      same_type <- same_type & block_rep[al$block] == block_rep[b1]
    }
    cand_q <- which(same_type)
    cand_q <- cand_q[!pres[cbind(al$block[cand_q], a)] & !pres[b1, al$entry[cand_q]]]
    if (length(cand_q) == 0) next
    if (length(cand_q) > 25L) cand_q <- cand_q[sample.int(length(cand_q), 25L)]

    m1 <- block_entries[[l]][[b1]]
    m1 <- m1[m1 != a]
    delta <- vapply(cand_q, function(q) {
      cc <- al$entry[q]
      m2 <- block_entries[[l]][[al$block[q]]]
      m2 <- m2[m2 != cc]
      lose <- m1[!m1 %in% m2]
      gain <- m2[!m2 %in% m1]
      ca_g <- conc[a, gain]
      ca_l <- conc[a, lose]
      cc_l <- conc[cc, lose]
      cc_g <- conc[cc, gain]
      sum(pair_cost(ca_g + 1L) - pair_cost(ca_g)) +
        sum(pair_cost(ca_l - 1L) - pair_cost(ca_l)) +
        sum(pair_cost(cc_l + 1L) - pair_cost(cc_l)) +
        sum(pair_cost(cc_g - 1L) - pair_cost(cc_g))
    }, numeric(1))

    best <- which(delta == min(delta))
    if (min(delta) > 0 || (min(delta) == 0 && stats::runif(1) < 0.5)) next
    order_q <- cand_q[best[sample.int(length(best))]]

    q <- NA_integer_
    for (q_try in order_q) {
      if (al$is_prep[p]) {
        trial_block <- al$block
        trial_block[c(p, q_try)] <- trial_block[c(q_try, p)]
        comps <- prep_block_components(al$entry, trial_block, length(setup$block_sizes[[l]])) - 1L
        if (comps > al$extra_components) next
      }
      q <- q_try
      break
    }
    if (is.na(q)) next

    cc <- al$entry[q]
    b2 <- al$block[q]
    m2 <- block_entries[[l]][[b2]]
    m2 <- m2[m2 != cc]
    conc[a, m1] <- conc[a, m1] - 1L
    conc[m1, a] <- conc[a, m1]
    conc[cc, m2] <- conc[cc, m2] - 1L
    conc[m2, cc] <- conc[cc, m2]
    conc[a, m2] <- conc[a, m2] + 1L
    conc[m2, a] <- conc[a, m2]
    conc[cc, m1] <- conc[cc, m1] + 1L
    conc[m1, cc] <- conc[cc, m1]

    pres[b1, a] <- FALSE
    pres[b2, cc] <- FALSE
    pres[b1, cc] <- TRUE
    pres[b2, a] <- TRUE
    presence[[l]] <- pres
    block_entries[[l]][[b1]] <- c(m1, cc)
    block_entries[[l]][[b2]] <- c(m2, a)
    al$block[c(p, q)] <- al$block[c(q, p)]
    if (al$is_prep[p]) {
      al$extra_components <- prep_block_components(al$entry, al$block, length(setup$block_sizes[[l]])) - 1L
    }
    allocs[[l]] <- al

    if (min(delta) < 0) {
      candidate$concurrence <- conc
      candidate$allocations <- allocs
      score <- prep_candidate_score(candidate, setup)
      if (score < best_score) {
        best_score <- score
        best_state <- list(conc = conc, allocs = allocs)
        stale <- 0L
      }
    }
  }

  candidate$concurrence <- best_state$conc
  candidate$allocations <- best_state$allocs
  candidate$score <- best_score
  candidate
}

# Randomized search: greedy candidates, swap improvement of the best ones,
# and selection by average efficiency factor among optimal candidates.
# Candidates are generated on integer entry ids only; the plot-level design
# is built once, for the chosen candidate. The score has a known lower bound
# (no entry repeats in a block, connected locations, no pair of entries
# sharing a block more than once), so the search stops once
# `n_keep` candidates reach it.
prep_search_design <- function(setup, n_candidates, n_keep, n_treatments) {
  bound <- setup$score_lower_bound
  best <- NULL
  optimal <- list()
  top <- list()
  needs_improvement <- FALSE
  extra <- 0L
  n_evaluated <- 0L

  consider <- function(candidate) {
    if (is.null(best) || candidate$score < best$score) best <<- candidate
    if (candidate$score <= bound) optimal[[length(optimal) + 1L]] <<- candidate
  }
  improve_top <- function() {
    for (candidate in top) {
      improved <- prep_improve_candidate(candidate, setup)
      consider(improved)
      if (improved$score <= bound) {
        needs_improvement <<- TRUE
        break
      }
    }
    top <<- list()
  }

  for (i in seq_len(n_candidates)) {
    n_evaluated <- i
    candidate <- prep_make_candidate(setup)
    if (needs_improvement && candidate$score > bound) {
      candidate <- prep_improve_candidate(candidate, setup)
    }
    consider(candidate)
    if (length(optimal) >= n_keep) break
    if (length(optimal) > 0) {
      # Once an optimum exists, look a little further for more efficient ones.
      extra <- extra + 1L
      # Improved candidates are expensive: look at fewer of them.
      if (extra >= (if (needs_improvement) min(n_keep, 2L) else 3L * n_keep)) break
      next
    }
    top <- c(top, list(candidate))
    if (length(top) > 2L) {
      top <- top[order(vapply(top, `[[`, numeric(1), "score"))[1:2]]
    }
    if (i == 50L) improve_top()
  }
  if (length(optimal) == 0 && length(top) > 0) improve_top()

  chosen <- best
  if (length(optimal) > 1L) {
    eff <- vapply(optimal, prep_candidate_efficiency, numeric(1), v = n_treatments)
    if (!all(is.na(eff))) chosen <- optimal[[which.max(eff)]]
  }
  list(candidate = chosen, n_evaluated = n_evaluated)
}

# Moves plots of entries that the greedy step had to put twice in a block to
# another block (swapping with a plot of the same type), then recounts the
# concurrences of the candidate.
prep_repair_duplicates <- function(candidate, setup) {
  allocs <- candidate$allocations
  for (l in seq_along(allocs)) {
    al <- allocs[[l]]
    block_rep <- setup$block_rep[[l]]
    dup <- which(duplicated(paste(al$block, al$entry)))
    for (p in dup) {
      e <- al$entry[p]
      b1 <- al$block[p]
      ok <- which(al$is_prep == al$is_prep[p] & al$block != b1 &
                    (!al$is_prep[p] | !setup$resolvable[l] | block_rep[al$block] == block_rep[b1]))
      ok <- ok[vapply(ok, function(q) {
        !any(al$entry[al$block == al$block[q]] == e) && !any(al$entry[al$block == b1] == al$entry[q])
      }, logical(1))]
      if (length(ok) == 0) next
      q <- ok[sample.int(length(ok), 1L)]
      al$block[c(p, q)] <- al$block[c(q, p)]
    }
    al$extra_components <- prep_block_components(al$entry, al$block, length(setup$block_sizes[[l]])) - 1L
    allocs[[l]] <- al
  }
  n <- setup$n_treatments
  presence <- do.call(rbind, lapply(seq_along(allocs), function(l) {
    al <- allocs[[l]]
    m <- matrix(0, length(setup$block_sizes[[l]]), n)
    m[cbind(al$block, al$entry)] <- 1
    m
  }))
  conc <- crossprod(presence)
  diag(conc) <- 0
  storage.mode(conc) <- "integer"
  candidate$allocations <- allocs
  candidate$concurrence <- conc
  candidate$same_entry_same_block <- sum(vapply(allocs, function(al) sum(duplicated(paste(al$block, al$entry))), integer(1)))
  candidate$score <- prep_candidate_score(candidate, setup)
  candidate
}

# Plot-level design, concurrence matrix, and summary for the chosen candidate.
prep_build_candidate_design <- function(candidate, inputs, setup, serpentine) {
  treatments <- inputs$treatments
  locations <- inputs$locations

  groups <- stats::setNames(
    lapply(candidate$groups, function(ids) treatments[ids]),
    locations
  )
  repeat_map <- rep("not_repeated", length(treatments))
  for (l in seq_along(locations)) {
    repeat_map[candidate$groups[[l]]] <- locations[l]
  }

  rows <- lapply(seq_along(locations), function(l) {
    loc <- locations[l]
    alloc <- candidate$allocations[[l]]
    block_sizes <- setup$block_sizes[[l]]
    # Plots within a block keep the order in which they were allocated.
    o <- order(alloc$block, seq_along(alloc$block))
    entry <- alloc$entry[o]
    block <- alloc$block[o]
    repeated_at <- repeat_map[entry]
    data.frame(
      all_entries = treatments[entry],
      occurrence = alloc$occurrence[o],
      repeated_at = repeated_at,
      is_extra_repeat = as.integer(alloc$occurrence[o] > 1L),
      entry_type = ifelse(
        repeated_at == loc,
        "p_rep_entry_at_this_location",
        "single_entry_at_this_location"
      ),
      location = loc,
      replicate = setup$block_rep[[l]][block],
      block = block,
      block_uid = paste0(loc, "_B", sprintf("%03d", block)),
      block_size = block_sizes[block],
      plot_in_block = sequence(tabulate(block, nbins = length(block_sizes))),
      stringsAsFactors = FALSE
    )
  })

  design <- do.call(rbind, rows)
  design <- prep_add_field_coordinates(
    design,
    locations,
    inputs$field_cols,
    inputs$rows_per_block,
    inputs$blocks_per_row,
    serpentine
  )
  design <- prep_finalize_design(design, resolvable = isTRUE(inputs$resolvable))
  concurrence <- augmented_prep_concurrence(design, treatments)

  list(
    design = design,
    summary = prep_score_design(design, concurrence, resolvable = isTRUE(inputs$resolvable)),
    concurrence_matrix = concurrence,
    repeated_groups = groups
  )
}

prep_add_field_coordinates <- function(design, locations, field_cols,
                                       rows_per_block, blocks_per_row,
                                       serpentine) {
  design <- design[order(match(design$location, locations), design$block, design$plot_in_block), ]
  design$row <- NA_integer_
  design$col <- NA_integer_

  for (loc in locations) {
    idx <- which(design$location == loc)
    if (!is.na(rows_per_block[[loc]]) && !is.na(blocks_per_row[[loc]])) {
      design <- prep_add_block_layout_coordinates(
        design = design,
        idx = idx,
        rows_per_block = rows_per_block[[loc]],
        blocks_per_row = blocks_per_row[[loc]],
        serpentine = serpentine
      )
    } else if (!is.null(field_cols)) {
      plot_number <- seq_along(idx)
      n_cols <- field_cols[[loc]]
      row <- ceiling(plot_number / n_cols)
      col <- ((plot_number - 1L) %% n_cols) + 1L
      if (isTRUE(serpentine)) {
        col <- ifelse(row %% 2 == 0, n_cols - col + 1L, col)
      }
      design$row[idx] <- row
      design$col[idx] <- col
    }
  }

  rownames(design) <- NULL
  design
}

prep_add_block_layout_coordinates <- function(design, idx, rows_per_block,
                                              blocks_per_row, serpentine) {
  location_design <- design[idx, , drop = FALSE]
  block_ids <- unique(location_design$block)
  block_sizes <- vapply(
    block_ids,
    function(block) sum(location_design$block == block),
    integer(1)
  )
  block_width <- ceiling(max(block_sizes) / rows_per_block)

  for (block_index in seq_along(block_ids)) {
    block <- block_ids[block_index]
    block_idx <- idx[design$block[idx] == block]
    plot_number <- seq_along(block_idx)
    block_band <- (block_index - 1L) %/% blocks_per_row
    block_column <- (block_index - 1L) %% blocks_per_row
    local_row <- ceiling(plot_number / block_width)
    local_col <- ((plot_number - 1L) %% block_width) + 1L

    if (isTRUE(serpentine)) {
      local_col <- ifelse(local_row %% 2 == 0, block_width - local_col + 1L, local_col)
    }

    design$row[block_idx] <- (block_band * rows_per_block) + local_row
    design$col[block_idx] <- (block_column * block_width) + local_col
  }

  design
}

prep_finalize_design <- function(design, resolvable = FALSE) {
  design |>
    dplyr::arrange(.data$location, .data$block, .data$plot_in_block) |>
    dplyr::mutate(
      plots = dplyr::row_number(),
      rep = if (resolvable) .data$replicate else .data$occurrence,
      is_control = 0L
    ) |>
    dplyr::select(dplyr::all_of(c(
      "plots",
      "location",
      "row",
      "col",
      "replicate",
      "block",
      "block_uid",
      "block_size",
      "plot_in_block",
      "all_entries",
      "occurrence",
      "entry_type",
      "is_extra_repeat",
      "repeated_at",
      "rep",
      "is_control"
    ))) |>
    tibble::as_tibble()
}

prep_score_design <- function(design, concurrence, resolvable = FALSE) {
  upper <- concurrence[upper.tri(concurrence)]
  block_entries <- split(design$all_entries, interaction(design$location, design$block, drop = TRUE))
  same_entry_same_block <- sum(vapply(block_entries, function(x) sum(duplicated(x)), numeric(1)))
  replication_counts <- table(design$all_entries)
  max_conc <- if (length(upper) == 0) 0 else max(upper)
  n_pairs_gt1 <- sum(upper > 1)
  total_excess <- sum(pmax(upper - 1, 0))
  components <- vapply(split(design, design$location), function(x) {
    prep_block_components(as.character(x$all_entries), x$block, length(unique(x$block)))
  }, integer(1))
  score <- same_entry_same_block * 1e12 +
    sum(components - 1L) * 1e10 +
    sum(pmax(upper - 1, 0)^2) * 1e7 +
    max_conc * 1000 +
    sum((as.numeric(replication_counts) - mean(as.numeric(replication_counts)))^2)

  tibble::tibble(
    score = score,
    max_pairwise_concurrence = max_conc,
    n_pairs_with_concurrence_gt_1 = n_pairs_gt1,
    total_excess_concurrence = total_excess,
    same_entry_same_block = same_entry_same_block,
    n_disconnected_locations = sum(components > 1L),
    efficiency_factor = prep_efficiency_factor(design),
    resolvable = resolvable,
    min_replication = min(replication_counts),
    max_replication = max(replication_counts),
    mean_replication = mean(replication_counts),
    n_locations = dplyr::n_distinct(design$location),
    n_entries = dplyr::n_distinct(design$all_entries),
    n_plots = nrow(design),
    n_blocks = dplyr::n_distinct(design$block_uid)
  )
}

# Average efficiency factor E = (v - 1) / sum(1 / e_i) of the block design,
# where e_i are the canonical efficiency factors (eigenvalues of
# R^-1/2 C R^-1/2, with C = R - N K^-1 N'), for the model with entries and
# blocks only (Williams, Piepho and Whitaker, 2011, eq. 4). Returns 0 for a
# disconnected design and NA when there are more than `max_entries` entries.
prep_efficiency_factor <- function(design, max_entries = 2000L) {
  entries <- factor(design$all_entries)
  prep_efficiency_ids(
    as.integer(entries),
    paste(design$location, design$block),
    nlevels(entries),
    max_entries = max_entries
  )
}

prep_balanced_counts <- function(total, n_groups) {
  base <- total %/% n_groups
  rem <- total %% n_groups
  out <- rep(base, n_groups)
  if (rem > 0) out[seq_len(rem)] <- out[seq_len(rem)] + 1L
  out
}

prep_balanced_block_sizes <- function(n_plots, max_block_size, n_blocks = NULL) {
  if (n_plots <= 0) stop("`n_plots` must be positive.", call. = FALSE)
  if (max_block_size <= 1) stop("`block_size` must be greater than 1.", call. = FALSE)

  if (is.null(n_blocks)) {
    n_blocks <- ceiling(n_plots / max_block_size)
  } else {
    n_blocks <- arc_positive_integer(n_blocks, "n_blocks_per_location")
  }

  if (n_blocks > n_plots) {
    stop("`n_blocks_per_location` cannot exceed the number of plots at a location.", call. = FALSE)
  }

  out <- prep_balanced_counts(n_plots, n_blocks)
  if (max(out) > max_block_size) {
    stop(
      "`n_blocks_per_location` is too small for `block_size`; ",
      "at least one block would contain ", max(out), " plots.",
      call. = FALSE
    )
  }

  out
}

