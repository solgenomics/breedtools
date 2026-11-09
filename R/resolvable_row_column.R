# ============================================================================
# Resolvable (spatial) row-column designs
# Williams, E. R., John, J. A. and Whitaker, D. (2006). Construction of
# resolvable spatial row-column designs. Biometrics 62, 103-108.
# ============================================================================

#' Generate a resolvable (spatial) row-column design
#'
#' Builds a resolvable row-column design for `v <= k * s` treatments in `r`
#' replicates, each a `k x s` array of rows and columns that contains every
#' treatment once (Williams, John and Whitaker, 2006). The design is found by
#' a search over swaps of treatments within replicates:
#' \enumerate{
#'   \item The average efficiency factor of the resolvable row-column design
#'     (rows and columns within replicates), together with the efficiency
#'     factors of the component row design and column design, is maximized.
#'   \item With `spatial = TRUE`, the spatial average efficiency factor of the
#'     two-dimensional linear variance (LV + RC) model is then maximized while
#'     keeping the three efficiency factors of step 1 within
#'     `efficiency_tolerance` of their optimized values. This spreads
#'     treatments over edge and interior plots and over row and column
#'     neighbours, without needing prior values for spatial parameters.
#' }
#' When replicates are contiguous, `latinize = TRUE` builds a latinized
#' design: no treatment appears twice in the same long column (replicates
#' stacked vertically) or long row (replicates side by side), within every
#' group of `latin_t` adjacent replicates (t-latinized designs; John and
#' Williams, 1998).
#'
#' When there are fewer treatments than plots per replicate, the remaining
#' plots are fillers (`all_entries = "FILLER"`, `is_filler = 1`), placed by the
#' search like an extra treatment.
#'
#' @param treatments Character vector of treatment names (`v <= k * s`).
#' @param n_reps Number of replicates (`r`).
#' @param rows_per_rep Number of rows in each replicate (`k`).
#' @param cols_per_rep Number of columns in each replicate (`s`).
#' @param spatial Logical. If `TRUE` (default), optimize the spatial
#'   efficiency factor in a second stage.
#' @param latinize Logical. If `TRUE`, build a latinized design for
#'   contiguous replicates.
#' @param replicate_layout How replicates are placed in the field:
#'   `"vertical"` (stacked, giving long columns) or `"horizontal"` (side by
#'   side, giving long rows).
#' @param latin_t Optional size of the groups of adjacent replicates within
#'   which latinization holds (t-latinized design). Defaults to all
#'   replicates.
#' @param n_starts Number of random starting designs.
#' @param n_iter Number of swap proposals per stage and start. Defaults to
#'   `min(20000, max(3000, 60 * v * r))`.
#' @param efficiency_tolerance Largest loss in the row-column, row and column
#'   efficiency factors accepted in the spatial stage.
#' @param plot_type Plot numbering order: `"serpentine"` or `"cartesian"`.
#' @param seed Optional random seed.
#' @param output_file Optional tab-delimited output file for the design.
#' @param trial_name Optional Breedbase trial name. When given, a
#'   `plot_name` column with the plot names used by [save_breedbase()] is
#'   added to the design (see [add_plot_names()]), so the names can be
#'   checked with [breedbase_check()] and changed before saving.
#' @param location Breedbase location of the trial, for the `{location}`
#'   field of the plot names. Used only with `trial_name`.
#' @param plot_name_pattern Pattern for the plot names, as in
#'   [add_plot_names()]. Used only with `trial_name`.
#' @param start_design Optional existing design to start from (for example
#'   one produced by CycDesigN), as a data frame with `replicate`,
#'   `row_in_rep`, `col_in_rep` and `all_entries` columns and complete
#'   replicates. The search then improves this design; use `n_iter = 1` and
#'   `spatial = FALSE` to keep it essentially as it is and only lay it out.
#'
#' @return An object of class `resolvable_row_column_design`: a list with
#' `input`, `summary` (efficiency factors, spatial upper bound, latinization
#' check), `design` (one row per plot with `replicate`, `row_in_rep`,
#' `col_in_rep`, field `row` and `col`, `long_row` / `long_col`,
#' `all_entries`, `is_filler`), `breedbase_design`, and `efficiency` (see
#' [resolvable_row_column_efficiency()]).
#'
#' @references
#' Williams, E. R., John, J. A. and Whitaker, D. (2006). Construction of
#' resolvable spatial row-column designs. Biometrics 62, 103-108.
#' \doi{10.1111/j.1541-0420.2005.00393.x}
#'
#' John, J. A. and Williams, E. R. (1998). t-Latinized designs. Australian
#' and New Zealand Journal of Statistics 40, 111-118.
#'
#' @seealso [resolvable_row_column_efficiency()],
#'   [analyze_resolvable_row_column()], [plot_resolvable_row_column_design()]
#'
#' @examples
#' # Example 2 of Williams, John and Whitaker (2006): 42 treatments in three
#' # replicates of 6 rows x 7 columns.
#' design <- resolvable_row_column_design(
#'   treatments = paste0("T", 1:42),
#'   n_reps = 3,
#'   rows_per_rep = 6,
#'   cols_per_rep = 7,
#'   n_starts = 1,
#'   n_iter = 2000,
#'   seed = 1
#' )
#' design$summary
#'
#' @export
resolvable_row_column_design <- function(
    treatments,
    n_reps,
    rows_per_rep,
    cols_per_rep,
    spatial = TRUE,
    latinize = FALSE,
    replicate_layout = c("vertical", "horizontal"),
    latin_t = NULL,
    n_starts = 3,
    n_iter = NULL,
    efficiency_tolerance = 0.002,
    plot_type = c("serpentine", "cartesian"),
    seed = NULL,
    output_file = NULL,
    start_design = NULL,
    trial_name = NULL,
    location = NULL,
    plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}"
) {
  replicate_layout <- match.arg(replicate_layout)
  plot_type <- match.arg(plot_type)
  treatments <- as.character(treatments)
  treatments <- treatments[!is.na(treatments) & treatments != ""]
  if (anyDuplicated(treatments)) stop("`treatments` must be unique.", call. = FALSE)
  r <- arc_positive_integer(n_reps, "n_reps")
  k <- arc_positive_integer(rows_per_rep, "rows_per_rep")
  s <- arc_positive_integer(cols_per_rep, "cols_per_rep")
  n_starts <- arc_positive_integer(n_starts, "n_starts")
  v_real <- length(treatments)
  v <- k * s
  if (v_real < 2) stop("At least two treatments are required.", call. = FALSE)
  if (v_real > v) {
    stop("There are more treatments (", v_real, ") than plots per replicate (", v, ").", call. = FALSE)
  }
  if (r < 2) stop("`n_reps` must be at least 2.", call. = FALSE)
  if (k < 2 || s < 2) stop("Replicates need at least 2 rows and 2 columns.", call. = FALSE)
  latin_dim <- if (replicate_layout == "vertical") "col" else "row"
  t_size <- if (is.null(latin_t)) r else arc_positive_integer(latin_t, "latin_t")
  if (isTRUE(latinize)) {
    n_lines <- if (latin_dim == "col") s else k
    if (min(t_size, r) > n_lines) {
      stop(
        "A latinized design needs at most ", n_lines, " replicates per latin group ",
        "(the number of ", if (latin_dim == "col") "columns" else "rows",
        " per replicate). Reduce `latin_t`.",
        call. = FALSE
      )
    }
  }
  if (is.null(n_iter)) n_iter <- min(20000L, max(3000L, 60L * v * r))
  n_iter <- arc_positive_integer(n_iter, "n_iter")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  setup <- rrc_setup(v, r, k, s)
  latin <- if (isTRUE(latinize)) list(dim = latin_dim, group = (seq_len(r) - 1L) %/% t_size + 1L) else NULL

  start_pos <- NULL
  labels <- NULL
  if (!is.null(start_design)) {
    sp <- rrc_pos_from_design(start_design, treatments, r, k, s)
    start_pos <- sp$pos
    labels <- sp$labels
    if (!is.null(latin) && rrc_latin_conflicts(start_pos, setup, latin) > 0) {
      stop("`start_design` is not latinized; use `latinize = FALSE`.", call. = FALSE)
    }
    n_starts <- 1L
  }

  best <- NULL
  for (start in seq_len(n_starts)) {
    pos <- if (is.null(start_pos)) rrc_initial(v, r, latin, setup) else start_pos
    res <- if (n_iter <= 1L) {
      list(pos = pos, eff = rrc_all_efficiencies(pos, setup))
    } else {
      rrc_search(pos, setup, latin, n_iter, spatial, efficiency_tolerance)
    }
    if (is.null(best) || rrc_better(res$eff, best$eff, spatial)) best <- res
  }

  design <- rrc_build_design(best$pos, setup, treatments, replicate_layout, plot_type, labels)
  efficiency <- resolvable_row_column_efficiency(design)
  efficiency[names(best$eff)] <- best$eff
  latin_ok <- rrc_latin_ok(design, latin_dim, if (is.null(latin)) rep(1L, r) else latin$group)

  summary <- tibble::tibble(
    n_treatments = v_real,
    n_fillers = v - v_real,
    n_reps = r,
    rows_per_rep = k,
    cols_per_rep = s,
    n_plots = nrow(design),
    efficiency_rc = efficiency$efficiency_rc,
    efficiency_rows = efficiency$efficiency_rows,
    efficiency_cols = efficiency$efficiency_cols,
    efficiency_spatial = efficiency$efficiency_spatial,
    upper_bound_spatial = efficiency$upper_bound_spatial,
    spatial = isTRUE(spatial),
    latinized = isTRUE(latinize),
    latin_ok = if (isTRUE(latinize)) latin_ok else NA,
    replicate_layout = replicate_layout
  )

  breedbase <- design |>
    dplyr::transmute(
      plots = .data$plots,
      block = (.data$replicate - 1L) * k + .data$row_in_rep,
      all_entries = .data$all_entries,
      rep = .data$replicate,
      is_control = 0L,
      row_number = .data$row,
      col_number = .data$col
    )

  result <- list(
    input = list(
      treatments = treatments, n_reps = r, rows_per_rep = k, cols_per_rep = s,
      spatial = isTRUE(spatial), latinize = isTRUE(latinize),
      replicate_layout = replicate_layout, latin_t = t_size,
      n_starts = n_starts, n_iter = n_iter, efficiency_tolerance = efficiency_tolerance,
      plot_type = plot_type, seed = seed
    ),
    summary = summary,
    design = design,
    breedbase_design = breedbase,
    efficiency = efficiency
  )
  class(result) <- c("resolvable_row_column_design", class(result))
  if (!is.null(trial_name)) result <- add_plot_names(result, trial_name, location, plot_name_pattern)
  if (!is.null(output_file)) {
    utils::write.table(result$design, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  result
}

#' Efficiency of a resolvable row-column design
#'
#' Computes the design criteria of Williams, John and Whitaker (2006) for a
#' resolvable row-column design with complete replicates: the average
#' efficiency factors of the row-column design, the row design and the
#' column design (rows and columns nested within replicates), the spatial
#' average efficiency factor of the two-dimensional linear variance model,
#' \deqn{A = rI - X'(I_r \otimes W_R \otimes W_C)X,}
#' with its upper bound \eqn{U = (v - t)/(v - 1)}, \eqn{t = tr(W_R) tr(W_C)},
#' and treatment neighbour and edge statistics. Optionally, it evaluates the
#' average pairwise variance criterion \eqn{z(X)} under an assumed
#' autoregressive (AR) model, which is used to compare designs (Section 4 of
#' the paper).
#'
#' @param design A design object from [resolvable_row_column_design()] or a
#'   data frame with replicate, row, column and treatment columns. Each
#'   replicate must be complete (`k * s` plots, each treatment once; filler
#'   plots count as separate treatments).
#' @param ar Optional named list of AR model parameters
#'   `eta_r`, `eta_c`, `psi`, `rho_r`, `rho_c` (Eq. 4 of the paper).
#' @param genotype_col,replicate_col,row_col,col_col Column names.
#'
#' @return A list with `efficiency_rc`, `efficiency_rows`,
#'   `efficiency_cols`, `efficiency_spatial`, `upper_bound_spatial`,
#'   `ar_criterion` (`NA` unless `ar` is given), and `neighbours`, a tibble
#'   with each treatment's number of distinct row and column neighbours and
#'   the number of replicates in which it is on an edge or in a corner.
#'
#' @examples
#' design <- resolvable_row_column_design(paste0("T", 1:20), 2, 4, 5,
#'                                        n_starts = 1, n_iter = 500, seed = 1)
#' eff <- resolvable_row_column_efficiency(design,
#'   ar = list(eta_r = 35, eta_c = 60, psi = 10, rho_r = 0.5, rho_c = 0.5))
#' eff[c("efficiency_rc", "efficiency_spatial", "ar_criterion")]
#'
#' @export
resolvable_row_column_efficiency <- function(design, ar = NULL,
                                             genotype_col = "all_entries",
                                             replicate_col = "replicate",
                                             row_col = "row_in_rep",
                                             col_col = "col_in_rep") {
  if (is.list(design) && !is.data.frame(design) && !is.null(design$design)) design <- design$design
  d <- as.data.frame(design, stringsAsFactors = FALSE)
  missing_cols <- setdiff(c(genotype_col, replicate_col, row_col, col_col), names(d))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  rep_id <- as.integer(factor(d[[replicate_col]]))
  row <- as.integer(d[[row_col]])
  col <- as.integer(d[[col_col]])
  k <- max(row)
  s <- max(col)
  r <- max(rep_id)
  trt <- as.character(d[[genotype_col]])
  filler <- if ("is_filler" %in% names(d)) d$is_filler == 1 else rep(FALSE, nrow(d))
  # Each filler plot is its own pseudo-treatment, matched across replicates
  # by its order within the replicate.
  trt[filler] <- paste0(".filler", stats::ave(seq_along(trt), rep_id, filler, FUN = seq_along))[filler]
  v <- k * s
  levels_t <- unique(trt)
  if (nrow(d) != r * v || length(levels_t) != v ||
      any(table(factor(trt, levels = levels_t), rep_id) != 1)) {
    stop(
      "Efficiency factors need complete replicates: every replicate must have ",
      "k x s plots with each treatment once.",
      call. = FALSE
    )
  }
  setup <- rrc_setup(v, r, k, s)
  pos <- matrix(0L, r, v)
  pos[cbind(rep_id, match(trt, levels_t))] <- (col - 1L) * k + row
  eff <- rrc_all_efficiencies(pos, setup)

  ar_criterion <- NA_real_
  if (!is.null(ar)) {
    ar_criterion <- rrc_ar_criterion(pos, setup, ar)
  }

  # Neighbours and edge positions.
  grid <- array(NA_integer_, c(r, k, s))
  grid[cbind(rep_id, row, col)] <- match(trt, levels_t)
  row_nb <- col_nb <- vector("list", v)
  for (q in seq_len(r)) {
    g <- grid[q, , ]
    for (i in seq_len(k)) for (j in seq_len(s)) {
      t0 <- g[i, j]
      if (j > 1) row_nb[[t0]] <- c(row_nb[[t0]], g[i, j - 1])
      if (j < s) row_nb[[t0]] <- c(row_nb[[t0]], g[i, j + 1])
      if (i > 1) col_nb[[t0]] <- c(col_nb[[t0]], g[i - 1, j])
      if (i < k) col_nb[[t0]] <- c(col_nb[[t0]], g[i + 1, j])
    }
  }
  t_index <- match(trt, levels_t)
  edge <- row %in% c(1L, k) | col %in% c(1L, s)
  corner <- row %in% c(1L, k) & col %in% c(1L, s)
  neighbours <- tibble::tibble(
    treatment = levels_t,
    row_neighbours = vapply(row_nb, function(x) length(unique(x)), integer(1)),
    col_neighbours = vapply(col_nb, function(x) length(unique(x)), integer(1)),
    neighbours = vapply(seq_len(v), function(t0) length(unique(c(row_nb[[t0]], col_nb[[t0]]))), integer(1)),
    reps_on_edge = as.integer(tapply(edge, factor(t_index, levels = seq_len(v)), sum)),
    reps_in_corner = as.integer(tapply(corner, factor(t_index, levels = seq_len(v)), sum))
  )
  neighbours <- neighbours[!grepl("^\\.filler", neighbours$treatment), , drop = FALSE]

  c(eff, list(ar_criterion = ar_criterion, neighbours = neighbours))
}

#' Plot a resolvable row-column design
#'
#' Draws the field with replicate boundaries.
#'
#' @param design An object from [resolvable_row_column_design()] or its
#'   `design` data frame.
#' @param label Cell label: `"all_entries"`, `"replicate"`, or `"none"`.
#' @param fill Cell colouring: `"replicate"` or `"is_filler"`.
#' @param main Optional title.
#' @param ... Passed to [plot_augmented_row_column_design()].
#'
#' @return Invisibly, the plotted data.
#'
#' @examples
#' design <- resolvable_row_column_design(paste0("T", 1:20), 2, 4, 5,
#'                                        n_starts = 1, n_iter = 300, seed = 1)
#' plot_resolvable_row_column_design(design)
#'
#' @export
plot_resolvable_row_column_design <- function(design, label = c("all_entries", "replicate", "none"),
                                              fill = c("replicate", "is_filler"), main = NULL, ...) {
  label <- match.arg(label)
  fill <- match.arg(fill)
  if (is.list(design) && !is.data.frame(design) && !is.null(design$design)) design <- design$design
  d <- as.data.frame(design, stringsAsFactors = FALSE)
  d$block <- d$replicate
  d$type <- ifelse(d$is_filler == 1, "filler", "treatment")
  d$is_control <- d$is_filler
  plot_label <- if (label == "replicate") "block" else label
  plot_fill <- if (fill == "replicate") "block" else "is_control"
  plot_augmented_row_column_design(
    d, label = plot_label, fill = plot_fill,
    main = if (is.null(main)) "Resolvable row-column design" else main, ...
  )
}

# ---- Internals: criteria ------------------------------------------------------

# Weight matrix of the LV + RC model (Williams, John and Whitaker, 2006):
# element (i, j), i <= j, is 6 / (k (k^2 - 1)) * (sum_{h<i} h + sum_{h<=k-j} h).
rrc_weight_matrix <- function(k) {
  cst <- 6 / (k * (k^2 - 1))
  i <- matrix(seq_len(k), k, k)
  j <- t(i)
  lo <- pmin(i, j)
  hi <- pmax(i, j)
  cst * (lo * (lo - 1) / 2 + (k - hi) * (k - hi + 1) / 2)
}

# Within-replicate plot matrices B (plots ordered column by column:
# p = (col - 1) k + row) for the four criteria. Each criterion matrix is
# C = r I - sum_q B[pos_q(a), pos_q(b)].
rrc_setup <- function(v, r, k, s) {
  PR <- kronecker(matrix(1 / s, s, s), diag(k))
  PC <- kronecker(diag(s), matrix(1 / k, k, k))
  K <- matrix(1 / (k * s), k * s, k * s)
  WR <- rrc_weight_matrix(s)
  WC <- rrc_weight_matrix(k)
  list(
    v = v, r = r, k = k, s = s,
    B = list(rc = PR + PC - K, rows = PR, cols = PC, spatial = kronecker(WR, WC)),
    upper_bound_spatial = (v - sum(diag(WR)) * sum(diag(WC))) / (v - 1)
  )
}

rrc_G <- function(pos, B) {
  G <- 0
  for (q in seq_len(nrow(pos))) G <- G + B[pos[q, ], pos[q, ]]
  G
}

rrc_inverse_state <- function(pos, setup, name) {
  C <- setup$r * diag(setup$v) - rrc_G(pos, setup$B[[name]])
  M <- tryCatch(solve(C + 1 / setup$v), error = function(e) NULL)
  if (is.null(M)) return(NULL)
  list(M = M, trace = sum(diag(M)))
}

rrc_eff_from_trace <- function(trace, setup) {
  tr <- trace - 1
  if (!is.finite(tr) || tr <= 0) return(0)
  (setup$v - 1) / (setup$r * tr)
}

rrc_all_efficiencies <- function(pos, setup) {
  out <- vapply(c("rc", "rows", "cols", "spatial"), function(nm) {
    st <- rrc_inverse_state(pos, setup, nm)
    if (is.null(st)) 0 else rrc_eff_from_trace(st$trace, setup)
  }, numeric(1))
  list(
    efficiency_rc = out[["rc"]],
    efficiency_rows = out[["rows"]],
    efficiency_cols = out[["cols"]],
    efficiency_spatial = out[["spatial"]],
    upper_bound_spatial = setup$upper_bound_spatial
  )
}

# Woodbury update for swapping treatments a and b in replicate q. The change
# of G is supported on rows and columns a and b:
# dG = u d' + d u' + dd (e_a e_a' - e_b e_b'), u = e_a - e_b.
rrc_swap_update <- function(st, pos, q, a, b, B, setup) {
  pa <- pos[q, a]
  pb <- pos[q, b]
  d <- B[pb, pos[q, ]] - B[pa, pos[q, ]]
  d[c(a, b)] <- 0
  dd <- B[pb, pb] - B[pa, pa]
  v <- setup$v
  M <- st$M
  Mu <- M[, a] - M[, b]
  Md <- drop(M %*% d)
  if (abs(dd) > 1e-12) {
    MU <- cbind(Mu, Md, M[, a], M[, b])
    Kinv <- -rbind(c(0, 1, 0, 0), c(1, 0, 0, 0), c(0, 0, 1 / dd, 0), c(0, 0, 0, -1 / dd))
    UtMU <- rbind(
      c(Mu[a] - Mu[b], Md[a] - Md[b], M[a, a] - M[b, a], M[a, b] - M[b, b]),
      c(sum(d * Mu), sum(d * Md), Md[a], Md[b]),
      c(Mu[a], Md[a], M[a, a], M[a, b]),
      c(Mu[b], Md[b], M[b, a], M[b, b])
    )
  } else {
    MU <- cbind(Mu, Md)
    Kinv <- -rbind(c(0, 1), c(1, 0))
    UtMU <- rbind(
      c(Mu[a] - Mu[b], Md[a] - Md[b]),
      c(sum(d * Mu), sum(d * Md))
    )
  }
  # dC = -dG = U K U' with K = -[[0,1],[1,0]] (and the dd block), so the
  # inverse changes by -MU S^-1 (MU)' with S = K^-1 + U'MU.
  S <- Kinv + UtMU
  S_inv <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(S_inv)) return(NULL)
  delta_trace <- -sum((MU %*% S_inv) * MU)
  list(MU = MU, S_inv = S_inv, trace = st$trace + delta_trace)
}

rrc_apply_update <- function(st, upd) {
  list(M = st$M - upd$MU %*% upd$S_inv %*% t(upd$MU), trace = upd$trace)
}

rrc_better <- function(a, b, spatial) {
  if (spatial) {
    a$efficiency_spatial > b$efficiency_spatial + 1e-12
  } else {
    (a$efficiency_rc + 0.05 * (a$efficiency_rows + a$efficiency_cols)) >
      (b$efficiency_rc + 0.05 * (b$efficiency_rows + b$efficiency_cols)) + 1e-12
  }
}

# ---- Internals: search ----------------------------------------------------------

# Random starting design: pos[q, t] is the plot of treatment t in replicate
# q. For latinized designs, conflicts (a treatment twice in a long line
# within a latin group) are removed by swaps first.
rrc_initial <- function(v, r, latin, setup) {
  for (attempt in 1:20) {
    pos <- t(vapply(seq_len(r), function(q) sample.int(v), integer(v)))
    if (is.null(latin)) return(pos)
    pos <- rrc_fix_latin(pos, setup, latin)
    if (!is.null(pos)) return(pos)
  }
  stop("Could not construct a latinized starting design; reduce `latin_t`.", call. = FALSE)
}

rrc_line <- function(p, setup, dim) {
  if (dim == "col") (p - 1L) %/% setup$k + 1L else (p - 1L) %% setup$k + 1L
}

rrc_latin_conflicts <- function(pos, setup, latin) {
  lines <- matrix(rrc_line(pos, setup, latin$dim), nrow(pos))
  sum(vapply(unique(latin$group), function(g) {
    l <- lines[latin$group == g, , drop = FALSE]
    sum(apply(l, 2, function(x) sum(duplicated(x))))
  }, numeric(1)))
}

rrc_fix_latin <- function(pos, setup, latin, max_iter = 200000L) {
  v <- setup$v
  conflicts <- rrc_latin_conflicts(pos, setup, latin)
  iter <- 0L
  while (conflicts > 0 && iter < max_iter) {
    iter <- iter + 1L
    q <- sample.int(nrow(pos), 1L)
    ab <- sample.int(v, 2L)
    trial <- pos
    trial[q, ab] <- pos[q, rev(ab)]
    tc <- rrc_latin_conflicts(trial, setup, latin)
    if (tc <= conflicts) {
      pos <- trial
      conflicts <- tc
    }
  }
  if (conflicts > 0) NULL else pos
}

# Would swapping a and b in replicate q break latinization?
rrc_swap_keeps_latin <- function(pos, setup, latin, q, a, b) {
  if (is.null(latin)) return(TRUE)
  la <- rrc_line(pos[q, a], setup, latin$dim)
  lb <- rrc_line(pos[q, b], setup, latin$dim)
  if (la == lb) return(TRUE)
  others <- which(latin$group == latin$group[q] & seq_len(nrow(pos)) != q)
  if (length(others) == 0) return(TRUE)
  !any(rrc_line(pos[others, a], setup, latin$dim) == lb) &&
    !any(rrc_line(pos[others, b], setup, latin$dim) == la)
}

rrc_search <- function(pos, setup, latin, n_iter, spatial, tol, t0 = 2e-3, t1 = 1e-6) {
  names_c <- c("rc", "rows", "cols", "spatial")
  states <- lapply(names_c, function(nm) rrc_inverse_state(pos, setup, nm))
  names(states) <- names_c
  if (any(vapply(states, is.null, logical(1)))) {
    states <- lapply(names_c, function(nm) {
      C <- setup$r * diag(setup$v) - rrc_G(pos, setup$B[[nm]])
      M <- rrc_ginv(C + 1 / setup$v)
      list(M = M, trace = sum(diag(M)))
    })
    names(states) <- names_c
  }
  eff_of <- function(st) vapply(st, function(x) rrc_eff_from_trace(x$trace, setup), numeric(1))

  # Simulated annealing over swaps within replicates; the best feasible
  # design visited is kept.
  run_stage <- function(pos, states, objective, feasible, n_iter) {
    current <- eff_of(states)
    best_pos <- pos
    best_obj <- objective(current)
    accepted <- 0L
    cooling <- (t1 / t0)^(1 / max(1, n_iter - 1))
    temp <- t0
    for (it in seq_len(n_iter)) {
      temp <- temp * cooling
      q <- sample.int(setup$r, 1L)
      ab <- rrc_propose_swap(pos[q, ], setup)
      if (!rrc_swap_keeps_latin(pos, setup, latin, q, ab[1], ab[2])) next
      upd <- lapply(names_c, function(nm) rrc_swap_update(states[[nm]], pos, q, ab[1], ab[2], setup$B[[nm]], setup))
      if (any(vapply(upd, is.null, logical(1)))) next
      names(upd) <- names_c
      new_eff <- vapply(upd, function(x) rrc_eff_from_trace(x$trace, setup), numeric(1))
      if (!feasible(new_eff)) next
      delta <- objective(new_eff) - objective(current)
      if (delta >= -1e-12 || stats::runif(1) < exp(delta / temp)) {
        for (nm in names_c) states[[nm]] <- rrc_apply_update(states[[nm]], upd[[nm]])
        pos[q, ab] <- pos[q, rev(ab)]
        current <- new_eff
        accepted <- accepted + 1L
        if (accepted %% 400L == 0L) {
          states <- lapply(names_c, function(nm) rrc_inverse_state(pos, setup, nm))
          names(states) <- names_c
          current <- eff_of(states)
        }
        if (objective(current) > best_obj + 1e-12) {
          best_obj <- objective(current)
          best_pos <- pos
        }
      }
    }
    states <- lapply(names_c, function(nm) rrc_inverse_state(best_pos, setup, nm))
    names(states) <- names_c
    list(pos = best_pos, states = states, eff = eff_of(states))
  }

  # Stage 0: balance row and column concurrences (fast, combinatorial).
  # Stage 1: tabu search on the row-column efficiency factor.
  # Stage 2 (spatial): annealing on the spatial efficiency factor, keeping
  # the row-column, row and column efficiency factors within `tol`.
  pos <- rrc_concurrence_stage(pos, setup, latin, n_iter = min(80000L, 10L * n_iter))
  pos <- rrc_tabu_stage(pos, setup, latin, n_steps = max(100L, n_iter %/% 20L))
  states <- lapply(names_c, function(nm) rrc_inverse_state(pos, setup, nm))
  names(states) <- names_c
  stage1 <- run_stage(
    pos, states,
    objective = function(e) e[["rc"]] + 0.05 * (e[["rows"]] + e[["cols"]]),
    feasible = function(e) TRUE,
    n_iter = n_iter %/% 4L
  )
  result <- stage1
  if (isTRUE(spatial)) {
    target <- stage1$eff
    result <- run_stage(
      stage1$pos, stage1$states,
      objective = function(e) e[["spatial"]],
      feasible = function(e) all(e[c("rc", "rows", "cols")] >= target[c("rc", "rows", "cols")] - tol),
      n_iter = n_iter
    )
  }
  list(pos = result$pos, eff = rrc_all_efficiencies(result$pos, setup))
}

# Stage 0: simulated annealing on the concurrence penalty
# sum over pairs of R^2 + C^2 + (R + C)^2, with R and C the numbers of rows
# and columns in which a pair of treatments meets.
rrc_concurrence_stage <- function(pos, setup, latin, n_iter, t0 = 2, t1 = 0.01) {
  v <- setup$v
  k <- setup$k
  r <- setup$r
  rowof <- (pos - 1L) %% k
  colof <- (pos - 1L) %/% k
  R <- matrix(0L, v, v)
  C <- matrix(0L, v, v)
  for (q in seq_len(r)) {
    R <- R + outer(rowof[q, ], rowof[q, ], "==")
    C <- C + outer(colof[q, ], colof[q, ], "==")
  }
  diag(R) <- 0L
  diag(C) <- 0L
  pen <- function(Rv, Cv) Rv^2 + Cv^2 + (Rv + Cv)^2
  cooling <- (t1 / t0)^(1 / max(1, n_iter - 1))
  temp <- t0
  score <- 0
  best <- pos
  best_score <- 0
  for (it in seq_len(n_iter)) {
    temp <- temp * cooling
    q <- sample.int(r, 1L)
    ab <- rrc_propose_swap(pos[q, ], setup)
    a <- ab[1]
    b <- ab[2]
    if (!rrc_swap_keeps_latin(pos, setup, latin, q, a, b)) next
    ra <- rowof[q, a]; rb <- rowof[q, b]; ca <- colof[q, a]; cb <- colof[q, b]
    others <- seq_len(v) != a & seq_len(v) != b
    mRa <- which(rowof[q, ] == ra & others)
    mRb <- which(rowof[q, ] == rb & others)
    mCa <- which(colof[q, ] == ca & others)
    mCb <- which(colof[q, ] == cb & others)
    idx <- unique(c(mRa, mRb, mCa, mCb))
    dRa <- dCa <- numeric(v)
    if (ra != rb) {
      dRa[mRa] <- -1
      dRa[mRb] <- 1
    }
    if (ca != cb) {
      dCa[mCa] <- -1
      dCa[mCb] <- 1
    }
    # b moves the opposite way.
    delta <- sum(pen(R[a, idx] + dRa[idx], C[a, idx] + dCa[idx]) - pen(R[a, idx], C[a, idx])) +
      sum(pen(R[b, idx] - dRa[idx], C[b, idx] - dCa[idx]) - pen(R[b, idx], C[b, idx]))
    if (delta <= 0 || stats::runif(1) < exp(-delta / temp)) {
      R[a, ] <- R[a, ] + dRa
      R[, a] <- R[a, ]
      C[a, ] <- C[a, ] + dCa
      C[, a] <- C[a, ]
      R[b, ] <- R[b, ] - dRa
      R[, b] <- R[b, ]
      C[b, ] <- C[b, ] - dCa
      C[, b] <- C[b, ]
      pos[q, ab] <- pos[q, rev(ab)]
      rowof[q, ab] <- rowof[q, rev(ab)]
      colof[q, ab] <- colof[q, rev(ab)]
      score <- score + delta
      if (score < best_score) {
        best_score <- score
        best <- pos
      }
    }
  }
  best
}

# Stage 1: tabu search on the row-column efficiency factor. Each step takes
# the best non-tabu swap within a row or within a column (the full
# neighbourhood for small designs, a random sample of 300 swaps otherwise).
rrc_tabu_stage <- function(pos, setup, latin, n_steps, tenure = 20L, max_candidates = 300L) {
  k <- setup$k
  v <- setup$v
  st <- rrc_inverse_state(pos, setup, "rc")
  if (is.null(st)) return(pos)
  best <- pos
  best_e <- rrc_eff_from_trace(st$trace, setup)
  tabu <- character(0)
  for (step in seq_len(n_steps)) {
    cand <- list()
    for (q in seq_len(setup$r)) {
      rowof <- (pos[q, ] - 1L) %% k
      colof <- (pos[q, ] - 1L) %/% k
      for (line in list(rowof, colof)) {
        for (g in unique(line)) {
          m <- which(line == g)
          if (length(m) < 2) next
          pr <- utils::combn(m, 2)
          cand[[length(cand) + 1L]] <- cbind(q, t(pr))
        }
      }
    }
    cand <- do.call(rbind, cand)
    if (nrow(cand) > max_candidates) cand <- cand[sample.int(nrow(cand), max_candidates), , drop = FALSE]
    best_val <- -Inf
    chosen <- NULL
    for (i in seq_len(nrow(cand))) {
      q <- cand[i, 1]
      a <- cand[i, 2]
      b <- cand[i, 3]
      if (!rrc_swap_keeps_latin(pos, setup, latin, q, a, b)) next
      upd <- rrc_swap_update(st, pos, q, a, b, setup$B$rc, setup)
      if (is.null(upd)) next
      val <- rrc_eff_from_trace(upd$trace, setup)
      key <- paste(q, a, b)
      if (key %in% tabu && val <= best_e) next
      if (val > best_val) {
        best_val <- val
        chosen <- list(q = q, ab = c(a, b), upd = upd, key = key)
      }
    }
    if (is.null(chosen)) break
    st <- rrc_apply_update(st, chosen$upd)
    pos[chosen$q, chosen$ab] <- pos[chosen$q, rev(chosen$ab)]
    tabu <- c(utils::tail(tabu, tenure - 1L), chosen$key)
    if (best_val > best_e + 1e-12) {
      best_e <- best_val
      best <- pos
    }
    if (step %% 100L == 0L) st <- rrc_inverse_state(pos, setup, "rc")
  }
  best
}

rrc_ginv <- function(x) {
  e <- eigen((x + t(x)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-10 * max(e$values)
  e$vectors[, keep, drop = FALSE] %*% (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
}

# Swap proposals: two treatments in the same row (changes columns only), in
# the same column (changes rows only), or anywhere in the replicate.
rrc_propose_swap <- function(p, setup) {
  k <- setup$k
  u <- stats::runif(1)
  a <- sample.int(setup$v, 1L)
  row_a <- (p[a] - 1L) %% k
  col_a <- (p[a] - 1L) %/% k
  cand <- if (u < 0.4) {
    which((p - 1L) %% k == row_a)
  } else if (u < 0.8) {
    which((p - 1L) %/% k == col_a)
  } else {
    seq_len(setup$v)
  }
  cand <- cand[cand != a]
  if (length(cand) == 0) cand <- setdiff(seq_len(setup$v), a)
  c(a, cand[sample.int(length(cand), 1L)])
}

rrc_build_design <- function(pos, setup, treatments, layout, plot_type, labels = NULL) {
  k <- setup$k
  s <- setup$s
  r <- setup$r
  v <- setup$v
  if (is.null(labels)) {
    labels <- c(treatments, rep("FILLER", v - length(treatments)))
    labels <- labels[sample.int(v)]
  }
  rows <- lapply(seq_len(r), function(q) {
    p <- pos[q, ]
    data.frame(
      replicate = q,
      row_in_rep = (p - 1L) %% k + 1L,
      col_in_rep = (p - 1L) %/% k + 1L,
      all_entries = labels,
      stringsAsFactors = FALSE
    )
  })
  d <- do.call(rbind, rows)
  if (layout == "vertical") {
    d$row <- (d$replicate - 1L) * k + d$row_in_rep
    d$col <- d$col_in_rep
  } else {
    d$row <- d$row_in_rep
    d$col <- (d$replicate - 1L) * s + d$col_in_rep
  }
  d$long_row <- d$row
  d$long_col <- d$col
  d$is_filler <- as.integer(d$all_entries == "FILLER")
  d$type <- ifelse(d$is_filler == 1L, "filler", "treatment")
  d <- arc_make_plot_order(d, plot_type = plot_type)
  d$rep <- d$replicate
  d$is_control <- 0L
  tibble::as_tibble(d[, c("plots", "replicate", "row_in_rep", "col_in_rep", "row", "col",
                          "long_row", "long_col", "all_entries", "type", "is_filler",
                          "rep", "is_control")])
}

# Positions (pos[q, t]) and treatment labels of an existing design.
rrc_pos_from_design <- function(d, treatments, r, k, s) {
  if (is.list(d) && !is.data.frame(d) && !is.null(d$design)) d <- d$design
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  need <- c("replicate", "row_in_rep", "col_in_rep", "all_entries")
  if (!all(need %in% names(d))) {
    stop("`start_design` needs columns ", paste(need, collapse = ", "), ".", call. = FALSE)
  }
  rep_id <- as.integer(factor(d$replicate))
  if (max(rep_id) != r || max(d$row_in_rep) != k || max(d$col_in_rep) != s || nrow(d) != r * k * s) {
    stop("`start_design` does not match `n_reps`, `rows_per_rep` and `cols_per_rep`.", call. = FALSE)
  }
  trt <- as.character(d$all_entries)
  real <- trt[trt %in% treatments]
  if (!setequal(unique(real), treatments)) {
    stop("`start_design` must contain every treatment in `treatments`.", call. = FALSE)
  }
  filler <- !trt %in% treatments
  trt[filler] <- paste0(".filler", stats::ave(seq_along(trt), rep_id, filler, FUN = seq_along))[filler]
  labels_id <- unique(trt)
  if (length(labels_id) != k * s || any(table(factor(trt, levels = labels_id), rep_id) != 1)) {
    stop("`start_design` must have complete replicates with each treatment once.", call. = FALSE)
  }
  pos <- matrix(0L, r, k * s)
  pos[cbind(rep_id, match(trt, labels_id))] <- (as.integer(d$col_in_rep) - 1L) * k + as.integer(d$row_in_rep)
  labels <- ifelse(grepl("^\\.filler", labels_id), "FILLER", labels_id)
  list(pos = pos, labels = labels)
}

rrc_latin_ok <- function(design, dim, group) {
  d <- design[design$is_filler == 0, , drop = FALSE]
  line <- if (dim == "col") d$col_in_rep else d$row_in_rep
  g <- group[d$replicate]
  !any(duplicated(paste(g, d$all_entries, line)))
}

# Average pairwise variance criterion under the AR model (Eq. 4):
# z(X) = trace((X'V*+X)^+) / (v - 1), V* = (I - P_G) V (I - P_G).
rrc_ar_criterion <- function(pos, setup, ar) {
  need <- c("eta_r", "eta_c", "psi", "rho_r", "rho_c")
  if (!all(need %in% names(ar))) {
    stop("`ar` must contain ", paste(need, collapse = ", "), ".", call. = FALSE)
  }
  k <- setup$k
  s <- setup$s
  r <- setup$r
  v <- setup$v
  n <- r * v
  DC <- ar$rho_c^abs(outer(seq_len(s), seq_len(s), "-"))
  DR <- ar$rho_r^abs(outer(seq_len(k), seq_len(k), "-"))
  within <- diag(k * s) + ar$eta_r * setup$B$rows + ar$eta_c * setup$B$cols + ar$psi * kronecker(DC, DR)
  V <- kronecker(diag(r), within)
  PG <- matrix(1 / n, n, n)
  Vs <- (diag(n) - PG) %*% V %*% (diag(n) - PG)
  X <- matrix(0, n, v)
  for (q in seq_len(r)) X[cbind((q - 1L) * v + pos[q, ], seq_len(v))] <- 1
  Vp <- rrc_ginv(Vs)
  info <- t(X) %*% Vp %*% X
  sum(diag(rrc_ginv(info))) / (v - 1)
}
