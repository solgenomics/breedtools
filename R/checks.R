# ============================================================================
# Simple checks for users: record [OK] / [FAIL] lines, standard design checks,
# and the accuracy of adjusted means in simulations.
# ============================================================================

.check_env <- new.env(parent = emptyenv())
.check_env$log <- NULL

#' Record and print a check
#'
#' Prints `[OK]` or `[FAIL]` followed by a description, and records the
#' result so that [check_summary()] can report how many checks passed. Use it
#' to write your own checks while exploring a design or an analysis.
#'
#' @param what Description of the check, e.g. `"every entry planted once"`.
#' @param ok Logical. `TRUE` when the check passes; anything else (including
#'   `NA`) counts as a failure.
#' @param section Optional label for grouping checks, e.g. `"p-rep"`.
#'   Defaults to the label set with [check_section()].
#' @param note Logical. If `TRUE`, the check is informational: it prints
#'   `[NOTE]` instead of `[FAIL]` when `ok` is not `TRUE` and is not counted
#'   as a failure.
#'
#' @return Invisibly, `TRUE` if the check passed.
#'
#' @seealso [check_section()], [check_summary()], [check_reset()], [check_design()]
#'
#' @examples
#' check_reset()
#' check_result("two plus two is four", 2 + 2 == 4)
#' check_result("a deliberate failure", FALSE, section = "demo")
#' check_summary()
#' check_reset()
#'
#' @export
check_result <- function(what, ok, section = NULL, note = FALSE) {
  if (is.null(section)) section <- if (is.null(.check_env$section)) "" else .check_env$section
  passed <- isTRUE(ok)
  status <- if (passed) "OK" else if (isTRUE(note)) "NOTE" else "FAIL"
  bt_say(sprintf("  %-7s%s\n", paste0("[", status, "]"), what))
  row <- data.frame(section = as.character(section), check = as.character(what),
                    status = status, stringsAsFactors = FALSE)
  .check_env$log <- rbind(.check_env$log, row)
  invisible(passed)
}

#' Summary of the recorded checks
#'
#' Reports what was checked and what is correct: how many of the checks
#' recorded with [check_result()] and [check_design()] passed, broken down by
#' section, followed by the checks that failed and the notes.
#'
#' @param details Logical. If `TRUE`, also list every check with its status.
#' @param print Logical. If `TRUE` (default), print the report.
#'
#' @return Invisibly, a tibble of all recorded checks with columns `section`,
#'   `check` and `status` (`"OK"`, `"FAIL"` or `"NOTE"`).
#'
#' @examples
#' check_reset()
#' check_result("one equals one", 1 == 1, section = "arithmetic")
#' check_result("two is even", 2 %% 2 == 0, section = "arithmetic")
#' check_summary()
#' check_summary(details = TRUE)
#' check_reset()
#'
#' @export
check_summary <- function(details = FALSE, print = TRUE) {
  log <- .check_env$log
  if (is.null(log)) {
    log <- data.frame(section = character(), check = character(), status = character())
  }
  log <- tibble::as_tibble(log)
  if (!isTRUE(print)) return(invisible(log))

  if (nrow(log) == 0) {
    cat("No checks recorded yet. Use check_design() or check_result().\n")
    return(invisible(log))
  }
  n_ok <- sum(log$status == "OK")
  n_fail <- sum(log$status == "FAIL")
  n_note <- sum(log$status == "NOTE")
  cat(sprintf("\n%d of %d checks passed%s\n", n_ok, n_ok + n_fail,
              if (n_note > 0) sprintf(" (plus %d notes)", n_note) else ""))

  sections <- ifelse(nzchar(log$section), log$section, "(no section)")
  by_section <- data.frame(
    section = unique(sections),
    checks = as.integer(tapply(log$status != "NOTE", factor(sections, unique(sections)), sum)),
    passed = as.integer(tapply(log$status == "OK", factor(sections, unique(sections)), sum)),
    failed = as.integer(tapply(log$status == "FAIL", factor(sections, unique(sections)), sum)),
    notes = as.integer(tapply(log$status == "NOTE", factor(sections, unique(sections)), sum)),
    stringsAsFactors = FALSE
  )
  cat("\nBy section:\n")
  print(by_section, row.names = FALSE, right = FALSE)

  if (n_fail > 0) {
    cat("\nFailed checks (not correct):\n")
    print(as.data.frame(log[log$status == "FAIL", c("section", "check")]), row.names = FALSE, right = FALSE)
  } else {
    cat("\nAll checks are correct.\n")
  }
  if (n_note > 0) {
    cat("\nNotes (informational, not failures):\n")
    print(as.data.frame(log[log$status == "NOTE", c("section", "check")]), row.names = FALSE, right = FALSE)
  }
  if (isTRUE(details)) {
    cat("\nEvery check:\n")
    print(as.data.frame(log[, c("status", "section", "check")]), row.names = FALSE, right = FALSE)
  }
  invisible(log)
}

#' Check that everything is ready before starting
#'
#' Run this at the start of a script. It starts a fresh list of checks (so
#' the final [check_summary()] only counts this run) and checks that you are
#' ready to design and analyse trials:
#' \itemize{
#'   \item breedtools is loaded, and which version;
#'   \item the optional packages 'lme4' (mixed models with random blocks,
#'     random entries or random genotype-by-location) and 'SpATS' (spatial
#'     models) are installed; a missing one is a note, not a failure, and the
#'     message says what it is needed for;
#'   \item the folder for your files exists (it is created if needed) and
#'     files can be written to it.
#' }
#'
#' @param out_dir Optional folder where the script will write field books
#'   and other files.
#'
#' @return Invisibly, `TRUE` if everything needed is in place.
#'
#' @seealso [check_summary()], [check_design()], [check_result()]
#'
#' @examples
#' check_setup(out_dir = file.path(tempdir(), "my_trials"))
#'
#' @export
check_setup <- function(out_dir = NULL) {
  .check_env$log <- NULL
  .check_env$section <- NULL
  check_section("0. Setup")
  version <- as.character(utils::packageVersion("breedtools"))
  ok <- check_result(sprintf("breedtools %s is loaded", version), TRUE)
  check_result("lme4 is installed (needed for mixed models: random blocks, random entries, random genotype x location)",
               requireNamespace("lme4", quietly = TRUE), note = TRUE)
  check_result("SpATS is installed (needed for spatial models)",
               requireNamespace("SpATS", quietly = TRUE), note = TRUE)
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    probe <- file.path(out_dir, ".breedtools_write_test")
    writable <- isTRUE(tryCatch(file.create(probe), error = function(e) FALSE, warning = function(w) FALSE))
    if (writable) unlink(probe)
    ok <- check_result(sprintf("files can be saved in %s", out_dir), dir.exists(out_dir) && writable) && ok
  }
  missing <- c(if (!requireNamespace("lme4", quietly = TRUE)) "lme4",
               if (!requireNamespace("SpATS", quietly = TRUE)) "SpATS")
  # Later checks are grouped by their own section or kind of design.
  .check_env$section <- NULL
  if (ok) {
    bt_say("Everything is ready to start a new design.\n")
    if (length(missing) > 0) {
      bt_say(sprintf("(Optional: install.packages(%s) to use the analyses that need %s.)\n",
                  paste0("c(", paste0('"', missing, '"', collapse = ", "), ")"),
                  paste(missing, collapse = " and ")))
    }
  } else {
    bt_say("Something is not ready: see the [FAIL] line above.\n")
  }
  invisible(ok)
}

#' Set the section label for the next checks
#'
#' Checks recorded after this call with [check_result()] and
#' [check_design()] are grouped under `name` in [check_summary()].
#'
#' @param name Section label, e.g. `"1. Row-column design"`. `NULL` removes
#'   the label.
#'
#' @return Invisibly, the previous label.
#'
#' @examples
#' check_reset()
#' check_section("arithmetic")
#' check_result("one equals one", 1 == 1)
#' check_summary()
#' check_reset()
#'
#' @export
check_section <- function(name) {
  old <- .check_env$section
  .check_env$section <- if (is.null(name)) NULL else as.character(name)
  if (!is.null(name)) bt_say("\n== ", name, " ==\n", sep = "")
  invisible(old)
}

#' Clear the recorded checks
#'
#' Removes the checks recorded so far with [check_result()] and
#' [check_design()], so that the next [check_summary()] only counts new
#' checks. It does not check or change any design or data; it only empties
#' the list. A message says how many earlier checks were cleared and how many
#' of them had passed.
#'
#' @return Invisibly, a tibble of the checks that were cleared.
#'
#' @examples
#' check_reset()
#' check_result("one equals one", 1 == 1)
#' check_reset()
#'
#' @export
check_reset <- function() {
  old <- check_summary(print = FALSE)
  if (nrow(old) == 0) {
    message("Ready to start: no earlier checks to clear. ",
            "The checks you run from now on will be counted together.")
  } else {
    n_ok <- sum(old$status == "OK")
    n_fail <- sum(old$status == "FAIL")
    n <- n_ok + n_fail
    result <- if (n == 1) {
      if (n_fail == 0) "1 check, which passed" else "1 check, which failed"
    } else if (n_fail == 0) {
      sprintf("%d checks, all passed", n)
    } else {
      sprintf("%d checks: %d passed, %d failed", n, n_ok, n_fail)
    }
    message(sprintf("Earlier checks cleared (%s). Ready to start new checks.", result))
  }
  .check_env$log <- NULL
  .check_env$section <- NULL
  invisible(old)
}

#' Accuracy of adjusted means against known values
#'
#' In a simulation the true genotype values are known. This function measures
#' how well an analysis recovered them: the correlation between the adjusted
#' means and the true values.
#'
#' @param fit An analysis object (from any `analyze_*()` function) or a data
#'   frame with `genotype` and `adjusted_mean` columns.
#' @param truth Named numeric vector of true genotype values; names are the
#'   genotype names.
#' @param genotypes Optional subset of genotypes to use (for example only the
#'   test entries, without the checks).
#' @param method Correlation method: `"pearson"` (default) or `"spearman"`
#'   (rank correlation, i.e. how well the ranking is recovered).
#'
#' @return A single number between -1 and 1.
#'
#' @examples
#' design <- augmented_prep_design(paste0("G", 1:30), c("A", "B"), block_size = 6, seed = 1)
#' trial <- design$design
#' set.seed(1)
#' truth <- setNames(rnorm(30, sd = 3), paste0("G", 1:30))
#' trial$yield <- 50 + truth[trial$all_entries] + rnorm(nrow(trial))
#' fit <- analyze_augmented_prep_design(trial, "yield")
#' prediction_accuracy(fit, truth)
#'
#' @export
prediction_accuracy <- function(fit, truth, genotypes = NULL, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  means <- if (is.list(fit) && !is.data.frame(fit) && !is.null(fit$adjusted_means)) fit$adjusted_means else fit
  if (!is.data.frame(means) || !all(c("genotype", "adjusted_mean") %in% names(means))) {
    stop("`fit` must be an analysis object or a data frame with `genotype` and `adjusted_mean`.", call. = FALSE)
  }
  if (is.null(names(truth))) stop("`truth` must be a named vector (names = genotypes).", call. = FALSE)
  g <- as.character(means$genotype)
  if (!is.null(genotypes)) g <- intersect(g, as.character(genotypes))
  g <- intersect(g, names(truth))
  if (length(g) < 3) stop("Fewer than 3 genotypes in common between `fit` and `truth`.", call. = FALSE)
  est <- means$adjusted_mean[match(g, means$genotype)]
  stats::cor(est, truth[g], use = "complete.obs", method = method)
}

#' Standard checks of a design
#'
#' Runs the structural checks that every design of its kind should pass and
#' prints one `[OK]` / `[FAIL]` line per check (`[NOTE]` for informational
#' ones). The results are recorded for [check_summary()].
#'
#' \describe{
#'   \item{Augmented row-column designs}{every plot used once, every entry
#'     planted once, every field row and column has a check, the requested
#'     number of check plots in each block, positive error df, all row,
#'     column and block effects estimable from the checks, and no check plots
#'     next to each other inside a block (a note for neighbours across block
#'     boundaries).}
#'   \item{Augmented p-rep designs}{every entry at every location, no entry
#'     twice in a block, block sizes within a location differ by at most one,
#'     every location connected, the pairwise concurrence (a note when pairs
#'     share more than one block), and the replicate structure of resolvable
#'     designs.}
#'   \item{Resolvable row-column designs}{every treatment once in every
#'     replicate, every plot used once, the spatial efficiency below its upper
#'     bound, and latinization when requested.}
#' }
#'
#' @param design A design object from [augmented_row_column_design()],
#'   [augmented_prep_design()] or [resolvable_row_column_design()].
#' @param section Optional label for the recorded checks. Defaults to the
#'   label set with [check_section()], or else the kind of design.
#'
#' @return Invisibly, `TRUE` if all (non-note) checks passed.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' check_reset()
#' check_design(rc)
#' check_summary()
#' check_reset()
#'
#' @export
check_design <- function(design, section = NULL) {
  UseMethod("check_design")
}

#' @export
check_design.default <- function(design, section = NULL) {
  stop("`design` must come from augmented_row_column_design(), augmented_prep_design() ",
       "or resolvable_row_column_design().", call. = FALSE)
}

#' @export
check_design.augmented_row_column_design <- function(design, section = NULL) {
  sec <- if (!is.null(section)) section else if (!is.null(.check_env$section)) .check_env$section else "augmented row-column design"
  d <- design$design
  inp <- design$input
  ck <- d[d$is_control == 1, ]
  en <- d[d$is_control == 0, ]
  R <- max(d$row)
  C <- max(d$col)
  bt_say("Checking the augmented row-column design\n")
  ok <- c(
    check_result("every plot used exactly once", !anyDuplicated(paste(d$row, d$col)) && nrow(d) == R * C, sec),
    check_result("every entry planted exactly once",
                 !anyDuplicated(en$all_entries) && setequal(en$all_entries, inp$treatments), sec),
    check_result("every field row has a check", setequal(unique(ck$row), seq_len(R)), sec),
    check_result("every field column has a check", setequal(unique(ck$col), seq_len(C)), sec),
    check_result("check plots per block as requested",
                 all(as.vector(table(factor(ck$block, levels = seq_along(inp$checks_per_block)))) ==
                       as.vector(t(inp$checks_per_block))), sec),
    if (design$summary$error_df == 0) {
      # Saturated with fixed effects but analyzable with random effects.
      check_result("error df > 0 (0): analyze with random_blocks = TRUE or entries = \"random\"",
                   FALSE, sec, note = TRUE) || TRUE
    } else {
      check_result(sprintf("error df > 0 (%d)", design$summary$error_df), design$summary$error_df > 0, sec)
    },
    check_result("row, column and block effects estimable from the checks", design$summary$d_criterion > 0, sec)
  )
  # Blocks with at least as many check plots as checks hold every check;
  # incomplete blocks hold each check at most once.
  pb <- as.vector(t(inp$checks_per_block))
  vc <- length(inp$controls)
  tab <- table(factor(ck$block, levels = seq_along(pb)), factor(ck$all_entries, levels = inp$controls))
  balanced <- all(vapply(seq_along(pb), function(b) {
    counts <- tab[b, ]
    all(counts >= floor(pb[b] / vc)) && all(counts <= ceiling(pb[b] / vc))
  }, logical(1)))
  ok <- c(ok, check_result("checks spread evenly in each block (every check in full blocks)", balanced, sec))
  # Rows (columns) within a row (column) group get floor or ceiling of the
  # group's check plots divided by its rows (columns).
  even_lines <- function(line, group, totals, sizes, n) {
    counts <- tabulate(line, n)
    lo <- floor(totals / sizes)[group]
    hi <- ceiling(totals / sizes)[group]
    all(counts >= lo & counts <= hi)
  }
  row_group <- rep(seq_along(inp$rows_per_rowgroup), inp$rows_per_rowgroup)
  col_group <- rep(seq_along(inp$cols_per_colgroup), inp$cols_per_colgroup)
  ok <- c(ok, check_result(
    "even numbers of check plots in the rows and columns of each group",
    even_lines(ck$row, row_group, rowSums(inp$checks_per_block), inp$rows_per_rowgroup, R) &&
      even_lines(ck$col, col_group, colSums(inp$checks_per_block), inp$cols_per_colgroup, C), sec))
  covered <- function(line, n) {
    tab <- table(factor(line, levels = seq_len(n)), factor(ck$all_entries, levels = inp$controls))
    all(tab >= floor(rowSums(tab) / vc))
  }
  ok <- c(ok, check_result("every check in each row and column with room for all checks",
                           covered(ck$row, R) && covered(ck$col, C), sec))
  # Neighbouring checks inside a block can be avoided unless a block has more
  # check plots than half its plots.
  avoidable <- all(inp$checks_per_block <= ceiling(outer(inp$rows_per_rowgroup, inp$cols_per_colgroup) / 2))
  n_inside <- design$summary$adjacent_checks_in_blocks
  ok <- c(ok, check_result(
    sprintf("no checks next to each other inside a block (%d pairs)", n_inside),
    n_inside == 0, sec, note = !avoidable
  ))
  check_result(sprintf("no checks next to each other across block boundaries (%d pairs)",
                       design$summary$adjacent_checks_across_blocks),
               design$summary$adjacent_checks_across_blocks == 0, sec, note = TRUE)
  if (identical(inp$method, "latin")) {
    ok <- c(ok, check_result("no check twice in a field row or column",
                             !any(table(ck$all_entries, ck$row) > 1) && !any(table(ck$all_entries, ck$col) > 1), sec))
  }
  invisible(all(ok))
}

#' @export
check_design.augmented_prep_design <- function(design, section = NULL) {
  sec <- if (!is.null(section)) section else if (!is.null(.check_env$section)) .check_env$section else "augmented p-rep design"
  d <- design$design
  inp <- design$input
  s <- design$summary
  bt_say("Checking the augmented p-rep design\n")
  sizes <- table(paste(d$location, d$block))
  by_loc <- split(as.integer(sizes), sub(" .*", "", names(sizes)))
  ok <- c(
    check_result("every entry planted at every location",
                 all(table(factor(d$all_entries, inp$treatments), factor(d$location, inp$locations)) >= 1), sec),
    check_result("no entry twice in the same block", s$same_entry_same_block == 0, sec),
    check_result("block sizes within a location differ by at most one",
                 all(vapply(by_loc, function(x) max(x) - min(x) <= 1, logical(1))), sec),
    check_result("every location connected (each can be analysed alone)", s$n_disconnected_locations == 0, sec)
  )
  check_result(sprintf("no pair of entries shares a block more than once (max %d)", s$max_pairwise_concurrence),
               s$max_pairwise_concurrence <= 1, sec, note = TRUE)
  if (isTRUE(inp$resolvable)) {
    pp <- d[d$entry_type == "p_rep_entry_at_this_location", ]
    ok <- c(ok, check_result("resolvable: each repeated entry once per replicate",
                             all(table(paste(pp$location, pp$all_entries), pp$replicate) <= 1), sec))
  }
  invisible(all(ok))
}

#' @export
check_design.resolvable_row_column_design <- function(design, section = NULL) {
  sec <- if (!is.null(section)) section else if (!is.null(.check_env$section)) .check_env$section else "resolvable row-column design"
  d <- design$design
  s <- design$summary
  real <- d[d$is_filler == 0, ]
  bt_say("Checking the resolvable row-column design\n")
  ok <- c(
    check_result("every treatment once in every replicate",
                 all(table(real$all_entries, real$replicate) == 1), sec),
    check_result("every plot used exactly once", !anyDuplicated(paste(d$row, d$col)), sec),
    check_result(sprintf("spatial efficiency %.4f below its upper bound %.4f", s$efficiency_spatial, s$upper_bound_spatial),
                 s$efficiency_spatial <= s$upper_bound_spatial + 1e-8, sec)
  )
  if (isTRUE(s$latinized)) {
    ok <- c(ok, check_result("latinized: no treatment twice in a long line", isTRUE(s$latin_ok), sec))
  }
  invisible(all(ok))
}
