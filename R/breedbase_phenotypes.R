# Phenotypic data to Breedbase: a template from the design and the upload of
# observations through BrAPI v2 (POST /observations).
#
# Breedbase stores an observation from three ids: the plot
# (observationUnitDbId), the trait (observationVariableDbId) and the value.
# The ids are looked up from the plot and trait names given by the user.

#' Phenotyping template with the Breedbase plot names of a design
#'
#' Makes a table with one row per plot, the plot names that
#' [save_breedbase()] gives (or gave) to the plots, the design columns
#' (accession, block, rep, row, column) and one empty column per trait. Fill
#' it in (in R or in a spreadsheet) and upload it with
#' [save_breedbase_phenotypes()].
#'
#' @param design A design from [augmented_row_column_design()],
#'   [augmented_prep_design()], [resolvable_row_column_design()],
#'   [as_breedbase_design()] or [prepare_breedbase_trial()].
#' @param trial_name The trial name used when saving the design. For p-rep
#'   designs saved as one trial per location, the location is added as in
#'   [save_breedbase()]. Can be omitted for a prepared Breedbase trial or a
#'   design with plot names (see [add_plot_names()]).
#' @param traits Trait names as in Breedbase, e.g. `"fresh root yield"` or
#'   `"fresh root yield|CO_334:0000013"`. Can be empty (no trait columns).
#' @param plot_name_pattern,plot_names The same plot naming used in
#'   [save_breedbase()]. By default the plot names stored in the design are
#'   used, when it has them.
#' @param location The Breedbase location of the trial, needed when the
#'   plot names are made here with the `{location}` field (the default
#'   pattern) and the design has no stored plot names.
#'
#' @return A tibble: `plot_name`, `trial_name`, `accession`, `plot_number`,
#'   `block`, `rep`, `row`, `col`, `is_check` and one column per trait.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' sheet <- phenotype_template(rc, "AYT_2026", location = "Ibadan",
#'                             traits = c("fresh root yield", "dry matter content percentage"))
#' head(sheet)
#'
#' @export
phenotype_template <- function(design, trial_name, traits = character(0),
                               plot_name_pattern = NULL, plot_names = NULL, location = NULL) {
  if (missing(trial_name)) trial_name <- NULL
  out <- bb_design_plot_table(design, trial_name, plot_name_pattern, plot_names, location = location)
  out$.design_location <- NULL
  out$.design_plots <- NULL
  for (t in unique(as.character(traits))) out[[t]] <- NA_real_
  out
}

# Breedbase plot names of a design with the design's own plot number (and
# location for p-rep designs), to match analysis results.
bb_design_plot_table <- function(design, trial_name, plot_name_pattern = NULL,
                                 plot_names = NULL, location = NULL) {
  if (is.null(trial_name) && inherits(design, "breedbase_trial")) trial_name <- design$trial$trial_name
  if (is.null(trial_name)) trial_name <- design$plot_name_trial
  studies <- bb_design_to_studies(design, trial_name, year = 2000, location = bb_template_location(design, location),
                                  description = "template", trial_type = "phenotyping_trial", design_type = NULL,
                                  plot_names = plot_names, plot_name_pattern = plot_name_pattern)
  rows <- lapply(studies, function(s) {
    level <- function(p, name) {
      rel <- Filter(function(r) identical(r$levelName, name), p$observationUnitPosition$observationLevelRelationships)
      if (length(rel)) as.integer(rel[[1]]$levelCode) else NA_integer_
    }
    tibble::tibble(
      plot_name = vapply(s$plots, `[[`, character(1), "observationUnitName"),
      trial_name = s$study$studyName,
      accession = vapply(s$plots, `[[`, character(1), "germplasmName"),
      plot_number = vapply(s$plots, function(p) as.integer(p$observationUnitPosition$observationLevel$levelCode), integer(1)),
      block = vapply(s$plots, level, integer(1), name = "block"),
      rep = vapply(s$plots, level, integer(1), name = "rep"),
      row = vapply(s$plots, function(p) as.integer(p$observationUnitPosition$positionCoordinateY %||% NA), integer(1)),
      col = vapply(s$plots, function(p) as.integer(p$observationUnitPosition$positionCoordinateX %||% NA), integer(1)),
      is_check = vapply(s$plots, function(p) identical(p$observationUnitPosition$entryType, "CHECK"), logical(1))
    )
  })
  out <- dplyr::bind_rows(rows)
  # The design's own plot number (and location) of each row, in the same
  # order save_breedbase() uses: per location, fillers left out, by plot.
  d <- as.data.frame(design$design, stringsAsFactors = FALSE)
  if ("is_filler" %in% names(d)) d <- d[d$is_filler != 1, , drop = FALSE]
  groups <- if (inherits(design, "augmented_prep_design")) {
    lapply(design$input$locations, function(l) d[d$location == l, , drop = FALSE])
  } else {
    list(d)
  }
  ordered <- dplyr::bind_rows(lapply(groups, function(x) x[order(x$plots), , drop = FALSE]))
  out$.design_location <- if ("location" %in% names(ordered)) as.character(ordered$location) else NA_character_
  out$.design_plots <- as.integer(ordered$plots)
  out
}

# The template only needs the location names of p-rep designs (to name the
# trials per location); other designs need a placeholder.
bb_template_location <- function(design, location = NULL) {
  if (!is.null(location)) return(location)
  if (inherits(design, "augmented_prep_design")) return(NULL)
  if (inherits(design, "breedbase_trial") && !is.null(design$trial$location)) return(design$trial$location)
  bb_no_location
}

# Placeholder location for templates of designs saved without a known
# location: a {location} field in the plot name pattern asks for it.
bb_no_location <- ".breedtools_no_location."


#' Save phenotypic data to Breedbase
#'
#' Uploads trait values of plots that are already in Breedbase (for example
#' saved with [save_breedbase()]) through BrAPI. Give the plot name and the
#' trait name with each value; the plots and traits are looked up in
#' Breedbase first, and every problem is listed before anything is written.
#'
#' The data can be wide, one row per plot and one column per trait (as made
#' by [phenotype_template()]):
#'
#' | plot_name                 | fresh root yield | dry matter content percentage |
#' |---------------------------|------------------|-------------------------------|
#' | AYT_2026_Ibadan_G31_1_1   | 25.4             | 31.2                          |
#'
#' or long, one row per value, with `trait` and `value` columns
#' (`trait = "trait", value = "value"`).
#'
#' Trait names are the Breedbase names, with or without the ontology id:
#' `"fresh root yield"` or `"fresh root yield|CO_334:0000013"` (the column
#' names of a Breedbase phenotype download). Values are checked by Breedbase
#' against the trait (e.g. numbers within the allowed range).
#'
#' @section Saving again (values already in Breedbase):
#' Before writing, the values already stored for these plots and traits are
#' read, and every value is sorted into new, unchanged (same value, not sent
#' again) or changed:
#' \describe{
#'   \item{`overwrite = FALSE` (default)}{If changed values exist, a message is
#'     printed and the function returns without an error; nothing is written.}
#'   \item{`overwrite = TRUE`}{Changed values replace the stored ones.}
#' }
#' Empty cells (`NA` or `""`) are never sent. When Breedbase has a value
#' for an empty cell (e.g. an outlier set to `NA` with [remove_outliers()]):
#' \describe{
#'   \item{`remove_empty = FALSE` (default)}{The stored value is kept (and
#'     reported).}
#'   \item{`remove_empty = TRUE`}{The stored value is deleted. Like the
#'     "remove values" option of Breedbase's upload page, this needs
#'     `overwrite = TRUE`, a curator account and access to the web pages
#'     (the site password, see [breedbase_login()]).}
#' }
#'
#' @param data A data frame with the plot names and the values.
#' @param con A connection from [breedbase_login()]; `NULL` with
#'   `dry_run = TRUE` returns the observations without checking them.
#' @param plot Name of the column with the plot names.
#' @param traits For wide data: the trait columns. By default every column
#'   except `plot`, `timestamp` and the design columns of
#'   [phenotype_template()] (`trial_name`, `accession`, `plot_number`,
#'   `block`, `rep`, `row`, `col`, `is_check`).
#' @param trait,value For long data: the names of the trait and value columns.
#' @param timestamp Optional name of a column with the date and time of each
#'   measurement (e.g. `"2026-08-15 10:30:00"`). Breedbase uses it to tell
#'   repeated measurements apart.
#' @param overwrite Logical. Replace values that Breedbase already has with
#'   different ones (see "Saving again").
#' @param remove_empty Logical. Delete stored values whose cell is empty
#'   (`NA` or `""`) in `data`; needs `overwrite = TRUE` (see "Saving again").
#' @param dry_run Logical. Check everything but write nothing.
#' @param confirm Logical. Ask before writing (default in interactive
#'   sessions).
#' @param refresh Logical. After saving, refresh the phenotype cache of
#'   Breedbase so the values show on the trial pages (see
#'   [breedbase_refresh()]).
#'
#' @inheritSection breedbase_check Slow servers
#'
#' @return Invisibly, a tibble with one row per trait: `trait`,
#'   `observation_variable_db_id`, `values_sent`, `values_found` (read back
#'   from Breedbase after saving), and how the values were handled:
#'   `new`, `replaced`, `unchanged` and `removed`. Returns invisible `NULL`
#'   when saving is declined, including changed values with `overwrite = FALSE`.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' sheet <- phenotype_template(rc, "AYT_2026", location = "Ibadan", traits = "fresh root yield")
#' sheet$`fresh root yield` <- round(rnorm(nrow(sheet), 25, 4), 1)
#'
#' # What would be sent (no connection needed):
#' obs <- save_breedbase_phenotypes(sheet, dry_run = TRUE)
#' head(obs)
#'
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' save_breedbase_phenotypes(sheet, con, dry_run = TRUE)   # check first
#' save_breedbase_phenotypes(sheet, con)
#'
#' # Saving a corrected sheet again: replace changed values, and delete the
#' # stored values of cells now empty (e.g. outliers set to NA).
#' save_breedbase_phenotypes(sheet, con, overwrite = TRUE, remove_empty = TRUE)
#'
#' # Long data: one row per value
#' long <- data.frame(plot_name = sheet$plot_name, trait = "fresh root yield",
#'                    value = sheet$`fresh root yield`)
#' save_breedbase_phenotypes(long, con, trait = "trait", value = "value")
#' }
#'
#' @export
save_breedbase_phenotypes <- function(data, con = NULL, plot = "plot_name", traits = NULL,
                                      trait = NULL, value = NULL, timestamp = NULL,
                                      overwrite = FALSE, remove_empty = FALSE,
                                      dry_run = FALSE, confirm = interactive(), refresh = TRUE) {
  bb_save_phenotypes(data, con, plot, traits, trait, value, timestamp, dry_run, confirm, refresh = refresh,
                     overwrite = overwrite, remove_empty = remove_empty)
}

# `trait_ids`: trait ids already known, named by trait (e.g. composed traits,
# which the BrAPI variable search may not list).
bb_save_phenotypes <- function(data, con, plot = "plot_name", traits = NULL, trait = NULL, value = NULL,
                               timestamp = NULL, dry_run = FALSE, confirm = interactive(), trait_ids = NULL,
                               refresh = TRUE, overwrite = FALSE, remove_empty = FALSE) {
  if (isTRUE(remove_empty) && !isTRUE(overwrite)) {
    stop("`remove_empty = TRUE` needs `overwrite = TRUE` (as on Breedbase's upload page): deleting a stored ",
         "value replaces it.", call. = FALSE)
  }
  all_obs <- bb_phenotypes_long(data, plot, traits, trait, value, timestamp, keep_empty = TRUE)
  is_empty <- is.na(all_obs$value) | !nzchar(all_obs$value)
  obs <- all_obs[!is_empty, , drop = FALSE]
  empty <- all_obs[is_empty, , drop = FALSE]
  if (!nrow(obs) && !isTRUE(remove_empty)) stop("There are no values to save (all are empty).", call. = FALSE)

  if (is.null(con)) {
    if (!isTRUE(dry_run)) stop("`con` is needed to save; use breedbase_login() first, or `dry_run = TRUE`.", call. = FALSE)
    message("Dry run without a connection: returning the observations that would be sent.")
    return(obs)
  }
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)

  cat("Checking phenotypes before saving to ", con$url, "\n", sep = "")
  cat("  [OK]   ", nrow(obs), " values of ", length(unique(all_obs$trait)), " trait",
      if (length(unique(all_obs$trait)) > 1) "s" else "", " on ", length(unique(obs$plot_name)), " plots",
      if (nrow(empty)) paste0(" (", nrow(empty), " empty cells)"), "\n", sep = "")
  web_ok <- TRUE
  if (isTRUE(remove_empty) && nrow(empty)) {
    web_ok <- !identical(bb_web_access(con), "password")
    if (!web_ok) {
      stop("`remove_empty = TRUE` deletes values through Breedbase's upload page, which this site closes with ",
           "a web password; log in with site_username = \"...\" (see ?breedbase_login).", call. = FALSE)
    }
  }

  # Plots and traits are both checked before stopping, so every problem is
  # listed at once.
  plot_names <- unique(all_obs$plot_name)
  plot_ids <- bb_search_names(con, "observationunits", "observationUnitNames", "observationUnitName", plot_names,
                              batch = getOption("breedtools.plot_batch", 200), label = "plots",
                              id_field = "observationUnitDbId")
  known_ids <- trait_ids
  found_vars <- bb_find_variables(con, setdiff(unique(all_obs$trait), names(known_ids)))
  trait_full <- c(attr(found_vars, "full_name"), stats::setNames(names(known_ids), names(known_ids)))
  trait_ids <- c(found_vars, known_ids)[unique(all_obs$trait)]
  problems <- character(0)
  missing_plots <- setdiff(plot_names, names(plot_ids))
  if (length(missing_plots)) {
    problems <- c(problems, paste0(
      length(missing_plots), " of ", length(plot_names), " plots are not in Breedbase: ",
      paste(utils::head(missing_plots, 10), collapse = ", "), if (length(missing_plots) > 10) ", ..." else "",
      ". Save the trial first (save_breedbase()) and use its plot names (phenotype_template())."))
  } else {
    cat("  [OK]   all ", length(plot_names), " plots are in Breedbase\n", sep = "")
  }
  missing_traits <- names(trait_ids)[is.na(trait_ids)]
  if (length(missing_traits)) {
    problems <- c(problems, paste0(
      "Trait", if (length(missing_traits) > 1) "s" else "", " not found in Breedbase: ",
      paste0("'", missing_traits, "'", collapse = ", "),
      ". Use the trait name as shown in Breedbase, e.g. \"fresh root yield\" or ",
      "\"fresh root yield|CO_334:0000013\"."))
  } else {
    for (t in names(trait_ids)) cat("  [OK]   trait '", t, "' (id ", trait_ids[[t]], ")\n", sep = "")
  }
  if (length(problems)) stop(paste(problems, collapse = "\n"), call. = FALSE)

  obs$observation_unit_db_id <- unname(plot_ids[obs$plot_name])
  obs$observation_variable_db_id <- unname(trait_ids[obs$trait])
  empty$observation_unit_db_id <- unname(plot_ids[empty$plot_name])
  empty$observation_variable_db_id <- unname(trait_ids[empty$trait])

  # Values already stored for these plots and traits: new, unchanged or
  # changed; and stored values whose cell is now empty.
  stored <- bb_stored_values(con, rbind(obs, empty))
  key_obs <- paste(obs$observation_unit_db_id, obs$observation_variable_db_id, sep = "\r")
  hit <- match(key_obs, stored$key)
  # Repeated measurements (with a timestamp) are always new values.
  hit[!is.na(obs$timestamp)] <- NA
  obs$observation_db_id <- stored$observation_db_id[hit]
  obs$stored_value <- stored$value[hit]
  obs$status <- ifelse(is.na(hit), "new", ifelse(bb_same_value(obs$value, obs$stored_value), "unchanged", "changed"))
  empty_hit <- match(paste(empty$observation_unit_db_id, empty$observation_variable_db_id, sep = "\r"), stored$key)
  to_remove <- empty[!is.na(empty_hit), , drop = FALSE]
  to_remove$stored_value <- stored$value[empty_hit[!is.na(empty_hit)]]

  n_new <- sum(obs$status == "new")
  n_same <- sum(obs$status == "unchanged")
  n_changed <- sum(obs$status == "changed")
  cat("  [OK]   ", n_new, " new values", if (n_same) paste0(", ", n_same, " already stored with the same value (not sent again)"),
      "\n", sep = "")
  if (n_changed) {
    if (!isTRUE(overwrite)) {
      cat("You have decided to not replace. Phenotypes are already there.\n")
      return(invisible(NULL))
    }
    ch <- obs[obs$status == "changed", , drop = FALSE]
    show <- utils::head(paste0(ch$plot_name, " / ", ch$trait, ": ", ch$stored_value, " -> ", ch$value), 10)
    cat("  [OK]   ", n_changed, " stored values will be replaced (overwrite = TRUE), e.g. ", show[1], "\n", sep = "")
  }
  if (nrow(to_remove)) {
    if (isTRUE(remove_empty)) {
      cat("  [OK]   ", nrow(to_remove), " stored values will be deleted (empty in `data`, remove_empty = TRUE), e.g. ",
          to_remove$plot_name[1], " / ", to_remove$trait[1], " = ", to_remove$stored_value[1], "\n", sep = "")
    } else {
      cat("  [NOTE] ", nrow(to_remove), " cells are empty in `data` but have a value in Breedbase; those values are ",
          "kept (remove_empty = TRUE deletes them)\n", sep = "")
      to_remove <- to_remove[0, , drop = FALSE]
    }
  }
  send_new <- obs[obs$status == "new", , drop = FALSE]
  send_put <- obs[obs$status == "changed", , drop = FALSE]
  as_record <- function(o, i) {
    r <- list(observationUnitDbId = o$observation_unit_db_id[i],
              observationVariableDbId = o$observation_variable_db_id[i],
              value = o$value[i])
    if (!is.na(o$timestamp[i])) r$observationTimeStamp <- o$timestamp[i]
    r
  }
  records <- lapply(seq_len(nrow(send_new)), function(i) as_record(send_new, i))

  plan <- sprintf("%d new value%s%s%s for %d trait%s on %d plots at %s", nrow(send_new),
                  if (nrow(send_new) == 1) "" else "s",
                  if (nrow(send_put)) sprintf(", %d replaced", nrow(send_put)) else "",
                  if (nrow(to_remove)) sprintf(", %d deleted", nrow(to_remove)) else "",
                  length(trait_ids), if (length(trait_ids) > 1) "s" else "", length(plot_names), con$url)
  summary_of <- function(found) {
    out <- bb_phenotype_summary(obs, found)
    count <- function(x, t) sum(x$trait == t)
    out$new <- vapply(out$trait, function(t) count(send_new, t), integer(1), USE.NAMES = FALSE)
    out$replaced <- vapply(out$trait, function(t) count(send_put, t), integer(1), USE.NAMES = FALSE)
    out$unchanged <- vapply(out$trait, function(t) count(obs[obs$status == "unchanged", ], t), integer(1),
                            USE.NAMES = FALSE)
    out$removed <- vapply(out$trait, function(t) count(to_remove, t), integer(1), USE.NAMES = FALSE)
    out
  }
  if (!nrow(send_new) && !nrow(send_put) && !nrow(to_remove)) {
    cat("Nothing to write: Breedbase already has these values.\n")
    return(invisible(summary_of(found = NULL)))
  }
  if (isTRUE(dry_run)) {
    cat("Dry run: everything is ready to save ", plan, ". Nothing was written.\n", sep = "")
    return(invisible(summary_of(found = NULL)))
  }
  if (isTRUE(confirm)) {
    answer <- readline(sprintf("Save %s? Type 'yes' to continue: ", plan))
    if (!identical(tolower(trimws(answer)), "yes")) {
      message("Nothing was saved.")
      return(invisible(NULL))
    }
  }

  batch <- getOption("breedtools.observation_batch", 500)
  groups <- split(seq_along(records), ceiling(seq_along(records) / batch))
  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    bb_write_checked(
      what = sprintf("save values %d-%d of %d", min(idx), max(idx), length(records)),
      write = function() {
        bb_post(con, "observations", records[idx], what = "save the phenotypes")
        "done"
      },
      find = function() {
        n <- bb_count_observations(con, send_new[idx, , drop = FALSE])
        if (n >= length(idx)) "done" else NA_character_
      }
    )
    if (length(groups) > 1) {
      cat(sprintf("  ...    batch %d of %d saved: %d/%d values\n", g, length(groups), max(idx), length(records)))
    }
  }
  # Changed values: BrAPI PUT /observations, keyed by the stored observation.
  if (nrow(send_put)) {
    for (idx in split(seq_len(nrow(send_put)), ceiling(seq_len(nrow(send_put)) / batch))) {
      payload <- stats::setNames(lapply(idx, function(i) as_record(send_put, i)), send_put$observation_db_id[idx])
      bb_put(con, "observations", payload, what = "replace stored phenotypes")
    }
    cat("  [OK]   ", nrow(send_put), " stored values replaced\n", sep = "")
  }
  # Empty cells with a stored value: Breedbase's upload with "overwrite" and
  # "remove values" (BrAPI cannot delete observations). The rows of these
  # plots are sent with all their values, so only the empty cells change.
  if (nrow(to_remove)) {
    bb_remove_values(con, all_obs[all_obs$plot_name %in% to_remove$plot_name, , drop = FALSE], trait_full)
    left <- bb_stored_values(con, to_remove)
    n_left <- sum(paste(to_remove$observation_unit_db_id, to_remove$observation_variable_db_id, sep = "\r") %in% left$key)
    cat("  ", if (n_left) "[FAIL]" else "[OK]  ", " ", nrow(to_remove) - n_left, " of ", nrow(to_remove),
        " stored values deleted\n", sep = "")
  }

  found <- bb_observations_found(con, obs)
  out <- summary_of(found)
  for (i in seq_len(nrow(out))) {
    status <- if (out$values_found[i] >= out$values_sent[i]) "[OK]  " else "[FAIL]"
    cat(sprintf("  %s trait '%s': %d of %d values stored (%d new, %d replaced, %d unchanged, %d deleted)\n",
                status, out$trait[i], out$values_found[i], out$values_sent[i], out$new[i], out$replaced[i],
                out$unchanged[i], out$removed[i]))
  }
  # Breedbase's BrAPI upload does not refresh the phenotype cache that the
  # trial pages and the search wizard read; refresh it so the values show.
  if (isTRUE(refresh)) breedbase_refresh(con, quiet = FALSE)
  invisible(out)
}

#' Refresh the Breedbase phenotype cache
#'
#' Values saved through BrAPI ([save_breedbase_phenotypes()],
#' [save_breedbase_first_stage()]) are in the database at once (the BrAPI
#' search finds them), but the trial pages and the search wizard of
#' Breedbase read a cache ("materialized views") that the BrAPI upload does
#' not refresh. This asks Breedbase to refresh it, like the "Update" button
#' of the search wizard. The save functions call it after saving. If it
#' cannot be reached, the values are still saved: they appear after the next
#' refresh of the site (usually overnight), or when you click "Update" on the
#' Search Wizard page.
#'
#' @param con A connection from [breedbase_login()].
#' @param views Which cache to refresh: `"phenotypes"` (default) or
#'   `"fullview"`.
#' @param quiet Logical. If `FALSE`, print what happened.
#'
#' @return Invisibly, `TRUE` if Breedbase accepted the refresh, else `FALSE`.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' breedbase_refresh(con)
#' }
#'
#' @export
breedbase_refresh <- function(con, views = c("phenotypes", "fullview"), quiet = FALSE) {
  views <- match.arg(views)
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  res <- tryCatch(bb_ajax_get(con, "ajax/breeder/refresh", what = "refresh the phenotype cache",
                              query = list(matviews = views)),
                  breedbase_error = function(e) e)
  ok <- !inherits(res, "breedbase_error")
  if (!isTRUE(quiet)) {
    if (ok) {
      cat("  [OK]   Breedbase is refreshing its phenotype cache: the values appear on the trial pages",
          "and in the search wizard in a few minutes\n")
    } else {
      cat("  [NOTE] The values are saved, but the phenotype cache could not be refreshed (",
          conditionMessage(res), ").\n         They appear on the trial pages after the next refresh of the site, ",
          "or click \"Update\" on the Search Wizard page.\n", sep = "")
    }
  }
  invisible(ok)
}

# Data (wide or long) to one row per value: plot_name, trait, value, timestamp.
bb_phenotypes_long <- function(data, plot, traits, trait, value, timestamp, keep_empty = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c(plot, trait, value, timestamp)
  absent <- setdiff(need, names(data))
  if (length(absent)) {
    stop("Column", if (length(absent) > 1) "s" else "", " not found in `data`: ",
         paste0("'", absent, "'", collapse = ", "), ". Columns are: ", paste(names(data), collapse = ", "), ".",
         call. = FALSE)
  }
  ts <- if (is.null(timestamp)) rep(NA_character_, nrow(data)) else bb_format_time(data[[timestamp]])
  if (!is.null(trait) || !is.null(value)) {
    if (is.null(trait) || is.null(value)) stop("Long data needs both `trait` and `value` columns.", call. = FALSE)
    long <- data.frame(plot_name = as.character(data[[plot]]), trait = as.character(data[[trait]]),
                       value = bb_value_text(data[[value]]), timestamp = ts, stringsAsFactors = FALSE)
  } else {
    design_cols <- c("trial_name", "accession", "plot_number", "block", "rep", "row", "col", "is_check")
    if (is.null(traits)) traits <- setdiff(names(data), c(plot, timestamp, design_cols))
    absent <- setdiff(traits, names(data))
    if (length(absent)) stop("Trait columns not found in `data`: ", paste(absent, collapse = ", "), call. = FALSE)
    if (!length(traits)) stop("No trait columns in `data`. Give them with `traits = c(...)`.", call. = FALSE)
    long <- do.call(rbind, lapply(traits, function(t) {
      data.frame(plot_name = as.character(data[[plot]]), trait = t, value = bb_value_text(data[[t]]),
                 timestamp = ts, stringsAsFactors = FALSE)
    }))
  }
  long$plot_name <- trimws(long$plot_name)
  long$trait <- trimws(long$trait)
  if (!isTRUE(keep_empty)) {
    long <- long[!is.na(long$value) & nzchar(long$value), , drop = FALSE]
    if (!nrow(long)) stop("There are no values to save (all are empty).", call. = FALSE)
  }
  if (anyNA(long$plot_name) || any(!nzchar(long$plot_name))) stop("Some rows have no plot name.", call. = FALSE)
  key <- paste(long$plot_name, long$trait, long$timestamp, sep = "\r")
  if (anyDuplicated(key)) {
    d <- long[duplicated(key), , drop = FALSE][1, ]
    stop("Plot '", d$plot_name, "' has more than one value of '", d$trait, "'",
         if (is.na(d$timestamp)) "; give a `timestamp` column to save repeated measurements" else " at the same time",
         ".", call. = FALSE)
  }
  rownames(long) <- NULL
  tibble::as_tibble(long)
}

bb_value_text <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  out <- as.character(x)
  out[is.na(x)] <- NA_character_
  trimws(out)
}

bb_format_time <- function(x) {
  if (inherits(x, "POSIXt")) return(format(x, "%Y-%m-%d %H:%M:%S"))
  if (inherits(x, "Date")) return(format(x, "%Y-%m-%d 00:00:00"))
  out <- trimws(as.character(x))
  out[out == ""] <- NA_character_
  out
}

# Trait ids by name. Breedbase searches the plain trait name and shows it as
# "name|ONTOLOGY:id"; both forms are accepted.
bb_find_variables <- function(con, traits) {
  if (!length(traits)) {
    out <- stats::setNames(character(0), character(0))
    attr(out, "full_name") <- out
    return(out)
  }
  base <- trimws(sub("\\|.*$", "", traits))
  data <- bb_search(con, "variables", list(observationVariableNames = as.list(unique(base)),
                                           pageSize = 10 * length(unique(base))))
  full <- vapply(data, function(x) as.character(x$observationVariableName %||% NA_character_), character(1))
  ids <- vapply(data, function(x) as.character(x$observationVariableDbId %||% NA_character_), character(1))
  short <- trimws(sub("\\|.*$", "", full))
  hits <- lapply(seq_along(traits), function(i) {
    if (grepl("|", traits[i], fixed = TRUE)) which(full == traits[i]) else which(short == base[i])
  })
  pick <- function(v) vapply(hits, function(h) if (length(h) == 1) v[h] else NA_character_, character(1))
  out <- stats::setNames(pick(ids), traits)
  # The full name "trait|ONTOLOGY:id" (what Breedbase analyses need).
  attr(out, "full_name") <- stats::setNames(pick(full), traits)
  out
}

# Stored observations of these plots and traits, searched a few plots at a
# time so each answer stays small.
bb_search_observations <- function(con, obs) {
  units <- unique(obs$observation_unit_db_id)
  vars <- unique(obs$observation_variable_db_id)
  per_search <- max(1L, floor(getOption("breedtools.observation_batch", 500) / length(vars)))
  out <- list()
  for (chunk in split(units, ceiling(seq_along(units) / per_search))) {
    out <- c(out, bb_search(con, "observations", list(
      observationUnitDbIds = as.list(chunk),
      observationVariableDbIds = as.list(vars),
      pageSize = 2 * length(chunk) * length(vars)
    )))
  }
  out
}

# How many of these values are stored (same plot, trait and value).
bb_observations_found <- function(con, obs) {
  data <- bb_search_observations(con, obs)
  stored <- vapply(data, function(x) paste(x$observationUnitDbId, x$observationVariableDbId, x$value, sep = "\r"),
                   character(1))
  paste(obs$observation_unit_db_id, obs$observation_variable_db_id, obs$value, sep = "\r") %in% stored
}

bb_count_observations <- function(con, obs) sum(bb_observations_found(con, obs))

# Values stored for these plots and traits: key (unit, variable), value and
# observationDbId (the latest one when a plot has several).
bb_stored_values <- function(con, obs) {
  empty <- tibble::tibble(key = character(), value = character(), observation_db_id = character())
  if (!nrow(obs)) return(empty)
  data <- bb_search_observations(con, obs)
  if (!length(data)) return(empty)
  field <- function(f) vapply(data, function(x) as.character(x[[f]] %||% NA_character_), character(1))
  out <- tibble::tibble(key = paste(field("observationUnitDbId"), field("observationVariableDbId"), sep = "\r"),
                        value = field("value"), observation_db_id = field("observationDbId"))
  out <- out[!is.na(out$value) & nzchar(out$value), , drop = FALSE]
  out[!duplicated(out$key, fromLast = TRUE), , drop = FALSE]
}

# Same value: equal as numbers when both are numbers, else as text.
bb_same_value <- function(new, old) {
  a <- suppressWarnings(as.numeric(new))
  b <- suppressWarnings(as.numeric(old))
  num <- !is.na(a) & !is.na(b)
  out <- !is.na(old) & trimws(new) == trimws(old)
  out[num] <- abs(a[num] - b[num]) <= 1e-8 * pmax(1, abs(b[num]))
  out
}

# Delete stored values: the phenotype upload of the website with "overwrite"
# and "remove values", one row per plot; empty cells are deleted, the other
# cells are written with the values given.
bb_remove_values <- function(con, long, trait_full) {
  plots <- unique(long$plot_name)
  traits <- unique(long$trait)
  f <- data.frame(observationunit_name = plots, stringsAsFactors = FALSE, check.names = FALSE)
  for (t in traits) {
    part <- long[long$trait == t, , drop = FALSE]
    v <- part$value[match(plots, part$plot_name)]
    nm <- unname(trait_full[t])
    f[[if (is.na(nm)) t else nm]] <- ifelse(is.na(v), "", v)
  }
  path <- tempfile("breedtools_remove_", fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(f, path, row.names = FALSE, na = "")
  fields <- list(upload_spreadsheet_phenotype_file_format = "simple",
                 upload_spreadsheet_phenotype_data_level = "plots",
                 phenotype_upload_overwrite_values = "1",
                 phenotype_upload_remove_values = "1")
  bb_phenotype_upload(con, "verify", fields, path)
  bb_phenotype_upload(con, "store", fields, path)
  invisible(TRUE)
}

bb_phenotype_summary <- function(obs, found) {
  traits <- unique(obs$trait)
  tibble::tibble(
    trait = traits,
    observation_variable_db_id = vapply(traits, function(t) obs$observation_variable_db_id[obs$trait == t][1],
                                        character(1), USE.NAMES = FALSE),
    values_sent = vapply(traits, function(t) sum(obs$trait == t), integer(1), USE.NAMES = FALSE),
    values_found = if (is.null(found)) NA_integer_ else
      vapply(traits, function(t) sum(found[obs$trait == t]), integer(1), USE.NAMES = FALSE)
  )
}
