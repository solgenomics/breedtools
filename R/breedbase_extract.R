#' Read phenotypic data of trials from Breedbase
#'
#' Reads, through BrAPI, the plots of the given trials (plot name and
#' number, accession, rep, block, row, column, check) and their values of
#' the given traits, and returns one table with one row per plot and one
#' column per trait, ready for [plot_trait()], [remove_outliers()],
#' [analyze_two_stage()] and [analyze_single_stage()]. Repeated
#' measurements of a trait on a plot are averaged (by Breedbase).
#'
#' @param con A connection from [breedbase_login()].
#' @param trials Trial names, as in Breedbase.
#' @param traits Trait names, as in Breedbase, with or without the ontology
#'   id (e.g. `"fresh root yield"` or `"fresh root yield|CO_334:0000013"`).
#'
#' @return A tibble: `trial_name`, `location`, `year`, `plot_name`, `plots`
#'   (plot number), `accession`, `rep`, `block`, `row`, `col`, `is_check`
#'   and one numeric column per trait (named as given in `traits`). The
#'   trait names in Breedbase are in `attr(result, "trait_names")`.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassava-test.sgn.cornell.edu", "jdoe", site_username = "web_user_name")
#' pheno <- breedbase_phenotypes(con, trials = c("FALL_2026_Ibadan", "FALL_2026_Ikenne"),
#'                               traits = c("fresh root yield", "dry matter content percentage"))
#' plot_trait(pheno, "fresh root yield", by = "trial_name")
#' }
#'
#' @export
breedbase_phenotypes <- function(con, trials, traits) {
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  trials <- unique(trimws(as.character(trials)))
  traits <- unique(trimws(as.character(traits)))
  if (!length(trials) || any(!nzchar(trials))) stop("Give the trial names in `trials`.", call. = FALSE)
  if (!length(traits) || any(!nzchar(traits))) stop("Give the trait names in `traits`.", call. = FALSE)
  cat("Reading ", length(trials), " trial", if (length(trials) > 1) "s" else "", " from ", con$url, "\n", sep = "")

  ids <- bb_find_variables(con, traits)
  full <- attr(ids, "full_name")
  if (anyNA(ids)) {
    stop("Trait(s) not found in Breedbase: ", paste0("'", traits[is.na(ids)], "'", collapse = ", "),
         ". Use the trait name as shown in Breedbase.", call. = FALSE)
  }
  info <- lapply(trials, function(t) bb_find_study_info(con, t))
  missing <- trials[vapply(info, function(x) is.na(x$id), logical(1))]
  if (length(missing)) {
    stop("Trial(s) not found in Breedbase: ", paste0("'", missing, "'", collapse = ", "), ".", call. = FALSE)
  }

  out <- list()
  for (k in seq_along(trials)) {
    pl <- breedbase_trial_plots(con, trials[k])
    obs <- bb_search(con, "observations", list(
      studyDbIds = list(info[[k]]$id),
      observationVariableDbIds = as.list(unname(ids)),
      observationLevel = list("plot"),
      pageSize = max(1000L, 2L * nrow(pl) * length(traits))
    ))
    unit <- vapply(obs, function(o) as.character(o$observationUnitName %||% NA), character(1))
    var <- vapply(obs, function(o) as.character(o$observationVariableDbId %||% NA), character(1))
    val <- suppressWarnings(as.numeric(vapply(obs, function(o) as.character(o$value %||% NA), character(1))))
    tab <- tibble::tibble(trial_name = trials[k], location = info[[k]]$location, year = info[[k]]$year,
                          plot_name = pl$plot_name, plots = pl$plot_number, accession = pl$accession,
                          rep = pl$rep, block = pl$block, row = pl$row, col = pl$col,
                          is_check = pl$is_check)
    counts <- integer(0)
    for (j in seq_along(traits)) {
      here <- var == unname(ids[j]) & !is.na(val)
      v <- tapply(val[here], unit[here], mean)
      tab[[traits[j]]] <- unname(v[tab$plot_name])
      counts[j] <- sum(!is.na(tab[[traits[j]]]))
    }
    cat(sprintf("  [OK]   %s (%s): %d plots; %s\n", trials[k], info[[k]]$location %||% "?", nrow(tab),
                paste(sprintf("%s %d values", traits, counts), collapse = ", ")))
    out[[k]] <- tab
  }
  res <- dplyr::bind_rows(out)
  if (all(vapply(traits, function(t) all(is.na(res[[t]])), logical(1)))) {
    cat("  [NOTE] No values were found for these traits in these trials.\n")
  }
  attr(res, "trait_names") <- stats::setNames(unname(full), traits)
  res
}

# Id, location and year of a trial (BrAPI study) by name.
bb_find_study_info <- function(con, name) {
  hit <- list()
  for (nm in unique(c(name, bb_clean_name(name, quiet = TRUE)))) {
    data <- bb_search(con, "studies", list(studyNames = list(nm), pageSize = 10))
    hit <- Filter(function(x) identical(x$studyName, nm), data)
    if (length(hit)) break
  }
  if (!length(hit)) return(list(id = NA_character_, location = NA_character_, year = NA_character_))
  h <- hit[[1]]
  list(id = as.character(h$studyDbId %||% NA), location = as.character(h$locationName %||% NA),
       year = as.character(unlist(h$seasons)[1] %||% NA))
}
