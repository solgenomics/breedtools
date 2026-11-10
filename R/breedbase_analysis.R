# Analysis results to Breedbase.
#
# First stage (one value per field plot, e.g. spatially corrected plot
# values): saved as a trait on the plots of the field trial, through BrAPI
# observations (save_breedbase_phenotypes()).
#
# Second stage (one value per genotype, e.g. adjusted means across
# locations): saved as a Breedbase analysis (an "analysis" project), like the
# Breedbase mixed-models tool does, through POST /ajax/analysis/store/json
# (not BrAPI; it accepts the BrAPI login token as `sgn_session_id`). Each
# trait is saved composed with a statistic term, e.g. "fresh root
# yield|CO_334:0000013" + "Adjusted means from BLUEs using LMER R|SGNSTAT:0000036".

#' First-stage plot values from an analysis
#'
#' Turns analyses of a trial into one value per field plot, ready for
#' [save_breedbase_first_stage()]: the genotype's adjusted mean plus the plot's
#' residual. This removes the design effects (blocks, rows, columns,
#' replicates) and the spatial trend, and keeps the plot-to-plot variation, so
#' the result has as many plots as the field trial. With SpATS
#' (`spatial = TRUE` or `stage1 = "spatial"`) these are spatially corrected
#' plot values.
#'
#' @param fits One analysis, or a named list of analyses, one per trait (the
#'   names are the Breedbase trait names). Analyses from
#'   [analyze_augmented_row_column_design()], [analyze_augmented_prep_design()],
#'   [analyze_augmented_prep_two_stage()] (stage-1 fits per location) or
#'   [analyze_resolvable_row_column()].
#' @param design The design that was saved to Breedbase. Leave it `NULL`
#'   for data read from Breedbase ([breedbase_phenotypes()]): the plot names
#'   are then taken from the data.
#' @param trial_name The trial name used when saving the design (not needed
#'   without a design).
#' @param trait Trait name when `fits` is a single analysis.
#' @param plot_name_pattern,plot_names,location The plot naming used in
#'   [save_breedbase()]; `location` is the Breedbase location of the trial,
#'   needed for the `{location}` field when the design has no stored plot
#'   names.
#' @param breedbase_plots Optional plots of the saved trial, from
#'   [breedbase_trial_plots()]. The design is then checked against the saved
#'   trial (the same accession at every plot number), and the plot names are
#'   taken exactly as saved.
#'
#' @return A tibble with `plot_name`, `trial_name`, `accession`,
#'   `plot_number`, `block`, `rep`, `row`, `col`, `is_check` and one column
#'   per trait. Its attribute `"method"` is `"spatial"` when any analysis used
#'   SpATS, otherwise `"mixed"`.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' set.seed(1)
#' d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
#' fit <- analyze_augmented_row_column_design(d, "yield")
#' values <- first_stage_values(fit, rc, "AYT_2026_Ibadan", trait = "fresh root yield",
#'                             location = "Ibadan")
#' head(values)
#'
#' @export
first_stage_values <- function(fits, design = NULL, trial_name = NULL, trait = NULL,
                               plot_name_pattern = NULL, plot_names = NULL,
                               breedbase_plots = NULL, location = NULL) {
  if (bb_is_analysis(fits)) {
    if (is.null(trait) || length(trait) != 1 || !nzchar(trait)) {
      stop("Give the Breedbase `trait` name of this analysis, e.g. trait = \"fresh root yield\".", call. = FALSE)
    }
    fits <- stats::setNames(list(fits), trait)
  }
  if (!is.list(fits) || is.null(names(fits)) || any(!nzchar(names(fits)))) {
    stop("`fits` must be one analysis (with `trait`) or a list of analyses named by trait.", call. = FALSE)
  }
  for (nm in names(fits)) {
    if (!bb_is_analysis(fits[[nm]])) stop("`fits$", nm, "` is not an analysis from breedtools.", call. = FALSE)
  }
  # Data read from Breedbase (breedbase_phenotypes()): the plot names are in
  # the data, no design object is needed.
  if (is.null(design)) return(bb_first_stage_from_data(fits))

  table <- bb_design_plot_table(design, trial_name, plot_name_pattern, plot_names, location = location)
  is_prep <- inherits(design, "augmented_prep_design")
  key_table <- paste(if (is_prep) table$.design_location else "", table$.design_plots, sep = "\r")
  spatial <- FALSE
  for (nm in names(fits)) {
    plots <- bb_first_stage_plots(fits[[nm]])
    spatial <- spatial || isTRUE(attr(plots, "spatial"))
    if (!"plots" %in% names(plots)) stop("The analysis of '", nm, "' has no `plots` column to match.", call. = FALSE)
    key <- paste(if (is_prep) as.character(plots$location) else "", plots$plots, sep = "\r")
    unknown <- setdiff(key, key_table)
    if (length(unknown)) {
      stop("The analysis of '", nm, "' has ", length(unknown), " plots that are not in the design. ",
           "Was it run on this design?", call. = FALSE)
    }
    table[[nm]] <- plots$value[match(key_table, key)]
  }
  table$.design_location <- NULL
  table$.design_plots <- NULL
  if (!is.null(breedbase_plots)) {
    # The saved trial must be this design: same accession at each plot number.
    bp <- as.data.frame(breedbase_plots)
    i <- match(table$plot_number, bp$plot_number)
    differ <- is.na(i) | bp$accession[pmax(i, 1L)] != table$accession
    if (any(differ)) {
      stop(sum(differ), " of ", nrow(table), " plots differ from the trial saved in Breedbase (another accession ",
           "at the same plot number, e.g. plot ", table$plot_number[which(differ)[1]], ": '",
           table$accession[which(differ)[1]], "' here, '", if (is.na(i[which(differ)[1]])) "none" else
           bp$accession[i[which(differ)[1]]], "' in Breedbase). The design in R is not the one that was saved ",
           "(another seed or settings?).", call. = FALSE)
    }
    table$plot_name <- bp$plot_name[i]       # the names exactly as saved
  }
  attr(table, "method") <- if (spatial) "spatial" else "mixed"
  table
}

bb_first_stage_from_data <- function(fits) {
  out <- NULL
  spatial <- FALSE
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    parts <- if (!is.null(fit$stage1_fits) && is.null(fit$plot_results)) fit$stage1_fits else list(fit)
    pieces <- lapply(seq_along(parts), function(k) {
      pr <- parts[[k]]$plot_results
      if (is.null(pr) || !"plot_name" %in% names(pr)) {
        stop("The analysis of '", nm, "' has no plot names: give the `design` and `trial_name`, or analyse ",
             "data that have a `plot_name` column (e.g. from breedbase_phenotypes()).", call. = FALSE)
      }
      v <- bb_first_stage_plots(parts[[k]])
      spatial <<- spatial || isTRUE(attr(v, "spatial"))
      gcol <- parts[[k]]$input$genotype_col %||% "accession"
      tibble::tibble(plot_name = as.character(pr$plot_name),
                     trial_name = if (length(parts) > 1 || !is.null(names(parts))) names(parts)[k] else
                       if ("trial_name" %in% names(pr)) as.character(pr$trial_name) else NA_character_,
                     accession = as.character(pr[[gcol]]), value = v$value)
    })
    piece <- dplyr::bind_rows(pieces)
    piece <- piece[!is.na(piece$value), , drop = FALSE]
    names(piece)[names(piece) == "value"] <- nm
    out <- if (is.null(out)) piece else dplyr::full_join(out, piece, by = c("plot_name", "trial_name", "accession"))
  }
  attr(out, "method") <- if (spatial) "spatial" else "mixed"
  out
}

bb_is_analysis <- function(x) {
  inherits(x, c("augmented_row_column_analysis", "augmented_prep_analysis", "augmented_prep_two_stage",
                "resolvable_row_column_analysis", "two_stage_analysis", "stage1_analysis",
                "single_stage_analysis")) ||
    (is.list(x) && !is.null(x$plot_results) && !is.null(x$adjusted_means))
}

# Plot values (adjusted mean of the genotype + plot residual) of one analysis.
bb_first_stage_plots <- function(fit, location = NA_character_) {
  if (!is.null(fit$stage1_fits) && is.null(fit$plot_results)) {
    parts <- lapply(names(fit$stage1_fits), function(loc) bb_first_stage_plots(fit$stage1_fits[[loc]], loc))
    out <- dplyr::bind_rows(parts)
    attr(out, "spatial") <- any(vapply(parts, function(p) isTRUE(attr(p, "spatial")), logical(1))) ||
      ("stage1" %in% names(fit$diagnostics) && identical(as.character(fit$diagnostics$stage1[1]), "spatial"))
    return(out)
  }
  if (inherits(fit, "SpATS")) {
    # A SpATS stage-1 fit (analyze_augmented_prep_two_stage(stage1 = "spatial")).
    pred <- stats::predict(fit, which = "geno", predFixed = "marginal")
    means <- stats::setNames(pred$predicted.values, as.character(pred$geno))
    out <- tibble::tibble(plots = as.integer(fit$data$plots), location = location,
                          value = unname(means[as.character(fit$data$geno)]) + unname(fit$residuals))
    attr(out, "spatial") <- TRUE
    return(out)
  }
  pr <- fit$plot_results
  if (is.null(pr) || !all(c("residual") %in% names(pr))) {
    stop("This analysis has no plot results (residuals).", call. = FALSE)
  }
  gcol <- fit$input$genotype_col %||% "all_entries"
  if (!gcol %in% names(pr)) stop("The plot results have no genotype column '", gcol, "'.", call. = FALSE)
  am <- fit$adjusted_means
  means <- stats::setNames(am$adjusted_mean, as.character(am$genotype))
  out <- tibble::tibble(
    plots = as.integer(pr$plots),
    location = if ("location" %in% names(pr)) as.character(pr$location) else location,
    value = unname(means[as.character(pr[[gcol]])]) + pr$residual
  )
  diag_value <- function(name) {
    dg <- fit$diagnostics
    if (is.data.frame(dg) && name %in% names(dg)) as.character(dg[[name]][1]) else ""
  }
  attr(out, "spatial") <- isTRUE(fit$input$spatial) ||
    grepl("spats|spatial", paste(diag_value("model_type"), diag_value("stage1_model")), ignore.case = TRUE)
  out
}

#' Find or create the Breedbase trait for adjusted values
#'
#' First-stage values need their own trait, so raw and adjusted values stay
#' apart. Breedbase makes such traits by composing the measured trait with a
#' statistic term, e.g. "fresh root yield" + "Adjusted Means from Spatial
#' Correction using SpATS R" gives the composed trait
#' `"fresh root yield|Adjusted Means from Spatial Correction using SpATS R"`
#' (ontology `COMP`). This function finds that trait, or creates it with
#' `create = TRUE` through the same calls as the "Compose a new trait" page
#' of Breedbase.
#'
#' @param con A connection from [breedbase_login()].
#' @param traits Measured trait names, as in Breedbase.
#' @param statistic The statistic term, with or without its id, as listed by
#'   [breedbase_statistics()]. By default it is chosen from this Breedbase's
#'   own list of statistics: the spatial correction term (SpATS). It is
#'   always checked against the server.
#' @param create Logical. Create the traits that do not exist yet.
#'
#' @return A tibble with `trait`, `adjusted_trait`, `adjusted_trait_db_id`
#'   and `status` (`"exists"`, `"created"` or `"missing"`).
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' breedbase_adjusted_trait(con, "fresh root yield")                  # look only
#' breedbase_adjusted_trait(con, "fresh root yield", create = TRUE)   # create if needed
#' }
#'
#' @export
breedbase_adjusted_trait <- function(con, traits, statistic = NULL, create = FALSE) {
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  ids <- bb_find_variables(con, traits)
  full <- attr(ids, "full_name")
  if (anyNA(ids)) {
    stop("Trait(s) not found in Breedbase: ", paste0("'", traits[is.na(ids)], "'", collapse = ", "),
         ". Use the trait name as shown in Breedbase.", call. = FALSE)
  }
  base <- sub("\\|[^|]*$", "", full)          # trait name without its ontology id
  # A first guess for a quick BrAPI-only lookup of an existing trait: the given
  # statistic, or the name of Breedbase's spatial correction statistic. It is
  # used only if the server has a trait with exactly that name.
  guess <- bb_choose_statistic(NULL, statistic, kind = "spatial", offline = TRUE)$name
  out <- tibble::tibble(trait = traits, adjusted_trait = paste(base, guess, sep = "|"),
                        adjusted_trait_db_id = NA_character_, status = "missing")
  # An existing composed trait is found through BrAPI (Breedbase links it as
  # a variable), so no web page is needed when it exists.
  for (i in seq_along(traits)) {
    id <- bb_find_variable_exact(con, out$adjusted_trait[i])
    if (!is.na(id)) {
      out$adjusted_trait_db_id[i] <- id
      out$status[i] <- "exists"
    }
  }
  if (all(out$status == "exists")) return(out)
  # Otherwise the statistic comes from the server's own list of terms, and
  # the trait is found by its components (trait id + statistic id), so names
  # never decide whether a trait exists.
  stat <- bb_choose_statistic(con, statistic, kind = "spatial")
  out$adjusted_trait[out$status != "exists"] <- paste(base, stat$name, sep = "|")[out$status != "exists"]
  for (i in which(out$status != "exists")) {
    components <- c(ids[[i]], stat$id)
    id <- bb_composed_trait_id(con, components)
    if (is.na(id) && isTRUE(create)) {
      bb_ajax_post(con, "ajax/onto/store_composed_term",
                   list(new_trait_names = stats::setNames(list(paste(components, collapse = ",")),
                                                          out$adjusted_trait[i])),
                   what = paste0("create trait '", out$adjusted_trait[i], "'"), json = "new_trait_names")
      id <- bb_composed_trait_id(con, components)
      if (is.na(id)) stop("Breedbase did not create the trait '", out$adjusted_trait[i], "'.", call. = FALSE)
      out$status[i] <- "created"
    } else if (!is.na(id)) {
      out$status[i] <- "exists"
    }
    out$adjusted_trait_db_id[i] <- id
  }
  out
}

# Full name ("name|ONTOLOGY:id") of the variable named exactly `name`.
bb_find_variable_full <- function(con, name) {
  data <- bb_search(con, "variables", list(observationVariableNames = list(name), pageSize = 20))
  full <- vapply(data, function(x) as.character(x$observationVariableName %||% NA_character_), character(1))
  hit <- full[sub("\\|[^|]*$", "", full) == name | full == name]
  if (length(hit)) hit[1] else NA_character_
}

# Id of the variable whose name (without its "|ONTOLOGY:id") is exactly `name`.
bb_find_variable_exact <- function(con, name) {
  data <- bb_search(con, "variables", list(observationVariableNames = list(name), pageSize = 20))
  full <- vapply(data, function(x) as.character(x$observationVariableName %||% NA_character_), character(1))
  ids <- vapply(data, function(x) as.character(x$observationVariableDbId %||% NA_character_), character(1))
  hit <- which(sub("\\|[^|]*$", "", full) == name | full == name)
  if (length(hit)) ids[hit[1]] else NA_character_
}

#' List the statistic terms of Breedbase (SGNSTAT ontology)
#'
#' Breedbase saves analysis results as a trait combined with a statistic
#' term of its statistics ontology (SGNSTAT), e.g. "Adjusted Means from
#' Spatial Correction using SpATS R". Breedbase loads these terms from the
#' file `ontology/cxgn_statistics.obo` of its source code, which breedtools
#' includes (`system.file("extdata", "cxgn_statistics.obo", package =
#' "breedtools")`), so the terms can be listed and checked without a server.
#' With a connection, the terms of that Breedbase are listed and compared
#' with the ontology file. When Breedbase adds statistics, the bundled copy
#' is refreshed with the script `data-raw/update_statistics_ontology.R` of
#' the package source.
#'
#' @param con Optional connection from [breedbase_login()]. Without it, the
#'   terms of the ontology file are listed.
#' @param pattern Optional text the term must contain (any case).
#'
#' @return A tibble with `statistic` (the name), `accession` (e.g.
#'   `"SGNSTAT:0000054"`) and `definition` (from the ontology file); with a
#'   connection also `db_id` (the server's id) and `in_ontology`
#'   (`TRUE` when name and accession match the ontology file).
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' breedbase_statistics(con)
#' breedbase_statistics(con, "BLUE")
#' }
#' breedbase_statistics(pattern = "spatial")      # from the ontology file
#'
#' @export
breedbase_statistics <- function(con = NULL, pattern = "") {
  obo <- bb_obo_statistics()
  if (is.null(con)) {
    keep <- grepl(pattern, obo$statistic, ignore.case = TRUE, fixed = FALSE)
    return(obo[keep, , drop = FALSE])
  }
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  out <- bb_statistic_terms(con, pattern)
  i <- match(out$accession, obo$accession)
  out$definition <- obo$definition[i]
  out$in_ontology <- !is.na(i) & tolower(obo$statistic[i]) == tolower(out$statistic)
  out
}

# The SGNSTAT terms of the ontology file Breedbase loads (bundled copy).
bb_obo_statistics <- function() {
  path <- system.file("extdata", "cxgn_statistics.obo", package = "breedtools")
  if (!nzchar(path)) stop("The statistics ontology file is missing from breedtools.", call. = FALSE)
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  starts <- which(lines == "[Term]")
  ends <- c(starts[-1] - 1L, length(lines))
  field <- function(block, key) {
    hit <- grep(paste0("^", key, ":"), block, value = TRUE)
    if (length(hit)) trimws(sub(paste0("^", key, ":"), "", hit[1])) else NA_character_
  }
  terms <- lapply(seq_along(starts), function(k) {
    block <- lines[starts[k]:ends[k]]
    def <- field(block, "def")
    if (is.na(def)) def <- field(block, "definition")
    tibble::tibble(statistic = field(block, "name"), accession = field(block, "id"),
                   definition = gsub('^"|"\\s*\\[.*$|"$', "", def %||% NA_character_),
                   obsolete = identical(field(block, "is_obsolete"), "true"))
  })
  out <- dplyr::bind_rows(terms)
  out <- out[!out$obsolete & out$accession != "SGNSTAT:0000000", c("statistic", "accession", "definition")]
  out[order(out$statistic), , drop = FALSE]
}

# Check a statistic against the ontology file. Accepts "name", "name|ID" or
# "ID"; returns the ontology's name and accession, or stops with suggestions.
bb_check_statistic_ontology <- function(statistic) {
  obo <- bb_obo_statistics()
  statistic <- trimws(statistic)
  id_rx <- "[A-Za-z][A-Za-z0-9_]*:[0-9]+"
  has_id <- grepl(paste0(id_rx, "\\s*$"), statistic)
  id <- if (has_id) trimws(sub(paste0("^.*?(", id_rx, ")\\s*$"), "\\1", statistic)) else NA_character_
  name <- trimws(sub(paste0("\\|?\\s*", id_rx, "\\s*$"), "", statistic))
  if (!is.na(id)) {
    if (!startsWith(id, "SGNSTAT:")) {
      stop("'", id, "' is not a statistic: use a term of the Breedbase statistics ontology (SGNSTAT), ",
           "see breedbase_statistics().", call. = FALSE)
    }
    row <- obo[obo$accession == id, , drop = FALSE]
    if (!nrow(row)) stop("'", id, "' is not a term of the Breedbase statistics ontology (SGNSTAT).", call. = FALSE)
    if (nzchar(name) && tolower(name) != tolower(row$statistic)) {
      stop("'", id, "' is \"", row$statistic, "\" in the Breedbase statistics ontology, not \"", name,
           "\". Use the name and id of the same term.", call. = FALSE)
    }
    return(list(name = row$statistic, accession = row$accession))
  }
  row <- obo[tolower(obo$statistic) == tolower(name), , drop = FALSE]
  if (!nrow(row)) {
    close <- obo$statistic[agrepl(name, obo$statistic, ignore.case = TRUE, max.distance = 0.2)]
    if (!length(close)) {
      words <- strsplit(tolower(name), "[^a-z0-9]+")[[1]]
      words <- words[nchar(words) > 3]
      score <- vapply(obo$statistic, function(x) sum(vapply(words, grepl, logical(1), x = tolower(x), fixed = TRUE)),
                      numeric(1))
      close <- obo$statistic[order(-score)][seq_len(min(3, sum(score > 0)))]
    }
    stop("\"", name, "\" is not a term of the Breedbase statistics ontology (SGNSTAT)",
         if (length(close)) paste0(". Did you mean: ", paste0("\"", utils::head(close, 3), "\"", collapse = " or "), "?")
         else ". See breedbase_statistics().", call. = FALSE)
  }
  list(name = row$statistic[1], accession = row$accession[1])
}

bb_statistic_terms <- function(con, pattern = "", db = "SGNSTAT") {
  hits <- bb_ajax_get(con, "ajax/onto/match", what = "list the statistic terms",
                      query = list(db_name = db, term_name = pattern))
  rows <- if (is.null(names(hits))) hits else list()
  out <- tibble::tibble(
    statistic = vapply(rows, function(h) as.character(h$cvterm_name %||% NA_character_), character(1)),
    accession = vapply(rows, function(h) as.character(h$accession %||% NA_character_), character(1)),
    db_id = vapply(rows, function(h) as.character(h$cvterm_id %||% NA_character_), character(1))
  )
  out[order(out$statistic), , drop = FALSE]
}

# The statistic to use: the one given (checked on the server), or one chosen
# from the server's list by meaning. Returns id, name and accession.
bb_choose_statistic <- function(con, statistic = NULL, kind = c("spatial", "blue", "blup"), offline = FALSE) {
  kind <- match.arg(kind)
  if (!is.null(statistic)) {
    # Must be a term of the statistics ontology; then found on the server.
    onto <- bb_check_statistic_ontology(statistic)
    if (isTRUE(offline) || is.null(con)) return(list(id = NA_character_, name = onto$name, accession = onto$accession))
    return(bb_find_statistic(con, paste0(onto$name, "|", onto$accession)))
  }
  terms <- if (isTRUE(offline) || is.null(con)) {
    o <- bb_obo_statistics()
    tibble::tibble(statistic = o$statistic, accession = o$accession, db_id = NA_character_)
  } else {
    bb_statistic_terms(con)
  }
  wanted <- switch(kind,
    spatial = list(c("spatial", "spats"), c("spatial", "adjusted")),
    blue = list(c("adjusted means", "blue", "lmer"), c("adjusted means", "blue"), c("blue")),
    blup = list(c("adjusted means", "blup", "lmer"), c("adjusted means", "blup"), c("blup"))
  )
  for (words in wanted) {
    ok <- Reduce(`&`, lapply(words, function(w) grepl(w, terms$statistic, ignore.case = TRUE)))
    if (any(ok)) {
      hit <- terms[ok, , drop = FALSE]
      hit <- hit[order(nchar(hit$statistic)), , drop = FALSE][1, ]
      if (!(isTRUE(offline) || is.null(con))) {
        onto <- bb_obo_statistics()
        if (!any(onto$accession == hit$accession & tolower(onto$statistic) == tolower(hit$statistic))) {
          bt_say("  [NOTE] the statistic '", hit$statistic, "' (", hit$accession, ") of this Breedbase is not in the ",
              "statistics ontology file of Breedbase.\n", sep = "")
        }
      }
      return(list(id = hit$db_id, name = hit$statistic, accession = hit$accession))
    }
  }
  stop("No ", switch(kind, spatial = "spatial correction", blue = "BLUE", blup = "BLUP"),
       " statistic was found", if (isTRUE(offline) || is.null(con)) " in the statistics ontology" else " in this Breedbase", ". Choose one of its statistics with `statistic = ...`: ",
       paste(utils::head(terms$statistic, 30), collapse = "; "), call. = FALSE)
}

# The SGNSTAT term id of a statistic (name with or without "|SGNSTAT:id").
bb_find_statistic <- function(con, statistic) {
  name <- trimws(sub("\\|[^|]*$", "", statistic))
  db <- if (grepl("|", statistic, fixed = TRUE)) sub(":.*$", "", sub("^.*\\|", "", statistic)) else "SGNSTAT"
  hits <- bb_ajax_get(con, "ajax/onto/match", what = "find the statistic term",
                      query = list(db_name = db, term_name = name))
  rows <- if (is.null(names(hits))) hits else list()
  match <- Filter(function(h) identical(tolower(trimws(h$cvterm_name %||% "")), tolower(name)), rows)
  if (!length(match)) {
    stop("The statistic '", name, "' was not found in this Breedbase (ontology ", db, ").", call. = FALSE)
  }
  list(id = as.character(match[[1]]$cvterm_id), name = match[[1]]$cvterm_name,
       accession = as.character(match[[1]]$accession %||% NA_character_))
}

# Id of the composed trait made of exactly these components, or NA.
bb_composed_trait_id <- function(con, components) {
  body <- bb_ajax_get(con, "ajax/onto/get_trait_from_exact_components", what = "look up the composed trait",
                      query = list(`ids[]` = components), allow_error = TRUE)
  if (is.null(body$trait_id)) NA_character_ else as.character(body$trait_id)
}

#' Save first-stage results to Breedbase as a trait of the field plots
#'
#' Saves plot-level first-stage results (from [first_stage_values()]) on the
#' plots of the field trial, as values of a trait made for adjusted values,
#' so raw and adjusted values stay apart. By default that trait is the
#' measured trait composed with a statistic, e.g. "fresh root yield|Adjusted
#' Means from Spatial Correction using SpATS R": it is found, or created
#' (after you confirm) when it does not exist yet; see
#' [breedbase_adjusted_trait()]. Or name your own traits with `as_trait`.
#' The values are saved through BrAPI, like [save_breedbase_phenotypes()].
#'
#' @param values A table from [first_stage_values()] (or any table with a
#'   `plot_name` column and one numeric column per trait).
#' @param con A connection from [breedbase_login()]; `NULL` with
#'   `dry_run = TRUE` returns what would be sent.
#' @param as_trait Optional named vector: for each trait column of `values`,
#'   an existing Breedbase trait to save it under, e.g.
#'   `c("fresh root yield" = "fresh root yield adjusted|COMP:0000123")`.
#' @param statistic The statistic composed with each trait when `as_trait`
#'   is not given. By default chosen from this Breedbase's own statistics
#'   (the spatial correction term); see [breedbase_statistics()].
#' @param create_trait Logical. Create the composed traits that do not exist
#'   yet (you are asked first when `confirm = TRUE`).
#' @param dry_run,confirm,overwrite,remove_empty As in
#'   [save_breedbase_phenotypes()]: run it again after a new analysis with
#'   `overwrite = TRUE` to replace the stored values.
#'
#' @inheritSection breedbase_check Slow servers
#'
#' @return As [save_breedbase_phenotypes()]: one row per trait with the
#'   values sent and found.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' set.seed(1)
#' d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
#' fit <- analyze_augmented_row_column_design(d, "yield")
#' values <- first_stage_values(fit, rc, "AYT_2026_Ibadan", trait = "fresh root yield",
#'                             location = "Ibadan")
#' obs <- save_breedbase_first_stage(values, dry_run = TRUE)
#' head(obs)
#'
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' save_breedbase_first_stage(values, con, dry_run = TRUE)   # check; shows the trait to create
#' save_breedbase_first_stage(values, con)                   # creates the trait if needed, saves
#' }
#'
#' @export
save_breedbase_first_stage <- function(values, con = NULL, as_trait = NULL, statistic = NULL,
                                       create_trait = TRUE, dry_run = FALSE, confirm = interactive(),
                                       overwrite = FALSE, remove_empty = FALSE) {
  if (!is.data.frame(values) || !"plot_name" %in% names(values)) {
    stop("`values` must be a table with a `plot_name` column, as from first_stage_values().", call. = FALSE)
  }
  design_cols <- c("plot_name", "trial_name", "accession", "plot_number", "block", "rep", "row", "col", "is_check")
  traits <- setdiff(names(values), design_cols)
  if (!length(traits)) stop("No trait columns in `values`.", call. = FALSE)
  values <- as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  # Without a connection the statistic cannot be checked; the dry run shows
  # the usual Breedbase name.
  stat_name <- bb_choose_statistic(NULL, statistic, kind = "spatial", offline = TRUE)$name
  trait_ids <- NULL

  if (!is.null(as_trait)) {
    if (is.null(names(as_trait)) || any(!nzchar(names(as_trait)))) {
      stop("`as_trait` must be named by the trait columns of `values`.", call. = FALSE)
    }
    absent <- setdiff(names(as_trait), traits)
    if (length(absent)) {
      stop("Not trait columns of `values`: ", paste(absent, collapse = ", "), ". Trait columns are: ",
           paste(traits, collapse = ", "), ".", call. = FALSE)
    }
    same <- names(as_trait)[trimws(sub("\\|.*$", "", as_trait)) == trimws(sub("\\|.*$", "", names(as_trait)))]
    if (length(same)) {
      stop("First-stage values of '", same[1], "' would be saved as the measured trait itself, mixing raw and ",
           "adjusted values on the same plots. Choose a trait made for adjusted values.", call. = FALSE)
    }
    targets <- as_trait[traits[traits %in% names(as_trait)]]
  } else if (is.null(con)) {
    targets <- stats::setNames(paste(sub("\\|.*$", "", traits), stat_name, sep = "|"), traits)
  } else {
    if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
    bt_say("Finding the traits for the adjusted values in ", con$url, "\n", sep = "")
    adj <- breedbase_adjusted_trait(con, traits, statistic = statistic, create = FALSE)
    for (i in seq_len(nrow(adj))) {
      if (adj$status[i] == "exists") {
        bt_say("  [OK]   trait '", adj$adjusted_trait[i], "' (id ", adj$adjusted_trait_db_id[i], ")\n", sep = "")
      } else {
        bt_say("  [NEW]  trait '", adj$adjusted_trait[i], "' will be created\n", sep = "")
      }
    }
    missing <- adj$status == "missing"
    if (any(missing) && !isTRUE(dry_run)) {
      if (!isTRUE(create_trait)) {
        stop("The trait(s) for adjusted values do not exist: ",
             paste0("'", adj$adjusted_trait[missing], "'", collapse = ", "),
             ". Use create_trait = TRUE, or name existing traits with `as_trait`.", call. = FALSE)
      }
      if (isTRUE(confirm)) {
        answer <- readline(sprintf("Create %d new trait%s in Breedbase (%s)? Type 'yes' to continue: ",
                                   sum(missing), if (sum(missing) > 1) "s" else "",
                                   paste(adj$adjusted_trait[missing], collapse = ", ")))
        if (!identical(tolower(trimws(answer)), "yes")) {
          message("Nothing was created or saved.")
          return(invisible(NULL))
        }
      }
      adj <- breedbase_adjusted_trait(con, traits, statistic = statistic, create = TRUE)
      for (i in which(adj$status == "created")) {
        bt_say("  [OK]   trait '", adj$adjusted_trait[i], "' created (id ", adj$adjusted_trait_db_id[i], ")\n", sep = "")
      }
    }
    targets <- stats::setNames(adj$adjusted_trait, traits)
    trait_ids <- stats::setNames(ifelse(is.na(adj$adjusted_trait_db_id), "(new)", adj$adjusted_trait_db_id),
                                 adj$adjusted_trait)
  }

  out <- values[, c("plot_name", names(targets)), drop = FALSE]
  names(out) <- c("plot_name", unname(targets))
  bb_save_phenotypes(out, con, plot = "plot_name", traits = unname(targets), dry_run = dry_run,
                     confirm = confirm, trait_ids = trait_ids, overwrite = overwrite, remove_empty = remove_empty)
}

#' Read the plots of a trial saved in Breedbase
#'
#' Returns the plots of a Breedbase trial (through BrAPI): plot name,
#' number, accession, rep, block, row and column. Use it to check that a
#' design in R is the one saved in Breedbase (see `breedbase_plots` in
#' [first_stage_values()]), or to work with a trial made elsewhere.
#'
#' @param con A connection from [breedbase_login()].
#' @param trial_name The trial name in Breedbase.
#'
#' @return A tibble with `plot_name`, `plot_number`, `accession`, `rep`,
#'   `block`, `row`, `col` and `is_check`, ordered by plot number.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' plots <- breedbase_trial_plots(con, "AYT_2026_Ibadan")
#' }
#'
#' @export
breedbase_trial_plots <- function(con, trial_name) {
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  id <- bb_find_study(con, trial_name)
  if (is.na(id)) stop("Trial '", trial_name, "' was not found in Breedbase.", call. = FALSE)
  pl <- bb_get_all(con, "observationunits", list(studyDbId = id, observationUnitLevelName = "plot"))
  level <- function(p, name) {
    rel <- Filter(function(r) identical(r$levelName, name), p$observationUnitPosition$observationLevelRelationships %||% list())
    if (length(rel)) suppressWarnings(as.integer(rel[[1]]$levelCode)) else NA_integer_
  }
  num <- function(x) suppressWarnings(as.integer(x %||% NA))
  out <- tibble::tibble(
    plot_name = vapply(pl, function(p) as.character(p$observationUnitName %||% NA), character(1)),
    plot_number = vapply(pl, function(p) num(p$observationUnitPosition$observationLevel$levelCode), integer(1)),
    accession = vapply(pl, function(p) as.character(p$germplasmName %||% NA), character(1)),
    rep = vapply(pl, level, integer(1), name = "rep"),
    block = vapply(pl, level, integer(1), name = "block"),
    row = vapply(pl, function(p) num(p$observationUnitPosition$positionCoordinateY), integer(1)),
    col = vapply(pl, function(p) num(p$observationUnitPosition$positionCoordinateX), integer(1)),
    is_check = vapply(pl, function(p) identical(p$observationUnitPosition$entryType, "CHECK"), logical(1))
  )
  out[order(out$plot_number), , drop = FALSE]
}

#' Second-stage genotype values from analyses
#'
#' Collects the adjusted means of analyses (one per trait) into one table
#' with one row per accession, ready for [save_breedbase_analysis()].
#'
#' @param fits One analysis, or a named list of analyses, one per trait (the
#'   names are the Breedbase trait names). Any breedtools analysis with
#'   `adjusted_means`, e.g. [analyze_augmented_prep_two_stage()].
#' @param trait Trait name when `fits` is a single analysis.
#'
#' @return A tibble with `accession` and one column per trait. Its attribute
#'   `"estimate"` is `"blup"` when the genotypes were random, else `"blue"`.
#'
#' @examples
#' prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
#' d <- prep$design
#' set.seed(1)
#' d$yield <- rnorm(nrow(d), 25, 3)
#' fit <- analyze_augmented_prep_two_stage(d, "yield")
#' head(second_stage_values(fit, trait = "fresh root yield"))
#'
#' @export
second_stage_values <- function(fits, trait = NULL) {
  if (bb_is_analysis(fits) || (is.list(fits) && !is.null(fits$adjusted_means))) {
    if (is.null(trait) || length(trait) != 1 || !nzchar(trait)) {
      stop("Give the Breedbase `trait` name of this analysis, e.g. trait = \"fresh root yield\".", call. = FALSE)
    }
    fits <- stats::setNames(list(fits), trait)
  }
  if (!is.list(fits) || is.null(names(fits)) || any(!nzchar(names(fits)))) {
    stop("`fits` must be one analysis (with `trait`) or a list of analyses named by trait.", call. = FALSE)
  }
  out <- NULL
  random <- FALSE
  for (nm in names(fits)) {
    am <- fits[[nm]]$adjusted_means
    if (is.null(am) || !all(c("genotype", "adjusted_mean") %in% names(am))) {
      stop("`fits$", nm, "` has no adjusted means.", call. = FALSE)
    }
    random <- random || identical(fits[[nm]]$input$entries, "random")
    part <- tibble::tibble(accession = as.character(am$genotype), value = am$adjusted_mean)
    names(part)[2] <- nm
    out <- if (is.null(out)) part else dplyr::full_join(out, part, by = "accession")
  }
  attr(out, "estimate") <- if (random) "blup" else "blue"
  out
}

#' Save analysis results to Breedbase as an analysis
#'
#' Saves analysis results as a Breedbase analysis (a project of type
#' analysis), through the same upload as the "Upload analysis" form of
#' Breedbase:
#' \describe{
#'   \item{`level = "accession"`}{"Analysis of: accessions": one value per
#'     genotype, e.g. second-stage adjusted means from
#'     [second_stage_values()].}
#'   \item{`level = "plot"`}{"Analysis of: plots": one value per field plot,
#'     e.g. first-stage values from [first_stage_values()]. Breedbase finds
#'     the field trial(s) of the plots and gives the analysis one plot per
#'     field plot (named `"<analysis_name>_<field plot name>"`), with the
#'     same accession, rep, block, row and column.}
#' }
#' Each trait is saved together with a statistic term, as Breedbase does for
#' every analysis. The analysis then appears on the Analyses page and can be
#' used in later analyses. Every name (analysis, accessions or plots, traits,
#' statistic) is checked first, and all problems are listed before anything
#' is written. The field trials do not change.
#'
#' The form's fields map to the arguments: Analysis Name (`analysis_name`,
#' unique), Analysis Description (`description`), Analysis Protocol
#' (`protocol`), Analysis Of (`level`), Breeding Program (`program`),
#' Statistical Ontology Term (`statistic`) and the optional time terms
#' (`time_terms`). The model type is Breedbase's own one for uploaded
#' results (`uploaded_generic_analysis_model`).
#'
#' @param data A data frame with one row per accession (`level =
#'   "accession"`) or per plot (`level = "plot"`) and one numeric column per
#'   trait.
#' @param con A connection from [breedbase_login()]; `NULL` with
#'   `dry_run = TRUE` returns the file that would be uploaded.
#' @param analysis_name Name of the new analysis (must not exist yet). Spaces
#'   and special characters are replaced with `"_"` (only letters, digits,
#'   `_`, `-` and `.` are kept).
#' @param program Breeding program of the analysis.
#' @param description Description of the analysis (required).
#' @param year Year of the analysis (default: this year).
#' @param level `"accession"` or `"plot"` (Breedbase's "Analysis of").
#' @param traits Trait columns (default: all except the name and design
#'   columns). Their names are the Breedbase trait names.
#' @param accession,plot Names of the accession and plot-name columns.
#' @param estimate For the default statistic: `"blue"` or `"blup"` (adjusted
#'   means of that kind), or `"spatial"` (spatial correction, the default for
#'   `level = "plot"`). By default taken from the data
#'   ([second_stage_values()], [first_stage_values()]).
#' @param statistic Optional statistic term: its name, `"name|SGNSTAT:id"`
#'   or the id, as listed by [breedbase_statistics()]. It must be a term of
#'   the Breedbase statistics ontology (checked against the ontology file,
#'   with suggestions for misspelled names, and then on the server). By
#'   default chosen by `estimate`.
#' @param protocol The model, e.g. a formula: "SpATS(yield ~ genotype +
#'   PSANOVA(col, row))". Default: the description.
#' @param time_terms Optional time terms per trait, as in the form's "Compose
#'   Trait Info With Time Term", e.g. `c("fresh root yield|CO_334:0000013" =
#'   "day 300|TIME:0000300")`.
#' @param upload_file Where to write the upload file when the site's web
#'   pages are closed by a web password (see Details). Default: the analysis
#'   name + `"_upload.csv"` in the working folder.
#' @param dry_run Logical. Check everything but write nothing.
#' @param confirm Logical. Ask before writing (default in interactive
#'   sessions).
#' @param folder Optional folder of the program to put the analysis in, e.g.
#'   the folder of its trials. Without a folder (the default, `NULL`) the
#'   analysis is in the breeding program itself. An existing folder is used
#'   as it is; a new one (special characters replaced with `"_"`) is created when
#'   `create_folder = TRUE`.
#' @param create_folder Logical. Create `folder` when it does not exist.
#'
#' @section Sites with a web password:
#' Breedbase creates analyses through its web pages; BrAPI cannot create an
#' analysis project. Test sites such as cassava-test.sgn.cornell.edu close
#' their web pages with an extra web password. With it (see
#' [breedbase_login()]), the analysis is uploaded directly. Without it,
#' everything is checked through BrAPI and the file of the "Upload analysis"
#' form is written, with the value of every form field printed: upload it in
#' the browser and Breedbase creates the same analysis project.
#'
#' @section The stored trait names:
#' Breedbase stores the values of an analysis under the trait composed with
#' the statistic, "trait|statistic|COMP:id", e.g. "fresh root yield|Adjusted
#' Means from Spatial Correction using SpATS R|COMP:0000379", never under the
#' plain trait. Before saving (and in a dry run) each trait is listed with
#' the composed name it will get, and whether that composed trait already
#' exists or Breedbase will create it. To use another statistic, list them
#' with [breedbase_statistics()] and give one with `statistic`; run with
#' `dry_run = TRUE` to see the names again.
#'
#' @return Invisibly, a tibble with the analysis name and id, the number of
#'   rows (accessions or plots) sent and read back, and the analysis page; or,
#'   when the web pages are closed, a list with the upload `file` and the
#'   `form` values. A dry run returns the `file`, the `form` and `traits`:
#'   each trait with the composed name it will be stored under
#'   (`stored_as`) and `status` (`"exists"` or `"new"`).
#'
#' @examples
#' prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
#' d <- prep$design
#' set.seed(1)
#' d$yield <- rnorm(nrow(d), 25, 3)
#' fit <- analyze_augmented_prep_two_stage(d, "yield")
#' means <- second_stage_values(fit, trait = "fresh root yield")
#' sent <- save_breedbase_analysis(means, analysis_name = "PREP_2026_stage2", program = "IITA",
#'                                 description = "Two-stage analysis, weighted", dry_run = TRUE)
#' head(sent$file)
#'
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' # Second stage: one value per accession.
#' save_breedbase_analysis(means, con, analysis_name = "PREP_2026_stage2", program = "IITA",
#'                         description = "Two-stage analysis, weighted")
#' # First stage: one value per field plot ("Analysis of: plots").
#' stage1 <- first_stage_values(fit, prep, "PREP_2026", trait = "fresh root yield")
#' save_breedbase_analysis(stage1, con, analysis_name = "PREP_2026_stage1", program = "IITA",
#'                         description = "SpATS per location", level = "plot")
#' }
#'
#' @export
save_breedbase_analysis <- function(data, con = NULL, analysis_name, program, description,
                                    year = format(Sys.Date(), "%Y"), level = c("accession", "plot"),
                                    traits = NULL, accession = "accession", plot = "plot_name",
                                    estimate = NULL, statistic = NULL, protocol = NULL, time_terms = NULL,
                                    upload_file = NULL, dry_run = FALSE, confirm = interactive(),
                                    folder = NULL, create_folder = FALSE) {
  level <- match.arg(level)
  if (!is.null(folder) && (length(folder) != 1 || is.na(folder) || !nzchar(trimws(folder)))) {
    stop("`folder` must be one folder name, or NULL to keep the analysis in the program.", call. = FALSE)
  }
  if (missing(analysis_name) || length(analysis_name) != 1 || is.na(analysis_name) || !nzchar(trimws(analysis_name))) {
    stop("Give the `analysis_name` of the new analysis, e.g. \"PREP_2026_stage2\".", call. = FALSE)
  }
  analysis_name <- bb_clean_name(analysis_name, "analysis name")
  if (missing(program) || length(program) != 1 || !nzchar(program)) stop("Give the breeding `program`.", call. = FALSE)
  if (missing(description) || length(description) != 1 || is.na(description) || !nzchar(trimws(description))) {
    stop("Give a `description` of the analysis, e.g. the model.", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  estimate <- estimate %||% attr(data, "estimate") %||%
    if (level == "plot" && !identical(attr(data, "method"), "mixed")) "spatial" else "blue"
  estimate <- match.arg(estimate, c("blue", "blup", "spatial"))
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  key <- if (level == "plot") plot else accession
  if (!key %in% names(data)) {
    stop("Column '", key, "' not found in `data` (needed for level = \"", level, "\"). Columns are: ",
         paste(names(data), collapse = ", "), ".", call. = FALSE)
  }
  design_cols <- c("plot_name", "trial_name", "accession", "plot_number", "block", "rep", "row", "col", "is_check")
  if (is.null(traits)) traits <- setdiff(names(data), c(key, accession, plot, design_cols))
  if (!length(traits)) stop("No trait columns in `data`. Give them with `traits = c(...)`.", call. = FALSE)
  absent <- setdiff(traits, names(data))
  if (length(absent)) stop("Trait columns not found in `data`: ", paste(absent, collapse = ", "), call. = FALSE)
  for (t in traits) if (!is.numeric(data[[t]])) stop("Trait column '", t, "' must be numeric.", call. = FALSE)
  stocks <- trimws(as.character(data[[key]]))
  if (anyNA(stocks) || any(!nzchar(stocks))) stop("Some rows have no ", key, ".", call. = FALSE)
  if (anyDuplicated(stocks)) {
    stop("Each ", if (level == "plot") "plot" else "accession", " must appear once; repeated: ",
         paste(utils::head(unique(stocks[duplicated(stocks)]), 5), collapse = ", "), call. = FALSE)
  }
  if (all(vapply(traits, function(t) all(is.na(data[[t]])), logical(1)))) {
    stop("There are no values to save (all are empty).", call. = FALSE)
  }
  protocol <- protocol %||% description

  # The file of Breedbase's analysis upload: "stock_names" then one column per
  # trait (full Breedbase names); empty cells for missing values.
  build_file <- function(trait_full) {
    f <- data.frame(stock_names = stocks, stringsAsFactors = FALSE, check.names = FALSE)
    for (t in traits) f[[trait_full[[t]]]] <- ifelse(is.na(data[[t]]), "", as.character(data[[t]]))
    f
  }
  build_form <- function(trait_full, statistic_full) {
    compose <- if (length(time_terms)) {
      paste(paste(names(time_terms), unname(time_terms), sep = "||||"), collapse = ",")
    } else {
      ""
    }
    list(
      upload_new_analysis_name = analysis_name,
      upload_new_analysis_description = description,
      upload_new_analysis_year = as.character(year),
      upload_new_analysis_protocol = protocol,
      upload_new_analysis_dataset_id = "",
      upload_new_analysis_statistical_ontology_term = statistic_full,
      upload_new_analysis_result_values_type = if (level == "plot") {
        "analysis_result_values_match_precomputed_design"
      } else {
        "analysis_result_values_match_accession_names"
      },
      upload_new_analysis_result_summary_string = paste0("method:breedtools ", gsub("[,:]", " ", protocol)),
      upload_new_analysis_result_trait_compose_info_string = compose,
      upload_new_analysis_model_name = paste(analysis_name, "model"),
      upload_new_analysis_model_description = description,
      upload_new_analysis_model_is_public = "yes",
      upload_new_analysis_model_language = "R",
      upload_new_analysis_model_properties_string = paste0("protocol:", gsub("[,:]", " ", protocol)),
      upload_new_analysis_model_application_name = "breedtools",
      upload_new_analysis_model_application_version = as.character(utils::packageVersion("breedtools"))
    )
  }

  if (is.null(con)) {
    if (!isTRUE(dry_run)) stop("`con` is needed to save; use breedbase_login() first, or `dry_run = TRUE`.", call. = FALSE)
    message("Dry run without a connection: returning the upload that would be sent (trait names and the ",
            "statistic are completed from Breedbase when connected).")
    onto <- bb_choose_statistic(NULL, statistic, kind = estimate, offline = TRUE)
    return(list(file = build_file(stats::setNames(as.list(traits), traits)),
                form = build_form(stats::setNames(as.list(traits), traits), paste0(onto$name, "|", onto$accession))))
  }
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)

  bt_say("Checking the analysis before saving to ", con$url, "\n", sep = "")
  # Analyses are created by Breedbase's web pages (BrAPI cannot create an
  # analysis project). If the site's web pages are closed by a web password
  # we do not have, everything is checked through BrAPI and the upload file
  # for the "Upload analysis" form of the website is written instead.
  web_ok <- !identical(bb_web_access(con), "password")
  prog <- bb_find_program(con, program)
  bt_say("  [OK]   breeding program '", prog$program, "'\n", sep = "")
  problems <- character(0)
  existing_analysis <- if (web_ok) bb_find_analysis(con, analysis_name) else NA_character_
  if (!is.na(existing_analysis)) {
    problems <- c(problems, paste0(
      "An analysis named '", analysis_name, "' already exists (", con$url, "/analyses/", existing_analysis,
      "). It may be left from an earlier save that failed: open it and check its values. ",
      "Delete it on that page to save again with this name, add values with add_breedbase_analysis_values(), ",
      "or choose another `analysis_name`."))
  } else if (!is.na(bb_find_study(con, analysis_name))) {
    problems <- c(problems, paste0("A trial named '", analysis_name,
                                   "' already exists. Choose another `analysis_name`."))
  } else {
    bt_say("  [OK]   analysis name '", analysis_name, "' is free", if (!web_ok) " (as a trial; analyses are checked by the website)", "\n", sep = "")
  }
  folder_id <- NA_character_
  if (!is.null(folder)) {
    folders <- breedbase_folders(con, prog$program)
    if (!folder %in% folders$folder) folder <- bb_clean_name(folder, "folder name")
    folder_id <- bb_folder_id(con, prog$program, folder)
    if (!is.na(folder_id)) {
      bt_say("  [OK]   folder '", folder, "' exists: the analysis will be put inside it (",
          con$url, "/folder/", folder_id, ")\n", sep = "")
    } else if (isTRUE(create_folder)) {
      bt_say("  [NEW]  folder '", folder, "' will be created\n", sep = "")
    } else {
      problems <- c(problems, paste0("Folder '", folder, "' was not found in program '", prog$program,
                                     "'. Existing folders: ", paste(utils::head(folders$folder, 20), collapse = ", "),
                                     ". Use create_folder = TRUE to create it."))
    }
  }
  if (level == "plot") {
    missing_stocks <- setdiff(stocks, bb_existing_plots(con, stocks))
    what <- "plots"
  } else {
    missing_stocks <- setdiff(stocks, bb_existing_germplasm(con, stocks))
    what <- "accessions"
  }
  if (length(missing_stocks)) {
    problems <- c(problems, paste0(length(missing_stocks), " of ", length(stocks), " ", what,
                                   " are not in Breedbase: ", paste(utils::head(missing_stocks, 10), collapse = ", "),
                                   if (level == "plot") ". Use the plot names of the saved trials (first_stage_values())." else ""))
  } else {
    bt_say("  [OK]   all ", length(stocks), " ", what, " are in Breedbase\n", sep = "")
  }
  full <- attr(bb_find_variables(con, traits), "full_name")
  if (anyNA(full)) {
    problems <- c(problems, paste0("Trait(s) not found in Breedbase: ",
                                   paste0("'", traits[is.na(full)], "'", collapse = ", "),
                                   ". Use the trait name as shown in Breedbase."))
  } else {
    for (t in traits) bt_say("  [OK]   trait '", full[[t]], "'\n", sep = "")
  }
  if (length(problems)) stop(paste(problems, collapse = "\n"), call. = FALSE)
  if (!web_ok) {
    onto <- bb_choose_statistic(NULL, statistic, kind = estimate, offline = TRUE)
    statistic_full <- paste0(onto$name, "|", onto$accession)
    bt_say("  [OK]   statistic '", statistic_full, "' (checked in the Breedbase statistics ontology)\n", sep = "")
    bb_show_composed(con, traits, full, onto, web = FALSE)
    if (!is.null(folder)) bt_say("  [NOTE] after uploading, move the analysis to folder '", folder, "' on the website\n", sep = "")
    return(invisible(bb_analysis_upload_file(
      build_file(as.list(full)), analysis_name = analysis_name, description = description, protocol = protocol,
      level = level, program = prog$program, year = year, statistic = statistic_full,
      time_terms = time_terms, upload_file = upload_file, con = con
    )))
  }
  stat <- bb_choose_statistic(con, statistic, kind = estimate)
  statistic_full <- paste0(stat$name, "|", stat$accession)
  bt_say("  [OK]   statistic '", statistic_full, "' (from this Breedbase)\n", sep = "")
  stored_as <- bb_show_composed(con, traits, full, stat, web = TRUE)
  bt_say("  [OK]   analysis of ", what, ": ", length(stocks), " rows, ", length(traits), " trait",
      if (length(traits) > 1) "s" else "", "\n", sep = "")

  file <- build_file(as.list(full))
  form <- build_form(as.list(full), statistic_full)
  form$upload_new_analysis_breeding_program_id <- prog$program_db_id
  plan <- sprintf("analysis '%s' of %d %s in program '%s'%s at %s", analysis_name, length(stocks), what,
                  prog$program, if (is.null(folder)) "" else paste0(", folder '", folder, "',"), con$url)
  if (isTRUE(dry_run)) {
    bt_say("Dry run: everything is ready to save ", plan, ". Nothing was written.\n", sep = "")
    return(invisible(list(file = file, form = form, traits = stored_as)))
  }
  if (isTRUE(confirm)) {
    answer <- readline(sprintf("Save %s? Type 'yes' to continue: ", plan))
    if (!identical(tolower(trimws(answer)), "yes")) {
      message("Nothing was saved.")
      return(invisible(NULL))
    }
  }

  path <- tempfile("breedtools_analysis_", fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(file, path, row.names = FALSE, na = "")
  analysis_id <- bb_write_checked(
    what = paste0("save analysis '", analysis_name, "'"),
    write = function() {
      res <- bb_ajax_upload(con, "ajax/analysis/store/spreadsheet", form, file_field = "upload_new_analysis_file",
                            file = path, what = paste0("save analysis '", analysis_name, "'"))
      as.character(res$analysis_id %||% NA)
    },
    find = function() bb_find_analysis(con, analysis_name),
    attempts = 1, check_on_error = TRUE
  )
  server_error <- attr(analysis_id, "server_error")
  analysis_id <- as.character(analysis_id)
  if (!is.null(folder)) {
    if (is.na(folder_id)) folder_id <- bb_create_folder(con, prog, folder, saved = "the analysis")
    bb_move_to_folder(con, analysis_id, folder_id, what = paste0("put analysis '", analysis_name,
                                                                "' in folder '", folder, "'"))
    bt_say("  [OK]   analysis put in folder '", folder, "' (", con$url, "/folder/", folder_id, ")\n", sep = "")
  }
  n_back <- bb_analysis_row_count(con, analysis_id)
  page <- paste0(con$url, "/analyses/", analysis_id)
  status <- if (!is.na(n_back) && n_back >= length(stocks)) "[OK]  " else "[CHECK]"
  bt_say(sprintf("  %s analysis '%s' saved with %s of %d %s\n", status, analysis_name,
              if (is.na(n_back)) "?" else as.character(n_back), length(stocks), what))
  bt_say("Open the analysis in Breedbase: ", page, "\n", sep = "")
  if (!identical(status, "[OK]  ")) {
    warning("Analysis '", analysis_name, "' was created, but ",
            if (is.na(n_back)) "its values could not be read back" else paste(n_back, "of", length(stocks), what, "have values"),
            if (!is.null(server_error)) paste0(" (", server_error, ")") else "",
            ". Open ", page, " to check it. If values are missing, delete the analysis on that page and save again, ",
            "or add them with add_breedbase_analysis_values().", call. = FALSE)
  }
  invisible(tibble::tibble(analysis_name = analysis_name, analysis_db_id = analysis_id, level = level,
                           rows_sent = length(stocks), rows_found = n_back, analysis_page = page))
}

# Without access to the web pages: write the file of the "Upload analysis"
# form (checked through BrAPI) and print what to fill in the form.
bb_analysis_upload_file <- function(file, analysis_name, description, protocol, level, program, year, statistic,
                                    time_terms, upload_file, con) {
  path <- upload_file %||% file.path(getwd(), paste0(gsub("[^A-Za-z0-9_.-]+", "_", analysis_name), "_upload.csv"))
  utils::write.csv(file, path, row.names = FALSE, na = "")
  stat_hint <- statistic
  bt_say("  [NOTE] The web pages of ", con$url, " need their own web password, so the analysis cannot be sent from R.\n",
      "         Everything was checked through BrAPI. Upload this file in your browser instead:\n\n", sep = "")
  bt_say("  File:   ", normalizePath(path), "\n", sep = "")
  bt_say("  Where:  ", con$url, " > Analyze > Analyses > \"Upload analysis\" (or open ", con$url, "/analyses)\n\n",
      sep = "")
  rows <- c(
    "Analysis Name" = analysis_name,
    "Analysis Description" = description,
    "Analysis Protocol" = protocol,
    "Analysis Of" = if (level == "plot") "plots (results per field plot)" else "accessions",
    "Breeding Program" = program,
    "Year" = as.character(year),
    "Statistical Ontology Term (SGNSTAT)" = stat_hint,
    "Compose Trait Info With Time Term" = if (length(time_terms)) paste(paste(names(time_terms), unname(time_terms), sep = "||||"), collapse = ",") else "(leave empty)",
    "File" = basename(path)
  )
  for (i in seq_along(rows)) bt_say(sprintf("  %-38s %s\n", paste0(names(rows)[i], ":"), rows[[i]]))
  bt_say("\n  The file has \"stock_names\" (", if (level == "plot") "field plot names" else "accession names",
      ") and one column per trait with its full Breedbase name.\n", sep = "")
  list(file = path, form = as.list(rows))
}

#' Add results to an analysis that is already in Breedbase
#'
#' Adds values (for example a second trait, or the same trait with another
#' statistic) to the plots of an analysis saved before, e.g. with
#' [save_breedbase_analysis()]. Breedbase names the plots of an analysis
#' `"<analysis_name>_<field plot name>"` ("Analysis of: plots") or
#' `"<analysis_name>_<accession>"` ("Analysis of: accessions"); the values
#' are uploaded to those plots with the phenotype upload of the website
#' (simple spreadsheet), which first verifies the file and then stores it.
#'
#' Values already present for a plot and trait are not replaced unless
#' `overwrite = TRUE`: the verification lists them and the upload stops.
#'
#' @param data A data frame with the field plot names (`level = "plot"`) or
#'   the accession names (`level = "accession"`) and one numeric column per
#'   trait, e.g. from [first_stage_values()] or [second_stage_values()].
#' @param con A connection from [breedbase_login()] (the site's web pages
#'   must be open; see [breedbase_login()]).
#' @param analysis_name Name of the existing analysis.
#' @param level `"plot"` or `"accession"`, as the analysis was saved.
#' @param traits Trait columns (default: all except the name and design
#'   columns). Their names are Breedbase trait names, including composed
#'   traits such as `"dry matter content percentage|Adjusted Means from
#'   Spatial Correction using SpATS R"` (see [breedbase_adjusted_trait()]).
#' @param plot,accession Names of the plot-name and accession columns.
#' @param statistic Optional statistic (see [breedbase_statistics()]). Values
#'   are always stored under "trait|statistic" (a composed trait, created if
#'   needed), never under the plain trait, as Breedbase does for analyses.
#'   By default the statistic the analysis already uses (from its stored
#'   traits); give another to keep, e.g., single-stage means next to
#'   two-stage means in the same analysis.
#' @param estimate When the analysis has no statistic yet: `"blue"`,
#'   `"blup"` or `"spatial"`, as in [save_breedbase_analysis()] (by default
#'   from the data).
#' @param overwrite Logical. Replace values already stored.
#' @param dry_run Logical. Only verify (Breedbase checks the file), store
#'   nothing.
#' @param confirm Logical. Ask before storing (default in interactive
#'   sessions).
#'
#' @return Invisibly, a list with the analysis name, the number of plots and
#'   values sent, and Breedbase's messages.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassava-test.sgn.cornell.edu", "jdoe",
#'                        site_username = "web_user_name")
#' # dm1: first-stage values of dry matter, one row per field plot
#' add_breedbase_analysis_values(dm1, con, analysis_name = "FALL_2026_stage1", level = "plot")
#' }
#'
#' @export
add_breedbase_analysis_values <- function(data, con, analysis_name, level = c("plot", "accession"), traits = NULL,
                                          plot = "plot_name", accession = "accession", statistic = NULL,
                                          estimate = NULL, overwrite = FALSE, dry_run = FALSE,
                                          confirm = interactive()) {
  level <- match.arg(level)
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  if (missing(analysis_name) || length(analysis_name) != 1 || !nzchar(analysis_name)) {
    stop("Give the `analysis_name` of the existing analysis.", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  key <- if (level == "plot") plot else accession
  if (!key %in% names(data)) stop("Column '", key, "' not found in `data`.", call. = FALSE)
  design_cols <- c("plot_name", "trial_name", "accession", "plot_number", "block", "rep", "row", "col", "is_check")
  if (is.null(traits)) traits <- setdiff(names(data), c(key, accession, plot, design_cols))
  if (!length(traits)) stop("No trait columns in `data`.", call. = FALSE)
  for (t in traits) if (!is.numeric(data[[t]])) stop("Trait column '", t, "' must be numeric.", call. = FALSE)
  stocks <- trimws(as.character(data[[key]]))
  if (anyDuplicated(stocks)) stop("Each ", key, " must appear once.", call. = FALSE)

  bt_say("Adding values to analysis '", analysis_name, "' at ", con$url, "\n", sep = "")
  if (identical(bb_web_access(con), "password")) {
    stop("This site's web pages need their own web password; log in with site_username = \"...\" ",
         "(see ?breedbase_login).", call. = FALSE)
  }
  # Analyses saved by breedtools have no spaces in their names; one made on
  # the website may have them, so the name as given is tried as well.
  id <- bb_find_analysis(con, bb_clean_name(analysis_name, quiet = TRUE))
  if (!is.na(id)) {
    analysis_name <- bb_clean_name(analysis_name, quiet = TRUE)
  } else {
    id <- bb_find_analysis(con, analysis_name)
  }
  if (is.na(id)) stop("Analysis '", analysis_name, "' was not found among your analyses.", call. = FALSE)
  bt_say("  [OK]   analysis '", analysis_name, "' (", con$url, "/analyses/", id, ")\n", sep = "")
  # Analysis values are never stored under a plain trait: each trait is
  # stored as "trait|statistic|COMP:id" (a composed trait, created if needed),
  # as Breedbase does when it saves an analysis. By default the statistic is
  # the one the analysis already uses, so the project stays consistent;
  # otherwise it is chosen as in save_breedbase_analysis().
  if (is.null(statistic)) {
    statistic <- bb_analysis_statistic(con, id)
    if (!is.null(statistic)) {
      bt_say("  [OK]   statistic '", statistic, "' (the one of this analysis)\n", sep = "")
    } else {
      estimate <- estimate %||% attr(data, "estimate") %||%
        if (level == "plot" && !identical(attr(data, "method"), "mixed")) "spatial" else "blue"
      estimate <- match.arg(estimate, c("blue", "blup", "spatial"))
      stat <- bb_choose_statistic(con, NULL, kind = estimate)
      statistic <- paste0(stat$name, "|", stat$accession)
      bt_say("  [OK]   statistic '", statistic, "'\n", sep = "")
    }
  }
  plain <- traits[!grepl("\\|COMP:[0-9]+$", traits)]        # composed traits given are kept
  if (length(plain)) {
    adj <- breedbase_adjusted_trait(con, plain, statistic = statistic, create = !isTRUE(dry_run))
    if (any(adj$status == "missing")) {
      new <- adj$adjusted_trait[adj$status == "missing"]
      if (!isTRUE(dry_run)) {
        stop("Could not create the trait(s) ", paste0("'", new, "'", collapse = ", "), ".", call. = FALSE)
      }
      # A dry run creates nothing, so Breedbase cannot verify values of a
      # trait that does not exist yet.
      for (n in new) bt_say("  [NEW]  trait '", n, "' will be created when storing\n", sep = "")
      bt_say("Dry run: values of new traits are verified by Breedbase when storing. Nothing was stored.\n")
      return(invisible(list(analysis_name = analysis_name, plots = length(stocks), values = NA_integer_,
                            new_traits = new)))
    }
    composed <- vapply(adj$adjusted_trait, function(n) bb_find_variable_full(con, n), character(1))
    for (i in seq_along(plain)) bt_say("  [OK]   ", plain[i], " -> '", composed[i], "' (", adj$status[i], ")\n", sep = "")
    names(data)[match(plain, names(data))] <- unname(composed)
    traits[match(plain, traits)] <- unname(composed)
  }
  # Trait names as Breedbase knows them (composed traits by their exact name).
  full <- vapply(traits, function(t) {
    f <- attr(bb_find_variables(con, t), "full_name")[[1]]
    if (is.na(f)) {
      cid <- bb_find_variable_exact(con, sub("\\|[A-Za-z][A-Za-z0-9_]*:[0-9]+$", "", t))
      if (!is.na(cid)) f <- t
    }
    f
  }, character(1))
  if (anyNA(full)) {
    stop("Trait(s) not found in Breedbase: ", paste0("'", traits[is.na(full)], "'", collapse = ", "),
         ". For first-stage values use a trait for adjusted values (breedbase_adjusted_trait()).", call. = FALSE)
  }
  for (t in traits) bt_say("  [OK]   trait '", full[[t]], "'\n", sep = "")
  units <- paste0(analysis_name, "_", stocks)
  f <- data.frame(observationunit_name = units, stringsAsFactors = FALSE, check.names = FALSE)
  for (t in traits) f[[full[[t]]]] <- ifelse(is.na(data[[t]]), "", as.character(data[[t]]))
  keep <- rowSums(f[, -1, drop = FALSE] != "") > 0
  f <- f[keep, , drop = FALSE]
  n_values <- sum(f[, -1, drop = FALSE] != "")
  path <- tempfile("breedtools_analysis_values_", fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(f, path, row.names = FALSE, na = "")
  fields <- list(upload_spreadsheet_phenotype_file_format = "simple",
                 upload_spreadsheet_phenotype_data_level = "plots",
                 phenotype_upload_overwrite_values = if (isTRUE(overwrite)) "1" else "")
  if (!isTRUE(overwrite)) fields$phenotype_upload_overwrite_values <- NULL

  # 1. Breedbase verifies the file (plot and trait names, existing values).
  verified <- bb_phenotype_upload(con, "verify", fields, path)
  if (length(verified$warning) && !isTRUE(overwrite)) {
    stop("Breedbase found values already stored for these plots and traits (", length(verified$warning),
         " messages), e.g. ", gsub("<[^>]+>", "", verified$warning[[1]]),
         ". Use another trait, or overwrite = TRUE to replace them.", call. = FALSE)
  }
  bt_say("  [OK]   Breedbase verified ", nrow(f), " analysis plots and ", n_values, " values\n", sep = "")
  plan <- sprintf("%d values on %d plots of analysis '%s'", n_values, nrow(f), analysis_name)
  if (isTRUE(dry_run)) {
    bt_say("Dry run: ready to add ", plan, ". Nothing was stored.\n", sep = "")
    return(invisible(list(analysis_name = analysis_name, plots = nrow(f), values = n_values, verified = verified)))
  }
  if (isTRUE(confirm)) {
    answer <- readline(sprintf("Add %s? Type 'yes' to continue: ", plan))
    if (!identical(tolower(trimws(answer)), "yes")) {
      message("Nothing was stored.")
      return(invisible(NULL))
    }
  }
  # 2. Breedbase stores the values (and refreshes its phenotype cache).
  stored <- bb_phenotype_upload(con, "store", fields, path)
  bt_say("  [OK]   added ", plan, "\n", sep = "")
  bt_say("Open the analysis in Breedbase: ", con$url, "/analyses/", id, "\n", sep = "")
  invisible(list(analysis_name = analysis_name, analysis_db_id = id, plots = nrow(f), values = n_values,
                 messages = unlist(stored$success)))
}

# The phenotype spreadsheet upload of the website: "verify" or "store".
# Breedbase answers {success: [...], warning: [...], error: [...]}.
bb_phenotype_upload <- function(con, step, fields, file) {
  if (!requireNamespace("curl", quietly = TRUE)) stop("Uploading needs the 'curl' package.", call. = FALSE)
  parts <- c(lapply(fields, as.character),
             list(upload_spreadsheet_phenotype_file_input = curl::form_file(file, type = "text/csv")))
  req <- do.call(httr2::req_body_multipart,
                 c(list(bb_ajax_request(con, paste0("ajax/phenotype/upload_", step, "/spreadsheet"))), parts))
  body <- bb_perform(req, what = paste(step, "the values"))
  errors <- unlist(body$error)
  if (length(errors)) {
    stop("Breedbase could not ", step, " the values: ", paste(gsub("<[^>]+>", "", errors), collapse = "; "),
         call. = FALSE)
  }
  body
}

# Upload a form with a file to a Breedbase web page (multipart, like the
# browser), with the login as sgn_session_id. Errors come back as
# {"error": ...} or {"error_messages": [...]} with HTTP 200.
bb_ajax_upload <- function(con, path, fields, file_field, file, what) {
  bb_require_httr2()
  if (!requireNamespace("curl", quietly = TRUE)) stop("Uploading needs the 'curl' package.", call. = FALSE)
  parts <- c(lapply(fields, as.character), list(sgn_session_id = con$token))
  parts[[file_field]] <- curl::form_file(file, type = "text/csv")
  req <- do.call(httr2::req_body_multipart, c(list(bb_ajax_request(con, path)), parts))
  body <- bb_perform(req, what = what)
  errors <- c(unlist(body$error), unlist(body$error_messages))
  if (length(errors)) bb_stop(paste0("Could not ", what, ": ", paste(errors, collapse = "; ")), code = 400L)
  body
}

# POST to a Breedbase (non-BrAPI) AJAX endpoint as a form, like the web page
# does: lists and objects are sent as JSON text, and the login token as
# sgn_session_id. Errors come back as {"error": "..."} with HTTP 200.
bb_ajax_post <- function(con, path, params, what, json = character(0)) {
  bb_require_httr2()
  form <- stats::setNames(lapply(names(params), function(nm) {
    v <- params[[nm]]
    if (nm %in% json) {
      # Character vectors stay JSON arrays even with one element.
      if (is.list(v)) as.character(jsonlite::toJSON(v, auto_unbox = TRUE, null = "null", digits = NA))
      else as.character(jsonlite::toJSON(as.character(v)))
    } else {
      as.character(v)
    }
  }), names(params))
  form$sgn_session_id <- con$token
  req <- bb_ajax_request(con, path)
  req <- do.call(httr2::req_body_form, c(list(req), form))
  body <- bb_perform(req, what = what)
  if (!is.null(body$error)) {
    bb_stop(paste0("Could not ", what, ": ", paste(unlist(body$error), collapse = "; ")), code = 400L)
  }
  body
}

# A request to a Breedbase web (non-BrAPI) page: the login as the session
# cookie of the website, and the site's web password when it has one.
bb_ajax_request <- function(con, path) {
  req <- httr2::request(paste0(con$url, "/", path))
  req <- httr2::req_headers(req, Cookie = paste0("sgn_session_id=", con$token))
  if (!is.null(con$site_auth)) req <- httr2::req_auth_basic(req, con$site_auth$username, con$site_auth$password)
  req
}

bb_ajax_get <- function(con, path, what, query = list(), allow_error = FALSE) {
  bb_require_httr2()
  req <- bb_ajax_request(con, path)
  req <- do.call(httr2::req_url_query, c(list(req), query, list(sgn_session_id = con$token, .multi = "explode")))
  body <- bb_perform(req, what = what)
  if (!isTRUE(allow_error) && is.list(body) && !is.null(names(body)) && !is.null(body$error)) {
    bb_stop(paste0("Could not ", what, ": ", paste(unlist(body$error), collapse = "; ")), code = 400L)
  }
  body
}

# Put a project (trial or analysis) in a folder, as the website's "Move to
# folder" does: only the old folder link is replaced, the breeding program
# is kept.
bb_move_to_folder <- function(con, project_id, folder_id, what) {
  bb_ajax_get(con, paste0("ajax/folder/", project_id, "/associate/parent/", folder_id), what = what)
  invisible(TRUE)
}

# Id of the user's analysis with this name (from the Analyses list), or NA.
bb_find_analysis <- function(con, name) {
  body <- tryCatch(bb_ajax_get(con, "ajax/analyses/by_user", what = "list your analyses"),
                   breedbase_error = function(e) NULL)
  if (is.null(body)) return(NA_character_)
  for (row in body$data) {
    link <- as.character(row[[1]])
    label <- sub("^.*>([^<]*)</a>.*$", "\\1", link)
    if (identical(label, name)) return(sub("^.*/analyses/([0-9]+).*$", "\\1", link))
  }
  NA_character_
}

# Entries of a saved analysis, read back from Breedbase.
bb_analysis_row_count <- function(con, analysis_id) {
  body <- tryCatch(bb_ajax_get(con, paste0("ajax/analysis/", analysis_id, "/retrieve"),
                               what = "read the saved analysis"),
                   breedbase_error = function(e) NULL)
  if (is.null(body) || is.null(body$data)) return(NA_integer_)
  max(length(body$data) - 1L, 0L)
}

#' Plot and accession names an analysis will have, before saving
#'
#' Shows, without contacting Breedbase, the names under which
#' [save_breedbase_analysis()] (or [add_breedbase_analysis_values()]) stores
#' each row of `data`: Breedbase gives every analysis its own plots, named
#' `"<analysis_name>_<field plot name>"` ("Analysis of: plots", with the
#' accession of the field plot) or `"<analysis_name>_<accession>"`
#' ("Analysis of: accessions"). Use it to see what will be saved, and to
#' compare with the plots of a saved analysis
#' ([breedbase_analysis_plots()]).
#'
#' @param data A data frame with one row per field plot (`level = "plot"`,
#'   e.g. from [first_stage_values()]) or per accession (`level =
#'   "accession"`, e.g. from [second_stage_values()]) and one column per
#'   trait.
#' @param analysis_name The analysis name, as given to
#'   [save_breedbase_analysis()] (special characters are replaced with `"_"`
#'   the same way).
#' @param level `"plot"` or `"accession"`.
#' @param plot,accession Names of the plot-name and accession columns.
#' @param traits Trait columns (default: all except the name and design
#'   columns).
#'
#' @return A tibble with `analysis_plot_name` (the name in the analysis),
#'   `plot_name` (the field plot; for accession analyses, the accession),
#'   `accession`, and the trait values.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' set.seed(1)
#' d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
#' fit <- analyze_augmented_row_column_design(d, "yield")
#' values <- first_stage_values(fit, rc, "AYT_2026", trait = "fresh root yield", location = "Ibadan")
#' head(analysis_plot_names(values, "AYT_2026_stage1", level = "plot"))
#'
#' @export
analysis_plot_names <- function(data, analysis_name, level = c("plot", "accession"), plot = "plot_name",
                                accession = "accession", traits = NULL) {
  level <- match.arg(level)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (missing(analysis_name) || length(analysis_name) != 1 || is.na(analysis_name) || !nzchar(trimws(analysis_name))) {
    stop("Give the `analysis_name`.", call. = FALSE)
  }
  analysis_name <- bb_clean_name(analysis_name, "analysis name", quiet = TRUE)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  key <- if (level == "plot") plot else accession
  if (!key %in% names(data)) stop("Column '", key, "' not found in `data`.", call. = FALSE)
  design_cols <- c("plot_name", "trial_name", "accession", "plot_number", "block", "rep", "row", "col", "is_check")
  if (is.null(traits)) traits <- setdiff(names(data), c(key, accession, plot, design_cols))
  stocks <- trimws(as.character(data[[key]]))
  out <- tibble::tibble(
    analysis_plot_name = paste0(analysis_name, "_", stocks),
    plot_name = stocks,
    accession = if (accession %in% names(data)) as.character(data[[accession]]) else NA_character_
  )
  for (t in traits) out[[t]] <- data[[t]]
  out
}

#' Plots of an analysis saved in Breedbase
#'
#' Reads back the plots of an existing analysis (from
#' [save_breedbase_analysis()] or made on the website): the analysis plot
#' names, the field plot (or accession) each one stands for, its accession
#' and design, and the values of every trait stored in the analysis. Use it
#' to see what an analysis holds before adding values with
#' [add_breedbase_analysis_values()], e.g. a second trait.
#'
#' Analysis plots are not part of BrAPI, so they are read through the web
#' pages of Breedbase (the site's web password is needed when it has one;
#' see [breedbase_login()]).
#'
#' Breedbase lists the plots of an analysis through their values: a plot
#' without any value (e.g. a lost plot, or an outlier removed before the
#' analysis) is not listed, although it is part of the analysis.
#'
#' @param con A connection from [breedbase_login()].
#' @param analysis_name Name of the analysis.
#'
#' @return A tibble with `analysis_plot_name`, `plot_name` (the field plot,
#'   or the accession for accession analyses: the name to use in the `plot`
#'   or `accession` column of [add_breedbase_analysis_values()]),
#'   `accession`, `plot_number`, `rep`, `block`, `row`, `col` and one column
#'   per trait, ordered by plot number. The analysis id is in the attribute
#'   `"analysis_db_id"`.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassava-test.sgn.cornell.edu", "jdoe",
#'                        site_username = "web_user_name")
#' plots <- breedbase_analysis_plots(con, "AYT_2026_stage1")
#' head(plots)
#' }
#'
#' @export
breedbase_analysis_plots <- function(con, analysis_name) {
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  if (identical(bb_web_access(con), "password")) {
    stop("This site's web pages need their own web password; log in with site_username = \"...\" ",
         "(see ?breedbase_login).", call. = FALSE)
  }
  id <- bb_find_analysis(con, bb_clean_name(analysis_name, quiet = TRUE))
  if (!is.na(id)) {
    analysis_name <- bb_clean_name(analysis_name, quiet = TRUE)
  } else {
    id <- bb_find_analysis(con, analysis_name)
  }
  if (is.na(id)) stop("Analysis '", analysis_name, "' was not found among your analyses.", call. = FALSE)
  trait_ids <- bb_analysis_traits(con, id)
  if (!length(trait_ids)) {
    info <- bb_ajax_get(con, paste0("ajax/analysis/", id, "/retrieve"), what = "read the analysis")
    for (t in unique(as.character(unlist(info$traits)))) {
      tid <- unname(bb_find_variables(con, t))
      if (is.na(tid)) tid <- bb_find_variable_exact(con, sub("\\|[A-Za-z][A-Za-z0-9_]*:[0-9]+$", "", t))
      if (!is.na(tid)) trait_ids[[t]] <- tid
    }
  }
  if (!length(trait_ids)) stop("Analysis '", analysis_name, "' has no trait values.", call. = FALSE)
  traits <- names(trait_ids)
  rows <- list()
  for (t in traits) {
    tid <- trait_ids[[t]]
    body <- bb_ajax_get(con, paste0("ajax/breeders/trial/", id, "/trait_phenotypes"),
                        query = list(display = "analysis_instance", trait = tid),
                        what = paste0("read the values of '", t, "' in the analysis"))
    m <- body$data
    if (length(m) < 2) next
    header <- vapply(m[[1]], function(x) as.character(x %||% NA), character(1))
    col <- function(name) match(name, header)
    tcol <- which(header == t)
    if (!length(tcol)) tcol <- length(header)            # the trait is the last column
    cell <- function(r, j) {
      v <- if (is.na(j) || j > length(r)) NULL else unlist(r[[j]])
      if (!length(v) || is.na(v[1]) || !nzchar(as.character(v[1]))) NA_character_ else as.character(v[1])
    }
    for (r in m[-1]) {
      unit <- cell(r, col("observationUnitName"))
      if (is.null(rows[[unit]])) {
        rows[[unit]] <- list(analysis_plot_name = unit, accession = cell(r, col("germplasmName")),
                             plot_number = cell(r, col("plotNumber")), rep = cell(r, col("replicate")),
                             block = cell(r, col("blockNumber")), row = cell(r, col("rowNumber")),
                             col = cell(r, col("colNumber")))
      }
      rows[[unit]][[t]] <- cell(r, tcol[1])
    }
  }
  if (!length(rows)) stop("No plots with values were found in analysis '", analysis_name, "'.", call. = FALSE)
  num <- function(x) suppressWarnings(as.numeric(x))
  get <- function(field) vapply(rows, function(x) x[[field]] %||% NA_character_, character(1), USE.NAMES = FALSE)
  units <- get("analysis_plot_name")
  prefix <- paste0(analysis_name, "_")
  out <- tibble::tibble(
    analysis_plot_name = units,
    plot_name = ifelse(startsWith(units, prefix), substring(units, nchar(prefix) + 1), units),
    accession = get("accession"),
    plot_number = as.integer(num(get("plot_number"))),
    rep = as.integer(num(get("rep"))),
    block = as.integer(num(get("block"))),
    row = as.integer(num(get("row"))),
    col = as.integer(num(get("col")))
  )
  for (t in traits) out[[t]] <- num(get(t))
  out <- out[order(out$plot_number, out$analysis_plot_name), , drop = FALSE]
  attr(out, "analysis_db_id") <- id
  out
}

# Traits stored in an analysis, as the trial page lists them ([trait id,
# trait name, ...] per trait): ids named by full trait name.
bb_analysis_traits <- function(con, id) {
  assayed <- bb_ajax_get(con, paste0("ajax/breeders/trial/", id, "/traits_assayed"),
                         what = "list the traits of the analysis")$traits_assayed
  while (length(assayed) == 1 && is.list(assayed[[1]]) && length(assayed[[1]]) && is.list(assayed[[1]][[1]])) {
    assayed <- assayed[[1]]                          # Breedbase wraps the list once more
  }
  stats::setNames(
    vapply(assayed, function(x) as.character(unlist(x[1])), character(1)),
    vapply(assayed, function(x) as.character(unlist(x[2])), character(1)))
}

# The statistic of an analysis, from its composed traits
# ("trait|statistic|COMP:id"); NULL when it has none.
bb_analysis_statistic <- function(con, id) {
  traits <- names(bb_analysis_traits(con, id))
  composed <- grep("\\|COMP:[0-9]+$", traits, value = TRUE)
  stats_found <- unique(sub("^[^|]*\\|", "", sub("\\|COMP:[0-9]+$", "", composed)))
  stats_found <- stats_found[nzchar(stats_found)]
  ok <- vapply(stats_found, function(st) {
    !inherits(try(bb_check_statistic_ontology(st), silent = TRUE), "try-error")
  }, logical(1))
  stats_found <- stats_found[ok]
  if (!length(stats_found)) return(NULL)
  if (length(stats_found) > 1) {
    bt_say("  [NOTE] the analysis has traits with ", length(stats_found), " statistics; using '", stats_found[1],
        "' (choose another with `statistic = ...`)\n", sep = "")
  }
  stats_found[1]
}

# Show (and return) the composed name each trait will be stored under in an
# analysis: "trait|statistic|COMP:id" when it exists, "trait|statistic" when
# Breedbase will create it. `full` are the traits' full names, `stat` the
# statistic (name, accession, and the server id when the web pages are open).
bb_show_composed <- function(con, traits, full, stat, web = TRUE) {
  base <- sub("\\|[^|]*$", "", unname(full[traits]))
  composed <- paste(base, stat$name, sep = "|")
  ids <- if (web) unname(bb_find_variables(con, traits)) else rep(NA_character_, length(traits))
  stored <- vapply(seq_along(traits), function(i) {
    hit <- if (web && !is.na(ids[i]) && !is.na(stat$id %||% NA)) {
      cid <- bb_composed_trait_id(con, c(ids[i], stat$id))
      if (is.na(cid)) NA_character_ else bb_find_variable_full(con, composed[i])
    } else {
      bb_find_variable_full(con, composed[i])
    }
    if (is.na(hit)) NA_character_ else hit
  }, character(1))
  status <- ifelse(is.na(stored), "new", "exists")
  stored[is.na(stored)] <- composed[is.na(stored)]
  for (i in seq_along(traits)) {
    bt_say("  [OK]   '", traits[i], "' will be stored as '", stored[i], "'",
        if (status[i] == "new") " (new composed trait: Breedbase creates it)" else "", "\n", sep = "")
  }
  bt_say("         Another statistic? See breedbase_statistics(con) and give statistic = \"...\"\n")
  tibble::tibble(trait = traits, stored_as = stored, status = status)
}
