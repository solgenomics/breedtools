# ============================================================================
# Save designs to Breedbase through BrAPI v2
#
# Breedbase vocabulary vs BrAPI:
#   Breedbase breeding program = BrAPI program
#   Breedbase folder           = BrAPI trial
#   Breedbase trial            = BrAPI study
#   Breedbase plot             = BrAPI observation unit (level "plot")
# ============================================================================

#' Log in to a Breedbase server
#'
#' Logs in through the BrAPI token call and returns a connection used by the
#' other `breedbase_*()` functions and by [save_breedbase()]. Only the access
#' token is kept; the password is not stored.
#'
#' Never write your password in a script. Leave `password` empty: it is read
#' from the environment variable `BREEDBASE_PASSWORD` if set (for example in
#' your `.Renviron` file), and otherwise you are asked for it in a pop-up.
#'
#' @param url Address of the Breedbase site, e.g.
#'   `"https://cassavabase.org"`.
#' @param username Your Breedbase username.
#' @param password Optional password; see Details.
#' @param site_username,site_password Only for sites with an extra web
#'   password in front of the website (test servers such as
#'   cassava-test.sgn.cornell.edu ask for it in a browser pop-up, "Restricted
#'   Test Site"): that user name and password. They are used for the non-BrAPI
#'   calls (trait composer, analyses); BrAPI does not need them. Usually not
#'   needed: when a site asks for a web password, the login first tries your
#'   Breedbase user name and password there, and uses them if they work.
#'   Leave `site_password` empty to read it from `BREEDBASE_SITE_PASSWORD` or
#'   to be asked in a pop-up.
#'
#' @return A `breedbase_connection` object.
#'
#' @seealso [save_breedbase()], [breedbase_programs()], [breedbase_folders()],
#'   [breedbase_locations()]
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", username = "jdoe")
#' con
#' }
#'
#' @export
breedbase_login <- function(url, username, password = NULL, site_username = NULL, site_password = NULL) {
  bb_require_httr2()
  url <- sub("/+$", "", as.character(url))
  if (!grepl("^https?://", url)) url <- paste0("https://", url)
  if (is.null(password) || !nzchar(password)) password <- Sys.getenv("BREEDBASE_PASSWORD")
  if (!nzchar(password)) {
    if (!interactive()) {
      stop("No password: set the BREEDBASE_PASSWORD environment variable or pass `password`.", call. = FALSE)
    }
    password <- if (requireNamespace("askpass", quietly = TRUE)) {
      askpass::askpass(sprintf("Breedbase password for %s at %s: ", username, url))
    } else {
      readline(sprintf("Breedbase password for %s: ", username))
    }
  }
  if (is.null(password) || !nzchar(password)) stop("No password given.", call. = FALSE)
  site_auth <- NULL
  if (!is.null(site_username) && nzchar(site_username)) {
    if (is.null(site_password) || !nzchar(site_password)) site_password <- Sys.getenv("BREEDBASE_SITE_PASSWORD")
    if (!nzchar(site_password)) {
      if (!interactive()) stop("No site password: set BREEDBASE_SITE_PASSWORD or pass `site_password`.", call. = FALSE)
      site_password <- if (requireNamespace("askpass", quietly = TRUE)) {
        askpass::askpass(sprintf("Web (site) password for %s at %s: ", site_username, url))
      } else {
        readline(sprintf("Web (site) password for %s: ", site_username))
      }
    }
    site_auth <- list(username = site_username, password = site_password)
  }

  # Form-encoded body: newer Breedbase versions reject a JSON object here
  # ("JSON array body required"); every version accepts the form.
  resp <- bb_perform(
    httr2::req_body_form(
      httr2::request(paste0(url, "/brapi/v2/token")),
      username = username, password = password, grant_type = "password"
    ),
    what = "log in"
  )
  token <- resp$access_token
  if (!is.character(token) || length(token) != 1 || is.na(token) || !nzchar(token)) {
    stop("Login failed for user '", username, "' at ", url, ": ", bb_status_message(resp),
         call. = FALSE)
  }
  con <- structure(
    list(url = url, username = username, token = token,
         user_display_name = resp$userDisplayName %||% username,
         logged_in_at = Sys.time(), site_auth = site_auth),
    class = "breedbase_connection"
  )
  message("Logged in to ", url, " as ", con$user_display_name, ".")
  # Some sites (test servers) protect their web pages with a web password;
  # BrAPI does not need it, but the trait composer, analyses and cache
  # refresh do. Often it is the Breedbase login itself: try that.
  if (is.null(site_auth)) {
    web <- bb_web_access(con)
    if (identical(web, "password")) {
      with_login <- con
      with_login$site_auth <- list(username = username, password = password)
      if (identical(bb_web_access(with_login), "ok")) {
        con <- with_login
        message("This site's web pages ask for a password: your Breedbase login works there and will be used.")
      } else if (interactive()) {
        # Ask now, so the analyses can be uploaded directly later.
        message("This site's web pages ask for their own web password (the browser pop-up \"",
                bb_web_realm(con), "\"), not your Breedbase login.")
        site_username <- readline("Web user name of the site (leave empty to skip): ")
        if (nzchar(trimws(site_username))) {
          pw <- Sys.getenv("BREEDBASE_SITE_PASSWORD")
          if (!nzchar(pw)) {
            pw <- if (requireNamespace("askpass", quietly = TRUE)) {
              askpass::askpass(sprintf("Web password for %s at %s: ", site_username, url))
            } else {
              readline("Web password: ")
            }
          }
          with_site <- con
          with_site$site_auth <- list(username = trimws(site_username), password = pw)
          if (identical(bb_web_access(with_site), "ok")) {
            con <- with_site
            message("Web password accepted: analyses will be uploaded directly.")
          } else {
            message("That web password was not accepted. Trials and phenotypes still work; analyses will be ",
                    "written as files for the \"Upload analysis\" form of the website.")
          }
        } else {
          message("Skipped. Trials and phenotypes work; analyses will be written as files for the ",
                  "\"Upload analysis\" form of the website.")
        }
      } else {
        message("Note: this site's web pages ask for their own web password (not your Breedbase login). ",
                "Trials and phenotypes work without it; analyses will be written as files for the ",
                "\"Upload analysis\" form of the website, or log in with site_username = \"...\".")
      }
    }
  }
  con
}

# The name the site's web password pop-up shows.
bb_web_realm <- function(con) {
  req <- httr2::req_error(httr2::request(paste0(con$url, "/ajax/breeder/check_status")), is_error = function(resp) FALSE)
  resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
  basic <- if (is.null(resp)) "" else httr2::resp_header(resp, "WWW-Authenticate") %||% ""
  if (grepl("realm", basic)) sub('^.*realm="?([^"]*)"?.*$', "\\1", basic) else "web password"
}

# Can the web (non-BrAPI) pages be reached? "ok", "password" (a web password
# is asked) or "unknown" (no answer).
bb_web_access <- function(con) {
  req <- httr2::req_error(bb_ajax_request(con, "ajax/breeder/check_status"), is_error = function(resp) FALSE)
  req <- httr2::req_timeout(httr2::req_user_agent(req, "breedtools (R package)"), 60)
  resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
  if (is.null(resp)) return("unknown")
  basic <- httr2::resp_header(resp, "WWW-Authenticate") %||% ""
  if (httr2::resp_status(resp) == 401 && grepl("^Basic", basic)) "password" else "ok"
}

#' @export
print.breedbase_connection <- function(x, ...) {
  cat("Breedbase connection\n")
  cat("  site: ", x$url, "\n", sep = "")
  cat("  user: ", x$user_display_name, " (", x$username, ")\n", sep = "")
  cat("  logged in at: ", format(x$logged_in_at), "\n", sep = "")
  if (!is.null(x$site_auth)) cat("  site (web) password: set for ", x$site_auth$username, "\n", sep = "")
  invisible(x)
}

#' List breeding programs, folders and locations in Breedbase
#'
#' Use these to find the exact names to give to [save_breedbase()].
#'
#' @param con A connection from [breedbase_login()].
#' @param program For `breedbase_folders()`: the breeding program name.
#'
#' @return A tibble with names and database ids.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' breedbase_programs(con)
#' breedbase_folders(con, program = "IITA")
#' breedbase_locations(con)
#' }
#'
#' @name breedbase_lists
NULL

#' @rdname breedbase_lists
#' @export
breedbase_programs <- function(con) {
  data <- bb_get_all(con, "programs")
  tibble::tibble(
    program = vapply(data, function(x) x$programName %||% NA_character_, character(1)),
    program_db_id = vapply(data, function(x) as.character(x$programDbId %||% NA), character(1))
  )
}

#' @rdname breedbase_lists
#' @export
breedbase_folders <- function(con, program) {
  prog <- bb_find_program(con, program)
  data <- bb_get_all(con, "trials", list(programDbId = prog$program_db_id))
  tibble::tibble(
    folder = vapply(data, function(x) x$trialName %||% NA_character_, character(1)),
    folder_db_id = vapply(data, function(x) as.character(x$trialDbId %||% NA), character(1)),
    program = prog$program
  )
}

#' @rdname breedbase_lists
#' @export
breedbase_locations <- function(con) {
  data <- bb_get_all(con, "locations")
  tibble::tibble(
    location = vapply(data, function(x) x$locationName %||% NA_character_, character(1)),
    location_db_id = vapply(data, function(x) as.character(x$locationDbId %||% NA), character(1)),
    country = vapply(data, function(x) x$countryName %||% NA_character_, character(1))
  )
}

#' Check names in Breedbase before building a design
#'
#' Run this first, before designing the trial: it checks that the location(s)
#' and every entry and check are already in Breedbase (and, when given, the
#' breeding program, the folder and that the trial name is free). All
#' problems are reported together, with close matches as suggestions
#' (for example `"ibadan"` -> `"Ibadan"`), so they can be fixed in one go.
#' [save_breedbase()] repeats the same checks before writing.
#'
#' @param con A connection from [breedbase_login()].
#' @param accessions Entry and check names: a character vector, a design
#'   from [augmented_row_column_design()], [augmented_prep_design()],
#'   [resolvable_row_column_design()] or [as_breedbase_design()], or a list
#'   of such designs (e.g. one per location).
#' @param location Breedbase location name(s).
#' @param program,folder Optional breeding program and folder (the folder is
#'   optional for saving too: without it the trial goes directly in the
#'   program).
#' @param trial_name Optional trial name(s), e.g. one per location. A name
#'   already used by a trial (with plots) is reported as
#'   `"TAKEN"`: change it before saving.
#' @param plot_names Optional plot names to check: a character vector, a
#'   design, or a list of designs. By default the plot names stored in the
#'   design(s) given in `accessions` are checked (see [add_plot_names()]).
#'   A design without stored plot names gets the default names of
#'   [save_breedbase()] for the check when it can be matched to one trial
#'   name and location (one design and one `trial_name`, or a list of
#'   designs and one trial name and location per design). Plot names must be unique in the whole
#'   database: a name already used is reported as `"TAKEN"`, a name given to
#'   two plots as `"REPEATED"`. To save with other names, change the
#'   `plot_name` column of the design or use [add_plot_names()].
#'
#' @section Slow servers:
#' Accessions are searched 25 names at a time, because Breedbase looks up
#' the pedigree of each one. If the server times out (HTTP 504, 502, 503 or
#' 429), the batch is split in half and tried again, so a busy server is
#' slower but still works. A line is printed after each batch with the names
#' done and the time.
#' Accessions found are remembered until R is restarted, so
#' [save_breedbase()] after [breedbase_check()] does not search them again.
#'
#' When creating the folder, the trial or its plots times out, the save
#' does not stop: the server may still finish the job, so it checks for the
#' result every 20 seconds for 3 minutes, and tries once more if it is not
#' there (Breedbase refuses names that already exist, so this cannot make
#' duplicates). If a save is cut off anyway, run [save_breedbase()] again:
#' an existing folder is used, and a trial left without plots gets its plots.
#' Settings: `options(breedtools.write_wait = 30, breedtools.write_checks = 10)`. Settings: `options(breedtools.germplasm_batch =
#' 10)` for smaller batches, `options(breedtools.breedbase_timeout = 600)` to
#' wait longer (seconds) for each answer.
#'
#' @return Invisibly, a tibble with one row per name checked: `what`, `name`,
#'   `status` (`"OK"`, `"MISSING"`, `"NEW"`, `"TAKEN"`, `"EMPTY"` or `"REPEATED"`) and
#'   `suggestion`. Its attribute `"ok"` is `TRUE` when nothing needs fixing.
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' entries <- c(paste0("IITA-TMS-IBA", 1:108))
#' checks <- c("TMEB419", "IITA-TMS-IBA30572", "TMEB693", "IITA-TMS-IBA000070")
#' breedbase_check(con, accessions = c(entries, checks), location = "Ibadan",
#'                 program = "IITA", folder = "AYT 2026", trial_name = "AYT_2026_Ibadan")
#' }
#'
#' @export
breedbase_check <- function(con, accessions = NULL, location = NULL, program = NULL, folder = NULL,
                            trial_name = NULL, plot_names = NULL) {
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  rows <- list()
  add <- function(what, name, status, suggestion = NA_character_) {
    rows[[length(rows) + 1]] <<- tibble::tibble(what = what, name = name, status = status, suggestion = suggestion)
  }
  cat("Checking names in ", con$url, "\n", sep = "")
  if (!is.null(trial_name)) {
    trial_name <- bb_clean_name(trial_name, "trial name")
  }

  if (!is.null(program)) {
    progs <- breedbase_programs(con)
    if (program %in% progs$program) {
      add("program", program, "OK")
      if (!is.null(folder)) {
        folders <- breedbase_folders(con, program)$folder
        if (!folder %in% folders) folder <- bb_clean_name(folder, "folder name")
        if (folder %in% folders) add("folder", folder, "OK")
        else add("folder", folder, "NEW", bb_suggest(folder, folders, "it will be created if you save with create_folder = TRUE"))
      }
    } else {
      add("program", program, "MISSING", bb_suggest(program, progs$program))
    }
  }
  if (!is.null(location)) {
    loc <- bb_check_locations(con, unique(unname(location)))
    for (i in seq_len(nrow(loc))) add("location", loc$name[i], loc$status[i], loc$suggestion[i])
  }
  for (d in unique(trial_name[duplicated(trial_name)])) add("trial", d, "REPEATED", "each trial needs its own name")
  for (tn in unique(trial_name)) {
    trial_id <- bb_find_study(con, tn)
    status <- if (is.na(trial_id)) "OK" else if (bb_count_plots(con, trial_id) > 0) "TAKEN" else "EMPTY"
    add("trial", tn, status, switch(status,
      TAKEN = paste0("a trial with this name is already in Breedbase (", con$url, "/breeders/trial/", trial_id,
                     "): choose another trial_name"),
      EMPTY = "exists without plots (an earlier save was cut off): save_breedbase() will add the plots",
      NA_character_))
  }
  if (!is.null(accessions)) {
    names_in <- bb_accession_names(accessions)
    for (d in unique(names_in[duplicated(names_in)])) add("accession", d, "REPEATED", "each entry should appear once")
    acc <- bb_check_accessions(con, unique(names_in))
    rows[[length(rows) + 1]] <- tibble::tibble(what = "accession", name = acc$name, status = acc$status,
                                               suggestion = acc$suggestion)
  }
  # Plot names: given, or from the design(s) given as accessions (a plain
  # vector of accession names has no plot names).
  plot_source <- if (!is.null(plot_names)) plot_names else
    if (is.character(accessions) || is.factor(accessions)) NULL else accessions
  plots_in <- bb_check_plot_input(plot_source, trial_name, location)
  if (length(plots_in$names)) {
    pn <- plots_in$names
    for (d in unique(pn[duplicated(pn)])) add("plot", d, "REPEATED", "plot names must be unique: rename one of them")
    taken <- bb_existing_plots(con, unique(pn))
    rows[[length(rows) + 1]] <- tibble::tibble(
      what = "plot", name = unique(pn), status = ifelse(unique(pn) %in% taken, "TAKEN", "OK"),
      suggestion = ifelse(unique(pn) %in% taken,
                          "a plot with this name is already in Breedbase: change the plot_name column or use add_plot_names()",
                          NA_character_))
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    tibble::tibble(what = character(), name = character(), status = character(), suggestion = character())
  bb_print_check(out, plots_in$unnamed)
  attr(out, "ok") <- all(out$status %in% c("OK", "NEW", "EMPTY"))
  invisible(out)
}

# Plot names to check: a character vector, or the plot names stored in a
# design (or a list of designs). A design without plot names gets the
# default names when it matches one trial name; otherwise it is counted as
# unnamed. Filler plots are left out.
bb_check_plot_input <- function(x, trial_name = NULL, location = NULL) {
  none <- list(names = character(0), unnamed = 0L)
  if (is.null(x)) return(none)
  if (is.character(x) || is.factor(x)) {
    x <- trimws(as.character(x))
    return(list(names = bb_clean_name(x[!is.na(x) & nzchar(x)], "plot names"), unnamed = 0L))
  }
  designs <- if (bb_is_design(x)) list(x) else if (is.list(x) && length(x) && all(vapply(x, bb_is_design, logical(1)))) x else NULL
  if (is.null(designs)) stop("`plot_names` must be a character vector, a design or a list of designs.", call. = FALSE)
  out <- character(0)
  unnamed <- 0L
  for (i in seq_along(designs)) {
    dz <- designs[[i]]
    if (!"plot_name" %in% names(dz$design)) {
      loc <- if (inherits(dz, "augmented_prep_design")) NULL else
        if (length(location) == length(designs)) unname(location[i]) else
        if (length(location) == 1) unname(location) else NULL
      named <- if (length(trial_name) == length(designs)) {
        tryCatch(suppressMessages(add_plot_names(dz, trial_name[i], location = loc)), error = function(e) NULL)
      }
      if (!is.null(named)) {
        dz <- named
      } else {
        unnamed <- unnamed + 1L
        next
      }
    }
    d <- as.data.frame(dz$design, stringsAsFactors = FALSE)
    if ("is_filler" %in% names(d)) d <- d[d$is_filler != 1, , drop = FALSE]
    pn <- trimws(as.character(d$plot_name))
    if (anyNA(pn) || any(!nzchar(pn))) {
      stop("Some plots of the design have no plot name: fill the `plot_name` column or use add_plot_names().",
           call. = FALSE)
    }
    out <- c(out, bb_clean_name(pn, "plot names"))
  }
  list(names = out, unnamed = unnamed)
}

bb_is_design <- function(x) {
  inherits(x, c("augmented_row_column_design", "augmented_prep_design", "resolvable_row_column_design",
                "breedbase_design"))
}

# Accession names from a vector or a design (filler plots are left out).
bb_accession_names <- function(x) {
  if (is.character(x) || is.factor(x)) return(trimws(as.character(x)))
  if (is.list(x) && is.null(x$design) && length(x) && all(vapply(x, bb_is_design, logical(1)))) {
    return(unique(unlist(lapply(x, bb_accession_names))))
  }
  if (is.list(x) && !is.null(x$design)) {
    d <- as.data.frame(x$design, stringsAsFactors = FALSE)
    if ("is_filler" %in% names(d)) d <- d[d$is_filler != 1, , drop = FALSE]
    return(unique(as.character(d$all_entries)))
  }
  stop("`accessions` must be a character vector or a design from breedtools.", call. = FALSE)
}

# Close matches: same name ignoring case and spaces first, then small typos.
bb_suggest <- function(name, choices, otherwise = NA_character_) {
  choices <- unique(choices[!is.na(choices)])
  key <- function(z) gsub("[[:space:]_-]+", "", tolower(z))
  hit <- choices[key(choices) == key(name)]
  if (!length(hit)) hit <- choices[utils::adist(tolower(name), tolower(choices))[1, ] <= max(1, nchar(name) %/% 5)]
  if (length(hit)) return(paste0("did you mean '", paste(utils::head(hit, 3), collapse = "' or '"), "'?"))
  otherwise
}

bb_check_locations <- function(con, names, locs = breedbase_locations(con)) {
  found <- names %in% locs$location
  tibble::tibble(
    name = names,
    status = ifelse(found, "OK", "MISSING"),
    location_db_id = vapply(names, function(n) {
      hit <- locs$location_db_id[locs$location == n]
      if (length(hit) == 1) hit else NA_character_
    }, character(1), USE.NAMES = FALSE),
    suggestion = vapply(seq_along(names), function(i) {
      if (found[i]) NA_character_ else bb_suggest(names[i], locs$location, "see breedbase_locations(con)")
    }, character(1))
  )
}

# Exact search first; missing names are searched again in upper, lower and
# trimmed spellings to suggest the name stored in Breedbase.
bb_check_accessions <- function(con, names) {
  found <- bb_existing_germplasm(con, names)
  missing <- setdiff(names, found)
  suggestion <- stats::setNames(rep(NA_character_, length(names)), names)
  if (length(missing)) {
    variants <- unique(c(toupper(missing), tolower(missing), gsub("[[:space:]]+", "", missing)))
    variants <- setdiff(variants, names)
    stored <- if (length(variants)) bb_existing_germplasm(con, variants) else character(0)
    for (m in missing) {
      suggestion[[m]] <- if (length(stored)) bb_suggest(m, stored, "add it to Breedbase first") else
        "add it to Breedbase first"
    }
  }
  tibble::tibble(name = names, status = ifelse(names %in% found, "OK", "MISSING"),
                 suggestion = unname(suggestion))
}

bb_print_check <- function(out, unnamed = 0L) {
  for (w in unique(out$what)) {
    x <- out[out$what == w, , drop = FALSE]
    bad <- x[!x$status %in% c("OK", "NEW", "EMPTY"), , drop = FALSE]
    new <- x[x$status == "NEW", , drop = FALSE]
    if (w %in% c("accession", "plot") && nrow(x) > 3) {
      n_ok <- sum(x$status == "OK")
      n_names <- sum(x$status != "REPEATED")
      tag <- if (nrow(bad)) "[FIX]" else "[OK] "
      what_ok <- if (w == "accession") " accessions are in Breedbase" else
        paste0(" plot names are free, e.g. '", x$name[x$status == "OK"][1], "'")
      if (w == "plot" && n_ok == 0) what_ok <- " plot names are free"
      cat("  ", tag, "  ", n_ok, " of ", n_names, what_ok, "\n", sep = "")
    } else {
      ok_note <- switch(w, folder = " exists: the trial will be saved inside it", trial = " is free",
                        plot = " is free", "")
      for (i in seq_len(nrow(x))) if (x$status[i] == "OK") cat("  [OK]   ", w, " '", x$name[i], "'", ok_note, "\n", sep = "")
    }
    for (i in which(x$status == "EMPTY")) cat("  [NOTE] ", w, " '", x$name[i], "' ", x$suggestion[i], "\n", sep = "")
    for (i in seq_len(nrow(new))) {
      cat("  [NEW]  ", w, " '", new$name[i], "' does not exist yet",
          if (!is.na(new$suggestion[i])) paste0(": ", new$suggestion[i]), "\n", sep = "")
    }
    shown <- utils::head(seq_len(nrow(bad)), 20)
    for (i in shown) {
      status_text <- c(MISSING = "not found", TAKEN = "already exists", REPEATED = "is repeated")[[bad$status[i]]]
      cat("  [FIX]  ", w, " '", bad$name[i], "' ", status_text,
          if (!is.na(bad$suggestion[i])) paste0(": ", bad$suggestion[i]), "\n", sep = "")
    }
    if (nrow(bad) > length(shown)) cat("         ... and ", nrow(bad) - length(shown), " more\n", sep = "")
  }
  if (unnamed > 0) {
    cat("  [NOTE] ", unnamed, " design", if (unnamed > 1) "s have" else " has",
        " no plot names yet, so their plot names were not checked: give one design per trial ",
        "(with its trial_name and location), or add them with add_plot_names(design, trial_name, location)\n",
        sep = "")
  }
  # What was checked, e.g. "1 program, 3 locations, 112 accessions, 162 plot names".
  labels <- c(program = "program", folder = "folder", location = "location", accession = "accession",
              trial = "trial name", plot = "plot name")
  counts <- vapply(names(labels), function(w) length(unique(out$name[out$what == w])), integer(1))
  checked <- paste0(counts, " ", labels, ifelse(counts > 1, "s", ""))[counts > 0]
  if (any(out$what == "folder" & out$status == "NEW")) checked[names(labels)[counts > 0] == "folder"] <- "1 new folder"
  if (unnamed > 0) checked <- c(checked, "plot names NOT checked")
  checked <- paste(checked, collapse = ", ")
  bad <- out[!out$status %in% c("OK", "NEW", "EMPTY"), , drop = FALSE]
  cat("\n")
  if (nrow(bad) == 0) {
    cat("READY: everything checked can be saved to Breedbase (", checked, ").\n", sep = "")
    cat("Next: save_breedbase(..., dry_run = TRUE) for a last check, then without dry_run to save.\n")
  } else {
    per <- table(factor(bad$what, levels = names(labels)))
    per <- per[per > 0]
    what_bad <- paste0(per, " ", labels[names(per)], ifelse(per > 1, "s", ""), collapse = ", ")
    cat("NOT READY: fix ", nrow(bad), " problem", if (nrow(bad) > 1) "s" else "", " (", what_bad,
        ", marked [FIX] above) and run breedbase_check() again.\n", sep = "")
    cat("Checked: ", checked, ".\n", sep = "")
  }
}

#' Save a design to Breedbase
#'
#' Creates the trial and its plots in a Breedbase database through BrAPI, in
#' a breeding program and, optionally, a folder of your choice. Works for designs from
#' [augmented_row_column_design()], [augmented_prep_design()],
#' [resolvable_row_column_design()], and agricolae designs converted with
#' [as_breedbase_design()] or prepared with [prepare_breedbase_trial()].
#'
#' Before anything is written, every name is checked: the breeding program,
#' the folder if given (created when `create_folder = TRUE`), the location(s), that
#' every entry and check already exists as an accession, and that the trial
#' name is free. In an interactive session you are then asked to confirm.
#' After saving, the plots are read back from Breedbase to verify the upload.
#'
#' Each plot is saved with its name, plot number, row and column, rep, block
#' and whether it is a check, and its GPS coordinates when the design has
#' them (see [add_plot_coordinates()]; they are sent just after the plots
#' and read back to check). Plot names stored in the design (a `plot_name`
#' column, made by the design functions when given a `trial_name` or by
#' [add_plot_names()]) are used as they are, so they can be reviewed and
#' changed before saving. Otherwise they are made automatically from the
#' trial name, accession, block and rep (see `plot_name_pattern`), or you can
#' give your own with `plot_names`. Plot names must be unique in the whole
#' database, which is checked before saving. A multi-location p-rep design becomes one Breedbase
#' trial per location, named `"<trial_name>_<location>"`. Filler plots of
#' resolvable row-column designs are not saved.
#'
#' Run first with `dry_run = TRUE`: all the checks are made, but nothing is
#' written. Without a connection, `dry_run = TRUE` returns the data that would
#' be sent, without contacting any server.
#'
#' @param design A design object, including the result of [as_breedbase_design()]
#'   or [prepare_breedbase_trial()]. For a prepared trial, recorded trial
#'   name, year, location, program, description, folder and trial type are
#'   used when their arguments are omitted. Explicit arguments override them.
#' @param con A connection from [breedbase_login()].
#' @param program Breeding program name.
#' @param folder Optional folder inside the program. Without a folder (the
#'   default, `NULL`) the trial is saved directly in the breeding program. An
#'   existing folder is used as it is; a new one is created when
#'   `create_folder = TRUE`.
#' @param trial_name Name of the new trial (must not exist yet). In the
#'   names saved in Breedbase (trial, plot and new folder names, and the
#'   location part of p-rep trial names) spaces and special characters such
#'   as `% $ ( ) { } [ ] \\ / < > *` are replaced with `"_"` (only letters,
#'   digits, `_`, `-` and `.` are kept), e.g. `"Alcobaca (BA)"` gives
#'   `"Alcobaca_BA"`. Accession and location names are kept as they are, to
#'   match Breedbase.
#' @param year Year of the trial.
#' @param location Breedbase location name. For p-rep designs, a named vector
#'   mapping the design's locations to Breedbase locations, e.g.
#'   `c(North = "Ibadan", South = "Ubiaja")`; by default the design's
#'   location names are used.
#' @param description Description of the trial (required), e.g. its
#'   purpose, material and season: `"Advanced yield trial of 108 lines and 4
#'   checks, Ibadan 2026"`. It is shown on the trial page in Breedbase.
#' @param trial_type Breedbase trial type, e.g. `"Preliminary Yield Trial"`.
#'   Types the database does not know are saved as a miscellaneous trial.
#' @param design_type Breedbase design type. Defaults to `"Augmented"` for
#'   augmented row-column designs, `"p-rep"` for p-rep designs and `"Alpha"`
#'   (rows as incomplete blocks within replicates) for resolvable row-column
#'   designs. Converted agricolae designs use the type selected by
#'   [as_breedbase_design()].
#' @param plot_names Optional plot names given by you: either a character
#'   vector with one name per plot, in the row order of `design$design`, or
#'   the name of a column of `design$design` that holds them. Names must be
#'   unique.
#' @param plot_name_pattern Pattern for automatic plot names, used when
#'   `plot_names` is `NULL`. With `NULL` (the default) the `plot_name`
#'   column of the design is used when it has one, and otherwise the pattern
#'   `"{trial}_{location}_{accession}_{block}_{rep}"`, which gives names like
#'   `"AYT_2026_Ibadan_LINE_007_3_1"`. Available fields: `{trial}` (the
#'   Breedbase trial name), `{location}` (the Breedbase location, with spaces
#'   and special characters replaced by `"_"`), `{accession}`, `{block}`,
#'   `{rep}`, `{plot}` (plot number), `{row}` and `{col}`. If a pattern gives the same name to
#'   two plots (the same accession twice in one block and rep), the plot
#'   number is added to those names.
#' @param create_folder Logical. Create `folder` if it does not exist.
#' @param dry_run Logical. Make all checks but write nothing.
#' @inheritSection breedbase_check Slow servers
#' @param confirm Logical. Ask before writing (default in interactive
#'   sessions).
#'
#' @return Invisibly, a tibble with one row per Breedbase trial: location,
#'   trial name, trial id, number of plots sent and found after saving,
#'   number of plots with GPS coordinates sent and found, and the address of
#'   the trial page. With `dry_run = TRUE` and no connection,
#'   a list with the trial and plot data that would be sent.
#'
#' @seealso [breedbase_login()], [breedbase_programs()]
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' # What would be sent (no server needed):
#' payload <- save_breedbase(rc, trial_name = "AYT_2026_Ithaca", year = 2026,
#'                           location = "Ithaca",
#'                           description = "Advanced yield trial of 108 lines and 4 checks",
#'                           dry_run = TRUE)
#' str(payload$studies[[1]], max.level = 1)
#' payload$plots[[1]][[1]]
#'
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
#'                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
#'                description = "Advanced yield trial of 108 lines and 4 checks, Ibadan 2026",
#'                create_folder = TRUE, dry_run = TRUE)   # check first
#' save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
#'                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
#'                description = "Advanced yield trial of 108 lines and 4 checks, Ibadan 2026")
#' }
#'
#' @export
save_breedbase <- function(design, con = NULL, program = NULL, folder = NULL, trial_name, year,
                           location = NULL, description, trial_type = "phenotyping_trial",
                           design_type = NULL, plot_names = NULL,
                           plot_name_pattern = NULL,
                           create_folder = FALSE, dry_run = FALSE, confirm = interactive()) {
  if (inherits(design, "breedbase_trial")) {
    metadata <- design$trial
    if (missing(trial_name)) trial_name <- metadata$trial_name
    if (missing(year)) year <- metadata$year
    if (missing(location)) location <- metadata$location
    if (missing(program)) program <- metadata$program
    if (missing(description)) description <- metadata$description
    if (missing(folder)) folder <- metadata$folder
    if (missing(trial_type)) trial_type <- metadata$trial_type
  }
  if (missing(description) || is.null(description) || length(description) != 1 ||
      is.na(description) || !nzchar(trimws(description))) {
    stop("Give a `description` of the trial, e.g. description = \"Advanced yield trial of ",
         "108 lines and 4 checks, Ibadan 2026\". It is shown on the trial page in Breedbase.",
         call. = FALSE)
  }
  studies <- bb_design_to_studies(design, trial_name, year, location, description, trial_type, design_type,
                                  plot_names, plot_name_pattern)

  if (is.null(con)) {
    if (!isTRUE(dry_run)) stop("`con` is needed to save; use breedbase_login() first, or `dry_run = TRUE`.", call. = FALSE)
    message("Dry run without a connection: returning the data that would be sent.")
    return(list(
      studies = lapply(studies, `[[`, "study"),
      plots = lapply(studies, `[[`, "plots")
    ))
  }
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  if (is.null(program)) stop("Give the breeding `program`.", call. = FALSE)
  if (!is.null(folder) && (length(folder) != 1 || is.na(folder) || !nzchar(trimws(folder)))) {
    stop("`folder` must be one folder name, or NULL to save the trial directly in the program.", call. = FALSE)
  }

  cat("Checking before saving to ", con$url, "\n", sep = "")
  prog <- bb_find_program(con, program)
  cat("  [OK]   breeding program '", prog$program, "'\n", sep = "")

  if (is.null(folder)) {
    # No folder: the trial's parent is the breeding program itself.
    folder_id <- prog$program_db_id
    cat("  [OK]   no folder given: the trial will be saved directly in program '", prog$program, "'\n", sep = "")
  } else {
    folders <- breedbase_folders(con, prog$program)
    # An existing folder is used as it is, even with spaces in its name; a
    # new folder is created with spaces replaced by "_".
    if (!folder %in% folders$folder) folder <- bb_clean_name(folder, "folder name")
    folder_row <- folders[folders$folder == folder, , drop = FALSE]
    folder_id <- if (nrow(folder_row) == 1) folder_row$folder_db_id else NA_character_
  }
  if (is.na(folder_id)) {
    if (!isTRUE(create_folder)) {
      stop("Folder '", folder, "' was not found in program '", prog$program, "'. Existing folders: ",
           paste(utils::head(folders$folder, 20), collapse = ", "),
           ". Use create_folder = TRUE to create it.", call. = FALSE)
    }
    cat("  [NEW]  folder '", folder, "' will be created\n", sep = "")
  } else if (!is.null(folder)) {
    cat("  [OK]   folder '", folder, "' exists: the trial will be saved inside it (",
        con$url, "/folder/", folder_id, ")\n", sep = "")
  }

  # Locations and accessions are checked together, so every problem is
  # reported at once.
  loc_names <- unique(vapply(studies, `[[`, character(1), "breedbase_location"))
  loc_check <- bb_check_locations(con, loc_names)
  germplasm <- unique(unlist(lapply(studies, function(s) vapply(s$plots, `[[`, character(1), "germplasmName"))))
  acc_check <- bb_check_accessions(con, germplasm)
  problems <- character(0)
  for (i in which(loc_check$status != "OK")) {
    problems <- c(problems, paste0("Location '", loc_check$name[i], "' was not found in Breedbase",
                                   if (!is.na(loc_check$suggestion[i])) paste0(" (", loc_check$suggestion[i], ")"),
                                   ". Map design locations with `location = c(<design name> = \"<Breedbase name>\")`."))
  }
  bad_acc <- acc_check[acc_check$status != "OK", , drop = FALSE]
  if (nrow(bad_acc) > 0) {
    shown <- utils::head(seq_len(nrow(bad_acc)), 15)
    problems <- c(problems, paste0(
      nrow(bad_acc), " of ", length(germplasm), " entries/checks are not accessions in Breedbase: ",
      paste0(bad_acc$name[shown], ifelse(grepl("^did you mean", bad_acc$suggestion[shown]),
                                         paste0(" (", bad_acc$suggestion[shown], ")"), ""), collapse = ", "),
      if (nrow(bad_acc) > 15) ", ..." else "",
      ". Add them to Breedbase first (or fix the names). Tip: run breedbase_check() before designing."))
  }
  if (length(problems)) stop(paste(problems, collapse = "\n"), call. = FALSE)
  # Breedbase needs the database id of each accession to save the plots.
  germplasm_ids <- stats::setNames(bb_germplasm_ids(con, germplasm), germplasm)
  if (anyNA(germplasm_ids)) {
    stop("Breedbase did not return the id of accession(s) ",
         paste(utils::head(names(germplasm_ids)[is.na(germplasm_ids)], 10), collapse = ", "), ".", call. = FALSE)
  }
  for (i in seq_along(studies)) {
    studies[[i]]$study$locationDbId <- loc_check$location_db_id[loc_check$name == studies[[i]]$breedbase_location]
  }
  cat("  [OK]   location", if (length(loc_names) > 1) "s" else "", " ", paste0("'", loc_names, "'", collapse = ", "),
      "\n", sep = "")
  cat("  [OK]   all ", length(germplasm), " entries and checks exist as accessions\n", sep = "")

  # A trial that exists without plots was left by an earlier save that was
  # cut off (e.g. a server timeout): the plots are added to it.
  for (i in seq_along(studies)) {
    name <- studies[[i]]$study$studyName
    existing_id <- bb_find_study(con, name)
    if (!is.na(existing_id)) {
      if (bb_count_plots(con, existing_id) > 0) {
        stop("A trial named '", name, "' already exists in Breedbase. Choose another `trial_name`.", call. = FALSE)
      }
      studies[[i]]$existing_id <- existing_id
      cat("  [OK]   trial '", name, "' exists without plots (an earlier save was cut off): ",
          "its plots will be added\n", sep = "")
    }
  }
  new_names <- vapply(Filter(function(s) is.null(s$existing_id), studies), function(s) s$study$studyName, character(1))
  if (length(new_names)) {
    cat("  [OK]   trial name", if (length(new_names) > 1) "s are" else " is", " free: ",
        paste(new_names, collapse = ", "), "\n", sep = "")
  }

  plot_names_all <- unlist(lapply(studies, function(s) vapply(s$plots, `[[`, character(1), "observationUnitName")))
  taken <- bb_existing_plots(con, plot_names_all)
  if (length(taken) > 0) {
    stop(length(taken), " plot names already exist in Breedbase (plot names must be unique in the database), e.g. ",
         paste(utils::head(taken, 5), collapse = ", "), ". Use another `trial_name` or `plot_name_pattern`, or change the plot names with ",
         "add_plot_names() or in the design's `plot_name` column.",
         call. = FALSE)
  }
  cat("  [OK]   all ", length(plot_names_all), " plot names are new, e.g. ", plot_names_all[1], "\n", sep = "")
  n_geo <- sum(vapply(studies, function(s) length(bb_study_geo(s)), integer(1)))
  if (n_geo) {
    cat("  [OK]   GPS coordinates of ", n_geo, " plots will be saved after the plots\n", sep = "")
  }

  n_plots <- sum(vapply(studies, function(s) length(s$plots), integer(1)))
  where <- if (is.null(folder)) sprintf("directly in program '%s'", prog$program) else
    sprintf("in folder '%s' of program '%s'", folder, prog$program)
  plan <- sprintf("%d trial%s with %d plots %s at %s",
                  length(studies), if (length(studies) > 1) "s" else "", n_plots, where, con$url)
  if (isTRUE(dry_run)) {
    cat("Dry run: everything is ready to save ", plan, ". Nothing was written.\n", sep = "")
    return(invisible(tibble::tibble(
      location = vapply(studies, `[[`, character(1), "breedbase_location"),
      trial_name = vapply(studies, function(s) s$study$studyName, character(1)),
      trial_db_id = NA_character_,
      plots_sent = vapply(studies, function(s) length(s$plots), integer(1)),
      plots_found = NA_integer_,
      plots_with_gps = vapply(studies, function(s) length(bb_study_geo(s)), integer(1)),
      gps_found = NA_integer_,
      trial_page = NA_character_
    )))
  }
  if (isTRUE(confirm)) {
    answer <- readline(sprintf("Save %s? Type 'yes' to continue: ", plan))
    if (!identical(tolower(trimws(answer)), "yes")) {
      message("Nothing was saved.")
      return(invisible(NULL))
    }
  }

  if (is.na(folder_id)) folder_id <- bb_create_folder(con, prog, folder, saved = "the trial")

  out <- vector("list", length(studies))
  for (i in seq_along(studies)) {
    s <- studies[[i]]
    s$study$trialDbId <- folder_id
    study_id <- s$existing_id %||% bb_write_checked(
      what = paste0("create trial '", s$study$studyName, "'"),
      write = function() {
        res <- bb_post(con, "studies", list(s$study), what = paste0("create trial '", s$study$studyName, "'"))
        as.character(res$result$data[[1]]$studyDbId %||% NA)
      },
      find = function() bb_find_study(con, s$study$studyName)
    )
    plots <- lapply(s$plots, function(p) {
      p$studyDbId <- study_id
      p$germplasmDbId <- germplasm_ids[[p$germplasmName]]
      p
    })
    # Breedbase saves all the plots of a trial together, or none.
    bb_write_checked(
      what = paste0("save the ", length(plots), " plots of '", s$study$studyName, "'"),
      write = function() {
        bb_post(con, "observationunits", plots, what = paste0("save the plots of '", s$study$studyName, "'"))
        study_id
      },
      find = function() if (bb_count_plots(con, study_id) >= length(plots)) study_id else NA_character_
    )
    n_found <- bb_count_plots(con, study_id)
    status <- if (n_found == length(plots)) "[OK]  " else "[FAIL]"
    cat(sprintf("  %s trial '%s' saved with %d of %d plots\n", status, s$study$studyName, n_found, length(plots)))
    geo <- bb_study_geo(s)
    gps_found <- if (length(geo)) bb_save_plot_geo(con, study_id, geo, s$study$studyName) else NA_integer_
    out[[i]] <- tibble::tibble(
      location = s$breedbase_location,
      trial_name = s$study$studyName,
      trial_db_id = study_id,
      plots_sent = length(plots),
      plots_found = n_found,
      plots_with_gps = length(geo),
      gps_found = gps_found,
      trial_page = paste0(con$url, "/breeders/trial/", study_id)
    )
  }
  out <- dplyr::bind_rows(out)
  cat("Open the trial in Breedbase: ", paste(out$trial_page, collapse = ", "), "\n", sep = "")
  invisible(out)
}

# ---- Internals --------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

bb_require_httr2 <- function() {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Connecting to Breedbase needs the 'httr2' package: install.packages(\"httr2\").", call. = FALSE)
  }
}

bb_status_message <- function(body) {
  st <- body$metadata$status
  if (is.null(st) || length(st) == 0) return("no message from the server")
  msgs <- vapply(st, function(x) as.character(x$message %||% ""), character(1))
  msgs <- msgs[nzchar(msgs) & !grepl("^BrAPI base call found|^Loading CXGN|^Search", msgs)]
  if (length(msgs) == 0) "no message from the server" else paste(unique(msgs), collapse = "; ")
}

# Perform a request and return the parsed JSON body; turn HTTP and BrAPI
# errors into readable messages.
bb_perform <- function(req, what) {
  req <- httr2::req_error(req, is_error = function(resp) FALSE)
  req <- httr2::req_user_agent(req, "breedtools (R package)")
  req <- httr2::req_timeout(req, getOption("breedtools.breedbase_timeout", 300))
  resp <- tryCatch(httr2::req_perform(req), error = function(e) {
    bb_stop(paste0("Could not ", what, ": ", conditionMessage(e)), code = NA_integer_,
            busy = grepl("timeout|timed out", conditionMessage(e), ignore.case = TRUE))
  })
  body <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
  code <- httr2::resp_status(resp)
  if (code >= 400) {
    detail <- if (is.null(body)) httr2::resp_status_desc(resp) else bb_status_message(body)
    hint <- if (code %in% c(401, 403)) " (check your username, password and permissions; you may need to log in again)" else ""
    basic <- httr2::resp_header(resp, "WWW-Authenticate") %||% ""
    if (code == 401 && grepl("^Basic", basic)) {
      realm <- sub('^.*realm="?([^"]*)"?.*$', "\\1", basic)
      detail <- paste0("the website asks for its own web password (\"", realm, "\")")
      hint <- paste0(". This site has an extra password in front of the website (BrAPI does not need it). ",
                     "Log in again with breedbase_login(url, username, site_username = \"...\") ",
                     "and the web user name and password of the site")
    }
    bb_stop(paste0("Could not ", what, " (HTTP ", code, "): ", detail, hint), code = code,
            busy = code %in% c(429L, 502L, 503L, 504L))
  }
  if (is.null(body)) stop("Could not ", what, ": the server did not return JSON. Is the address right?", call. = FALSE)
  body
}

# Errors carry the HTTP code and whether the server was just busy or slow
# (timeouts, 502/503/504), so long searches can be retried in smaller parts.
bb_stop <- function(message, code, busy = FALSE) {
  stop(structure(class = c("breedbase_error", "error", "condition"),
                 list(message = message, call = NULL, code = code, busy = busy)))
}

bb_request <- function(con, path) {
  req <- httr2::request(paste0(con$url, "/brapi/v2/", path))
  httr2::req_headers(req, Authorization = paste("Bearer", con$token))
}

# GET all pages of a BrAPI list call.
bb_get_all <- function(con, path, query = list(), page_size = 1000) {
  bb_require_httr2()
  out <- list()
  page <- 0
  repeat {
    req <- do.call(httr2::req_url_query, c(list(bb_request(con, path)), query, list(page = page, pageSize = page_size)))
    body <- bb_perform(req, what = paste("read", path))
    out <- c(out, body$result$data %||% list())
    total_pages <- body$metadata$pagination$totalPages %||% 1
    page <- page + 1
    if (page >= total_pages) break
  }
  out
}

bb_post <- function(con, path, payload, what) {
  bb_require_httr2()
  req <- httr2::req_body_json(bb_request(con, path), payload, auto_unbox = TRUE, null = "null")
  bb_perform(httr2::req_method(req, "POST"), what = what)
}

bb_put <- function(con, path, payload, what) {
  bb_require_httr2()
  req <- httr2::req_body_json(bb_request(con, path), payload, auto_unbox = TRUE, null = "null")
  bb_perform(httr2::req_method(req, "PUT"), what = what)
}

bb_find_program <- function(con, program) {
  progs <- breedbase_programs(con)
  hit <- progs[progs$program == program, , drop = FALSE]
  if (nrow(hit) != 1) {
    stop("Breeding program '", program, "' was not found. Programs available: ",
         paste(progs$program, collapse = ", "), call. = FALSE)
  }
  hit
}

# Accessions found in Breedbase during this R session, by server, so the
# same names are not searched again (e.g. breedbase_check() then
# save_breedbase()). Only names that were found are remembered.
# The cache holds germplasmDbId values named by accession name: Breedbase
# needs the id (not only the name) to save plots.
bb_found_cache <- new.env(parent = emptyenv())

bb_cached_ids <- function(con) bb_found_cache[[con$url]] %||% stats::setNames(character(0), character(0))

bb_existing_germplasm <- function(con, names) {
  names <- unique(names)
  known <- intersect(names, names(bb_cached_ids(con)))
  to_search <- setdiff(names, known)
  if (length(known) && length(to_search)) {
    cat("  ...    ", length(known), " of ", length(names),
        " accessions were already found in this session; searching the other ", length(to_search), "\n", sep = "")
  } else if (length(known)) {
    cat("  ...    all ", length(known), " accessions were already found in this session (no new search)\n", sep = "")
  }
  # Breedbase builds the pedigree of every accession it finds, so germplasm
  # searches are slow: send small batches.
  found <- bb_search_names(con, "germplasm", "germplasmNames", "germplasmName", to_search,
                           batch = getOption("breedtools.germplasm_batch", 25), label = "accessions",
                           id_field = "germplasmDbId")
  found <- found[!is.na(found) & nzchar(found)]
  cache <- bb_cached_ids(con)
  cache[names(found)] <- found
  assign(con$url, cache, envir = bb_found_cache)
  intersect(names, c(known, names(found)))
}

# germplasmDbId of each accession name (NA when not found).
bb_germplasm_ids <- function(con, names) {
  ids <- bb_cached_ids(con)
  missing <- setdiff(names, names(ids))
  if (length(missing)) bb_existing_germplasm(con, missing)
  unname(bb_cached_ids(con)[names])
}

bb_existing_plots <- function(con, names) {
  bb_search_names(con, "observationunits", "observationUnitNames", "observationUnitName", names,
                  batch = getOption("breedtools.plot_batch", 200), label = "plot names")
}

# Progress of long searches: one line per finished batch, so it shows the
# same way in the console, RStudio, source(), Rscript and logs:
#   ...    batch 2 of 5 done:  50/112 accessions (45%)  0:42
bb_progress <- function(total, label) {
  start <- Sys.time()
  elapsed <- function() {
    secs <- as.integer(round(as.numeric(difftime(Sys.time(), start, units = "secs"))))
    sprintf("%d:%02d", secs %/% 60, secs %% 60)
  }
  say <- function(text) {
    cat(text, "\n", sep = "")
    utils::flush.console()
  }
  list(
    start = function(n_batches, batch) {
      say(sprintf("  ...    searching %d %s in Breedbase, %d batch%s of up to %d", total, label, n_batches,
                  if (n_batches > 1) "es" else "", batch))
    },
    update = function(batch_no, n_batches, done) {
      say(sprintf("  ...    batch %d of %d done: %*d/%d %s (%d%%)  %s", batch_no, n_batches, nchar(total), done,
                  total, label, floor(100 * done / total), elapsed()))
    },
    note = function(text) say(paste0("  ...    ", text, "  ", elapsed())),
    finish = function(batches, retries) {
      extra <- if (retries > 0) sprintf(", %d slow batch%s split and retried", retries, if (retries > 1) "es" else "") else ""
      say(sprintf("  [OK]   searched %d %s in %d batch%s (%s%s)", total, label, batches,
                  if (batches > 1) "es" else "", elapsed(), extra))
    }
  )
}

# Search names in batches. A batch that times out (or gets 502/503/504/429)
# is split in half and tried again, down to one name, with a short pause.
# With `id_field`, the result is the ids of the names found, named by name.
bb_search_names <- function(con, entity, field, name_field, names, batch, label, id_field = NULL) {
  names <- unique(names)
  if (!length(names)) return(if (is.null(id_field)) character(0) else stats::setNames(character(0), character(0)))
  queue <- split(names, ceiling(seq_along(names) / batch))
  found <- character(0)
  found_ids <- character(0)
  pause <- getOption("breedtools.retry_pause", 1)
  show <- length(names) > batch
  bar <- if (show) bb_progress(length(names), label) else NULL
  done <- 0L
  batches <- 0L
  retries <- 0L
  if (show) bar$start(length(queue), batch)
  while (length(queue)) {
    chunk <- queue[[1]]
    queue <- queue[-1]
    data <- tryCatch(
      bb_search(con, entity, stats::setNames(list(as.list(chunk), length(chunk) * 2), c(field, "pageSize"))),
      breedbase_error = function(e) {
        if (!isTRUE(e$busy)) stop(e)
        retries <<- retries + 1L
        if (length(chunk) == 1) {
          if (show) bar$note(sprintf("server slow on '%s': waiting and trying once more", chunk))
          Sys.sleep(5 * pause)
          return(tryCatch(bb_search(con, entity, stats::setNames(list(list(chunk), 2), c(field, "pageSize"))),
                          breedbase_error = function(e2) {
                            stop(conditionMessage(e2), ". The server is too slow or busy at the moment: ",
                                 "try again later, or ask its administrators.", call. = FALSE)
                          }))
        }
        if (show) bar$note(sprintf("server slow: splitting a batch of %d into two of %d and %d",
                                   length(chunk), ceiling(length(chunk) / 2), length(chunk) %/% 2))
        Sys.sleep(pause)
        half <- ceiling(length(chunk) / 2)
        queue <<- c(list(chunk[seq_len(half)], chunk[-seq_len(half)]), queue)
        NULL
      }
    )
    if (is.null(data)) next
    found <- c(found, vapply(data, function(x) as.character(x[[name_field]] %||% NA_character_), character(1)))
    if (!is.null(id_field)) {
      found_ids <- c(found_ids, vapply(data, function(x) as.character(x[[id_field]] %||% NA_character_), character(1)))
    }
    done <- done + length(chunk)
    batches <- batches + 1L
    if (show) bar$update(batches, batches + length(queue), done)
  }
  if (show) bar$finish(batches, retries)
  if (!is.null(id_field)) {
    keep <- found %in% names & !duplicated(found)
    return(stats::setNames(found_ids[keep], found[keep]))
  }
  intersect(names, found)
}

# Run a write. If the server times out (or is busy), the write may still have
# been done, so look for the result for a while before trying once more.
# Breedbase refuses folder, trial and plot names that already exist, so a
# second try cannot make duplicates.
bb_write_checked <- function(what, write, find, attempts = 2, exists_ok = FALSE, exists_hint = NULL) {
  wait <- getOption("breedtools.write_wait", 20)
  checks <- getOption("breedtools.write_checks", 9)
  for (attempt in seq_len(attempts)) {
    res <- tryCatch(write(), breedbase_error = function(e) e)
    if (!inherits(res, "breedbase_error")) {
      if (is.na(res)) stop("Breedbase did not return an id when trying to ", what, ".", call. = FALSE)
      return(res)
    }
    if (isTRUE(exists_ok) && identical(as.integer(res$code), 409L)) {
      # Already there (e.g. made by an earlier save that timed out): use it.
      id <- find()
      if (is.na(id)) stop("Could not ", what, ": ", exists_hint, call. = FALSE)
      return(structure(id, existed = TRUE))
    }
    if (!isTRUE(res$busy)) stop(res)
    why <- if (is.na(res$code)) "no answer in time" else paste("HTTP", res$code)
    cat("  ...    the server was slow to ", what, " (", why, "); checking whether it was done\n", sep = "")
    for (k in seq_len(checks)) {
      Sys.sleep(wait)
      id <- tryCatch(find(), breedbase_error = function(e) NA_character_)
      if (!is.na(id)) {
        cat("  ...    confirmed: done\n")
        return(id)
      }
      cat(sprintf("  ...    not there yet (check %d of %d)\n", k, checks))
    }
    if (attempt < attempts) cat("  ...    not done: trying again\n")
  }
  stop("Could not ", what, ": the server kept timing out. Nothing is lost: run save_breedbase() ",
       "again later; what was already saved is found and the save continues.", call. = FALSE)
}

# A trial is found by its name as given or, for names with spaces, by the
# name breedtools saves (spaces replaced with "_").
bb_find_study <- function(con, name) {
  for (nm in unique(c(name, bb_clean_name(name, quiet = TRUE)))) {
    data <- bb_search(con, "studies", list(studyNames = list(nm), pageSize = 10))
    hit <- Filter(function(x) identical(x$studyName, nm), data)
    if (length(hit)) return(as.character(hit[[1]]$studyDbId %||% NA))
  }
  NA_character_
}

bb_count_plots <- function(con, study_id) {
  length(bb_get_all(con, "observationunits", list(studyDbId = study_id, observationUnitLevelName = "plot")))
}

bb_study_exists <- function(con, name) {
  data <- bb_search(con, "studies", list(studyNames = list(name), pageSize = 10))
  any(vapply(data, function(x) identical(x$studyName, name), logical(1)))
}

# BrAPI search: POST /search/<entity>; the results come back directly or via
# a searchResultsDbId to GET.
bb_search <- function(con, entity, body) {
  res <- bb_post(con, paste0("search/", entity), body, what = paste("search", entity))
  if (!is.null(res$result$searchResultsDbId)) {
    return(bb_get_all(con, paste0("search/", entity, "/", res$result$searchResultsDbId)))
  }
  res$result$data %||% list()
}

# Create a folder (a BrAPI trial) in the program and return its id. A folder
# that already exists (e.g. made by an earlier save that timed out) is used.
bb_create_folder <- function(con, prog, folder, saved = "the trial") {
  folder_id <- bb_write_checked(
    what = paste0("create folder '", folder, "'"),
    write = function() {
      res <- bb_post(con, "trials", list(list(trialName = folder, programDbId = prog$program_db_id,
                                              trialDescription = paste("Created by breedtools"))),
                     what = paste0("create folder '", folder, "'"))
      as.character(res$result$data[[1]]$trialDbId %||% NA)
    },
    find = function() bb_folder_id(con, prog$program, folder),
    exists_ok = TRUE,
    exists_hint = paste0("a trial or a folder of another program already uses the name '", folder,
                         "'. Choose another `folder`.")
  )
  if (isTRUE(attr(folder_id, "existed"))) {
    cat("  [OK]   folder '", folder, "' already exists: ", saved, " is saved inside it\n", sep = "")
  } else {
    cat("  [OK]   folder '", folder, "' created\n", sep = "")
  }
  as.character(folder_id)
}

# Id of a folder of the program, or NA.
bb_folder_id <- function(con, program, folder) {
  f <- breedbase_folders(con, program)
  id <- f$folder_db_id[f$folder == folder]
  if (length(id) == 1) id else NA_character_
}

# Names saved in Breedbase (trials, plots, folders, analyses) contain only
# letters (accented letters too), digits, "_", "-" and ".": spaces and other
# characters (% $ ( ) { } [ ] \\ / < > * ...) become "_", repeated "_" are
# merged and "_" at the start or end are dropped, e.g. "Alcobaca (BA)" ->
# "Alcobaca_BA". Accession and location names are not changed, since they
# must match the names already in Breedbase.
bb_clean_name <- function(x, what = "name", quiet = FALSE) {
  out <- gsub("[^\\p{L}\\p{N}_.-]+", "_", trimws(as.character(x)), perl = TRUE)
  out <- gsub("_{2,}", "_", out)
  out <- gsub("^_+|_+$", "", out)
  if (any(!is.na(out) & !nzchar(out))) {
    stop("The ", what, " '", x[!is.na(out) & !nzchar(out)][1], "' has no letters or digits.", call. = FALSE)
  }
  changed <- !is.na(x) & out != x
  if (!quiet && any(changed)) {
    message("Spaces and special characters replaced with '_' in the ", what,
            if (sum(changed) > 1) paste0(" (", sum(changed), ")"),
            ": '", x[changed][1], "' -> '", out[changed][1], "'")
  }
  out
}

bb_default_plot_pattern <- "{trial}_{location}_{accession}_{block}_{rep}"

#' Add Breedbase plot names to a design
#'
#' Adds a `plot_name` column to `design$design` with the names the plots
#' will have in Breedbase, so they can be seen, checked with
#' [breedbase_check()] and changed before [save_breedbase()]. The design
#' functions call this when given a `trial_name`. Run it again to make the
#' names for another trial name or pattern, or edit the column by hand: any
#' unique names are accepted. [save_breedbase()], [phenotype_template()] and
#' [first_stage_values()] use the stored names.
#'
#' @param design A design from [augmented_row_column_design()],
#'   [augmented_prep_design()], [resolvable_row_column_design()],
#'   [as_breedbase_design()] or [prepare_breedbase_trial()].
#' @param trial_name The Breedbase trial name. For a multi-location p-rep
#'   design, the location is added as in [save_breedbase()] (one trial per
#'   location). Can be omitted for a prepared Breedbase trial.
#' @param location The Breedbase location of the trial, used by the
#'   `{location}` field of the pattern. Not needed for p-rep designs (their
#'   own locations are used, or give a named vector mapping them to
#'   Breedbase locations, as in [save_breedbase()]) or for a prepared
#'   Breedbase trial (its recorded location is used).
#' @param plot_name_pattern Pattern for the names, as in [save_breedbase()]:
#'   fields `{trial}`, `{location}`, `{accession}`, `{block}`, `{rep}`,
#'   `{plot}`, `{row}` and `{col}`.
#' @param plot_names Optional names given by you instead of a pattern: one
#'   per row of `design$design`, or the name of a column holding them.
#'
#' @return The design with a `plot_name` column (`NA` for filler plots,
#'   which are not saved) and the trial name in `plot_name_trial`.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' rc <- add_plot_names(rc, "AYT_2026", location = "Ibadan")
#' head(rc$design$plot_name)
#' # Other names, e.g. with the field position:
#' rc <- add_plot_names(rc, "AYT_2026_Ibadan", plot_name_pattern = "{trial}_{plot}_r{row}c{col}")
#' head(rc$design$plot_name)
#'
#' @export
add_plot_names <- function(design, trial_name = NULL, location = NULL,
                           plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}", plot_names = NULL) {
  if (is.null(trial_name) && inherits(design, "breedbase_trial")) trial_name <- design$trial$trial_name
  if (is.null(trial_name) || length(trial_name) != 1 || is.na(trial_name) || !nzchar(trimws(trial_name))) {
    stop("Give one `trial_name` for the plot names.", call. = FALSE)
  }
  if (!bb_is_design(design)) {
    stop("`design` must come from augmented_row_column_design(), augmented_prep_design(), ",
         "resolvable_row_column_design() or as_breedbase_design().", call. = FALSE)
  }
  if (is.null(plot_names) && is.null(plot_name_pattern)) plot_name_pattern <- bb_default_plot_pattern
  trial_name <- bb_clean_name(trial_name, "trial name")
  tab <- bb_design_plot_table(design, trial_name, plot_name_pattern, plot_names, location = location)
  d <- design$design
  is_prep <- inherits(design, "augmented_prep_design")
  key_d <- paste(if (is_prep) as.character(d$location) else "", d$plots, sep = "\r")
  key_t <- paste(if (is_prep) tab$.design_location else "", tab$.design_plots, sep = "\r")
  d$plot_name <- tab$plot_name[match(key_d, key_t)]
  design$design <- d
  design$plot_name_trial <- trial_name
  design
}

# Convert a design into BrAPI study and observation unit records (without
# database ids, which are filled in when saving).
# Plot names from a pattern such as "{trial}_{location}_{accession}_{block}_{rep}".
# Names that come out the same get the plot number added.
bb_plot_names <- function(pattern, trial, accession, block, rep, plot, row, col, location = NA) {
  fields <- list(trial = trial, location = location, accession = accession, block = block, rep = rep,
                 plot = plot, row = row, col = col)
  used <- regmatches(pattern, gregexpr("\\{[a-z]+\\}", pattern))[[1]]
  unknown <- setdiff(gsub("[{}]", "", used), names(fields))
  if (length(unknown) > 0) {
    stop("Unknown field(s) in `plot_name_pattern`: ", paste0("{", unknown, "}", collapse = ", "),
         ". Use {trial}, {location}, {accession}, {block}, {rep}, {plot}, {row} or {col}.", call. = FALSE)
  }
  if ("{location}" %in% used) {
    if (anyNA(location) || identical(location, bb_no_location)) {
      stop("`plot_name_pattern` uses {location}: give the Breedbase `location` of the trial, ",
           "or use a pattern without {location}.", call. = FALSE)
    }
    fields$location <- bb_clean_name(location, "location in plot names", quiet = TRUE)
  }
  n <- length(accession)
  out <- rep(pattern, n)
  for (f in names(fields)) {
    out <- vapply(seq_len(n), function(i) gsub(paste0("{", f, "}"), as.character(rep_len(fields[[f]], n)[i]),
                                               out[i], fixed = TRUE), character(1))
  }
  dup <- out %in% out[duplicated(out)]
  if (any(dup)) {
    message(sum(dup), " plots would have the same name with this pattern ",
            "(same accession in the same block and rep); the plot number was added to their names.")
    out[dup] <- paste0(out[dup], "_", plot[dup])
  }
  out
}

bb_design_to_studies <- function(design, trial_name, year, location, description, trial_type, design_type,
                                 plot_names = NULL, plot_name_pattern = NULL) {
  if (missing(trial_name) || is.null(trial_name) || !nzchar(trial_name)) stop("Give a `trial_name`.", call. = FALSE)
  if (missing(year)) stop("Give the `year` of the trial.", call. = FALSE)
  trial_name <- bb_clean_name(trial_name, "trial name")
  kind <- if (inherits(design, "augmented_row_column_design")) {
    "augmented"
  } else if (inherits(design, "augmented_prep_design")) {
    "prep"
  } else if (inherits(design, "resolvable_row_column_design")) {
    "rrc"
  } else if (inherits(design, "breedbase_design")) {
    "imported"
  } else {
    stop("`design` must come from augmented_row_column_design(), augmented_prep_design() or ",
         "resolvable_row_column_design(), or be converted with as_breedbase_design().", call. = FALSE)
  }
  if (is.null(design_type)) {
    design_type <- switch(kind, augmented = "Augmented", prep = "p-rep", rrc = "Alpha", imported = design$design_type)
  }
  supported <- bb_design_types()
  if (length(design_type) != 1L || is.na(design_type) || !design_type %in% supported) {
    stop("`design_type` must be one of the types Breedbase accepts: ", paste(supported, collapse = ", "), ".",
         call. = FALSE)
  }
  d <- as.data.frame(design$design, stringsAsFactors = FALSE)
  # Plot names stored in the design are used unless other names or a
  # pattern are given.
  if (is.null(plot_names) && is.null(plot_name_pattern)) {
    if ("plot_name" %in% names(d)) {
      plot_names <- "plot_name"
      named_for <- design$plot_name_trial
      if (!is.null(named_for) && !identical(named_for, trial_name)) {
        message("The plot names in the design were made for trial '", named_for, "', not '", trial_name,
                "'; they are used as they are. Use add_plot_names() to make them again.")
      }
    } else {
      plot_name_pattern <- bb_default_plot_pattern
    }
  }
  # User-given plot names: a vector in row order, or a column of the design.
  if (!is.null(plot_names)) {
    if (length(plot_names) == 1 && plot_names %in% names(d)) {
      d$.plot_name <- as.character(d[[plot_names]])
    } else if (length(plot_names) == nrow(d)) {
      d$.plot_name <- as.character(plot_names)
    } else {
      stop("`plot_names` must have one name per plot (", nrow(d), "), or be the name of a column of the design.",
           call. = FALSE)
    }
    keep <- if ("is_filler" %in% names(d)) d$is_filler != 1 else rep(TRUE, nrow(d))
    if (anyNA(d$.plot_name[keep]) || any(!nzchar(d$.plot_name[keep]))) stop("`plot_names` has empty names.", call. = FALSE)
    if (anyDuplicated(d$.plot_name[keep])) {
      stop("`plot_names` must be unique; repeated: ",
           paste(utils::head(unique(d$.plot_name[keep][duplicated(d$.plot_name[keep])]), 5), collapse = ", "),
           call. = FALSE)
    }
  }

  if (kind == "prep") {
    design_locations <- design$input$locations
    if (is.null(location)) {
      location <- stats::setNames(design_locations, design_locations)
    } else if (is.null(names(location))) {
      if (length(location) != length(design_locations)) {
        stop("Give one Breedbase location per design location (", paste(design_locations, collapse = ", "),
             "), as a named vector.", call. = FALSE)
      }
      names(location) <- design_locations
    }
    groups <- lapply(design_locations, function(l) list(
      design_location = l,
      breedbase_location = unname(location[[l]]),
      data = d[d$location == l, , drop = FALSE],
      name = if (length(design_locations) > 1) paste(trial_name, bb_clean_name(l), sep = "_") else trial_name
    ))
  } else {
    if (is.null(location) || length(location) != 1) stop("Give the Breedbase `location` of the trial.", call. = FALSE)
    groups <- list(list(design_location = NA, breedbase_location = unname(location), data = d, name = trial_name))
  }

  lapply(groups, function(g) {
    x <- g$data
    if ("is_filler" %in% names(x)) x <- x[x$is_filler != 1, , drop = FALSE]
    x <- x[order(x$plots), , drop = FALSE]
    plot_number <- if (kind == "prep") seq_len(nrow(x)) else x$plots
    rep_no <- switch(kind,
      augmented = x$rep,
      prep = if ("replicate" %in% names(x)) x$replicate else rep(1L, nrow(x)),
      rrc = x$replicate,
      imported = x$rep
    )
    block_no <- switch(kind,
      augmented = x$block,
      prep = x$block,
      rrc = (x$replicate - 1L) * max(x$row_in_rep) + x$row_in_rep,
      imported = x$block
    )
    is_check <- if ("is_control" %in% names(x)) x$is_control == 1 else rep(FALSE, nrow(x))
    has_coords <- all(c("row", "col") %in% names(x)) && !anyNA(x$row) && !anyNA(x$col)
    unit_names <- if (!is.null(x$.plot_name)) {
      x$.plot_name
    } else {
      bb_plot_names(plot_name_pattern, trial = g$name, location = g$breedbase_location %||% NA,
                    accession = as.character(x$all_entries),
                    block = block_no, rep = rep_no, plot = plot_number,
                    row = if (has_coords) x$row else NA, col = if (has_coords) x$col else NA)
    }
    unit_names <- bb_clean_name(unit_names, "plot names")
    if (anyDuplicated(unit_names)) {
      stop("Plot names are repeated after replacing spaces with '_': ",
           paste(utils::head(unique(unit_names[duplicated(unit_names)]), 5), collapse = ", "), ".", call. = FALSE)
    }
    plots <- lapply(seq_len(nrow(x)), function(i) {
      position <- list(
        observationLevel = list(levelName = "plot", levelCode = as.character(plot_number[i])),
        observationLevelRelationships = list(
          list(levelName = "rep", levelCode = as.character(rep_no[i])),
          list(levelName = "block", levelCode = as.character(block_no[i]))
        ),
        entryType = if (is_check[i]) "CHECK" else "TEST"
      )
      if (has_coords) {
        position$positionCoordinateX <- as.character(x$col[i])
        position$positionCoordinateXType <- "GRID_COL"
        position$positionCoordinateY <- as.character(x$row[i])
        position$positionCoordinateYType <- "GRID_ROW"
      }
      # GPS (see add_plot_coordinates()). Breedbase ignores it when the plots
      # are created; save_breedbase() sends it after, with PUT.
      geo <- bb_plot_geojson(x[i, , drop = FALSE])
      if (!is.null(geo)) position$geoCoordinates <- geo
      list(
        observationUnitName = unit_names[i],
        germplasmName = as.character(x$all_entries[i]),
        observationUnitPosition = position
      )
    })
    generator <- if (kind == "imported") "agricolae" else "breedtools"
    study <- list(
      studyName = g$name,
      studyDescription = description,
      studyType = trial_type,
      seasons = list(as.character(year)),
      experimentalDesign = list(PUI = design_type, description = paste(design_type, "design generated by", generator)),
      additionalInfo = list(design_generated_by = generator)
    )
    if (kind == "imported") {
      study$experimentalDesign$description <- "Design imported from an agricolae field book"
      study$additionalInfo$design_converted_by <- "breedtools"
      if (!is.null(design$parameters$design)) {
        study$additionalInfo$agricolae_design <- design$parameters$design
      }
    }
    list(study = study, plots = plots, breedbase_location = g$breedbase_location)
  })
}
