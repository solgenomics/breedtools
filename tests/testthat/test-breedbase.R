# A small fake Breedbase server for testing the BrAPI client without a
# network. It keeps state (folders, studies, plots) and checks the same
# things the real server checks (see CXGN::BrAPI::v2 in solgenomics/sgn).
reset_found_cache <- function() rm(list = ls(breedtools:::bb_found_cache), envir = breedtools:::bb_found_cache)

fake_breedbase <- function(accessions = c(paste0("G", 1:200), paste0("C", 1:6), LETTERS[1:4]),
                           existing_studies = character(0), password = "secret", max_search = Inf,
                           slow_writes = character(0), lost_writes = character(0), site_auth = NULL,
                           statistics = NULL) {
  reset_found_cache()
  state <- new.env()
  state$moved <- list()
  state$folders <- list(list(trialDbId = "500", trialName = "Old folder", programDbId = "10"))
  state$studies <- lapply(existing_studies, function(n) list(studyDbId = "900", studyName = n))
  state$plots <- list()
  state$observations <- list()
  state$variables <- list(list(observationVariableDbId = "70001", observationVariableName = "fresh root yield|CO_334:0000013"),
                          list(observationVariableDbId = "70002", observationVariableName = "dry matter content percentage|CO_334:0000092"),
                          list(observationVariableDbId = "70003", observationVariableName = "fresh root yield spatially adjusted|COMP:0000123"))
  state$analyses <- list()
  # SGNSTAT terms and composed traits (by component ids), like the chado cvterm tables.
  state$statistics <- list(list(cvterm_id = 80054, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Adjusted Means from Spatial Correction using SpATS R",
                                accession = "SGNSTAT:0000054"),
                           list(cvterm_id = 80036, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Adjusted means from BLUEs using LMER R", accession = "SGNSTAT:0000036"),
                           list(cvterm_id = 80034, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Adjusted means from BLUPs using LMER R", accession = "SGNSTAT:0000034"),
                           list(cvterm_id = 80037, cv_name = "SGNStatistics_ontology",
                                cvterm_name = "Phenotypic BLUEs using LMER R", accession = "SGNSTAT:0000037"))
  if (!is.null(statistics)) state$statistics <- statistics
  state$composed <- list()
  state$requests <- list()
  state$next_id <- 1000
  # slow_writes: the write is done but the gateway answers 504 (once).
  # lost_writes: the gateway answers 504 and the write is not done (once).
  state$slow <- slow_writes
  state$lost <- lost_writes
  ok <- function(data = NULL, extra = list()) {
    httr2::response_json(200, body = c(list(metadata = list(pagination = list(totalPages = 1), status = list()),
                                     result = list(data = data)), extra))
  }
  err <- function(code, msg) {
    httr2::response_json(code, body = list(metadata = list(status = list(list(message = msg, messageType = "ERROR"))),
                                    result = list()))
  }
  mock <- function(req) {
    url <- httr2::url_parse(req$url)
    path <- sub("^/brapi/v2/", "", url$path)
    method <- req$method %||% if (!is.null(req$body)) "POST" else "GET"
    body <- req$body$data
    auth <- httr2::req_get_headers(req, redacted = "reveal")$Authorization
    state$requests[[length(state$requests) + 1]] <- list(method = method, path = path, body = body, auth = auth)
    if (path == "token") {
      # The login is a form: its values are marked AsIs.
      if (identical(as.character(unclass(body$password)), password)) {
        return(httr2::response_json(200, body = list(access_token = "TOKEN123", userDisplayName = "Jane Doe")))
      }
      return(httr2::response_json(200, body = list(access_token = NULL,
                                            metadata = list(status = list(list(message = "Incorrect password"))))))
    }
    # Breedbase AJAX (not BrAPI): login token as sgn_session_id; errors as
    # {"error": ...} with HTTP 200, like CXGN's controllers.
    if (grepl("^/?ajax/", path)) {
      path <- sub("^/", "", path)
      # A test site with a web password (nginx basic auth) in front of the website.
      if (!is.null(site_auth)) {
        want <- paste("Basic", jsonlite::base64_enc(paste0(site_auth[1], ":", site_auth[2])))
        got <- httr2::req_get_headers(req, redacted = "reveal")$Authorization %||% ""
        if (!identical(got, want)) {
          return(httr2::response(401, headers = list(`WWW-Authenticate` = 'Basic realm="Restricted Test Site"',
                                                     `Content-Type` = "text/html"),
                                 body = charToRaw("<html>401 Authorization Required</html>")))
        }
      }
      form <- if (is.null(req$body)) {
        list()
      } else if (identical(req$body$type, "multipart")) {
        lapply(req$body$data, function(x) if (inherits(x, "form_file")) x else as.character(x))
      } else {
        lapply(req$body$data, function(x) utils::URLdecode(as.character(x)))
      }
      cookie <- httr2::req_get_headers(req, redacted = "reveal")$Cookie %||% ""
      session <- form$sgn_session_id %||% url$query$sgn_session_id %||% sub("^.*sgn_session_id=([^;]*).*$", "\\1", cookie)
      if (!identical(session, "TOKEN123")) return(httr2::response_json(200, body = list(error = "You must be logged in to do this!")))
      if (path == "ajax/breeder/check_status") {
        return(httr2::response_json(200, body = list(timestamp = "2026-11-04")))
      }
      if (path == "ajax/breeder/refresh") {
        state$refreshed <- c(state$refreshed, url$query$matviews)
        return(httr2::response_json(200, body = list(message = "refreshing")))
      }
      if (path == "ajax/onto/match") {
        q <- url$query
        hits <- Filter(function(s) grepl(tolower(q$term_name), tolower(s$cvterm_name), fixed = TRUE) &&
                         startsWith(s$accession, paste0(q$db_name, ":")), state$statistics)
        return(httr2::response_json(200, body = hits))
      }
      if (path == "ajax/onto/get_trait_from_exact_components") {
        ids <- regmatches(req$url, gregexpr("ids(%5B%5D|\\[\\])=[0-9]+", req$url))[[1]]
        ids <- sort(sub("^.*=", "", ids))
        hit <- Filter(function(cp) identical(sort(cp$components), ids), state$composed)
        if (!length(hit)) return(httr2::response_json(200, body = list(error = "No exact matches found.")))
        return(httr2::response_json(200, body = list(trait_id = as.integer(hit[[1]]$id))))
      }
      if (path == "ajax/onto/store_composed_term") {
        new <- jsonlite::fromJSON(form$new_trait_names, simplifyVector = FALSE)
        for (nm in names(new)) {
          comps <- strsplit(new[[nm]], ",")[[1]]
          if (length(comps) < 2) return(httr2::response_json(200, body = list(error = "less than 2 components")))
          id <- as.character(state$next_id <- state$next_id + 1)
          state$composed[[length(state$composed) + 1]] <- list(id = id, name = nm, components = comps)
        }
        return(httr2::response_json(200, body = list(success = "Saved new trait", names = as.list(names(new)))))
      }
      if (path %in% c("ajax/phenotype/upload_verify/spreadsheet", "ajax/phenotype/upload_store/spreadsheet")) {
        # Like Breedbase: units may be analysis plots ("<analysis>_<stock>"), traits must exist,
        # existing values are warned about unless overwritten.
        f <- utils::read.csv(form$upload_spreadsheet_phenotype_file_input$path, check.names = FALSE, colClasses = "character")
        units <- unlist(lapply(state$analyses, function(a) paste0(a$name, "_", names(a$values))))
        bad <- setdiff(f$observationunit_name, units)
        if (length(bad)) return(httr2::response_json(200, body = list(error = list(paste("Plots not in database:", paste(bad, collapse = ","))), success = list())))
        known <- c(vapply(state$variables, `[[`, "", "observationVariableName"),
                   vapply(state$composed, function(cp) cp$name, ""),
                   vapply(state$composed, function(cp) paste0(cp$name, "|COMP:000", cp$id), ""))
        badt <- setdiff(names(f)[-1], known)
        if (length(badt)) return(httr2::response_json(200, body = list(error = list(paste("Traits not in database:", paste(badt, collapse = ","))), success = list())))
        key <- unlist(lapply(names(f)[-1], function(tr) paste(f$observationunit_name[nzchar(f[[tr]])], tr)))
        existing <- intersect(key, names(state$analysis_values))
        warn <- if (length(existing) && is.null(form$phenotype_upload_overwrite_values)) as.list(paste(existing, "already has a <strong>different value</strong>")) else list()
        if (grepl("store", path)) {
          if (length(warn)) return(httr2::response_json(200, body = list(error = list("values exist"), success = list())))
          for (tr in names(f)[-1]) for (i in which(nzchar(f[[tr]]))) state$analysis_values[[paste(f$observationunit_name[i], tr)]] <- f[[tr]][i]
        }
        return(httr2::response_json(200, body = list(success = list("File data verified."), warning = warn, error = list())))
      }
      if (grepl("^ajax/folder/[0-9]+/associate/parent/[0-9]+$", path)) {
        ids <- regmatches(path, gregexpr("[0-9]+", path))[[1]]
        state$moved[[ids[1]]] <- ids[2]
        return(httr2::response_json(200, body = list(success = 1)))
      }
      if (path == "ajax/analysis/store/spreadsheet") {
        # Like Breedbase: name unique; "stock_names" then traits; stocks are
        # plots or accessions depending on "Analysis of"; traits must exist.
        nm <- form$upload_new_analysis_name
        if (nm %in% c(vapply(state$analyses, `[[`, "", "name"), vapply(state$studies, `[[`, "", "studyName"))) {
          return(httr2::response_json(200, body = list(error = paste("An analysis with name", nm, "already exists"))))
        }
        f <- utils::read.csv(form$upload_new_analysis_file$path, check.names = FALSE, colClasses = "character")
        if (names(f)[1] != "stock_names") return(httr2::response_json(200, body = list(error_messages = list("Header row must contain stock_names"))))
        known <- c(vapply(state$variables, `[[`, "", "observationVariableName"),
                   vapply(state$composed, function(cp) paste0(cp$name, "|COMP:000", cp$id), ""))
        bad_traits <- setdiff(names(f)[-1], known)
        if (length(bad_traits)) return(httr2::response_json(200, body = list(error_messages = list(paste("The following traits are not in the database:", paste(bad_traits, collapse = ","))))))
        plots_level <- form$upload_new_analysis_result_values_type == "analysis_result_values_match_precomputed_design"
        valid <- if (plots_level) vapply(state$plots, `[[`, "", "observationUnitName") else accessions
        bad_stocks <- setdiff(f$stock_names, valid)
        if (length(bad_stocks)) return(httr2::response_json(200, body = list(error_messages = list(paste("The following stocks are not in the database:", paste(bad_stocks, collapse = ","))))))
        # Like Breedbase: each trait is stored composed with the statistic of
        # the analysis, "trait|statistic|COMP:id" (created when new).
        stat <- form$upload_new_analysis_statistical_ontology_term
        stat_name <- sub("\\|[^|]*$", "", stat)
        stat_id <- Filter(function(x) identical(x$accession, sub("^.*\\|", "", stat)), state$statistics)
        stored <- f
        for (j in seq_along(names(f))[-1]) {
          t <- names(f)[j]
          if (grepl("\\|COMP:[0-9]+$", t)) next
          cname <- paste(sub("\\|[^|]*$", "", t), stat_name, sep = "|")
          hit <- Filter(function(cp) identical(cp$name, cname), state$composed)
          if (!length(hit)) {
            vid <- Filter(function(v) identical(v$observationVariableName, t), state$variables)[[1]]$observationVariableDbId
            cid <- as.character(state$next_id <- state$next_id + 1)
            state$composed[[length(state$composed) + 1]] <- list(
              id = cid, name = cname, components = c(vid, if (length(stat_id)) as.character(stat_id[[1]]$cvterm_id)))
            hit <- list(state$composed[[length(state$composed)]])
          }
          names(stored)[j] <- paste0(cname, "|COMP:000", hit[[1]]$id)
        }
        id <- as.character(state$next_id <- state$next_id + 1)
        state$analyses[[length(state$analyses) + 1]] <- list(id = id, name = nm, form = form, file = f, stored = stored,
                                                             values = stats::setNames(as.list(f$stock_names), f$stock_names))
        return(httr2::response_json(200, body = list(success = 1, analysis_id = as.integer(id), model_id = 5)))
      }
      if (path == "ajax/analysis/store/json") {
        if (form$analysis_name %in% c(vapply(state$analyses, `[[`, "", "name"), vapply(state$studies, `[[`, "", "studyName"))) {
          return(httr2::response_json(200, body = list(error = paste("An analysis with name", form$analysis_name, "already exists"))))
        }
        values <- jsonlite::fromJSON(form$analysis_result_values, simplifyVector = FALSE)
        traits <- unlist(jsonlite::fromJSON(form$analysis_trait_names))
        known <- vapply(state$variables, `[[`, "", "observationVariableName")
        if (!all(traits %in% known)) return(httr2::response_json(200, body = list(error = "trait not found")))
        id <- as.character(state$next_id <- state$next_id + 1)
        state$analyses[[length(state$analyses) + 1]] <- list(id = id, name = form$analysis_name, form = form,
                                                             values = values)
        return(httr2::response_json(200, body = list(success = 1, analysis_id = as.integer(id), model_id = 5)))
      }
      if (path == "ajax/analyses/by_user") {
        rows <- lapply(state$analyses, function(a) list(paste0('<a href="/analyses/', a$id, '">', a$name, "</a>"), "d"))
        return(httr2::response_json(200, body = list(data = rows)))
      }
      if (grepl("^ajax/analysis/[0-9]+/retrieve$", path)) {
        id <- sub("^ajax/analysis/([0-9]+)/retrieve$", "\\1", path)
        a <- Filter(function(a) a$id == id, state$analyses)[[1]]
        rows <- c(list(list("header")), lapply(names(a$values), function(n) list(n)))
        # Traits: those of the upload file and those added later.
        added <- unique(sub("^[^ ]+ ", "", grep(paste0("^", a$name, "_"), names(state$analysis_values), value = TRUE)))
        traits <- unique(c(if (!is.null(a$stored)) names(a$stored)[-1], added))
        # Like cassava-test: the summary's trait list can be empty.
        return(httr2::response_json(200, body = list(analysis_name = a$name, data = rows, traits = list())))
      }
      if (grepl("^ajax/breeders/trial/[0-9]+/traits_assayed$", path)) {
        # Like Breedbase: [[ [trait id, "name|ONTO:id", components, count], ... ]] (wrapped once more).
        id <- sub("^ajax/breeders/trial/([0-9]+)/traits_assayed$", "\\1", path)
        a <- Filter(function(a) a$id == id, state$analyses)[[1]]
        added <- unique(sub("^[^ ]+ ", "", grep(paste0("^", a$name, "_"), names(state$analysis_values), value = TRUE)))
        traits <- unique(c(if (!is.null(a$stored)) names(a$stored)[-1], added))
        all_vars <- c(lapply(state$variables, function(v) list(id = v$observationVariableDbId, name = v$observationVariableName)),
                      lapply(state$composed, function(cp) list(id = cp$id, name = paste0(cp$name, "|COMP:000", cp$id))))
        rows <- lapply(traits, function(t) {
          v <- Filter(function(v) identical(v$name, t), all_vars)[[1]]
          list(as.integer(v$id), v$name, NULL, 10L)
        })
        return(httr2::response_json(200, body = list(traits_assayed = list(rows))))
      }
      if (grepl("^ajax/breeders/trial/[0-9]+/trait_phenotypes$", path)) {
        # Like Breedbase: the phenotype matrix of one trait (only the units
        # with a value), with display = analysis_instance for analyses.
        id <- sub("^ajax/breeders/trial/([0-9]+)/trait_phenotypes$", "\\1", path)
        a <- Filter(function(a) a$id == id, state$analyses)[[1]]
        q <- url$query
        if (!identical(q$display, "analysis_instance")) return(httr2::response_json(200, body = list(data = list())))
        all_vars <- c(lapply(state$variables, function(v) list(id = v$observationVariableDbId, name = v$observationVariableName)),
                      lapply(state$composed, function(cp) list(id = cp$id, name = paste0(cp$name, "|COMP:000", cp$id))))
        tname <- Filter(function(v) identical(as.character(v$id), as.character(q$trait)), all_vars)[[1]]$name
        header <- list("studyName", "germplasmName", "observationUnitName", "replicate", "blockNumber", "plotNumber",
                       "rowNumber", "colNumber", tname)
        rows <- list(header)
        for (stock in a$stored$stock_names) {
          unit <- paste0(a$name, "_", stock)
          v <- if (tname %in% names(a$stored)) a$stored[[tname]][a$stored$stock_names == stock] else ""
          v2 <- state$analysis_values[[paste(unit, tname)]]
          if (!is.null(v2)) v <- v2
          if (!nzchar(v)) next
          fp <- Filter(function(u) u$observationUnitName == stock, state$plots)
          pos <- if (length(fp)) fp[[1]]$observationUnitPosition else NULL
          lvl <- function(n) {
            r <- Filter(function(x) x$levelName == n, pos$observationLevelRelationships %||% list())
            if (length(r)) r[[1]]$levelCode else NULL
          }
          rows[[length(rows) + 1]] <- list(a$name, if (length(fp)) fp[[1]]$germplasmName else stock, unit, lvl("rep"),
                                           lvl("block"), pos$observationLevel$levelCode, pos$positionCoordinateY,
                                           pos$positionCoordinateX, v)
        }
        return(httr2::response_json(200, body = list(status = "success", data = rows)))
      }
      return(httr2::response_json(404, body = list(error = "unknown ajax call")))
    }
    if (!identical(auth, "Bearer TOKEN123")) return(err(401, "You must be logged in"))
    if (path == "programs") return(ok(list(list(programDbId = "10", programName = "IITA"),
                                         list(programDbId = "11", programName = "NRCRI"))))
    if (path == "locations") return(ok(list(list(locationDbId = "1", locationName = "Ibadan"),
                                          list(locationDbId = "2", locationName = "Ubiaja"),
                                          list(locationDbId = "3", locationName = "Ikenne"))))
    if (path == "trials" && method == "GET") {
      return(ok(Filter(function(f) f$programDbId == url$query$programDbId, state$folders)))
    }
    if (path == "trials" && method == "POST") {
      taken <- c(vapply(state$folders, `[[`, "", "trialName"), vapply(state$studies, `[[`, "", "studyName"))
      if (body[[1]]$trialName %in% taken) {
        return(err(409, sprintf("A trial or study with the name %s already exists in the database.", body[[1]]$trialName)))
      }
      id <- as.character(state$next_id <- state$next_id + 1)
      state$folders[[length(state$folders) + 1]] <- list(trialDbId = id, trialName = body[[1]]$trialName,
                                                         programDbId = body[[1]]$programDbId)
      return(ok(list(list(trialDbId = id))))
    }
    if (path == "search/germplasm") {
      # A slow server: too many names in one search ends in a gateway timeout.
      if (length(unlist(body$germplasmNames)) > max_search) {
        return(httr2::response(504, headers = list(`Content-Type` = "text/html"), body = charToRaw("<html>Gateway Timeout</html>")))
      }
      hits <- intersect(unlist(body$germplasmNames), accessions)
      return(ok(lapply(hits, function(n) list(germplasmDbId = as.character(5000 + match(n, accessions)),
                                              germplasmName = n))))
    }
    if (path == "search/observationunits") {
      hits <- Filter(function(u) u$observationUnitName %in% unlist(body$observationUnitNames), state$plots)
      return(ok(hits))
    }
    if (path == "search/variables") {
      wanted <- unlist(body$observationVariableNames)
      hits <- Filter(function(v) sub("\\|.*$", "", v$observationVariableName) %in% wanted, state$variables)
      composed <- Filter(function(cp) cp$name %in% wanted, state$composed)
      hits <- c(hits, lapply(composed, function(cp) list(observationVariableDbId = cp$id,
                                                          observationVariableName = paste0(cp$name, "|COMP:000", cp$id))))
      return(ok(hits))
    }
    if (path == "search/observations") {
      unit_study <- vapply(state$plots, function(u) as.character(u$studyDbId), "")
      unit_name <- vapply(state$plots, `[[`, "", "observationUnitName")
      names(unit_study) <- names(unit_name) <- vapply(state$plots, `[[`, "", "observationUnitDbId")
      hits <- Filter(function(o) {
        (is.null(body$observationUnitDbIds) || o$observationUnitDbId %in% unlist(body$observationUnitDbIds)) &&
          (is.null(body$studyDbIds) || unit_study[[o$observationUnitDbId]] %in% unlist(body$studyDbIds)) &&
          o$observationVariableDbId %in% unlist(body$observationVariableDbIds)
      }, state$observations)
      hits <- lapply(hits, function(o) c(o, list(observationUnitName = unit_name[[o$observationUnitDbId]])))
      return(ok(hits))
    }
    if (path == "observations" && method == "POST") {
      # Like Breedbase: the three fields are required, the plot and trait must
      # exist, and a plot-trait-time combination is stored once.
      unit_ids <- vapply(state$plots, `[[`, "", "observationUnitDbId")
      var_ids <- c(vapply(state$variables, `[[`, "", "observationVariableDbId"),
                   vapply(state$composed, `[[`, "", "id"))
      seen <- vapply(state$observations, function(o) paste(o$observationUnitDbId, o$observationVariableDbId,
                                                            o$observationTimeStamp %||% ""), "")
      for (o in body) {
        if (is.null(o$observationUnitDbId) || is.null(o$observationVariableDbId) || is.null(o$value)) {
          return(err(400, "An observation object is missing a required field"))
        }
        if (!o$observationUnitDbId %in% unit_ids) return(err(400, "observationUnitDbId not found"))
        if (!o$observationVariableDbId %in% var_ids) return(err(400, "observationVariableDbId not found"))
        key <- paste(o$observationUnitDbId, o$observationVariableDbId, o$observationTimeStamp %||% "")
        if (key %in% seen) return(err(400, "The combination already exists in the database"))
        seen <- c(seen, key)
      }
      for (i in seq_along(body)) body[[i]]$observationDbId <- as.character(900000 + length(state$observations) + i)
      state$observations <- c(state$observations, unname(body))
      return(ok(list()))
    }
    if (path == "observations" && method == "PUT") {
      # Like Breedbase: {observationDbId: observation}, stored with overwrite.
      ids <- vapply(state$observations, function(o) o$observationDbId %||% "", "")
      for (id in names(body)) {
        i <- match(id, ids)
        if (is.na(i)) return(err(400, "observationDbId not found"))
        state$observations[[i]]$value <- body[[id]]$value
      }
      return(ok(list()))
    }
    if (path == "search/studies") {
      hits <- Filter(function(s) s$studyName %in% unlist(body$studyNames), state$studies)
      return(ok(hits))
    }
    if (path == "studies" && method == "POST") {
      s <- body[[1]]
      allowed <- c("CRD", "Alpha", "MAD", "Lattice", "Augmented", "RCBD", "p-rep", "splitplot", "stripplot",
                   "greenhouse", "Westcott", "Analysis")
      if (is.null(s$trialDbId) || is.null(s$locationDbId) || is.null(s$studyType) || is.null(s$experimentalDesign$PUI)) {
        return(err(400, "missing required field"))
      }
      if (!s$experimentalDesign$PUI %in% allowed) return(err(400, "Experimental Design not supported"))
      if (s$studyName %in% vapply(state$studies, `[[`, "", "studyName")) return(err(409, "Study exists"))
      id <- as.character(state$next_id <- state$next_id + 1)
      state$studies[[length(state$studies) + 1]] <- c(s, list(studyDbId = id))
      return(ok(list(list(studyDbId = id))))
    }
    if (path == "observationunits" && method == "POST") {
      nums <- vapply(body, function(u) u$observationUnitPosition$observationLevel$levelCode, "")
      if (anyDuplicated(paste(vapply(body, `[[`, "", "studyDbId"), nums))) return(err(422, "Plot number duplicated"))
      for (u in body) {
        # Like Breedbase: the germplasmDbId is required and must match.
        id <- u$germplasmDbId %||% ""
        if (!grepl("^[0-9]+$", id)) return(err(404, "Verify Germplasm Id."))
        if (!identical(accessions[as.integer(id) - 5000], u$germplasmName)) return(err(404, "germplasm not found"))
        state$next_id <- state$next_id + 1
        u$observationUnitDbId <- as.character(state$next_id)
        u$observationUnitPosition$geoCoordinates <- NULL   # like Breedbase: POST does not store the GPS
        state$plots[[length(state$plots) + 1]] <- u
      }
      return(ok(list()))
    }
    if (path == "observationunits" && method == "PUT") {
      # Like Breedbase: a map of observationUnitDbId to the unit; the level
      # must be plot; the GPS replaces the stored one.
      for (id in names(body)) {
        u <- body[[id]]
        if (!identical(u$observationUnitPosition$observationLevel$levelName, "plot")) {
          return(err(400, "Only 'plot', 'plant' or 'tissue_sample' allowed for observation level"))
        }
        k <- which(vapply(state$plots, `[[`, "", "observationUnitDbId") == id)
        if (!length(k)) return(err(404, "observation unit not found"))
        state$plots[[k]]$observationUnitPosition$geoCoordinates <- u$observationUnitPosition$geoCoordinates
      }
      return(ok(list()))
    }
    if (path == "observationunits" && method == "GET") {
      return(ok(Filter(function(u) u$studyDbId == url$query$studyDbId, state$plots)))
    }
    err(404, paste("unknown call", path))
  }
  timeout <- httr2::response(504, headers = list(`Content-Type` = "text/html"),
                             body = charToRaw("<html>Gateway Timeout</html>"))
  mock_with_timeouts <- function(req) {
    path <- sub("^/brapi/v2/", "", httr2::url_parse(req$url)$path)
    is_write <- !is.null(req$body) && !grepl("^search|^token", path)
    if (is_write && path %in% state$lost) {
      state$lost <- setdiff(state$lost, path)
      return(timeout)
    }
    res <- mock(req)
    if (is_write && path %in% state$slow) {
      state$slow <- setdiff(state$slow, path)
      return(timeout)
    }
    res
  }
  list(mock = mock_with_timeouts, state = state)
}

test_that("converted field books save as a new trial and round-trip through Breedbase", {
  skip_if_not_installed("httr2")
  # Small field book in agricolae's RCBD format, with an arbitrary treatment column name.
  book <- data.frame(plots = c(102, 101, 202, 201), block = c(1, 1, 2, 2),
                     genotype = c("G2", "C1", "C1", "G2"))
  converted <- as_breedbase_design(book, design_type = "RCBD", checks = "C1")
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    invisible(capture.output(checked <- save_breedbase(converted, con, program = "IITA",
      trial_name = "AG_2026", year = 2026, location = "Ibadan",
      description = "Trial generated with agricolae", dry_run = TRUE)))
    expect_length(fake$state$studies, 0)
    expect_length(fake$state$plots, 0)
    invisible(capture.output(saved <- save_breedbase(converted, con, program = "IITA",
      trial_name = "AG_2026", year = 2026, location = "Ibadan",
      description = "Trial generated with agricolae", confirm = FALSE)))
    expect_length(fake$state$studies, 1)
    expect_identical(fake$state$studies[[1]]$experimentalDesign$PUI, "RCBD")
    expect_equal(saved$plots_sent, 4)
    expect_equal(saved$plots_found, 4)
    expect_identical(fake$state$studies[[1]]$trialDbId, "10")
    readback <- breedbase_trial_plots(con, "AG_2026")
    expected <- phenotype_template(converted, "AG_2026", location = "Ibadan")
    expected$trial_name <- NULL
    expect_equal(readback[names(expected)], expected)
    expect_error(capture.output(save_breedbase(converted, con, program = "IITA",
      trial_name = "AG_2026", year = 2026, location = "Ibadan",
      description = "Trial generated with agricolae", confirm = FALSE)), "already exists")
  })
})

test_that("prepared agricolae trials save using recorded metadata and allow explicit overrides", {
  skip_if_not_installed("httr2")
  raw <- list(book = data.frame(plots = c(101, 102, 201, 202), block = c(1, 1, 2, 2),
                                entries = c("G1", "C1", "C1", "G1")),
              parameters = list(design = "rcbd"))
  trial <- prepare_breedbase_trial(raw, trial_name = "PREPARED_2026", year = 2026,
    location = "Ibadan", program = "IITA", folder = "Old folder",
    description = "Prepared RCBD", checks = "C1", prompt = FALSE)
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    invisible(capture.output(save_breedbase(trial, con, dry_run = TRUE)))
    expect_length(fake$state$studies, 0)
    invisible(capture.output(saved <- save_breedbase(trial, con, confirm = FALSE)))
    expect_equal(saved$plots_found, 4)
    expect_identical(fake$state$studies[[1]]$studyName, "PREPARED_2026")
    expect_identical(fake$state$studies[[1]]$trialDbId, "500")
    actual <- breedbase_trial_plots(con, "PREPARED_2026")
    expected <- phenotype_template(trial)
    expected$trial_name <- NULL
    expect_equal(actual[names(expected)], expected)
    invisible(capture.output(save_breedbase(trial, con, trial_name = "OVERRIDE_2027", year = 2027,
      folder = NULL, location = "Ubiaja", description = "Another field", confirm = FALSE)))
    expect_identical(fake$state$studies[[2]]$studyName, "OVERRIDE_2027")
    expect_identical(fake$state$studies[[2]]$trialDbId, "10")
    expect_identical(fake$state$studies[[2]]$locationDbId, "2")
    expect_identical(fake$state$studies[[2]]$seasons, list("2027"))
    expect_identical(fake$state$studies[[2]]$studyDescription, "Another field")
  })
})

test_that("dry run without a connection returns the BrAPI records", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  payload <- save_breedbase(rc, trial_name = "AYT_2026", year = 2026, location = "Ibadan", dry_run = TRUE,
  description = "Test trial")
  expect_length(payload$studies, 1)
  study <- payload$studies[[1]]
  expect_equal(study$studyName, "AYT_2026")
  expect_equal(study$experimentalDesign$PUI, "Augmented")
  expect_equal(study$seasons[[1]], "2026")
  plots <- payload$plots[[1]]
  expect_length(plots, 144)
  expect_equal(sum(vapply(plots, function(p) p$observationUnitPosition$entryType == "CHECK", logical(1))), 36)
  p1 <- plots[[1]]
  first <- rc$design[rc$design$plots == 1, ]
  expect_equal(p1$observationUnitName,
               paste("AYT_2026", "Ibadan", first$all_entries, first$block, first$rep, sep = "_"))
  names_all <- vapply(plots, `[[`, "", "observationUnitName")
  expect_false(anyDuplicated(names_all) > 0)
  expect_equal(p1$observationUnitPosition$observationLevel$levelName, "plot")
  expect_equal(p1$observationUnitPosition$positionCoordinateY, as.character(rc$design$row[rc$design$plots == 1]))
  expect_equal(vapply(p1$observationUnitPosition$observationLevelRelationships, `[[`, "", "levelName"), c("rep", "block"))
  expect_error(save_breedbase(rc, trial_name = "x", year = 2026, location = "Ibadan",
  description = "Test trial"), "breedbase_login")
  expect_error(save_breedbase(rc, trial_name = "x", year = 2026, location = "Ibadan", design_type = "RRC",
                              dry_run = TRUE,
                              description = "Test trial"), "Breedbase accepts")
})

test_that("p-rep designs become one trial per location; fillers are skipped", {
  prep <- augmented_prep_design(paste0("G", 1:42), c("North", "South"), block_size = 6, field_cols = 7, seed = 1)
  payload <- save_breedbase(prep, trial_name = "PREP_2026", year = 2026,
                            location = c(North = "Ibadan", South = "Ubiaja"), dry_run = TRUE,
                            description = "Test trial")
  expect_equal(vapply(payload$studies, `[[`, "", "studyName"), c("PREP_2026_North", "PREP_2026_South"))
  expect_equal(vapply(payload$plots, length, integer(1)), as.vector(table(prep$design$location)[c("North", "South")]))
  expect_equal(payload$studies[[1]]$experimentalDesign$PUI, "p-rep")
  nums <- vapply(payload$plots[[2]], function(p) p$observationUnitPosition$observationLevel$levelCode, "")
  expect_false(anyDuplicated(nums) > 0)

  rrc <- resolvable_row_column_design(paste0("P", 1:19), 2, 4, 5, n_starts = 1, n_iter = 200, seed = 1)
  rp <- save_breedbase(rrc, trial_name = "RRC_2026", year = 2026, location = "Ibadan", dry_run = TRUE,
  description = "Test trial")
  expect_length(rp$plots[[1]], 38)
  expect_false("FILLER" %in% vapply(rp$plots[[1]], `[[`, "", "germplasmName"))
  expect_equal(rp$studies[[1]]$experimentalDesign$PUI, "Alpha")
})

test_that("login, checks and saving work against a fake Breedbase", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase()
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    expect_error(breedbase_login("https://bb.example.org", "jdoe", password = "wrong"), "Login failed")
    con <- suppressMessages(breedbase_login("bb.example.org/", "jdoe", password = "secret"))
    expect_s3_class(con, "breedbase_connection")
    expect_equal(con$url, "https://bb.example.org")
    expect_false(any(grepl("secret", capture.output(print(con)))))
    expect_null(con$password)

    expect_equal(breedbase_programs(con)$program, c("IITA", "NRCRI"))
    expect_equal(breedbase_folders(con, "IITA")$folder, "Old folder")
    expect_error(breedbase_folders(con, "Nope"), "Programs available: IITA, NRCRI")

    # Dry run with a connection checks everything and writes nothing.
    n_before <- length(fake$state$requests)
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                                create_folder = TRUE, dry_run = TRUE,
                                                description = "Test trial"))
    expect_true(any(grepl("Nothing was written", out)))
    posts <- Filter(function(r) r$method == "POST" && !grepl("^search|token", r$path), fake$state$requests)
    expect_length(posts, 0)

    # Real save.
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                                create_folder = TRUE, confirm = FALSE,
                                                description = "Test trial"))
    expect_equal(res$plots_sent, 144)
    expect_equal(res$plots_found, 144)
    expect_match(res$trial_page, "^https://bb.example.org/breeders/trial/")
    expect_true("AYT_2026" %in% vapply(fake$state$folders, `[[`, "", "trialName"))
    study <- fake$state$studies[[1]]
    expect_equal(study$locationDbId, "1")
    expect_equal(study$trialDbId, fake$state$folders[[2]]$trialDbId)
    expect_true(all(vapply(Filter(function(r) r$path != "token" && !grepl("ajax/", r$path), fake$state$requests),
                           function(r) identical(r$auth, "Bearer TOKEN123"), logical(1))))

    # Saving again with the same name is refused before writing.
    expect_error(capture.output(save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                               trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                               confirm = FALSE,
                                               description = "Test trial")), "already exists")
  })
})

test_that("clear errors for missing accessions, folders, locations and logins", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(accessions = paste0("G", 1:100))
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "Ibadan", confirm = FALSE,
    description = "Test trial")),
                 "12 of 112 entries/checks are not accessions")
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Missing", "T1", 2026, "Ibadan", confirm = FALSE,
    description = "Test trial")),
                 "create_folder = TRUE")
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "Mars", confirm = FALSE,
    description = "Test trial")),
                 "Location 'Mars' was not found")
    bad <- con
    bad$token <- "expired"
    expect_error(breedbase_programs(bad), "HTTP 401.*log in again")
  })
  expect_false(any(grepl("secret", unlist(lapply(fake$state$requests, function(r) r$auth)))))
})

test_that("plot names: automatic pattern, own pattern, and user-given names", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  p <- save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                      plot_name_pattern = "{trial}-{plot}-r{row}c{col}",
                      description = "Test trial")
  expect_equal(p$plots[[1]][[1]]$observationUnitName, sprintf("T-1-r%dc%d", d$row[d$plots == 1], d$col[d$plots == 1]))
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                              plot_name_pattern = "{trial}_{plant}",
                              description = "Test trial"), "Unknown field")

  own <- paste0("MYPLOT-", seq_len(nrow(d)))
  p <- save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE, plot_names = own,
  description = "Test trial")
  sent <- vapply(p$plots[[1]], `[[`, "", "observationUnitName")
  expect_setequal(sent, own)
  d$my_names <- paste0("X", d$plots)
  rc2 <- rc
  rc2$design <- d
  p <- save_breedbase(rc2, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE, plot_names = "my_names",
  description = "Test trial")
  expect_equal(p$plots[[1]][[1]]$observationUnitName, "X1")
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                              plot_names = rep("A", nrow(d)),
                              description = "Test trial"), "unique")
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE,
                              plot_names = "A",
                              description = "Test trial"), "one name per plot")

  # Extra check replication: the same check twice in a block gets the plot number added.
  extra <- augmented_row_column_design(paste0("E", 1:90), LETTERS[1:4], 12, 12, 4, 4,
                                       checks_per_block = 6, n_candidates = 2, seed = 1)
  expect_message(p <- save_breedbase(extra, trial_name = "X", year = 2026, location = "Ibadan", dry_run = TRUE,
  description = "Test trial"),
                 "plot number was added")
  sent <- vapply(p$plots[[1]], `[[`, "", "observationUnitName")
  expect_false(anyDuplicated(sent) > 0)
})

test_that("plot names already in the database are refused before writing", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase()
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "Ibadan", confirm = FALSE,
    description = "Test trial"))
    expect_error(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T2", 2026, "Ibadan",
                                               plot_names = vapply(fake$state$plots, `[[`, "", "observationUnitName")[
                                                 order(vapply(fake$state$plots, function(u) as.integer(u$observationUnitPosition$observationLevel$levelCode), 1L))],
                                               confirm = FALSE,
                                               description = "Test trial")),
                 "plot names already exist")
  })
})

test_that("a trial description is required and is sent to Breedbase", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2, seed = 1)
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", dry_run = TRUE),
               "Give a `description` of the trial")
  expect_error(save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan", description = "  ",
                              dry_run = TRUE), "Give a `description`")
  p <- save_breedbase(rc, trial_name = "T", year = 2026, location = "Ibadan",
                      description = "Yield trial of 108 lines, Ibadan 2026", dry_run = TRUE)
  expect_equal(p$studies[[1]]$studyDescription, "Yield trial of 108 lines, Ibadan 2026")
})

test_that("breedbase_check reports every problem at once, with suggestions", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(accessions = c(paste0("G", 1:108), paste0("C", 1:4), "TMEB419"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))

    out <- capture.output(res <- breedbase_check(con, accessions = c(paste0("G", 1:100), "tmeb419", "X9", "G1"),
                                                 location = "ibadan", program = "IITA", folder = "New one",
                                                 trial_name = "T1"))
    expect_false(attr(res, "ok"))
    expect_equal(res$status[res$what == "location"], "MISSING")
    expect_match(res$suggestion[res$what == "location"], "did you mean 'Ibadan'")
    expect_equal(res$status[res$what == "folder"], "NEW")
    expect_equal(res$status[res$what == "trial"], "OK")
    acc <- res[res$what == "accession", ]
    expect_match(acc$suggestion[acc$name == "tmeb419"], "did you mean 'TMEB419'")
    expect_match(acc$suggestion[acc$name == "X9"], "add it to Breedbase first")
    expect_equal(acc$status[acc$name == "G1"], c("REPEATED", "OK"))
    expect_true(any(grepl("100 of 102 accessions are in Breedbase", out)))
    expect_false(any(res$what == "plot"))                    # accession names are not plot names
    expect_true(any(grepl("^NOT READY: fix 4 problems \\(1 location, 3 accessions, marked", out)))
    expect_true(any(grepl("^Checked: 1 program, 1 new folder, 1 location, 102 accessions, 1 trial name\\.", out)))

    # A design can be given directly; all good -> ok.
    rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                      n_candidates = 5, seed = 1)
    out <- capture.output(res <- breedbase_check(con, rc, location = "Ibadan"))
    expect_true(attr(res, "ok"))
    expect_true(any(grepl("plot names NOT checked", out)))    # no trial name: plot names unknown
    named <- add_plot_names(rc, "T1", location = "Ibadan")
    out <- capture.output(res <- breedbase_check(con, named, location = "Ibadan", trial_name = "T1"))
    expect_true(any(grepl("^READY: everything checked can be saved .*1 trial name, 144 plot names\\)", out)))
    out <- capture.output(res <- breedbase_check(con, c(paste0("G", 1:108), paste0("C", 1:4)), location = "Ibadan",
                                                 plot_names = c("P1", "P2", "P1")))
    expect_equal(res$status[res$what == "plot" & res$name == "P1"], c("REPEATED", "OK"))
    expect_true(any(grepl("^NOT READY: fix 1 problem \\(1 plot name", out)))
    # No write calls were made by the checks.
    posts <- Filter(function(r) r$method == "POST" && !grepl("^search|token", r$path), fake$state$requests)
    expect_length(posts, 0)
  })
})

test_that("save_breedbase lists location and accession problems together", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(accessions = paste0("G", 1:100))
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    err <- tryCatch(capture.output(save_breedbase(rc, con, "IITA", "Old folder", "T1", 2026, "ibadan",
                                                  confirm = FALSE, description = "Test trial")),
                    error = conditionMessage)
    expect_match(err, "Location 'ibadan' was not found.*did you mean 'Ibadan'")
    expect_match(err, "12 of 112 entries/checks are not accessions")
  })
})

test_that("slow servers: searches that time out are split and retried", {
  skip_if_not_installed("httr2")
  old_options <- options(breedtools.retry_pause = 0)
  on.exit(options(old_options), add = TRUE)
  fake <- fake_breedbase(max_search = 12)
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                                trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                                create_folder = TRUE, dry_run = TRUE, description = "Test trial"))
    expect_true(any(grepl("all 112 entries and checks exist as accessions", out)))
    expect_true(any(grepl("searched 112 accessions in [0-9]+ batches .*split and retried", out)))
    # Found accessions are remembered: a second check searches nothing.
    n_search <- length(Filter(function(r) r$path == "search/germplasm", fake$state$requests))
    out2 <- capture.output(breedbase_check(con, accessions = rc))
    expect_true(any(grepl("all 112 accessions were already found in this session", out2)))
    expect_equal(length(Filter(function(r) r$path == "search/germplasm", fake$state$requests)), n_search)
    sizes <- vapply(Filter(function(r) r$path == "search/germplasm", fake$state$requests),
                    function(r) length(unlist(r$body$germplasmNames)), integer(1))
    expect_true(all(sizes <= 25))
    expect_true(any(sizes > 12))

    # A server that cannot answer even one name gives a clear message.
    dead <- fake_breedbase(max_search = 0)
    httr2::with_mocked_responses(dead$mock, {
      con2 <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
      expect_error(capture.output(breedbase_check(con2, accessions = c("G1", "G2"))),
                   "HTTP 504.*too slow or busy")
    })
  })
})

test_that("timeouts while writing: the result is checked, retried, and a re-run continues", {
  skip_if_not_installed("httr2")
  old_options <- options(breedtools.write_wait = 0, breedtools.write_checks = 2, breedtools.retry_pause = 0)
  on.exit(options(old_options), add = TRUE)
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  save <- function(con) save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
                                       trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
                                       create_folder = TRUE, confirm = FALSE, description = "Test trial")
  n_studies <- function(f) length(f$state$studies)

  # The server finishes every write but the gateway times out: no duplicates.
  fake <- fake_breedbase(slow_writes = c("trials", "studies", "observationunits"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save(con))
  })
  expect_equal(res$plots_found, 144)
  expect_equal(sum(vapply(fake$state$folders, `[[`, "", "trialName") == "AYT_2026"), 1)
  expect_equal(n_studies(fake), 1)
  expect_length(fake$state$plots, 144)
  expect_true(any(grepl("slow to create folder 'AYT_2026'", out)))
  expect_true(any(grepl("confirmed: done", out)))

  # Writes that are lost are sent once more.
  fake <- fake_breedbase(lost_writes = c("trials", "studies", "observationunits"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save(con))
  })
  expect_equal(res$plots_found, 144)
  expect_equal(n_studies(fake), 1)
  expect_true(any(grepl("trying again", out)))

  # A save cut off after the trial was made: the re-run adds the plots.
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    # Lose the plots on every try.
    lose_all <- function(req) {
      if (grepl("/brapi/v2/observationunits$", req$url) && !is.null(req$body)) return(httr2::response(504))
      fake$mock(req)
    }
    httr2::with_mocked_responses(lose_all, {
      expect_error(capture.output(save(con)), "kept timing out.*run save_breedbase\\(\\) again")
    })
    expect_equal(n_studies(fake), 1)
    expect_length(fake$state$plots, 0)
    out <- capture.output(res <- save(con))
  })
  expect_true(any(grepl("exists without plots", out)))
  expect_equal(res$plots_found, 144)
  expect_equal(n_studies(fake), 1)
})

test_that("an existing trial is reported; an existing folder is used", {
  skip_if_not_installed("httr2")
  old_options <- options(breedtools.write_wait = 0, breedtools.write_checks = 1)
  on.exit(options(old_options), add = TRUE)
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", folder = "Old folder", trial_name = "T_2026",
                                  year = 2026, location = "Ibadan", confirm = FALSE, description = "Test"))
    # The folder existed: no new folder, the trial is inside it.
    expect_equal(sum(vapply(fake$state$folders, `[[`, "", "trialName") == "Old folder"), 1)
    expect_equal(fake$state$studies[[1]]$trialDbId, "500")

    out <- capture.output(res <- breedbase_check(con, accessions = rc, location = "Ibadan",
                                                 program = "IITA", folder = "Old folder", trial_name = "T_2026"))
    expect_false(attr(res, "ok"))
    expect_equal(res$status[res$what == "trial"], "TAKEN")
    expect_true(any(grepl("\\[FIX\\]  trial 'T_2026' already exists.*choose another trial_name", out)))
    expect_true(any(grepl("folder 'Old folder' exists: the trial will be saved inside it", out)))

    # A folder made by an earlier save that timed out (not listed yet when
    # checked): Breedbase answers 409 and the existing folder is used.
    fake$state$folders[[length(fake$state$folders) + 1]] <- list(trialDbId = "777", trialName = "AYT_2026",
                                                                 programDbId = "10")
    hide_once <- TRUE
    listing_hides_new <- function(req) {
      if (hide_once && grepl("/brapi/v2/trials\\?", req$url) && is.null(req$body)) {
        hide_once <<- FALSE
        return(httr2::response_json(200, body = list(metadata = list(pagination = list(totalPages = 1)),
                                                     result = list(data = list(fake$state$folders[[1]])))))
      }
      fake$mock(req)
    }
    httr2::with_mocked_responses(listing_hides_new, {
      out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", folder = "AYT_2026",
                                                  trial_name = "T2_2026", year = 2026, location = "Ibadan",
                                                  create_folder = TRUE, confirm = FALSE, description = "Test"))
    })
    expect_true(any(grepl("folder 'AYT_2026' already exists: the trial is saved inside it", out)))
    expect_equal(fake$state$studies[[2]]$trialDbId, "777")
    expect_equal(sum(vapply(fake$state$folders, `[[`, "", "trialName") == "AYT_2026"), 1)

    # A trial left without plots is a note, not a problem.
    fake$state$studies[[3]] <- list(studyDbId = "950", studyName = "EMPTY_2026")
    out <- capture.output(res <- breedbase_check(con, trial_name = "EMPTY_2026"))
    expect_true(attr(res, "ok"))
    expect_true(any(grepl("\\[NOTE\\] trial 'EMPTY_2026' exists without plots", out)))
  })
})

test_that("the folder is optional: the trial goes directly in the program", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- capture.output(res <- save_breedbase(rc, con, program = "IITA", trial_name = "NoFolder_2026",
                                                year = 2026, location = "Ibadan", confirm = FALSE,
                                                description = "Test"))
  })
  expect_true(any(grepl("no folder given: the trial will be saved directly in program 'IITA'", out)))
  expect_false(any(grepl("folder ''", out, fixed = TRUE)))
  expect_equal(fake$state$studies[[1]]$trialDbId, "10")   # the program's id
  expect_length(fake$state$folders, 1)                     # no folder created
  expect_equal(res$plots_found, 144)
  expect_error(save_breedbase(rc, con, program = "IITA", folder = "", trial_name = "X", year = 2026,
                              location = "Ibadan", description = "Test"), "or NULL")
})

test_that("phenotype_template gives the plot names that save_breedbase uses", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  sheet <- phenotype_template(rc, "AYT_2026", location = "Ibadan", traits = c("fresh root yield", "dry matter content percentage"))
  sent <- save_breedbase(rc, trial_name = "AYT_2026", year = 2026, location = "Ibadan", dry_run = TRUE,
                         description = "x")
  expect_setequal(sheet$plot_name, vapply(sent$plots[[1]], `[[`, "", "observationUnitName"))
  expect_equal(nrow(sheet), 144)
  expect_true(all(c("fresh root yield", "dry matter content percentage") %in% names(sheet)))
  expect_equal(sum(sheet$is_check), sum(rc$design$is_control))
  expect_false(anyNA(sheet$row))

  prep <- augmented_prep_design(paste0("E", 1:30), locations = c("North", "South"), seed = 1)
  ps <- phenotype_template(prep, "PREP_2026", traits = "fresh root yield")
  expect_setequal(unique(ps$trial_name), c("PREP_2026_North", "PREP_2026_South"))
})

test_that("phenotypes: dry run, wide and long data, checks and upload", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  sheet <- phenotype_template(rc, "AYT_2026", location = "Ibadan", traits = c("fresh root yield", "dry matter content percentage|CO_334:0000092"))
  set.seed(2)
  sheet$`fresh root yield` <- round(rnorm(144, 25, 4), 1)
  sheet$`dry matter content percentage|CO_334:0000092` <- round(rnorm(144, 32, 3), 1)
  sheet$`fresh root yield`[1:4] <- NA        # not measured: skipped

  obs <- suppressMessages(save_breedbase_phenotypes(sheet, dry_run = TRUE))
  expect_equal(nrow(obs), 2 * 144 - 4)
  expect_setequal(unique(obs$trait), c("fresh root yield", "dry matter content percentage|CO_334:0000092"))
  expect_error(save_breedbase_phenotypes(sheet, plot = "plot"), "Columns are")
  dup <- rbind(sheet[1:2, ], sheet[2, ])
  expect_error(save_breedbase_phenotypes(dup, dry_run = TRUE), "more than one value")

  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    # Plots not saved yet: listed together with a wrong trait.
    bad <- sheet
    names(bad)[names(bad) == "fresh root yield"] <- "root yield"
    err <- tryCatch(capture.output(save_breedbase_phenotypes(bad, con, confirm = FALSE)), error = conditionMessage)
    expect_match(err, "144 of 144 plots are not in Breedbase")
    expect_match(err, "not found in Breedbase: 'root yield'")

    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    n_before <- length(fake$state$requests)
    out <- capture.output(res <- save_breedbase_phenotypes(sheet, con, dry_run = TRUE))
    expect_true(any(grepl("Nothing was written", out)))
    expect_length(fake$state$observations, 0)

    old_options <- options(breedtools.observation_batch = 100)
    on.exit(options(old_options), add = TRUE)
    out <- capture.output(res <- save_breedbase_phenotypes(sheet, con, confirm = FALSE))
    expect_length(fake$state$observations, 2 * 144 - 4)
    expect_equal(res$values_found, res$values_sent)
    expect_equal(res$observation_variable_db_id, c("70001", "70002"))
    expect_true(any(grepl("batch 3 of 3 saved", out)))
    one <- fake$state$observations[[1]]
    expect_true(all(c("observationUnitDbId", "observationVariableDbId", "value") %in% names(one)))
    expect_type(one$value, "character")

    # Saving the same values again: nothing is sent (no duplicates).
    n_obs <- length(fake$state$observations)
    out <- capture.output(res <- save_breedbase_phenotypes(sheet, con, confirm = FALSE))
    expect_true(any(grepl("Nothing to write", out)))
    expect_length(fake$state$observations, n_obs)
    expect_equal(res$unchanged, c(140L, 144L))

    # Changed values: skipped without an error, so a loop can continue.
    changed <- sheet
    changed$`fresh root yield`[5:7] <- 99
    mixed <- changed
    mixed$`fresh root yield`[1] <- 42       # a new value alongside conflicts
    stored_before <- fake$state$observations
    n_requests <- length(fake$state$requests)
    completed <- 0L
    out <- capture.output(for (s in list(changed, mixed)) {
      result <- withVisible(save_breedbase_phenotypes(s, con, confirm = FALSE))
      completed <- completed + 1L
    })
    expect_equal(completed, 2L)
    expect_equal(sum(out == "You have decided to not replace. Phenotypes are already there."), 2L)
    expect_false(any(grepl(" -> ", out, fixed = TRUE)))
    expect_null(result$value)
    expect_false(result$visible)
    expect_identical(fake$state$observations, stored_before)
    requests <- fake$state$requests[seq.int(n_requests + 1L, length(fake$state$requests))]
    writes <- Filter(function(r) r$method %in% c("PUT", "DELETE") ||
                       (r$method == "POST" && r$path == "observations"), requests)
    expect_length(writes, 0L)

    # Changed values are still replaced when overwrite is enabled.
    out <- capture.output(res <- save_breedbase_phenotypes(changed, con, confirm = FALSE, overwrite = TRUE))
    expect_length(fake$state$observations, n_obs)            # replaced, not added
    expect_equal(res$replaced, c(3L, 0L))
    expect_equal(res$values_found, res$values_sent)
    puts <- Filter(function(r) r$method == "PUT", fake$state$requests)
    expect_length(puts, 1)
    expect_length(puts[[1]]$body, 3)

    # Empty cells with stored values: kept and reported; deleting needs overwrite.
    blank <- changed
    blank$`fresh root yield`[5] <- NA
    out <- capture.output(save_breedbase_phenotypes(blank, con, confirm = FALSE))
    expect_true(any(grepl("1 cells are empty in `data` but have a value in Breedbase", out)))
    expect_error(save_breedbase_phenotypes(blank, con, remove_empty = TRUE), "needs `overwrite = TRUE`")

    # Long data with time stamps: repeated measurements on other dates.
    long <- data.frame(plot = sheet$plot_name[1:10], trait = "fresh root yield", value = 20:29,
                       when = as.Date("2026-09-01"))
    out <- capture.output(res <- save_breedbase_phenotypes(long, con, plot = "plot", trait = "trait",
                                                           value = "value", timestamp = "when", confirm = FALSE))
    expect_equal(res$values_found, 10)
    expect_equal(fake$state$observations[[length(fake$state$observations)]]$observationTimeStamp,
                 "2026-09-01 00:00:00")
  })
})

test_that("first-stage plot values: one per field plot, design effects removed", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  d$dm <- rnorm(nrow(d), 32, 2)
  fit <- analyze_augmented_row_column_design(d, "yield")
  v <- first_stage_values(list("fresh root yield" = fit,
                               "dry matter content percentage" = analyze_augmented_row_column_design(d, "dm")),
                          rc, "AYT_2026", location = "Ibadan")
  expect_equal(nrow(v), 144)
  expect_setequal(v$plot_name, phenotype_template(rc, "AYT_2026", location = "Ibadan")$plot_name)
  # value = adjusted mean of the genotype + residual of the plot
  pr <- fit$plot_results
  am <- setNames(fit$adjusted_means$adjusted_mean, fit$adjusted_means$genotype)
  i <- match(v$plot_number, pr$plots)
  expect_equal(v$`fresh root yield`, unname(am[pr$all_entries[i]] + pr$residual[i]))
  # Row groups are removed (the trend within groups needs a spatial model).
  expect_lt(abs(coef(lm(v$`fresh root yield` ~ v$row))[2]), 0.2)
  expect_equal(attr(v, "method"), "mixed")
  if (requireNamespace("SpATS", quietly = TRUE)) {
    vs <- first_stage_values(analyze_augmented_row_column_design(d, "yield", spatial = TRUE), rc, "AYT_2026", location = "Ibadan",
                             trait = "fresh root yield")
    expect_lt(abs(coef(lm(vs$`fresh root yield` ~ vs$row))[2]), 0.15)   # most of the row trend (0.4) is removed
    expect_equal(attr(vs, "method"), "spatial")
  }
  expect_error(first_stage_values(fit, rc, "AYT_2026", location = "Ibadan"), "Give the Breedbase `trait`")

  # p-rep two-stage: stage-1 fits per location, one trial per location
  prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
  p <- prep$design
  p$y <- rnorm(nrow(p), 25, 3)
  v2 <- first_stage_values(analyze_augmented_prep_two_stage(p, "y"), prep, "PREP_2026", trait = "fresh root yield")
  expect_equal(nrow(v2), nrow(p))
  expect_false(anyNA(v2$`fresh root yield`))
  expect_setequal(unique(v2$trial_name), c("PREP_2026_North", "PREP_2026_South"))
})

test_that("resolvable row-column analyses give plot residuals", {
  r <- resolvable_row_column_design(paste0("G", 1:24), n_rep = 2, rows = 4, cols = 6, seed = 1)
  dd <- r$design
  set.seed(3)
  dd$y <- 20 + rnorm(nrow(dd)) + 0.5 * dd$row
  f <- analyze_resolvable_row_column(dd, "y", compare = FALSE)
  pr <- f$plot_results
  expect_equal(nrow(pr), nrow(dd))
  expect_equal(pr$fitted + pr$residual, pr$y)
  # X'e = 0: residuals sum to zero within each replicate and each treatment
  expect_true(all(abs(tapply(pr$residual, pr$replicate, sum)) < 1e-6))
  expect_true(all(abs(tapply(pr$residual, pr$all_entries, sum)) < 1e-6))
  v <- first_stage_values(f, r, "RRC_2026", location = "Ibadan", trait = "fresh root yield")
  expect_equal(nrow(v), nrow(dd))
  expect_false(anyNA(v$`fresh root yield`))
})

test_that("first stage saved as a trait of the field plots; second stage as an analysis", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  fit <- analyze_augmented_row_column_design(d, "yield")
  v <- first_stage_values(fit, rc, "AYT_2026", location = "Ibadan", trait = "fresh root yield")
  adjusted <- c("fresh root yield" = "fresh root yield spatially adjusted|COMP:0000123")

  # Without a connection the adjusted trait is "trait|statistic".
  obs0 <- suppressMessages(save_breedbase_first_stage(v, dry_run = TRUE))
  expect_equal(unique(obs0$trait), "fresh root yield|Adjusted Means from Spatial Correction using SpATS R")
  expect_error(save_breedbase_first_stage(v, as_trait = c("fresh root yield" = "fresh root yield|CO_334:0000013"),
                                          dry_run = TRUE), "mixing raw and adjusted")
  obs <- suppressMessages(save_breedbase_first_stage(v, as_trait = adjusted, dry_run = TRUE))
  expect_equal(nrow(obs), 144)
  expect_equal(unique(obs$trait), unname(adjusted))

  means <- second_stage_values(fit, trait = "fresh root yield")
  expect_equal(nrow(means), 112)
  sent <- suppressMessages(save_breedbase_analysis(means, analysis_name = "AYT_2026_means", program = "IITA",
                                                   description = "Row-column model", dry_run = TRUE))
  expect_equal(sent$form$upload_new_analysis_result_values_type, "analysis_result_values_match_accession_names")
  expect_equal(sent$form$upload_new_analysis_statistical_ontology_term,
               "Adjusted means from BLUEs using LMER R|SGNSTAT:0000036")
  expect_equal(names(sent$file), c("stock_names", "fresh root yield"))
  expect_equal(nrow(sent$file), 112)
  expect_error(save_breedbase_analysis(means, analysis_name = "A", program = "IITA", dry_run = TRUE),
               "Give a `description`")

  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    out <- capture.output(res <- save_breedbase_first_stage(v, con, as_trait = adjusted, confirm = FALSE))
    expect_equal(res$values_found, 144)
    expect_equal(unique(vapply(fake$state$observations, `[[`, "", "observationVariableDbId")), "70003")

    out <- capture.output(res <- save_breedbase_analysis(means, con, analysis_name = "AYT_2026_means",
                                                         program = "IITA", description = "Row-column model",
                                                         confirm = FALSE))
    expect_equal(res$rows_found, 112)
    expect_match(res$analysis_page, "/analyses/[0-9]+$")
    a <- fake$state$analyses[[1]]
    expect_equal(names(a$file), c("stock_names", "fresh root yield|CO_334:0000013"))   # full trait names
    expect_equal(a$form$upload_new_analysis_breeding_program_id, "10")
    expect_equal(a$form$sgn_session_id, "TOKEN123")
    expect_equal(a$form$upload_new_analysis_statistical_ontology_term,
                 "Adjusted means from BLUEs using LMER R|SGNSTAT:0000036")
    expect_equal(as.numeric(a$file[a$file$stock_names == "G1", 2]), means$`fresh root yield`[means$accession == "G1"],
                 tolerance = 1e-8)

    # In a folder: a missing folder stops before writing, unless created;
    # the analysis is then moved into it.
    expect_error(capture.output(save_breedbase_analysis(means, con, analysis_name = "AYT_2026_means_f",
                                                        program = "IITA", description = "x", confirm = FALSE,
                                                        folder = "Analyses 2026")),
                 "Folder 'Analyses_2026' was not found")
    expect_length(fake$state$analyses, 1)
    out <- capture.output(res_f <- save_breedbase_analysis(means, con, analysis_name = "AYT_2026_means_f",
                                                           program = "IITA", description = "x", confirm = FALSE,
                                                           folder = "Analyses 2026", create_folder = TRUE))
    fid <- vapply(fake$state$folders, `[[`, "", "trialDbId")[vapply(fake$state$folders, `[[`, "", "trialName") == "Analyses_2026"]
    expect_length(fid, 1)
    expect_equal(fake$state$moved[[res_f$analysis_db_id]], fid)
    expect_true(any(grepl("analysis put in folder 'Analyses_2026'", out)))

    # First stage as an analysis of the field plots ("Analysis of: plots").
    out <- capture.output(res1 <- save_breedbase_analysis(v, con, analysis_name = "AYT_2026_stage1", program = "IITA",
                                                          description = "Row-column model", level = "plot",
                                                          estimate = "blue", confirm = FALSE))
    a1 <- Filter(function(a) a$name == "AYT_2026_stage1", fake$state$analyses)[[1]]
    expect_equal(a1$form$upload_new_analysis_result_values_type, "analysis_result_values_match_precomputed_design")
    expect_setequal(a1$file$stock_names, v$plot_name)
    expect_equal(res1$rows_found, 144)
    bad_plots <- v
    bad_plots$plot_name[1] <- "NOT_A_PLOT"
    expect_error(capture.output(save_breedbase_analysis(bad_plots, con, analysis_name = "X1", program = "IITA",
                                                        description = "x", level = "plot", confirm = FALSE)),
                 "1 of 144 plots are not in Breedbase: NOT_A_PLOT")

    # The same name again is refused before writing; unknown traits are listed.
    expect_error(capture.output(save_breedbase_analysis(means, con, analysis_name = "AYT_2026_means",
                                                        program = "IITA", description = "x", confirm = FALSE)),
                 "already exists")
    bad <- means
    names(bad)[2] <- "root weight"
    expect_error(capture.output(save_breedbase_analysis(bad, con, analysis_name = "Other", program = "IITA",
                                                        description = "x", confirm = FALSE)),
                 "not found in Breedbase: 'root weight'")
  })
})

test_that("first stage: the adjusted trait is found or created (trait + statistic)", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  v <- first_stage_values(analyze_augmented_row_column_design(d, "yield"), rc, "AYT_2026", location = "Ibadan", trait = "fresh root yield")
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    a <- breedbase_adjusted_trait(con, "fresh root yield")
    expect_equal(a$status, "missing")
    expect_equal(a$adjusted_trait, "fresh root yield|Adjusted Means from Spatial Correction using SpATS R")

    # A dry run reports the new trait and creates nothing.
    out <- capture.output(save_breedbase_first_stage(v, con, dry_run = TRUE))
    expect_true(any(grepl("\\[NEW\\]  trait 'fresh root yield\\|Adjusted Means.*will be created", out)))
    expect_length(fake$state$composed, 0)
    expect_error(capture.output(save_breedbase_first_stage(v, con, create_trait = FALSE, confirm = FALSE)),
                 "do not exist")

    # Saving creates the composed trait (trait id + statistic id) and stores the values under it.
    out <- capture.output(res <- save_breedbase_first_stage(v, con, confirm = FALSE))
    expect_length(fake$state$composed, 1)
    expect_setequal(fake$state$composed[[1]]$components, c("70001", "80054"))
    expect_equal(res$values_found, 144)
    expect_equal(unique(vapply(fake$state$observations, `[[`, "", "observationVariableDbId")),
                 fake$state$composed[[1]]$id)

    # The next time it is found, not created again.
    a2 <- breedbase_adjusted_trait(con, "fresh root yield", create = TRUE)
    expect_equal(a2$status, "exists")
    expect_length(fake$state$composed, 1)
    # Another statistic gives another trait; an unknown statistic is a clear error.
    a3 <- breedbase_adjusted_trait(con, "fresh root yield", statistic = "Adjusted means from BLUEs using LMER R",
                                   create = TRUE)
    expect_equal(a3$status, "created")
    expect_error(breedbase_adjusted_trait(con, "fresh root yield", statistic = "Magic means"), "not a term of the Breedbase statistics ontology")
  })
})

test_that("sites with a web password: clear message, and site_username/site_password work", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(site_auth = c("tester", "webpw"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    expect_equal(breedbase_programs(con)$program, c("IITA", "NRCRI"))       # BrAPI needs no web password
    expect_error(breedbase_adjusted_trait(con, "fresh root yield"),
                 "own web password \\(\"Restricted Test Site\"\\).*site_username")
    con2 <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret",
                                             site_username = "tester", site_password = "webpw"))
    expect_true(any(grepl("site \\(web\\) password: set for tester", capture.output(print(con2)))))
    expect_false(any(grepl("webpw", capture.output(print(con2)))))
    a <- breedbase_adjusted_trait(con2, "fresh root yield", create = TRUE)
    expect_equal(a$status, "created")
    # BrAPI calls still send the login token, not the web password.
    brapi <- Filter(function(r) !grepl("^/?ajax|token", r$path), fake$state$requests)
    expect_true(all(vapply(brapi, function(r) identical(r$auth, "Bearer TOKEN123"), logical(1))))
  })
})

test_that("an existing adjusted trait is found through BrAPI, without the website", {
  skip_if_not_installed("httr2")
  fake <- fake_breedbase(site_auth = c("tester", "webpw"))
  fake$state$composed[[1]] <- list(id = "99001", name = "fresh root yield|Adjusted Means from Spatial Correction using SpATS R",
                                   components = c("70001", "80054"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    a <- breedbase_adjusted_trait(con, "fresh root yield")      # no web password given, still works
    expect_equal(a$status, "exists")
    expect_equal(a$adjusted_trait_db_id, "99001")
    paths <- vapply(fake$state$requests, `[[`, "", "path")
    expect_false(any(grepl("ajax", paths) & !grepl("check_status", paths)))   # only the login's web check
  })
})

test_that("after saving phenotypes the Breedbase phenotype cache is refreshed", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  sheet <- phenotype_template(rc, "AYT_2026", location = "Ibadan", traits = "fresh root yield")
  sheet$`fresh root yield` <- seq_len(nrow(sheet))
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    out <- capture.output(save_breedbase_phenotypes(sheet, con, confirm = FALSE))
    expect_equal(fake$state$refreshed, "phenotypes")
    expect_true(any(grepl("refreshing its phenotype cache", out)))
  })
  # Behind a web password the refresh cannot run: a note, not an error.
  fake2 <- fake_breedbase(site_auth = c("tester", "webpw"))
  httr2::with_mocked_responses(fake2$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    out <- capture.output(res <- save_breedbase_phenotypes(sheet, con, confirm = FALSE))
    expect_equal(res$values_found, 144)
    expect_true(any(grepl("\\[NOTE\\] The values are saved", out)))
    expect_true(any(grepl("Search Wizard", out)))
  })
})

test_that("statistics come from the server: other names and ids work, nothing is fixed", {
  skip_if_not_installed("httr2")
  # Another Breedbase: different ids, slightly different names, no LMER terms.
  other <- list(
    list(cvterm_id = 5, cv_name = "SGNStatistics_ontology", cvterm_name = "Spatially corrected values using SpATS",
         accession = "SGNSTAT:0000901"),
    list(cvterm_id = 6, cv_name = "SGNStatistics_ontology", cvterm_name = "Adjusted means from BLUEs",
         accession = "SGNSTAT:0000902"))
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  fit <- analyze_augmented_row_column_design(d, "yield")
  means <- second_stage_values(fit, trait = "fresh root yield")
  fake <- fake_breedbase(statistics = other)
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    s <- breedbase_statistics(con)
    expect_equal(s$accession, c("SGNSTAT:0000902", "SGNSTAT:0000901"))
    a <- breedbase_adjusted_trait(con, "fresh root yield", create = TRUE)
    expect_equal(a$adjusted_trait, "fresh root yield|Spatially corrected values using SpATS")
    expect_setequal(fake$state$composed[[1]]$components, c("70001", "5"))
    out <- capture.output(save_breedbase_analysis(means, con, analysis_name = "A1", program = "IITA",
                                                  description = "x", confirm = FALSE))
    expect_equal(fake$state$analyses[[1]]$form$upload_new_analysis_statistical_ontology_term,
                 "Adjusted means from BLUEs|SGNSTAT:0000902")
    # BLUPs asked, but this Breedbase has none: stop and list its statistics.
    expect_error(capture.output(save_breedbase_analysis(means, con, analysis_name = "A2", program = "IITA",
                                                        description = "x", estimate = "blup", confirm = FALSE)),
                 "No BLUP statistic.*Adjusted means from BLUEs; Spatially corrected values using SpATS")
    # A statistic given by name is checked against the server.
    expect_error(capture.output(save_breedbase_analysis(means, con, analysis_name = "A3", program = "IITA",
                                                        description = "x", statistic = "Magic means",
                                                        confirm = FALSE)), "not a term of the Breedbase statistics ontology")
  })
  # On the default fake server, the BLUP term is chosen for random genotypes.
  fake2 <- fake_breedbase()
  httr2::with_mocked_responses(fake2$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase_analysis(means, con, analysis_name = "A4", program = "IITA", description = "x",
                                           estimate = "blup", confirm = FALSE))
    expect_equal(fake2$state$analyses[[1]]$form$upload_new_analysis_statistical_ontology_term,
                 "Adjusted means from BLUPs using LMER R|SGNSTAT:0000034")
  })
})

test_that("a web password equal to the Breedbase login is found and used at login", {
  skip_if_not_installed("httr2")
  same <- fake_breedbase(site_auth = c("jdoe", "secret"))
  httr2::with_mocked_responses(same$mock, {
    expect_message(con <- breedbase_login("https://bb.example.org", "jdoe", password = "secret"),
                   "your Breedbase login works there")
    expect_equal(con$site_auth$username, "jdoe")
    expect_equal(breedbase_adjusted_trait(con, "fresh root yield", create = TRUE)$status, "created")
  })
  other <- fake_breedbase(site_auth = c("tester", "webpw"))
  httr2::with_mocked_responses(other$mock, {
    expect_message(con <- breedbase_login("https://bb.example.org", "jdoe", password = "secret"),
                   "own web password")
    expect_null(con$site_auth)
  })
  plain <- fake_breedbase()
  httr2::with_mocked_responses(plain$mock, {
    msgs <- character(0)
    con <- withCallingHandlers(breedbase_login("https://bb.example.org", "jdoe", password = "secret"),
                               message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
    expect_null(con$site_auth)
    expect_false(any(grepl("web", msgs)))
  })
})

test_that("web pages closed by a web password: the analysis upload file is written, checked via BrAPI", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  fit <- analyze_augmented_row_column_design(d, "yield")
  stage1 <- first_stage_values(fit, rc, "AYT_2026", location = "Ibadan", trait = "fresh root yield")
  means <- second_stage_values(fit, trait = "fresh root yield")
  fake <- fake_breedbase(site_auth = c("tester", "webpw"))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    f1 <- tempfile(fileext = ".csv")
    out <- capture.output(res <- save_breedbase_analysis(stage1, con, analysis_name = "AYT_2026_stage1",
                                                         program = "IITA", description = "SpATS", level = "plot",
                                                         upload_file = f1, confirm = FALSE))
    expect_true(file.exists(f1))
    up <- utils::read.csv(f1, check.names = FALSE)
    expect_equal(names(up), c("stock_names", "fresh root yield|CO_334:0000013"))   # full name from BrAPI
    expect_setequal(up$stock_names, stage1$plot_name)
    expect_true(any(grepl("Upload analysis", out)))
    expect_true(any(grepl("Analysis Of: +plots", out)))
    expect_equal(res$form[["Breeding Program"]], "IITA")
    expect_length(fake$state$analyses, 0)                    # nothing sent

    # Problems are still found through BrAPI before any file is written.
    bad <- means
    bad$accession[1] <- "NOT_AN_ACCESSION"
    expect_error(capture.output(save_breedbase_analysis(bad, con, analysis_name = "X", program = "IITA",
                                                        description = "x", upload_file = tempfile(), confirm = FALSE)),
                 "not in Breedbase: NOT_AN_ACCESSION")
  })
})

test_that("statistics are checked against the Breedbase statistics ontology (bundled obo)", {
  obo <- breedbase_statistics()
  expect_gt(nrow(obo), 40)
  expect_equal(obo$accession[obo$statistic == "Adjusted Means from Spatial Correction using SpATS R"], "SGNSTAT:0000054")
  expect_false("SGNSTAT:0000000" %in% obo$accession)                   # the root is not a statistic
  ok <- breedtools:::bb_check_statistic_ontology("adjusted means from BLUEs using LMER R")
  expect_equal(ok, list(name = "Adjusted means from BLUEs using LMER R", accession = "SGNSTAT:0000036"))
  expect_equal(breedtools:::bb_check_statistic_ontology("SGNSTAT:0000054")$name,
               "Adjusted Means from Spatial Correction using SpATS R")
  expect_error(breedtools:::bb_check_statistic_ontology("Adjusted means from BLUEs using LMER R|SGNSTAT:0000054"),
               "is \"Adjusted Means from Spatial Correction using SpATS R\" in the Breedbase statistics ontology")
  expect_error(breedtools:::bb_check_statistic_ontology("Adjusted means from BLUE using LMER"), "Did you mean")
  expect_error(breedtools:::bb_check_statistic_ontology("fresh root yield|CO_334:0000013"), "is not a statistic")
  # Offline choices come from the ontology.
  expect_equal(breedtools:::bb_choose_statistic(NULL, NULL, "blup", offline = TRUE)$accession, "SGNSTAT:0000034")
  # A dry run without a connection already carries the exact term.
  means <- data.frame(accession = c("G1", "G2"), `fresh root yield` = c(1, 2), check.names = FALSE)
  sent <- suppressMessages(save_breedbase_analysis(means, analysis_name = "A", program = "P", description = "d",
                                                   dry_run = TRUE))
  expect_equal(sent$form$upload_new_analysis_statistical_ontology_term,
               "Adjusted means from BLUEs using LMER R|SGNSTAT:0000036")
})

test_that("breedbase_trial_plots reads a saved trial; first_stage_values checks the design against it", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  other <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 2)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + rnorm(nrow(d))
  fit <- analyze_augmented_row_column_design(d, "yield")
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    fake$state$studies[[1]]$studyDbId -> sid
    plots <- breedbase_trial_plots(con, "AYT_2026")
    expect_equal(nrow(plots), 144)
    expect_equal(plots$accession, rc$design$all_entries[match(plots$plot_number, rc$design$plots)])
    v <- first_stage_values(fit, rc, "AYT_2026", location = "Ibadan", trait = "y", breedbase_plots = plots)
    expect_setequal(v$plot_name, plots$plot_name)
    expect_error(first_stage_values(fit, other, "AYT_2026", location = "Ibadan", trait = "y", breedbase_plots = plots),
                 "differ from the trial saved in Breedbase")
    expect_error(breedbase_trial_plots(con, "NOPE"), "was not found")
  })
})

test_that("values are added to an existing analysis (a second trait), without duplicating", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  d <- rc$design
  set.seed(1)
  d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d))
  d$dm <- 32 + rnorm(nrow(d))
  stage1 <- first_stage_values(analyze_augmented_row_column_design(d, "yield"), rc, "AYT_2026", location = "Ibadan", trait = "fresh root yield")
  dm1 <- first_stage_values(analyze_augmented_row_column_design(d, "dm"), rc, "AYT_2026", location = "Ibadan",
                            trait = "dry matter content percentage|CO_334:0000092")
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "AYT_2026", year = 2026,
                                  location = "Ibadan", confirm = FALSE, description = "x"))
    # Before saving: the names the analysis plots will have.
    preview <- analysis_plot_names(stage1, "S1", level = "plot")
    expect_equal(preview$analysis_plot_name, paste0("S1_", stage1$plot_name))
    expect_equal(preview$accession, stage1$accession)
    expect_equal(analysis_plot_names(stage1, "S 1")$analysis_plot_name[1], paste0("S_1_", stage1$plot_name[1]))
    # Before saving: the composed trait each trait will be stored under.
    out <- capture.output(dr <- save_breedbase_analysis(stage1, con, analysis_name = "S1", program = "IITA",
                                                        description = "x", level = "plot", dry_run = TRUE))
    expect_equal(dr$traits$status, "new")
    expect_match(dr$traits$stored_as, "^fresh root yield\\|[^|]+$")
    expect_true(any(grepl("'fresh root yield' will be stored as 'fresh root yield|", out, fixed = TRUE)))
    expect_true(any(grepl("breedbase_statistics(con)", out, fixed = TRUE)))
    # Another statistic, chosen from the list: another stored name.
    blup <- "Adjusted means from BLUPs using LMER R"
    expect_true(blup %in% breedbase_statistics(con)$statistic)
    out <- capture.output(dr2 <- save_breedbase_analysis(stage1, con, analysis_name = "S1", program = "IITA",
                                                         description = "x", level = "plot", dry_run = TRUE,
                                                         statistic = blup))
    expect_equal(dr2$traits$stored_as, paste0("fresh root yield|", blup))
    expect_error(capture.output(save_breedbase_analysis(stage1, con, analysis_name = "S1", program = "IITA",
                                                        description = "x", level = "plot", dry_run = TRUE,
                                                        statistic = "My own statistic")))
    capture.output(save_breedbase_analysis(stage1, con, analysis_name = "S1", program = "IITA", description = "x",
                                           level = "plot", confirm = FALSE))
    # The composed trait now exists: a dry run shows it with its COMP id.
    out <- capture.output(dr3 <- save_breedbase_analysis(stage1, con, analysis_name = "S1_again", program = "IITA",
                                                         description = "x", level = "plot", dry_run = TRUE))
    expect_equal(dr3$traits$status, "exists")
    expect_match(dr3$traits$stored_as, paste0("^", dr$traits$stored_as, "\\|COMP:[0-9]+$"))
    # Read back: the same plots and accessions, one trait so far.
    saved <- breedbase_analysis_plots(con, "S1")
    expect_equal(nrow(saved), 144)
    expect_setequal(saved$analysis_plot_name, preview$analysis_plot_name)
    expect_equal(saved$accession, preview$accession[match(saved$analysis_plot_name, preview$analysis_plot_name)])
    expect_equal(saved$plot_number, 1:144)
    # Breedbase stores the analysis values composed with its statistic.
    ytrait <- grep("^fresh root yield\\|[^|]+\\|COMP:[0-9]+$", names(saved), value = TRUE)
    expect_length(ytrait, 1)
    stat <- sub("^fresh root yield\\|([^|]+)\\|COMP:[0-9]+$", "\\1", ytrait)
    expect_false("fresh root yield|CO_334:0000013" %in% names(saved))
    expect_equal(saved[[ytrait]], stage1[["fresh root yield"]][match(saved$plot_name, stage1$plot_name)],
                 tolerance = 1e-6)
    expect_error(breedbase_analysis_plots(con, "NOPE"), "was not found")
    out <- capture.output(res <- add_breedbase_analysis_values(dm1, con, analysis_name = "S1", level = "plot",
                                                               confirm = FALSE))
    expect_equal(res$plots, 144)
    # The second trait too: same statistic as the analysis, never the plain trait.
    expect_true(any(grepl(paste0("statistic '", stat, "' (the one of this analysis)"), out, fixed = TRUE)))
    expect_false(any(grepl("CO_334:0000092", names(fake$state$analysis_values))))
    both <- breedbase_analysis_plots(con, "S1")
    dtrait <- grep(paste0("^dry matter content percentage\\|", stat, "\\|COMP:[0-9]+$"), names(both), value = TRUE)
    expect_length(dtrait, 1)
    expect_false("dry matter content percentage|CO_334:0000092" %in% names(both))
    expect_false(anyNA(both[[dtrait]]))
    expect_length(fake$state$analysis_values, 144)
    expect_true(all(startsWith(names(fake$state$analysis_values), "S1_AYT_2026_")))
    # The same values again: Breedbase reports them; nothing is replaced without overwrite.
    expect_error(capture.output(add_breedbase_analysis_values(dm1, con, analysis_name = "S1", level = "plot",
                                                              confirm = FALSE)), "already stored")
    expect_error(capture.output(add_breedbase_analysis_values(dm1, con, analysis_name = "NOPE", confirm = FALSE)),
                 "was not found")
    out <- capture.output(dr <- add_breedbase_analysis_values(dm1, con, analysis_name = "S1", overwrite = TRUE,
                                                              dry_run = TRUE))
    expect_true(any(grepl("Nothing was stored", out)))
  })
})

test_that("from Breedbase to analyses: read trials, outliers, first, second and single stage, store", {
  skip_if_not_installed("httr2")
  skip_if_not_installed("lme4")
  sites <- c("Ibadan", "Ikenne")
  designs <- lapply(1:2, function(i) augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                                                 n_candidates = 5, seed = i))
  set.seed(1)
  g <- stats::setNames(rnorm(112, sd = 3), c(paste0("G", 1:108), paste0("C", 1:4)))
  fake <- fake_breedbase(accessions = c(paste0("G", 1:108), paste0("C", 1:4)))
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    for (i in 1:2) {
      capture.output(save_breedbase(designs[[i]], con, program = "IITA", trial_name = paste0("T_", sites[i]),
                                    year = 2026, location = sites[i], confirm = FALSE, description = "x"))
      sheet <- phenotype_template(designs[[i]], paste0("T_", sites[i]), traits = "fresh root yield",
                                  location = sites[i])
      sheet$`fresh root yield` <- round(50 + g[sheet$accession] + 0.5 * sheet$row + rnorm(nrow(sheet)), 2)
      sheet$`fresh root yield`[5] <- 500                         # an outlier
      capture.output(save_breedbase_phenotypes(sheet, con, confirm = FALSE))
    }
    out <- capture.output(pheno <- breedbase_phenotypes(con, trials = paste0("T_", sites), traits = "fresh root yield"))
    expect_equal(nrow(pheno), 288)
    expect_equal(sum(!is.na(pheno$`fresh root yield`)), 288)
    expect_equal(attr(pheno, "trait_names")[["fresh root yield"]], "fresh root yield|CO_334:0000013")
    expect_true(all(c("plot_name", "accession", "rep", "block", "row", "col") %in% names(pheno)))
    expect_error(breedbase_phenotypes(con, "NOPE", "fresh root yield"), "Trial\\(s\\) not found")

    clean <- remove_outliers(pheno, "fresh root yield", by = "trial_name", quiet = TRUE)
    expect_true(500 %in% attr(clean, "outliers")$value)
    two <- analyze_two_stage(clean, "fresh root yield", stage1 = "blocks", genotype_col = "accession",
                             location_col = "trial_name")
    one <- analyze_single_stage(clean, "fresh root yield")
    n_geno <- length(unique(clean$accession[!is.na(clean$`fresh root yield`)]))   # genotypes with data left
    expect_equal(nrow(two$adjusted_means), n_geno)
    expect_equal(nrow(one$adjusted_means), n_geno)
    expect_gt(cor(two$adjusted_means$adjusted_mean, g[two$adjusted_means$genotype]), 0.8)
    m <- merge(two$adjusted_means, one$adjusted_means, by = "genotype")
    expect_gt(cor(m$adjusted_mean.x, m$adjusted_mean.y), 0.95)      # two-stage and single-stage agree
    if (requireNamespace("SpATS", quietly = TRUE)) {
      sp <- analyze_single_stage(clean, "fresh root yield", spatial = TRUE)
      expect_equal(sp$diagnostics$model, "spatial (SpATS)")
      expect_equal(nrow(sp$adjusted_means), n_geno)
      expect_false(any(grepl("f\\(", sp$variance_components$component)))
      expect_lt(abs(coef(lm(residual ~ row, sp$plot_results))[2]), 0.05)   # the row trend is in the surface
    }

    s1 <- first_stage_values(two, trait = "fresh root yield")
    expect_setequal(s1$plot_name, clean$plot_name[!is.na(clean$`fresh root yield`)])
    s2 <- second_stage_values(two, trait = "fresh root yield")
    ss <- second_stage_values(one, trait = "fresh root yield")

    capture.output(save_breedbase_analysis(s1, con, analysis_name = "MET_stage1", program = "IITA",
                                           description = "first stage", level = "plot", estimate = "blue",
                                           confirm = FALSE))
    capture.output(save_breedbase_analysis(s2, con, analysis_name = "MET_means", program = "IITA",
                                           description = "two-stage means", confirm = FALSE))
    # A second analysis in the same project: single-stage means under "trait|statistic".
    out <- capture.output(res <- add_breedbase_analysis_values(ss, con, analysis_name = "MET_means", level = "accession",
                                                               statistic = "Phenotypic BLUEs using LMER R",
                                                               confirm = FALSE))
    expect_equal(res$plots, n_geno)
    expect_true(all(grepl("^MET_means_.* fresh root yield\\|Phenotypic BLUEs using LMER R\\|COMP:", names(fake$state$analysis_values))))
  })
})

test_that("spaces in saved names become '_', accession and location names are kept", {
  prep <- augmented_prep_design(c(paste0("G", 1:40), "TME 419", "TMS 30572"), c("Alcobaca (BA)", "Itamaraju (BA)"),
                                block_size = 6, field_cols = 7, seed = 1)
  expect_message(
    payload <- save_breedbase(prep, trial_name = "EMBRAPA 2026", year = 2026,
                              location = c("Alcobaca (BA)" = "Alcobaca (BA)", "Itamaraju (BA)" = "Itamaraju (BA)"),
                              dry_run = TRUE, description = "Test trial"),
    "special characters replaced"
  )
  expect_equal(vapply(payload$studies, `[[`, "", "studyName"),
               c("EMBRAPA_2026_Alcobaca_BA", "EMBRAPA_2026_Itamaraju_BA"))
  units <- unlist(lapply(payload$plots, function(p) vapply(p, `[[`, "", "observationUnitName")))
  germ <- unlist(lapply(payload$plots, function(p) vapply(p, `[[`, "", "germplasmName")))
  expect_false(any(grepl(" ", units)))
  expect_true("TME 419" %in% germ)
  expect_true(any(grepl("TME_419", units)))
  expect_equal(bb_clean_name(c(" a  b ", "c"), quiet = TRUE), c("a_b", "c"))
  expect_equal(bb_clean_name(c("AYT 50% $x (BA) {1} [2] a\\b/c <d> e*f", "Alcoba\u00e7a-1.2"), quiet = TRUE),
               c("AYT_50_x_BA_1_2_a_b_c_d_e_f", "Alcoba\u00e7a-1.2"))
  expect_error(bb_clean_name("()", "trial name"), "no letters or digits")
})

test_that("plot names are made with the design, can be changed and are used when saving", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                    n_candidates = 5, seed = 1, trial_name = "AYT_2026", location = "Ibadan")
  first <- rc$design[rc$design$plots == 1, ]
  expect_equal(first$plot_name, paste("AYT_2026", "Ibadan", first$all_entries, first$block, first$rep, sep = "_"))
  expect_false(anyDuplicated(rc$design$plot_name) > 0)
  rc0 <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  expect_false("plot_name" %in% names(rc0$design))
  save <- function(d) suppressMessages(save_breedbase(d, trial_name = "AYT_2026", year = 2026, location = "Ibadan",
                                                      description = "x", dry_run = TRUE))
  expect_identical(save(rc), save(rc0))
  expect_identical(phenotype_template(rc)$plot_name, phenotype_template(rc0, "AYT_2026", location = "Ibadan")$plot_name)
  rc$design$plot_name[rc$design$plots == 1] <- "MY_OWN_NAME"
  expect_true("MY_OWN_NAME" %in% vapply(save(rc)$plots[[1]], `[[`, "", "observationUnitName"))
  expect_error(add_plot_names(rc0, "AYT_2026"), "uses \\{location\\}")
  rc2 <- add_plot_names(rc0, "AYT_2026", plot_name_pattern = "{trial}_{plot}")
  expect_equal(rc2$design$plot_name[rc2$design$plots == 7], "AYT_2026_7")
  prep <- augmented_prep_design(paste0("G", 1:60), c("North", "South"), block_size = 8, n_candidates = 5,
                                seed = 1, trial_name = "PREP")
  expect_true(all(startsWith(prep$design$plot_name[prep$design$location == "South"], "PREP_South_South_")))
  expect_identical(suppressMessages(save_breedbase(prep, trial_name = "PREP", year = 2026, description = "x",
                                                   dry_run = TRUE))$plots,
                   suppressMessages(save_breedbase(augmented_prep_design(paste0("G", 1:60), c("North", "South"),
                     block_size = 8, n_candidates = 5, seed = 1), trial_name = "PREP", year = 2026,
                     description = "x", dry_run = TRUE))$plots)
})

test_that("breedbase_check finds plot names from designs", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  named <- add_plot_names(rc, "A", location = "Ibadan")
  expect_length(bb_check_plot_input(named)$names, 144)
  expect_equal(bb_check_plot_input(rc)$unnamed, 1L)
  both <- bb_check_plot_input(list(rc, rc), trial_name = c("A", "B"), location = c("Ibadan", "Ubiaja"))
  expect_length(both$names, 288)
  expect_true(all(startsWith(both$names[145:288], "B_Ubiaja_")))
  expect_equal(bb_check_plot_input(rc, trial_name = c("A", "B"))$unnamed, 1L)
  expect_setequal(bb_accession_names(list(rc, rc)), c(paste0("G", 1:108), paste0("C", 1:4)))
})

test_that("plot coordinates: from the field corner, points and tables", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  g <- add_plot_coordinates(rc, origin = c(7.5, 3.9), plot_width = 3, plot_length = 5, row_gap = 1)
  d <- g$design
  p11 <- d[d$row == 1 & d$col == 1, ]
  expect_equal(p11$WGS84_bottom_left_y, 7.5)
  expect_equal(p11$WGS84_bottom_left_x, 3.9)
  # Distances in metres, with a local degree length.
  m_lat <- 110574; m_lon <- 110339
  p12 <- d[d$row == 1 & d$col == 2, ]
  p21 <- d[d$row == 2 & d$col == 1, ]
  expect_equal((p12$WGS84_bottom_left_x - 3.9) * m_lon, 3, tolerance = 0.01)      # columns go east
  expect_equal((p21$WGS84_bottom_left_y - 7.5) * m_lat, 6, tolerance = 0.01)      # rows go north (5 m + 1 m alley)
  expect_equal(p11$WGS84_top_right_y - p11$WGS84_bottom_left_y, 5 / m_lat, tolerance = 0.01)
  # Rotated: columns north-east.
  r <- add_plot_coordinates(rc, origin = c(7.5, 3.9), plot_width = 10, plot_length = 5,
                            col_bearing = 45, row_bearing = 315)$design
  q <- r[r$row == 1 & r$col == 2, ]
  expect_equal((q$WGS84_bottom_left_x - 3.9) * m_lon, 10 * sin(pi / 4), tolerance = 0.01)
  expect_equal((q$WGS84_bottom_left_y - 7.5) * m_lat, 10 * cos(pi / 4), tolerance = 0.01)
  expect_error(add_plot_coordinates(rc, origin = c(7.5, 3.9), plot_width = 3, plot_length = 5, row_bearing = 45),
               "right angles")
  expect_error(add_plot_coordinates(rc, origin = c(200, 3.9), plot_width = 3, plot_length = 5), "swapped")
  expect_error(add_plot_coordinates(rc), "Give one of")

  # GeoJSON in the payload, as Breedbase stores it: [lon, lat], closed ring.
  pay <- save_breedbase(g, trial_name = "T", year = 2026, location = "Ibadan", description = "x", dry_run = TRUE)
  geo <- pay$plots[[1]][[1]]$observationUnitPosition$geoCoordinates
  expect_equal(geo$geometry$type, "Polygon")
  ring <- geo$geometry$coordinates[[1]]
  expect_length(ring, 5)
  expect_equal(ring[[1]], ring[[5]])
  first <- g$design[g$design$plots == 1, ]
  expect_equal(ring[[1]], c(first$WGS84_bottom_left_x, first$WGS84_bottom_left_y))

  # Points from a table matched by plot number; Breedbase template names too.
  pts <- data.frame(plot_number = rc$design$plots, WGS84_y = 7.5 + rc$design$row / 1e4,
                    WGS84_x = 3.9 + rc$design$col / 1e4)
  pt <- add_plot_coordinates(rc, coordinates = pts[sample(nrow(pts)), ])
  expect_equal(pt$design$latitude, 7.5 + rc$design$row / 1e4)
  expect_false(any(bb_corner_cols %in% names(pt$design)))
  pay <- save_breedbase(pt, trial_name = "T", year = 2026, location = "Ibadan", description = "x", dry_run = TRUE)
  geo <- pay$plots[[1]][[1]]$observationUnitPosition$geoCoordinates
  expect_equal(geo$geometry$type, "Point")
  expect_equal(geo$geometry$coordinates, c(pt$design$longitude[pt$design$plots == 1], pt$design$latitude[pt$design$plots == 1]))
  expect_error(add_plot_coordinates(rc, coordinates = rbind(pts, transform(pts[1, ], plot_number = 999L))),
               "match no plot")
  # Vectors; a plot without coordinates is saved without them.
  lat <- 7.5 + rc$design$row / 1e4; lat[1] <- NA
  expect_message(v <- add_plot_coordinates(rc, latitude = lat, longitude = 3.9 + rc$design$col / 1e4), "1 plot has no")
  pay <- save_breedbase(v, trial_name = "T", year = 2026, location = "Ibadan", description = "x", dry_run = TRUE)
  has <- vapply(pay$plots[[1]], function(p) !is.null(p$observationUnitPosition$geoCoordinates), logical(1))
  expect_equal(sum(has), 143)

  # p-rep: one origin per location.
  bare <- augmented_prep_design(paste0("G", 1:60), c("North", "South"), block_size = 8, n_candidates = 5, seed = 1)
  expect_error(add_plot_coordinates(bare, origin = list(North = c(7.5, 3.9), South = c(6.1, 5.2)),
                                    plot_width = 3, plot_length = 5), "field layout")
  prep <- augmented_prep_design(paste0("G", 1:60), c("North", "South"), block_size = 8, n_candidates = 5, seed = 1,
                                rows_per_block = 2, blocks_per_row = 2)
  expect_false(anyNA(prep$design$row))
  expect_error(add_plot_coordinates(prep, origin = c(7.5, 3.9), plot_width = 3, plot_length = 5), "one c\\(latitude")
  pp <- add_plot_coordinates(prep, origin = list(North = c(7.5, 3.9), South = c(6.1, 5.2)),
                             plot_width = 3, plot_length = 5)
  expect_true(all(pp$design$latitude[pp$design$location == "South"] < 7))
})

test_that("plot coordinates are saved after the plots and read back", {
  skip_if_not_installed("httr2")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 1)
  g <- add_plot_coordinates(rc, origin = c(7.5, 3.9), plot_width = 3, plot_length = 5)
  fake <- fake_breedbase()
  httr2::with_mocked_responses(fake$mock, {
    con <- suppressMessages(breedbase_login("https://bb.example.org", "jdoe", password = "secret"))
    out <- NULL
    printed <- capture.output(out <- save_breedbase(g, con, program = "IITA", trial_name = "GPS_2026", year = 2026,
                                                    location = "Ibadan", description = "x", confirm = FALSE))
    expect_equal(out$plots_with_gps, 144)
    expect_equal(out$gps_found, 144)
    expect_true(any(grepl("GPS coordinates of trial 'GPS_2026' saved for 144 of 144", printed)))
    puts <- Filter(function(r) r$method == "PUT", fake$state$requests)
    expect_true(length(puts) >= 1)
    expect_null(puts[[1]]$body[[1]]$observationUnitPosition$observationLevel$levelCode)   # plot numbers untouched
    stored <- Filter(function(u) u$observationUnitPosition$observationLevel$levelCode == "1", fake$state$plots)[[1]]
    expect_equal(stored$observationUnitPosition$geoCoordinates$geometry$type, "Polygon")

    # A trial saved without GPS gets it later; plots matched by number,
    # accessions checked.
    capture.output(save_breedbase(rc, con, program = "IITA", trial_name = "LATER_2026", year = 2026,
                                  location = "Ibadan", description = "x", confirm = FALSE))
    capture.output(res <- save_breedbase_coordinates(g, con, trial_name = "LATER_2026", confirm = FALSE))
    expect_equal(res$gps_found, 144)
    other <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 5, seed = 2)
    other <- add_plot_coordinates(other, origin = c(7.5, 3.9), plot_width = 3, plot_length = 5)
    expect_error(capture.output(save_breedbase_coordinates(other, con, trial_name = "LATER_2026", confirm = FALSE)),
                 "another accession")
    expect_error(save_breedbase_coordinates(rc, con, trial_name = "LATER_2026"), "no plot coordinates")
  })
})
