# GPS coordinates of the plots: added to the design, saved to Breedbase with
# the plots.
#
# Breedbase keeps the GPS of a plot as a GeoJSON Feature (stock property
# "plot_geo_json"), a Polygon with the four corners of the plot or a Point,
# coordinates as [longitude, latitude] (WGS84), as its "Upload plot GPS"
# form does. BrAPI POST /observationunits does not store them, so they are
# sent after the plots are created, with PUT /observationunits.

bb_corner_cols <- c("WGS84_bottom_left_x", "WGS84_bottom_left_y", "WGS84_bottom_right_x", "WGS84_bottom_right_y",
                    "WGS84_top_right_x", "WGS84_top_right_y", "WGS84_top_left_x", "WGS84_top_left_y")

#' Add GPS coordinates to the plots of a design
#'
#' Adds the GPS position of every plot to `design$design`, so that
#' [save_breedbase()] saves it with the plots and Breedbase shows the trial
#' on its field map. Coordinates are WGS84 decimal degrees (as from a
#' phone or GPS). There are three ways to give them:
#'
#' * **From the field corner** (`origin`): give the GPS of the outer corner
#'   of the plot in row 1, column 1, the plot size and the directions in
#'   which rows and columns run. The four corners of every plot are
#'   calculated from its row and column.
#' * **One point per plot** (`latitude`, `longitude`): vectors with one value
#'   per row of `design$design`, or the names of columns of the design that
#'   hold them.
#' * **A table** (`coordinates`), e.g. read with [read_breedbase_file()]: one
#'   row per plot, matched by `plot_name` (when both have it) or by plot
#'   number (`plots`, `plot_number` or `plot`; for p-rep designs also
#'   `location`). It holds either `latitude` and `longitude` (or `WGS84_x`,
#'   `WGS84_y`) for a point, or the eight corner columns of the Breedbase
#'   GPS upload template (`WGS84_bottom_left_x`, `WGS84_bottom_left_y`, ...,
#'   `WGS84_top_left_y`, where `x` is the longitude and `y` the latitude)
#'   for a polygon.
#'
#' @param design A design from [augmented_row_column_design()],
#'   [augmented_prep_design()], [resolvable_row_column_design()],
#'   [as_breedbase_design()] or [prepare_breedbase_trial()].
#' @param coordinates Optional table of coordinates, as described above.
#' @param latitude,longitude Optional point per plot, as described above.
#' @param origin Optional `c(latitude, longitude)` of the outer corner of the
#'   plot in row 1, column 1. For a p-rep design with several locations, a
#'   list with one such pair per location, named by location.
#' @param plot_width,plot_length Plot size in metres: `plot_width` across
#'   columns (from column 1 towards column 2), `plot_length` across rows
#'   (from row 1 towards row 2). Used with `origin`.
#' @param col_gap,row_gap Space in metres between plots of neighbouring
#'   columns and rows (alleys). Used with `origin`.
#' @param col_bearing,row_bearing Compass direction (degrees clockwise from
#'   north) in which column numbers and row numbers increase; they must be
#'   at right angles. The default puts column 1 at the west, columns running
#'   east, and row 1 at the south, rows running north. Used with `origin`.
#'
#' @details
#' With `origin`, positions are calculated on a local flat approximation of
#' the Earth, accurate to centimetres over a field. The design needs `row`
#' and `col` columns (designs converted from agricolae need them in the
#' field book).
#'
#' The design gets `latitude` and `longitude` columns (the centre of the
#' plot), plus the eight corner columns for polygons. They can be edited.
#' [save_breedbase()] sends a polygon when all eight corner columns of a plot
#' are filled, and otherwise a point from `latitude` and `longitude`. Plots
#' without coordinates are saved without them. For a trial that is already
#' in Breedbase, use [save_breedbase_coordinates()].
#'
#' @return The design with the coordinate columns added.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' # Plots of 3 m x 5 m, 1 m alleys between rows, columns running east:
#' rc <- add_plot_coordinates(rc, origin = c(7.4931, 3.9035),
#'                            plot_width = 3, plot_length = 5, row_gap = 1)
#' head(rc$design[, c("plots", "row", "col", "latitude", "longitude")])
#'
#' # A field laid out at an angle: columns running north-east.
#' rc <- add_plot_coordinates(rc, origin = c(7.4931, 3.9035), plot_width = 3,
#'                            plot_length = 5, col_bearing = 45, row_bearing = 315)
#'
#' # One GPS point per plot, measured in the field:
#' pts <- data.frame(plots = rc$design$plots,
#'                   latitude = 7.4931 + rc$design$row * 5e-5,
#'                   longitude = 3.9035 + rc$design$col * 3e-5)
#' rc <- add_plot_coordinates(rc, coordinates = pts)
#'
#' @export
add_plot_coordinates <- function(design, coordinates = NULL, latitude = NULL, longitude = NULL,
                                 origin = NULL, plot_width = NULL, plot_length = NULL,
                                 col_gap = 0, row_gap = 0, col_bearing = 90, row_bearing = 0) {
  if (!bb_is_design(design)) {
    stop("`design` must come from augmented_row_column_design(), augmented_prep_design(), ",
         "resolvable_row_column_design() or as_breedbase_design().", call. = FALSE)
  }
  modes <- c(coordinates = !is.null(coordinates), points = !is.null(latitude) || !is.null(longitude),
             origin = !is.null(origin))
  if (sum(modes) != 1) {
    stop("Give one of: `origin` (with the plot size), `latitude` and `longitude`, or a `coordinates` table.",
         call. = FALSE)
  }
  d <- as.data.frame(design$design, stringsAsFactors = FALSE)
  n <- nrow(d)
  new <- data.frame(latitude = rep(NA_real_, n), longitude = rep(NA_real_, n))
  for (cc in bb_corner_cols) new[[cc]] <- NA_real_
  polygons <- FALSE

  if (modes[["points"]]) {
    if (is.null(latitude) || is.null(longitude)) stop("Give both `latitude` and `longitude`.", call. = FALSE)
    get <- function(v, what) {
      if (is.character(v) && length(v) == 1 && v %in% names(d)) v <- d[[v]]
      v <- suppressWarnings(as.numeric(v))
      if (length(v) != n) {
        stop("`", what, "` must have one value per plot (", n, ") or name a column of the design.", call. = FALSE)
      }
      v
    }
    new$latitude <- get(latitude, "latitude")
    new$longitude <- get(longitude, "longitude")
  } else if (modes[["coordinates"]]) {
    tab <- as.data.frame(coordinates, stringsAsFactors = FALSE)
    i <- bb_match_coordinate_rows(design, d, tab)
    num <- function(col) suppressWarnings(as.numeric(tab[[col]]))[i]
    if (all(bb_corner_cols %in% names(tab))) {
      polygons <- TRUE
      for (cc in bb_corner_cols) new[[cc]] <- num(cc)
    } else {
      lat_col <- intersect(c("latitude", "lat", "WGS84_y"), names(tab))[1]
      lon_col <- intersect(c("longitude", "lon", "long", "WGS84_x"), names(tab))[1]
      if (is.na(lat_col) || is.na(lon_col)) {
        stop("`coordinates` needs `latitude` and `longitude` columns (or `WGS84_y` and `WGS84_x`), ",
             "or the eight corner columns of the Breedbase GPS template (", bb_corner_cols[1], ", ...).",
             call. = FALSE)
      }
      new$latitude <- num(lat_col)
      new$longitude <- num(lon_col)
    }
  } else {
    polygons <- TRUE
    if (is.null(plot_width) || is.null(plot_length)) {
      stop("With `origin`, give the plot size in metres: `plot_width` and `plot_length`.", call. = FALSE)
    }
    for (a in list(list(plot_width, "plot_width"), list(plot_length, "plot_length"))) {
      if (!is.numeric(a[[1]]) || length(a[[1]]) != 1 || !is.finite(a[[1]]) || a[[1]] <= 0) {
        stop("`", a[[2]], "` must be one positive number (metres).", call. = FALSE)
      }
    }
    for (a in list(list(col_gap, "col_gap"), list(row_gap, "row_gap"))) {
      if (!is.numeric(a[[1]]) || length(a[[1]]) != 1 || !is.finite(a[[1]]) || a[[1]] < 0) {
        stop("`", a[[2]], "` must be one number of metres, 0 or more.", call. = FALSE)
      }
    }
    turn <- (row_bearing - col_bearing) %% 360
    if (min(abs(turn - 90), abs(turn - 270)) > 0.5) {
      stop("`row_bearing` and `col_bearing` must be at right angles (e.g. 0 and 90).", call. = FALSE)
    }
    if (!all(c("row", "col") %in% names(d)) || anyNA(d$row) || anyNA(d$col)) {
      stop("The design needs `row` and `col` for every plot to place the plots from `origin`",
           if (inherits(design, "augmented_prep_design")) paste0(
             ": make the p-rep design with a field layout (`field_cols`, or `rows_per_block` and ",
             "`blocks_per_row`)") else "",
           ". Or give one point per plot with `latitude` and `longitude`, or a `coordinates` table.", call. = FALSE)
    }
    is_prep <- inherits(design, "augmented_prep_design")
    locs <- if (is_prep) unique(as.character(d$location)) else NA
    for (l in locs) {
      o <- if (is.list(origin)) {
        if (is.na(l) || is.null(names(origin)) || !l %in% names(origin)) {
          stop("Give `origin` as list(", if (is.na(l)) "c(latitude, longitude)" else
            paste0(paste0(locs, " = c(lat, lon)"), collapse = ", "), ").", call. = FALSE)
        }
        origin[[l]]
      } else {
        if (is_prep && length(locs) > 1) {
          stop("This p-rep design has ", length(locs), " locations: give `origin` as a list with one ",
               "c(latitude, longitude) per location, named ", paste(locs, collapse = ", "), ".", call. = FALSE)
        }
        origin
      }
      o <- suppressWarnings(as.numeric(o))
      if (length(o) != 2 || anyNA(o)) stop("`origin` must be c(latitude, longitude).", call. = FALSE)
      bb_check_latlon(o[1], o[2], "origin")
      rows <- if (is.na(l)) seq_len(n) else which(d$location == l)
      corners <- bb_plot_corners(d$row[rows], d$col[rows], o[1], o[2], plot_width, plot_length,
                                 col_gap, row_gap, col_bearing, row_bearing)
      for (cc in bb_corner_cols) new[[cc]][rows] <- corners[[cc]]
    }
  }
  if (polygons) {
    new$latitude <- rowMeans(new[, bb_corner_cols[c(2, 4, 6, 8)]])
    new$longitude <- rowMeans(new[, bb_corner_cols[c(1, 3, 5, 7)]])
  }
  bb_check_latlon(new$latitude, new$longitude, "coordinates")
  if (all(is.na(new$latitude))) stop("No plot got coordinates: check the values given.", call. = FALSE)
  keep <- if ("is_filler" %in% names(d)) d$is_filler != 1 else rep(TRUE, n)
  missing <- sum(keep & is.na(new$latitude))
  if (missing) message(missing, " plot", if (missing > 1) "s have" else " has", " no coordinates.")

  out <- design$design
  for (cc in c("latitude", "longitude", bb_corner_cols)) out[[cc]] <- NULL
  out$latitude <- new$latitude
  out$longitude <- new$longitude
  if (polygons) for (cc in bb_corner_cols) out[[cc]] <- new[[cc]]
  design$design <- out
  design
}

# Rows of `tab` for each plot of the design: by plot name when both have it,
# otherwise by plot number (and location for p-rep designs).
bb_match_coordinate_rows <- function(design, d, tab) {
  if (!nrow(tab)) stop("`coordinates` has no rows.", call. = FALSE)
  if ("plot_name" %in% names(tab) && "plot_name" %in% names(d)) {
    key_t <- trimws(as.character(tab$plot_name))
    key_d <- as.character(d$plot_name)
    by <- "plot_name"
  } else {
    plot_col <- intersect(c("plots", "plot_number", "plot"), names(tab))[1]
    if (is.na(plot_col)) {
      stop("`coordinates` needs a `plot_name` column (matching the design's plot names) or a plot number ",
           "column (`plots`, `plot_number` or `plot`).", call. = FALSE)
    }
    is_prep <- inherits(design, "augmented_prep_design") && length(unique(d$location)) > 1
    if (is_prep && !"location" %in% names(tab)) {
      stop("Plot numbers repeat across the locations of this p-rep design: add a `location` column to ",
           "`coordinates`, or match by `plot_name`.", call. = FALSE)
    }
    key_t <- paste(if (is_prep) tab$location else "", suppressWarnings(as.integer(tab[[plot_col]])), sep = "\r")
    key_d <- paste(if (is_prep) d$location else "", d$plots, sep = "\r")
    by <- plot_col
  }
  if (anyDuplicated(key_t[!is.na(key_t)])) stop("`coordinates` has repeated plots (by ", by, ").", call. = FALSE)
  unknown <- setdiff(key_t, key_d)
  if (length(unknown)) {
    stop(length(unknown), " row", if (length(unknown) > 1) "s" else "", " of `coordinates` match no plot of the ",
         "design (by ", by, "), e.g. ", sub("^.*\r", "", unknown[1]), ".", call. = FALSE)
  }
  match(key_d, key_t)
}

bb_check_latlon <- function(lat, lon, what) {
  ok <- is.na(lat) | is.na(lon) | (abs(lat) <= 90 & abs(lon) <= 180)
  if (!all(ok)) {
    stop("The ", what, " are not latitude and longitude in decimal degrees (latitude -90 to 90, longitude ",
         "-180 to 180), e.g. ", lat[!ok][1], ", ", lon[!ok][1], ". Are they swapped, or in another system?",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Corners of each plot from its row and column. Local metres (x across
# columns, y across rows) are turned to east/north by the bearings, then to
# degrees with the length of a degree at the origin's latitude.
bb_plot_corners <- function(row, col, lat0, lon0, plot_width, plot_length, col_gap, row_gap,
                            col_bearing, row_bearing) {
  rad <- pi / 180
  phi <- lat0 * rad
  m_lat <- 111132.92 - 559.82 * cos(2 * phi) + 1.175 * cos(4 * phi) - 0.0023 * cos(6 * phi)
  m_lon <- 111412.84 * cos(phi) - 93.5 * cos(3 * phi) + 0.118 * cos(5 * phi)
  cb <- col_bearing * rad
  rb <- row_bearing * rad
  x0 <- (col - 1) * (plot_width + col_gap)
  y0 <- (row - 1) * (plot_length + row_gap)
  at <- function(x, y) {
    east <- x * sin(cb) + y * sin(rb)
    north <- x * cos(cb) + y * cos(rb)
    list(lon = lon0 + east / m_lon, lat = lat0 + north / m_lat)
  }
  bl <- at(x0, y0)
  br <- at(x0 + plot_width, y0)
  tr <- at(x0 + plot_width, y0 + plot_length)
  tl <- at(x0, y0 + plot_length)
  stats::setNames(list(bl$lon, bl$lat, br$lon, br$lat, tr$lon, tr$lat, tl$lon, tl$lat), bb_corner_cols)
}

# GeoJSON of one plot (a row of the design table), as Breedbase stores it;
# NULL without coordinates.
bb_plot_geojson <- function(x) {
  if (all(bb_corner_cols %in% names(x))) {
    v <- vapply(bb_corner_cols, function(cc) as.numeric(x[[cc]]), numeric(1))
    if (!anyNA(v)) {
      ring <- list(v[1:2], v[3:4], v[5:6], v[7:8], v[1:2])
      return(list(type = "Feature",
                  geometry = list(type = "Polygon", coordinates = list(lapply(ring, unname))),
                  properties = list(format = "WGS84")))
    }
  }
  if (all(c("latitude", "longitude") %in% names(x))) {
    lat <- suppressWarnings(as.numeric(x[["latitude"]]))
    lon <- suppressWarnings(as.numeric(x[["longitude"]]))
    if (length(lat) == 1 && !is.na(lat) && !is.na(lon)) {
      return(list(type = "Feature", geometry = list(type = "Point", coordinates = c(lon, lat)),
                  properties = list(format = "WGS84")))
    }
  }
  NULL
}

# Send the GPS of the plots of one saved trial (PUT /observationunits, a map
# of observationUnitDbId to the unit) and read them back. `geo` is a list of
# GeoJSON named by plot number; `accessions` (named by plot number) are
# checked against the saved plots first when given.
bb_save_plot_geo <- function(con, study_id, geo, trial, accessions = NULL) {
  units <- bb_get_all(con, "observationunits", list(studyDbId = study_id, observationUnitLevelName = "plot"))
  nums <- vapply(units, function(u) as.character(u$observationUnitPosition$observationLevel$levelCode %||% NA),
                 character(1))
  ids <- vapply(units, function(u) as.character(u$observationUnitDbId %||% NA), character(1))
  missing <- setdiff(names(geo), nums)
  if (length(missing)) {
    stop("Trial '", trial, "' has no plot number ", paste(utils::head(missing, 5), collapse = ", "),
         " in Breedbase: is it the trial of this design?", call. = FALSE)
  }
  if (!is.null(accessions)) {
    saved <- vapply(units, function(u) as.character(u$germplasmName %||% NA), character(1))[match(names(geo), nums)]
    differ <- which(saved != accessions[names(geo)])
    if (length(differ)) {
      p <- names(geo)[differ[1]]
      stop(length(differ), " plots of trial '", trial, "' have another accession in Breedbase than in the design ",
           "(e.g. plot ", p, ": '", accessions[[p]], "' here, '", saved[differ[1]], "' in Breedbase). ",
           "Is it the trial of this design?", call. = FALSE)
    }
  }
  unit_id <- ids[match(names(geo), nums)]
  batch <- getOption("breedtools.geo_batch", 200)
  parts <- split(seq_along(geo), ceiling(seq_along(geo) / batch))
  for (k in seq_along(parts)) {
    body <- stats::setNames(lapply(parts[[k]], function(j) list(
      observationUnitPosition = list(observationLevel = list(levelName = "plot"), geoCoordinates = geo[[j]])
    )), unit_id[parts[[k]]])
    what <- paste0("save the GPS of the plots of '", trial, "'")
    # Each PUT replaces the GPS of its plots, so a second try is safe.
    tryCatch(bb_put(con, "observationunits", body, what = what), breedbase_error = function(e) {
      if (!isTRUE(e$busy)) stop(e)
      Sys.sleep(getOption("breedtools.retry_pause", 1) * 5)
      bb_put(con, "observationunits", body, what = what)
    })
  }
  back <- bb_get_all(con, "observationunits", list(studyDbId = study_id, observationUnitLevelName = "plot"))
  has_geo <- vapply(back, function(u) length(u$observationUnitPosition$geoCoordinates) > 0, logical(1))
  back_nums <- vapply(back, function(u) as.character(u$observationUnitPosition$observationLevel$levelCode %||% NA),
                      character(1))
  n_found <- sum(has_geo & back_nums %in% names(geo))
  status <- if (n_found == length(geo)) "[OK]  " else "[FAIL]"
  bt_say(sprintf("  %s GPS coordinates of trial '%s' saved for %d of %d plots\n", status, trial, n_found, length(geo)))
  n_found
}

# GeoJSON of the plots of one study record, named by plot number.
bb_study_geo <- function(s) {
  geo <- lapply(s$plots, function(p) p$observationUnitPosition$geoCoordinates)
  names(geo) <- vapply(s$plots, function(p) p$observationUnitPosition$observationLevel$levelCode, character(1))
  Filter(Negate(is.null), geo)
}

#' Save plot GPS coordinates of a trial already in Breedbase
#'
#' Sends the plot coordinates of a design (see [add_plot_coordinates()]) to
#' a trial that is already saved in Breedbase, e.g. when the plots were
#' measured after the trial was saved. Plots are matched by plot number,
#' after checking that every plot has the same accession in Breedbase as in
#' the design. Coordinates already stored for these plots are replaced.
#' [save_breedbase()] does this by itself for a design that has coordinates.
#'
#' @param design A design with coordinates, from [add_plot_coordinates()].
#' @param con A connection from [breedbase_login()].
#' @param trial_name The trial name in Breedbase. For a multi-location p-rep
#'   design, the location is added as in [save_breedbase()]. Can be omitted
#'   for a prepared Breedbase trial or a design with plot names.
#' @param dry_run Logical. Check the trial and plots but write nothing.
#' @param confirm Logical. Ask before writing (default in interactive
#'   sessions).
#'
#' @return Invisibly, a tibble with one row per trial: `trial_name`,
#'   `plots_with_gps` (sent) and `gps_found` (read back from Breedbase).
#'
#' @examples
#' \dontrun{
#' con <- breedbase_login("https://cassavabase.org", "jdoe")
#' rc <- add_plot_coordinates(rc, origin = c(7.4931, 3.9035), plot_width = 3, plot_length = 5)
#' save_breedbase_coordinates(rc, con, trial_name = "AYT_2026_Ibadan")
#' }
#'
#' @export
save_breedbase_coordinates <- function(design, con, trial_name = NULL, dry_run = FALSE, confirm = interactive()) {
  if (!inherits(con, "breedbase_connection")) stop("`con` must come from breedbase_login().", call. = FALSE)
  if (is.null(trial_name) && inherits(design, "breedbase_trial")) trial_name <- design$trial$trial_name
  if (is.null(trial_name)) trial_name <- design$plot_name_trial
  if (is.null(trial_name)) stop("Give the `trial_name` of the saved trial.", call. = FALSE)
  # Plots are matched by number, so any plot names do here.
  studies <- bb_design_to_studies(design, trial_name, year = 2000, location = bb_template_location(design),
                                  description = "x", trial_type = "x", design_type = NULL,
                                  plot_name_pattern = "{trial}_{plot}")
  geos <- lapply(studies, bb_study_geo)
  if (!sum(lengths(geos))) stop("The design has no plot coordinates: use add_plot_coordinates() first.", call. = FALSE)
  bt_say("Checking before saving GPS coordinates to ", con$url, "\n", sep = "")
  ids <- character(length(studies))
  for (i in seq_along(studies)) {
    name <- studies[[i]]$study$studyName
    ids[i] <- bb_find_study(con, name)
    if (is.na(ids[i])) stop("Trial '", name, "' was not found in Breedbase.", call. = FALSE)
    bt_say("  [OK]   trial '", name, "': ", length(geos[[i]]), " plots with coordinates\n", sep = "")
  }
  out <- tibble::tibble(trial_name = vapply(studies, function(s) s$study$studyName, character(1)),
                        plots_with_gps = lengths(geos), gps_found = NA_integer_)
  if (isTRUE(dry_run)) {
    bt_say("Dry run: nothing was written.\n")
    return(invisible(out))
  }
  if (isTRUE(confirm)) {
    answer <- readline(sprintf("Save the GPS of %d plots (replacing any stored)? Type 'yes' to continue: ",
                               sum(out$plots_with_gps)))
    if (!identical(tolower(trimws(answer)), "yes")) {
      message("Nothing was saved.")
      return(invisible(NULL))
    }
  }
  for (i in seq_along(studies)) {
    s <- studies[[i]]
    acc <- stats::setNames(vapply(s$plots, `[[`, character(1), "germplasmName"),
                           vapply(s$plots, function(p) p$observationUnitPosition$observationLevel$levelCode, character(1)))
    out$gps_found[i] <- bb_save_plot_geo(con, ids[i], geos[[i]], s$study$studyName, accessions = acc)
  }
  invisible(out)
}
