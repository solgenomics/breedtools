# Looking at trait data before and after the first stage: plots, outliers,
# and raw versus adjusted values.

#' Find outliers in trait data (IQR rule)
#'
#' Flags values below \eqn{Q_1 - k \cdot IQR} or above
#' \eqn{Q_3 + k \cdot IQR}, where \eqn{Q_1}, \eqn{Q_3} are the quartiles and
#' \eqn{k} is `multiplier` (Tukey's fences; 1.5 by default, 3 for only the
#' extreme values). The rule is applied within each group of `by`, usually
#' each location, because locations have different means.
#'
#' @param data A data frame with one row per plot.
#' @param traits Names of the numeric trait columns.
#' @param multiplier The IQR multiplier \eqn{k}.
#' @param by Optional grouping column(s); by default `"location"` when the
#'   data have it.
#' @param id_cols Columns that identify the plot in the report (those
#'   present in `data`).
#'
#' @return A tibble with one row per outlier: the identifying columns, the
#'   group, `trait`, `value`, the fences `lower` and `upper`, and `side`
#'   (`"low"` or `"high"`).
#'
#' @examples
#' d <- data.frame(plots = 1:20, all_entries = paste0("G", 1:20),
#'                 yield = c(rnorm(18, 50, 2), 80, 10))
#' find_outliers(d, "yield")
#'
#' @export
find_outliers <- function(data, traits, multiplier = 1.5,
                          by = if ("location" %in% names(data)) "location" else NULL,
                          id_cols = c("location", "plot_name", "plots", "plot_number", "all_entries", "accession",
                                      "row", "col")) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.numeric(multiplier) || length(multiplier) != 1 || multiplier <= 0) {
    stop("`multiplier` must be one positive number, e.g. 1.5.", call. = FALSE)
  }
  absent <- setdiff(c(traits, by), names(data))
  if (length(absent)) stop("Column(s) not found: ", paste(absent, collapse = ", "), call. = FALSE)
  id_cols <- intersect(unique(c(by, id_cols)), names(data))
  group <- if (length(by)) interaction(data[by], drop = TRUE, lex.order = TRUE) else factor(rep("all", nrow(data)))
  rows <- list()
  for (t in traits) {
    x <- suppressWarnings(as.numeric(data[[t]]))
    for (g in levels(group)) {
      i <- which(group == g)
      v <- x[i]
      if (sum(!is.na(v)) < 4) next
      q <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      iqr <- q[2] - q[1]
      lower <- q[1] - multiplier * iqr
      upper <- q[2] + multiplier * iqr
      out <- i[!is.na(v) & (v < lower | v > upper)]
      if (!length(out)) next
      part <- tibble::as_tibble(data[out, id_cols, drop = FALSE])
      part$trait <- t
      part$value <- x[out]
      part$lower <- lower
      part$upper <- upper
      part$side <- ifelse(x[out] < lower, "low", "high")
      part$.row <- out
      rows[[length(rows) + 1]] <- part
    }
  }
  res <- if (length(rows)) dplyr::bind_rows(rows) else
    tibble::tibble(trait = character(), value = numeric(), lower = numeric(), upper = numeric(),
                   side = character(), .row = integer())
  attr(res, "multiplier") <- multiplier
  res
}

#' Remove outliers from trait data (IQR rule)
#'
#' Finds outliers with [find_outliers()] and sets those values to `NA`, so
#' the plots stay in the design but the values are left out of the
#' analysis. The plots and values removed are printed and returned with the
#' data.
#'
#' @inheritParams find_outliers
#' @param quiet Logical. If `FALSE` (default), print what was removed.
#'
#' @return `data` with the outlying values set to `NA`. The removed values
#'   are in `attr(result, "outliers")` (as from [find_outliers()]).
#'
#' @examples
#' d <- data.frame(plots = 1:20, all_entries = paste0("G", 1:20),
#'                 yield = c(rnorm(18, 50, 2), 80, 10))
#' clean <- remove_outliers(d, "yield")
#' attr(clean, "outliers")
#'
#' @export
remove_outliers <- function(data, traits, multiplier = 1.5,
                            by = if ("location" %in% names(data)) "location" else NULL,
                            id_cols = c("location", "plot_name", "plots", "plot_number", "all_entries",
                                        "accession", "row", "col"),
                            quiet = FALSE) {
  found <- find_outliers(data, traits, multiplier = multiplier, by = by, id_cols = id_cols)
  for (t in unique(found$trait)) data[[t]][found$.row[found$trait == t]] <- NA
  if (!isTRUE(quiet)) {
    n_values <- sum(vapply(traits, function(t) sum(!is.na(suppressWarnings(as.numeric(data[[t]])))), integer(1))) +
      nrow(found)
    if (!nrow(found)) {
      bt_say("No outliers (IQR x ", multiplier, ") among ", n_values, " values", if (length(by)) paste0(" (within each ", paste(by, collapse = ", "), ")"), ".\n", sep = "")
    } else {
      bt_say("Removed ", nrow(found), " of ", n_values, " values as outliers (outside Q1 - ", multiplier,
          " IQR and Q3 + ", multiplier, " IQR", if (length(by)) paste0(", within each ", paste(by, collapse = ", ")),
          "):\n", sep = "")
      show <- as.data.frame(found[, setdiff(names(found), ".row"), drop = FALSE])
      bt_show(show, row.names = FALSE, digits = 4)
    }
  }
  attr(data, "outliers") <- found[, setdiff(names(found), ".row"), drop = FALSE]
  data
}

#' Plot a trait: boxplot, distribution and field map
#'
#' Draws the trait per location (or other group) so gradients, odd plots and
#' outliers are easy to see:
#' \describe{
#'   \item{`"boxplot"`}{All locations side by side. Outliers (IQR rule, as in
#'     [find_outliers()], or the table given in `outliers`) are drawn in red
#'     and labelled with their plot number; the other plots in grey.}
#'   \item{`"histogram"`}{The distribution at each location.}
#'   \item{`"field"`}{A map of the field with the value of every plot, with a
#'     colour bar. The colour scale of each location spans the values inside
#'     the IQR fences, so a few extreme plots do not wash out the rest of the
#'     map; plots beyond the scale get the end colour and are circled in red.
#'     Missing plots are grey.}
#' }
#'
#' @param data A data frame with one row per plot.
#' @param trait Name of the numeric trait column.
#' @param by Optional grouping column; by default `"location"` when the data
#'   have it.
#' @param type `"both"` (histogram and field map), `"histogram"`, `"field"`,
#'   `"boxplot"` or `"all"` (boxplot, then histogram and field map).
#' @param row_col,col_col Field row and column columns (for the map).
#' @param outliers Optional table from [find_outliers()] or
#'   `attr(remove_outliers(...), "outliers")`: these plots are marked. By
#'   default the outliers are found in `data` with the IQR rule.
#' @param breaks Histogram breaks, as in [graphics::hist()].
#' @param multiplier IQR multiplier for the default outliers and the colour
#'   scale.
#' @param scale `"location"` (default: one colour scale per location, so
#'   each field shows its own gradient) or `"common"` (one scale for all
#'   locations, to compare their levels).
#' @param plot_col Column with the plot numbers used as labels.
#'
#' @return Invisibly, a tibble with a summary per group: plots, values,
#'   missing, outliers, mean, sd, min, max and CV (%).
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' d$yield <- 50 + 0.8 * d$row + rnorm(nrow(d), sd = 2)
#' d$yield[c(5, 40)] <- c(120, 10)
#' plot_trait(d, "yield", type = "boxplot")
#' plot_trait(d, "yield", type = "field")
#'
#' @export
plot_trait <- function(data, trait, by = if ("location" %in% names(data)) "location" else NULL,
                       type = c("both", "histogram", "field", "boxplot", "all"), row_col = "row", col_col = "col",
                       outliers = NULL, breaks = "Sturges", multiplier = 1.5, scale = c("location", "common"),
                       plot_col = intersect(c("plots", "plot_number"), names(data))[1]) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  absent <- setdiff(c(trait, by), names(data))
  if (length(absent)) stop("Column(s) not found: ", paste(absent, collapse = ", "), call. = FALSE)
  has_map <- all(c(row_col, col_col) %in% names(data))
  if (type == "field" && !has_map) stop("The field map needs `", row_col, "` and `", col_col, "` columns.", call. = FALSE)
  if (type == "both" && !has_map) type <- "histogram"
  if (type == "all" && !has_map) type <- "boxplot"
  x <- suppressWarnings(as.numeric(data[[trait]]))
  group <- if (length(by)) factor(as.character(data[[by]]), levels = unique(as.character(data[[by]]))) else
    factor(rep(trait, nrow(data)))
  groups <- levels(group)
  # Outliers as a logical per row: from the table given, else the IQR rule.
  is_out <- bb_outlier_rows(data, trait, x, group, outliers, by, row_col, col_col, plot_col, multiplier)
  labels <- if (!is.na(plot_col) && plot_col %in% names(data)) as.character(data[[plot_col]]) else
    as.character(seq_len(nrow(data)))
  common_zlim <- bb_robust_range(x[!is_out], multiplier)

  if (type %in% c("boxplot", "all")) {
    bb_trait_boxplot(x, group, is_out, labels, trait, by, multiplier)
    if (type == "all" && interactive()) grDevices::devAskNewPage(TRUE)
  }
  if (type %in% c("both", "histogram", "field", "all")) {
    per_group <- if (type %in% c("both", "all")) 2 else 1
    ncol <- if (per_group == 2) 2 else min(length(groups), 3)
    old <- graphics::par(mfrow = c(ceiling(length(groups) * per_group / ncol), ncol), mar = c(4, 4, 3, 6))
    on.exit(graphics::par(old), add = TRUE)
    if (type == "all" && interactive()) on.exit(grDevices::devAskNewPage(FALSE), add = TRUE)
    for (g in groups) {
      i <- which(group == g)
      title <- if (length(by)) paste(trait, "-", g) else trait
      if (type %in% c("both", "histogram", "all")) {
        graphics::hist(x[i], breaks = breaks, col = "steelblue", border = "white", main = title,
                       xlab = trait, ylab = "Plots")
        graphics::abline(v = mean(x[i], na.rm = TRUE), lty = 2)
        if (any(is_out[i])) graphics::rug(x[i][is_out[i]], col = "red", lwd = 2)
      }
      if (type %in% c("both", "field", "all")) {
        zlim <- if (scale == "common") common_zlim else bb_robust_range(x[i][!is_out[i]], multiplier)
        bb_field_image(data[[row_col]][i], data[[col_col]][i], x[i], main = paste(title, "(field)"),
                       zlim = zlim, mark = is_out[i], legend = trait)
      }
    }
  }
  summary <- dplyr::bind_rows(lapply(groups, function(g) {
    v <- x[group == g]
    n_out <- sum(is_out[group == g])
    tibble::tibble(group = g, plots = length(v), values = sum(!is.na(v)), missing = sum(is.na(v)),
                   outliers = n_out,
                   mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE),
                   min = suppressWarnings(min(v, na.rm = TRUE)), max = suppressWarnings(max(v, na.rm = TRUE)),
                   cv_percent = 100 * stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE))
  }))
  invisible(summary)
}

# Rows of `data` that are outliers of `trait`: from a find_outliers() table
# (matched by group and plot, or row/col), else by the IQR rule.
bb_outlier_rows <- function(data, trait, x, group, outliers, by, row_col, col_col, plot_col, multiplier) {
  if (is.null(outliers)) {
    found <- find_outliers(data, trait, multiplier = multiplier, by = by)
    out <- rep(FALSE, nrow(data))
    out[found$.row] <- TRUE
    return(out)
  }
  o <- outliers[as.character(outliers$trait) == trait, , drop = FALSE]
  if (!nrow(o)) return(rep(FALSE, nrow(data)))
  g_data <- if (length(by)) as.character(data[[by]]) else ""
  g_out <- if (length(by) && by %in% names(o)) as.character(o[[by]]) else ""
  if (!is.na(plot_col) && plot_col %in% names(o)) {
    key_data <- paste(g_data, data[[plot_col]], sep = "\r")
    key_out <- paste(g_out, o[[plot_col]], sep = "\r")
  } else if (all(c(row_col, col_col) %in% names(o))) {
    key_data <- paste(g_data, data[[row_col]], data[[col_col]], sep = "\r")
    key_out <- paste(g_out, o[[row_col]], o[[col_col]], sep = "\r")
  } else {
    stop("`outliers` needs a plot or row/col column to match the data.", call. = FALSE)
  }
  key_data %in% key_out
}

# Range of the values inside the IQR fences (the colour scale of a map).
bb_robust_range <- function(x, multiplier = 1.5) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  if (length(x) >= 4) {
    q <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
    x <- x[x >= q[1] - multiplier * diff(q) & x <= q[2] + multiplier * diff(q)]
  }
  range(x)
}

# Boxplot per group, one panel each with its own scale: the other plots grey
# (jittered), outliers red and labelled with the plot number. Outliers far
# beyond the box are drawn as triangles at the panel edge, with their value,
# so they do not squash the box.
bb_trait_boxplot <- function(x, group, is_out, labels, trait, by, multiplier = 1.5) {
  groups <- levels(group)
  old <- graphics::par(mfrow = c(1, length(groups)), mar = c(3, 4, 4, 1), oma = c(0, 0, 2, 0))
  on.exit(graphics::par(old), add = TRUE)
  for (g in groups) {
    i <- which(group == g)
    v <- x[i]
    out <- is_out[i] & !is.na(v)
    keep <- v[!out & !is.na(v)]
    if (!length(keep)) keep <- v[!is.na(v)]
    if (!length(keep)) {
      graphics::plot.new()
      graphics::title(main = g)
      next
    }
    r <- range(keep)
    pad <- 0.6 * max(diff(r), abs(r[1]) * 0.05, 1e-8)
    ylim <- c(r[1] - pad, r[2] + pad)
    graphics::boxplot(keep, outline = FALSE, col = "grey95", border = "grey40", ylim = ylim, ylab = trait,
                      main = g, xlim = c(0.4, 1.8))
    # Deterministic jitter (no change to the user's random numbers).
    jitter <- ((seq_along(keep) * 0.618034) %% 1 - 0.5) * 0.3
    graphics::points(1 + jitter, keep, pch = 16, cex = 0.7,
                     col = grDevices::adjustcolor("grey30", 0.5))
    if (any(out)) {
      vo <- v[out]
      lo <- labels[i][out]
      inside <- vo >= ylim[1] & vo <= ylim[2]
      if (any(inside)) {
        graphics::points(rep(1, sum(inside)), vo[inside], pch = 21, cex = 1.4, bg = "red", col = "darkred")
        graphics::text(rep(1.08, sum(inside)), vo[inside], lo[inside], adj = 0, cex = 0.75, col = "darkred")
      }
      for (side in c("high", "low")) {
        far <- if (side == "high") vo > ylim[2] else vo < ylim[1]
        if (!any(far)) next
        y <- if (side == "high") ylim[2] else ylim[1]
        ord <- order(if (side == "high") -vo[far] else vo[far])
        graphics::points(1, y, pch = if (side == "high") 24 else 25, cex = 1.6, bg = "red", col = "darkred")
        txt <- paste0("plot ", lo[far][ord], ": ", format(signif(vo[far][ord], 4)), collapse = "\n")
        graphics::text(1.12, y, txt, adj = c(0, if (side == "high") 1 else 0), cex = 0.7, col = "darkred")
      }
    }
    graphics::mtext(paste0(sum(out), " outlier", if (sum(out) == 1) "" else "s", " of ", sum(!is.na(v)), " plots"),
                    side = 3, line = 0.2, cex = 0.7, col = "darkred")
  }
  graphics::mtext(paste0(trait, ": outliers (red, labelled by plot; triangles = beyond the axis)"),
                  outer = TRUE, cex = 0.95, font = 2)
}

# Field map: one cell per plot, row 1 at the top; NA plots grey; values
# beyond `zlim` get the end colours; `mark` (logical per plot, or a table
# with row and col) circled in red; a colour bar on the right.
bb_field_image <- function(row, col, value, main, zlim = range(value, na.rm = TRUE), mark = NULL, legend = NULL) {
  row <- as.integer(row)
  col <- as.integer(col)
  nr <- max(row, na.rm = TRUE)
  nc <- max(col, na.rm = TRUE)
  pal <- grDevices::hcl.colors(50, "RdYlBu", rev = TRUE)
  ok_z <- all(is.finite(zlim)) && diff(zlim) > 0
  shown <- if (ok_z) pmin(pmax(value, zlim[1]), zlim[2]) else value
  m <- matrix(NA_real_, nr, nc)
  m[cbind(row, col)] <- shown
  graphics::plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(nr + 0.5, 0.5), xlab = "Column", ylab = "Row",
                 main = main, xaxs = "i", yaxs = "i")
  graphics::rect(0.5, nr + 0.5, nc + 0.5, 0.5, col = "grey80", border = NA)
  if (ok_z) {
    graphics::image(seq_len(nc), seq_len(nr), t(m), col = pal, zlim = zlim, add = TRUE)
  } else {
    graphics::image(seq_len(nc), seq_len(nr), t(m), col = pal[25], add = TRUE)
  }
  graphics::abline(h = seq_len(nr + 1) - 0.5, v = seq_len(nc + 1) - 0.5, col = "white", lwd = 0.5)
  if (is.logical(mark)) {
    if (any(mark, na.rm = TRUE)) {
      graphics::points(col[which(mark)], row[which(mark)], pch = 1, cex = 2.2, lwd = 2, col = "red")
    }
  } else if (!is.null(mark) && nrow(mark)) {
    graphics::points(mark$col, mark$row, pch = 1, cex = 2.2, lwd = 2, col = "red")
  }
  if (ok_z) bb_color_bar(zlim, pal, legend)
  graphics::mtext(if (is.null(mark)) "grey = missing" else "grey = missing   red circle = outlier", side = 3,
                  line = 0.2, cex = 0.65)
}

# Vertical colour bar just right of the current plot, with tick labels.
bb_color_bar <- function(zlim, pal, label = NULL) {
  usr <- graphics::par("usr")
  w <- diff(usr[1:2])
  x0 <- usr[2] + 0.04 * w
  x1 <- x0 + 0.05 * w
  y <- seq(usr[3], usr[4], length.out = length(pal) + 1)
  graphics::rect(x0, y[-length(y)], x1, y[-1], col = pal, border = NA, xpd = TRUE)
  graphics::rect(x0, usr[3], x1, usr[4], border = "grey40", xpd = TRUE)
  ticks <- pretty(zlim, n = 4)
  ticks <- ticks[ticks >= zlim[1] & ticks <= zlim[2]]
  ty <- usr[3] + (ticks - zlim[1]) / diff(zlim) * (usr[4] - usr[3])
  graphics::segments(x1, ty, x1 + 0.015 * w, ty, xpd = TRUE)
  graphics::text(x1 + 0.02 * w, ty, format(ticks), adj = 0, cex = 0.7, xpd = TRUE)
  if (!is.null(label)) graphics::text(x0, usr[4], label, adj = c(0, -0.6), cex = 0.7, xpd = TRUE)
}

#' Compare raw and first-stage adjusted values
#'
#' Puts each plot's raw value next to its first-stage value (the genotype's
#' adjusted mean plus the plot's residual, as in [first_stage_values()]),
#' summarises how much the first stage changed (per location: correlation of
#' raw and adjusted plot values, the spread before and after, and how much
#' the ranking of the genotypes changed) and draws, per location, the raw
#' and the adjusted field maps and a raw-versus-adjusted plot.
#'
#' @param fit An analysis (e.g. [analyze_augmented_row_column_design()] with
#'   `spatial = TRUE`) or a two-stage analysis ([analyze_two_stage()],
#'   [analyze_augmented_prep_two_stage()]), whose stage-1 fits are compared.
#' @param plot Logical. Draw the plots.
#' @param locations For two-stage analyses: the locations to draw (default
#'   all).
#'
#' @return Invisibly, a list with `plots` (one row per plot: location, plot,
#'   genotype, row, col, `raw`, `adjusted`, `change`) and `summary` (one row
#'   per location). The summary is also printed.
#'
#' @examples
#' rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
#'                                   12, 12, 4, 4, n_candidates = 5, seed = 1)
#' d <- rc$design
#' set.seed(1)
#' d$yield <- 50 + 0.8 * d$row + rnorm(nrow(d), sd = 2)
#' fit <- analyze_augmented_row_column_design(d, "yield")
#' cmp <- compare_first_stage(fit)
#'
#' @export
compare_first_stage <- function(fit, plot = TRUE, locations = NULL) {
  fits <- if (!is.null(fit$stage1_fits) && is.null(fit$plot_results)) fit$stage1_fits else list(fit)
  if (is.null(names(fits)) || length(fits) == 1 && identical(fits[[1]], fit)) names(fits) <- "trial"
  per_plot <- dplyr::bind_rows(lapply(names(fits), function(loc) bb_compare_one(fits[[loc]], loc)))
  summary <- dplyr::bind_rows(lapply(unique(per_plot$location), function(loc) {
    p <- per_plot[per_plot$location == loc & !is.na(per_plot$raw) & !is.na(per_plot$adjusted), , drop = FALSE]
    raw_g <- tapply(p$raw, p$genotype, mean)
    adj_g <- tapply(p$adjusted, p$genotype, mean)
    top <- max(1L, ceiling(0.1 * length(raw_g)))
    tibble::tibble(
      location = loc, plots = nrow(p),
      cor_raw_adjusted = stats::cor(p$raw, p$adjusted),
      sd_raw = stats::sd(p$raw), sd_adjusted = stats::sd(p$adjusted),
      genotype_rank_cor = stats::cor(raw_g, adj_g[names(raw_g)], method = "spearman"),
      top10pct_kept = mean(names(sort(-raw_g))[seq_len(top)] %in% names(sort(-adj_g))[seq_len(top)])
    )
  }))
  bt_say("Raw versus first-stage adjusted values:\n")
  bt_show(as.data.frame(summary), row.names = FALSE, digits = 3)
  bt_say("(genotype_rank_cor: Spearman correlation of genotype means; top10pct_kept: share of the raw top 10%",
      "still in the adjusted top 10%)\n")
  if (isTRUE(plot)) {
    draw <- if (is.null(locations)) unique(per_plot$location) else intersect(locations, unique(per_plot$location))
    has_map <- !all(is.na(per_plot$row)) && !all(is.na(per_plot$col))
    old <- graphics::par(mfrow = c(length(draw), if (has_map) 3 else 1), mar = c(4, 4, 3, 6))
    on.exit(graphics::par(old), add = TRUE)
    for (loc in draw) {
      p <- per_plot[per_plot$location == loc, , drop = FALSE]
      label <- if (loc == "trial") "" else paste0(" - ", loc)
      if (has_map) {
        zlim <- range(c(p$raw, p$adjusted), na.rm = TRUE)
        bb_field_image(p$row, p$col, p$raw, main = paste0("Raw", label), zlim = zlim)
        bb_field_image(p$row, p$col, p$adjusted, main = paste0("Adjusted (first stage)", label), zlim = zlim)
      }
      graphics::plot(p$raw, p$adjusted, pch = 19, col = grDevices::adjustcolor("steelblue", 0.6),
                     xlab = "Raw value", ylab = "Adjusted value", main = paste0("Raw vs adjusted", label))
      graphics::abline(0, 1, lty = 2)
    }
  }
  invisible(list(plots = per_plot, summary = summary))
}

bb_compare_one <- function(fit, loc) {
  values <- bb_first_stage_plots(fit, loc)
  if (inherits(fit, "SpATS")) {
    d <- fit$data
    return(tibble::tibble(location = loc, plots = as.integer(d$plots), genotype = as.character(d$geno),
                          row = as.numeric(d$row), col = as.numeric(d$col), raw = d$y,
                          adjusted = values$value, change = values$value - d$y))
  }
  pr <- fit$plot_results
  response <- fit$input$response
  gcol <- fit$input$genotype_col %||% "all_entries"
  rcol <- if ("row" %in% names(pr)) "row" else NA
  ccol <- if ("col" %in% names(pr)) "col" else NA
  raw <- if (!is.null(response) && response %in% names(pr)) suppressWarnings(as.numeric(pr[[response]])) else
    pr$fitted + pr$residual
  # A two-stage fit is grouped by its named stage-1 trials. Breedbase can
  # omit location labels, or several trials can share a physical location.
  location <- rep(loc, nrow(pr))
  location_col <- fit$input$location_col %||% "location"
  if (identical(loc, "trial") && location_col %in% names(pr)) {
    labels <- as.character(pr[[location_col]])
    present <- !is.na(labels) & nzchar(trimws(labels))
    location[present] <- labels[present]
  }
  tibble::tibble(
    location = location,
    plots = as.integer(pr$plots), genotype = as.character(pr[[gcol]]),
    row = if (is.na(rcol)) NA_real_ else as.numeric(pr[[rcol]]),
    col = if (is.na(ccol)) NA_real_ else as.numeric(pr[[ccol]]),
    raw = raw, adjusted = values$value, change = values$value - raw
  )
}
