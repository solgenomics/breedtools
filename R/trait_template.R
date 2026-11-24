#' Create a Field Book trait file
#'
#' Makes trait definitions for offline import into the Field Book Android
#' app. Use [phenotype_template()] with no trait columns for the companion
#' field CSV. No connection to Breedbase is needed to create either file.
#'
#' @param traits Character vector of unique trait names. Use names that
#'   Breedbase can resolve when uploading observations later.
#' @param format Field Book format: `"numeric"`, `"categorical"`, `"percent"`,
#'   `"date"`, `"boolean"`, `"text"`, `"photo"`, `"counter"`, `"location"`, or
#'   `"angle"`. A single value is repeated for all traits.
#' @param default_value,minimum,maximum,details Optional defaults, numeric
#'   bounds and descriptions, respectively. Each must have length one or
#'   the number of traits. `NA` becomes an empty cell. Units can be included
#'   in `details`.
#' @param categories For categorical traits, a character vector of
#'   slash-separated choices, or a list of character vectors, one per trait
#'   (e.g. `list(character(), c("low", "medium", "high"))`). A single
#'   value or list element is repeated. Other formats must leave this empty.
#' @param is_visible Logical; show the trait in Field Book's Collect screen.
#'   Length one or the number of traits.
#' @param output_file Optional path to a `.trt` file. Written as UTF-8 CSV
#'   with empty missing cells and no row names.
#'
#' @details
#' The columns follow Field Book's bundled
#' \href{https://github.com/PhenoApps/Field-Book/blob/main/app/src/main/assets/trait/trait_sample.trt}{sample trait file}.
#' Import the resulting `.trt` through Traits > Import/Export in the app.
#' See the \href{https://field-book.readthedocs.io/en/latest/traits.html}{Field Book trait documentation}.
#' These local definitions do not create ontology traits in Breedbase;
#' those must already exist before observations can be saved there.
#'
#' @return A tibble with `trait`, `format`, `defaultValue`, `minimum`,
#'   `maximum`, `details`, `categories`, `isVisible`, and `realPosition`.
#'   The table is also returned when `output_file` is supplied.
#' @seealso [phenotype_template()], [save_breedbase_phenotypes()]
#' @examples
#' traits <- trait_template(
#'   c("fresh root yield", "notes"), format = c("numeric", "text"),
#'   minimum = c(0, NA), details = c("Yield in t/ha", "Field observations")
#' )
#' trait_template("severity", format = "categorical",
#'                categories = list(c("low", "medium", "high")))
#' path <- tempfile(fileext = ".trt")
#' trait_template("fresh root yield", minimum = 0, output_file = path)
#' @export
trait_template <- function(traits, format = "numeric", default_value = "",
                           minimum = NA_real_, maximum = NA_real_, details = "",
                           categories = "", is_visible = TRUE, output_file = NULL) {
  if (!is.character(traits) || !length(traits) || anyNA(traits) ||
      any(!nzchar(trimws(traits))) || anyDuplicated(traits)) {
    stop("`traits` must contain unique, non-empty names without missing values.", call. = FALSE)
  }
  n <- length(traits)
  expand <- function(x, name) {
    if (!length(x) || !length(x) %in% c(1L, n)) {
      stop("`", name, "` must have length one or length(traits).", call. = FALSE)
    }
    rep(x, length.out = n)
  }
  chars <- function(x, name) {
    x <- as.character(expand(x, name))
    x[is.na(x)] <- ""
    x
  }
  format <- chars(format, "format")
  supported <- c("numeric", "categorical", "percent", "date", "boolean",
                 "text", "photo", "counter", "location", "angle")
  if (any(!format %in% supported)) {
    stop("Unsupported Field Book `format`. Choose: ", paste(supported, collapse = ", "), ".", call. = FALSE)
  }
  if (is.list(categories)) {
    categories <- vapply(categories, function(x) {
      if (!is.character(x) || anyNA(x) || any(!nzchar(trimws(x))) ||
          any(grepl("/", x, fixed = TRUE)) || anyDuplicated(x)) {
        stop("Category lists must contain unique, non-empty strings without '/' or missing values.", call. = FALSE)
      }
      paste(x, collapse = "/")
    }, character(1))
  }
  categories <- chars(categories, "categories")
  categorical <- format == "categorical"
  choices <- strsplit(categories, "/", fixed = TRUE)
  bad_choices <- vapply(choices, function(x) {
    !length(x) || any(!nzchar(trimws(x))) || anyDuplicated(x) > 0L
  }, logical(1)) | grepl("/$", categories)
  if (any(categorical & bad_choices) || any(!categorical & nzchar(categories))) {
    stop("Give slash-separated `categories` only for categorical traits, with unique, non-empty choices.", call. = FALSE)
  }
  minimum <- chars(minimum, "minimum")
  maximum <- chars(maximum, "maximum")
  lo <- suppressWarnings(as.numeric(minimum))
  hi <- suppressWarnings(as.numeric(maximum))
  if (any(nzchar(minimum) & !is.finite(lo)) || any(nzchar(maximum) & !is.finite(hi)) ||
      any(lo > hi, na.rm = TRUE)) {
    stop("`minimum` and `maximum` must be finite numbers with minimum <= maximum, or empty.", call. = FALSE)
  }
  is_visible <- expand(is_visible, "is_visible")
  if (!is.logical(is_visible) || anyNA(is_visible)) {
    stop("`is_visible` must be logical without missing values.", call. = FALSE)
  }
  out <- tibble::tibble(
    trait = traits, format = format,
    defaultValue = chars(default_value, "default_value"),
    minimum = minimum, maximum = maximum, details = chars(details, "details"),
    categories = categories, isVisible = ifelse(is_visible, "true", "false"),
    realPosition = seq_len(n)
  )
  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1L ||
        is.na(output_file) || !nzchar(output_file)) {
      stop("`output_file` must be a single non-empty path.", call. = FALSE)
    }
    utils::write.csv(out, output_file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
  out
}
