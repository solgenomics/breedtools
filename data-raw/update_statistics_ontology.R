# Update the Breedbase statistics ontology bundled with breedtools
# ---------------------------------------------------------------------------
# Breedbase stores analysis values under traits composed with a statistic of
# its statistics ontology (SGNSTAT), e.g.
#   "fresh root yield|Adjusted Means from Spatial Correction using SpATS R|COMP:0000379".
# breedtools checks statistics against a copy of the ontology file that
# Breedbase loads (ontology/cxgn_statistics.obo in solgenomics/sgn), kept in
# inst/extdata/cxgn_statistics.obo, so they can be listed and checked offline
# (breedbase_statistics(), statistic = "..." in the save functions).
#
# When Breedbase adds statistics, run this script from the package folder:
#   source("data-raw/update_statistics_ontology.R")
# It downloads the current file, checks that it can be read, lists the terms
# added and removed, and replaces the bundled copy. Then run the tests,
# reinstall the package and commit inst/extdata/cxgn_statistics.obo.

url  <- "https://raw.githubusercontent.com/solgenomics/sgn/master/ontology/cxgn_statistics.obo"
dest <- "inst/extdata/cxgn_statistics.obo"
if (!file.exists(dest)) stop("Run this script from the breedtools package folder.", call. = FALSE)

# The terms of an OBO file: id and name of every [Term].
read_terms <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  starts <- which(lines == "[Term]")
  ends <- c(starts[-1] - 1L, length(lines))
  field <- function(block, key) {
    hit <- grep(paste0("^", key, ":"), block, value = TRUE)
    if (length(hit)) trimws(sub(paste0("^", key, ":"), "", hit[1])) else NA_character_
  }
  terms <- lapply(seq_along(starts), function(k) {
    block <- lines[starts[k]:ends[k]]
    data.frame(id = field(block, "id"), name = field(block, "name"))
  })
  do.call(rbind, terms)
}

new_file <- tempfile(fileext = ".obo")
utils::download.file(url, new_file, quiet = TRUE, mode = "wb")

new <- read_terms(new_file)
old <- read_terms(dest)
if (!nrow(new) || anyNA(new$id) || !all(grepl("^SGNSTAT:[0-9]+$", new$id))) {
  stop("The downloaded file does not look like the SGNSTAT ontology; nothing was changed.", call. = FALSE)
}

added   <- new[!new$id %in% old$id, , drop = FALSE]
removed <- old[!old$id %in% new$id, , drop = FALSE]
renamed <- merge(old, new, by = "id", suffixes = c("_old", "_new"))
renamed <- renamed[renamed$name_old != renamed$name_new, , drop = FALSE]

cat("Statistics ontology: ", nrow(old), " terms bundled, ", nrow(new), " terms in Breedbase now\n", sep = "")
if (!nrow(added) && !nrow(removed) && !nrow(renamed) &&
    identical(readLines(new_file, warn = FALSE), readLines(dest, warn = FALSE))) {
  cat("Already up to date: nothing to do.\n")
} else {
  if (nrow(added))   { cat("\nAdded:\n");   print(added, row.names = FALSE) }
  if (nrow(removed)) { cat("\nRemoved:\n"); print(removed, row.names = FALSE) }
  if (nrow(renamed)) { cat("\nRenamed:\n"); print(renamed, row.names = FALSE) }
  file.copy(new_file, dest, overwrite = TRUE)
  cat("\nUpdated ", dest, " (source: ", url, ", ", format(Sys.Date()), ").\n", sep = "")
  cat("Next: devtools::test(), reinstall, and commit ", dest, ".\n", sep = "")
}
