test_that("trait files preserve Field Book's schema and CSV values", {
  path <- tempfile(fileext = ".trt")
  out <- trait_template(
    c("yield", "severity", "notes"),
    format = c("numeric", "categorical", "text"),
    minimum = c(0, NA, NA), maximum = c(100, NA, NA),
    details = c('Yield, in "t/ha"', "Sévérité", "Notes\nsecond line"),
    categories = list(character(), c("low", "medium", "high"), character()),
    is_visible = c(TRUE, TRUE, FALSE), output_file = path
  )
  expect_named(out, c("trait", "format", "defaultValue", "minimum", "maximum",
                      "details", "categories", "isVisible", "realPosition"))
  reread <- read.csv(path, check.names = FALSE, colClasses = "character",
                     fileEncoding = "UTF-8")
  expect_identical(reread$trait, out$trait)
  expect_identical(reread$details, out$details)
  expect_identical(reread$categories, c("", "low/medium/high", ""))
  expect_identical(reread$minimum, c("0", "", ""))
  expect_identical(reread$isVisible, c("true", "true", "false"))
  expect_identical(reread$realPosition, as.character(1:3))
  expect_identical(reread$defaultValue, rep("", 3))
})

test_that("trait definitions recycle scalar settings and reject invalid definitions", {
  out <- trait_template(c("yield", "height"), minimum = 0, default_value = NA)
  expect_identical(out$format, rep("numeric", 2))
  expect_identical(out$minimum, rep("0", 2))
  expect_identical(out$defaultValue, rep("", 2))
  for (traits in list(character(), c("x", "x"), " ", NA_character_)) {
    expect_error(trait_template(traits), "unique, non-empty")
  }
  expect_error(trait_template("x", format = "unknown"), "Unsupported")
  expect_error(trait_template(c("x", "y"), details = letters[1:3]), "length")
  expect_error(trait_template("x", minimum = 10, maximum = 1), "minimum <= maximum")
  expect_error(trait_template("x", minimum = "invalid"), "finite numbers")
  expect_error(trait_template("x", maximum = Inf), "finite numbers")
  expect_error(trait_template("x", is_visible = NA), "logical")
  expect_error(trait_template("x", format = "categorical"), "categories")
  expect_error(trait_template("x", categories = "a/b"), "categories")
  for (categories in c("a/a", "a//b", "a/", "/a")) {
    expect_error(trait_template("x", format = "categorical", categories = categories), "categories")
  }
  expect_error(trait_template("x", format = "categorical", categories = list("a/b")), "Category lists")
  expect_error(trait_template("x", output_file = NA_character_), "output_file")
})

test_that("offline field files keep plot identities through collection and upload preview", {
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                    12, 12, 4, 4, n_candidates = 5, seed = 1)
  field <- phenotype_template(rc, "AYT_2026", location = "Ibadan")
  path <- tempfile(fileext = ".csv")
  write.csv(field, path, row.names = FALSE, na = "")
  observed <- read.csv(path, check.names = FALSE, colClasses = "character")
  expect_identical(observed$plot_name, field$plot_name)
  expect_equal(anyDuplicated(observed$plot_name), 0L)
  observed[["fresh root yield"]] <- seq_len(nrow(observed)) / 10
  observed$collector <- "Extra app metadata"
  observed <- observed[rev(seq_len(nrow(observed))), ]
  obs <- save_breedbase_phenotypes(observed, traits = "fresh root yield", dry_run = TRUE)
  expect_equal(nrow(obs), nrow(field))
  expect_setequal(obs$plot_name, field$plot_name)
  expect_true(all(obs$trait == "fresh root yield"))
  expect_equal(as.numeric(obs$value[match(field$plot_name, obs$plot_name)]),
               seq_len(nrow(field)) / 10)
})
