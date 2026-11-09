agricolae_payload <- function(design, ...) {
  suppressMessages(save_breedbase(design, trial_name = "AG_2026", year = 2026,
    location = "Ibadan", description = "Design from agricolae", dry_run = TRUE, ...))
}

test_that("RCBD conversion keeps accessions, plot numbers and blocks through the save payload", {
  skip_if_not_installed("agricolae")
  withr::local_preserve_seed()
  accessions <- c("G1", "G2", "G3", "C1")
  original <- agricolae::design.rcbd(accessions, r = 3, seed = 2026)
  converted <- as_breedbase_design(original, checks = "C1")
  expect_s3_class(converted, "breedbase_design")
  expect_identical(converted$book, original$book)
  expect_identical(converted$parameters, original$parameters)
  expect_equal(converted$design$all_entries, as.character(original$book$accessions))
  expect_equal(converted$design$plots, original$book$plots)
  payload <- agricolae_payload(converted)
  expect_identical(payload$studies[[1]]$experimentalDesign$PUI, "RCBD")
  expect_identical(payload$studies[[1]]$additionalInfo$design_generated_by, "agricolae")
  expect_identical(payload$studies[[1]]$additionalInfo$agricolae_design, "rcbd")
  sheet <- phenotype_template(converted, "AG_2026", traits = "yield", location = "Ibadan")
  expect_equal(sheet$plot_number, original$book$plots)
  expect_equal(sheet$accession, as.character(original$book$accessions))
  expect_equal(sheet$block, as.integer(as.character(original$book$block)))
  expect_equal(sheet$rep, sheet$block)
  expect_equal(sheet$is_check, sheet$accession == "C1")
  expect_true(all(is.na(sheet$row)))
  expect_true(all(is.na(sheet$col)))
  expect_identical(sheet$plot_name, vapply(payload$plots[[1]], `[[`, "", "observationUnitName"))
  expect_true(all(is.na(sheet$yield)))
})

test_that("CRD, alpha, lattice and augmented blocks use their actual replication columns", {
  skip_if_not_installed("agricolae")
  withr::local_preserve_seed()
  accessions <- paste0("G", 1:12)
  crd <- agricolae::design.crd(accessions[1:3], r = c(1, 2, 3), seed = 1)
  converted <- as_breedbase_design(crd)
  expect_identical(converted$design_type, "CRD")
  expect_true(all(converted$design$block == 1L))
  expect_equal(converted$design$rep, crd$book$r)

  invisible(capture.output(alpha <- agricolae::design.alpha(accessions, k = 3, r = 2, seed = 1)))
  converted <- as_breedbase_design(alpha)
  expect_identical(converted$design_type, "Alpha")
  expect_equal(converted$design$rep, as.integer(as.character(alpha$book$replication)))
  expect_equal(converted$design$block, as.integer(as.character(alpha$book$block)))
  # agricolae's cols are within-block positions, not global field coordinates.
  expect_false("col" %in% names(converted$design))
  expect_length(agricolae_payload(converted)$plots[[1]], 24)

  invisible(capture.output(lattice <- agricolae::design.lattice(accessions[1:9], r = 2, seed = 1)))
  converted <- as_breedbase_design(lattice)
  expect_identical(converted$design_type, "Lattice")
  expect_equal(converted$design$rep, as.integer(as.character(lattice$book$r)))

  dau <- agricolae::design.dau(c("C1", "C2"), accessions, r = 3, seed = 1)
  converted <- as_breedbase_design(dau)
  expect_identical(converted$design_type, "Augmented")
  expect_equal(sum(converted$design$is_control), 6)
  expect_true(all(converted$design$rep == 1L))
  expect_equal(sum(as_breedbase_design(dau, checks = character(0))$design$is_control), 0)
})

test_that("bare and serpentine field books retain the plot to accession mapping", {
  skip_if_not_installed("agricolae")
  withr::local_preserve_seed()
  accessions <- paste0("G", 1:4)
  original <- agricolae::design.rcbd(accessions, r = 3, seed = 1)
  book <- agricolae::zigzag(original)
  expect_error(as_breedbase_design(book), "Give `design_type`")
  converted <- as_breedbase_design(book, design_type = "RCBD")
  own_names <- paste0("CUSTOM_", seq_len(nrow(book)))
  payload <- agricolae_payload(converted, plot_names = own_names)
  ordered <- order(book$plots)
  expect_identical(vapply(payload$plots[[1]], `[[`, "", "germplasmName"), as.character(book$accessions[ordered]))
  expect_identical(vapply(payload$plots[[1]], `[[`, "", "observationUnitName"), own_names[ordered])
  expect_equal(converted$design$rep, as.integer(as.character(book$block)))
})

test_that("custom column mapping works offline without an agricolae dependency", {
  book <- data.frame(id = c(102, 101), genotype = c("G2", "G1"),
    iblock = factor(c(7, 8)), repetition = factor(c(2, 1)),
    north = c(1, 1), east = c(2, 1))
  converted <- as_breedbase_design(book, design_type = "Alpha", plot = "id",
    block = "iblock", rep = "repetition", row = "north", col = "east")
  expect_equal(converted$design$all_entries, book$genotype)
  expect_equal(converted$design$block, c(7L, 8L))
  expect_equal(converted$design$rep, c(2L, 1L))
  payload <- agricolae_payload(converted)
  expect_equal(payload$plots[[1]][[1]]$observationUnitPosition$positionCoordinateX, "1")
  expect_equal(payload$plots[[1]][[2]]$observationUnitPosition$positionCoordinateX, "2")
})

test_that("multiple factors need an explicit accession column and split plots stay distinct", {
  skip_if_not_installed("agricolae")
  withr::local_preserve_seed()
  genotype <- paste0("G", 1:4)
  irrigation <- c("low", "high")
  original <- agricolae::design.split(genotype, irrigation, r = 3, seed = 1)
  expect_error(as_breedbase_design(original), "Give `accession`")
  expect_message(converted <- as_breedbase_design(original, accession = "genotype"), "unique plot numbers")
  expect_identical(converted$design_type, "splitplot")
  expect_equal(converted$design$plots, seq_len(nrow(original$book)))
  expect_identical(converted$book, original$book)
  sheet <- suppressMessages(phenotype_template(converted, "SPLIT_2026", location = "Ibadan"))
  expect_equal(nrow(sheet), 24)
  expect_equal(anyDuplicated(sheet$plot_name), 0)
  expect_equal(sheet$accession, as.character(original$book$genotype))
  expect_equal(sheet$rep, as.integer(as.character(original$book$block)))

  strip <- agricolae::design.strip(genotype, irrigation, r = 3, seed = 1)
  converted <- as_breedbase_design(strip, accession = "genotype")
  expect_identical(converted$design_type, "stripplot")
  expect_equal(converted$design$rep, as.integer(as.character(strip$book$block)))
  factorial <- agricolae::design.ab(c(2, 3), r = 3, seed = 1)
  expect_error(as_breedbase_design(factorial), "Give `accession`")
  expect_identical(as_breedbase_design(factorial, accession = "A")$design_type, "RCBD")
})

test_that("designs without a direct Breedbase type require a choice and keep coordinates", {
  skip_if_not_installed("agricolae")
  withr::local_preserve_seed()
  accessions <- paste0("G", 1:4)
  square <- agricolae::design.lsd(accessions, seed = 1)
  expect_error(as_breedbase_design(square), "no automatic Breedbase design type.*lsd")
  # Explicit metadata override; the source remains a Latin square.
  converted <- as_breedbase_design(square, design_type = "RCBD")
  expect_identical(converted$parameters$design, "lsd")
  expect_equal(converted$design$row, as.integer(as.character(square$book$row)))
  expect_equal(converted$design$col, as.integer(as.character(square$book$col)))
  expect_true(all(converted$design$rep == 1L))
  expect_length(agricolae_payload(converted)$plots[[1]], 16)
  expect_identical(agricolae_payload(converted)$studies[[1]]$additionalInfo$agricolae_design, "lsd")
})

test_that("malformed field books fail before any upload", {
  book <- data.frame(plots = c(101, 102), block = c(1, 1), genotype = c("G1", "G2"))
  convert <- function(x, ...) as_breedbase_design(x, design_type = "RCBD", ...)
  expect_error(convert(NULL), "non-empty")
  expect_error(convert(book[FALSE, ]), "non-empty")
  expect_error(convert(book, accession = "missing"), "must name one column")
  expect_error(convert(book, checks = "absent"), "absent from the field book")
  expect_error(convert(book, checks = NA), "missing or empty")
  expect_error(as_breedbase_design(book, design_type = "lsd"), "types Breedbase accepts")
  bad <- book
  bad$plots[2] <- bad$plots[1]
  expect_error(convert(bad), "must be unique")
  for (value in c(NA, Inf, -1, 0, 1.5)) {
    bad <- book
    bad$plots[1] <- value
    expect_error(convert(bad), "positive integers")
  }
  bad <- book
  bad$genotype[1] <- " "
  expect_error(convert(bad), "missing or empty")
  bad$genotype[1] <- NA
  expect_error(convert(bad), "missing or empty")
  bad <- book
  names(bad)[3] <- "plots"
  expect_error(convert(bad), "column names")
  bad <- book
  bad$row <- 1:2
  expect_error(convert(bad), "both `row` and `col`")
  bad$col <- c(1, NA)
  expect_error(convert(bad), "positive integers")
  bad$row <- bad$col <- c(1, 1)
  expect_error(convert(bad), "unique position")
})
