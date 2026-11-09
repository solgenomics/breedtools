trial_fixture <- function() {
  list(book = data.frame(plots = c(101L, 102L, 103L, 201L, 202L, 203L),
                         block = factor(c(1, 1, 1, 2, 2, 2)),
                         entries = c("G2", "C1", "G1", "G1", "G2", "C1")),
       parameters = list(design = "rcbd", seed = 2026))
}

prepare_trial_fixture <- function(design = trial_fixture(), ...) {
  args <- list(trial_name = "RCBD_2026", year = 2026, location = "Ibadan", program = "IITA",
               description = "RCBD with two complete blocks", prompt = FALSE)
  extra <- list(...)
  for (name in names(extra)) args[name] <- extra[name]
  do.call(prepare_breedbase_trial, c(list(design = design), args))
}

test_that("preparation converts the design and carries the trial information to existing functions", {
  trial <- prepare_trial_fixture(checks = "C1")
  expect_s3_class(trial, "breedbase_trial")
  expect_s3_class(trial, "breedbase_design")
  expect_identical(trial$book, trial_fixture()$book)
  expect_identical(trial$parameters, trial_fixture()$parameters)
  expect_identical(trial$design_type, "RCBD")
  expect_identical(trial$trial$year, 2026L)
  expect_identical(trial$trial$program, "IITA")
  expect_identical(trial$trial$location, "Ibadan")
  expect_null(trial$trial$folder)
  payload <- suppressMessages(save_breedbase(trial, dry_run = TRUE))
  expect_identical(payload$studies[[1]]$studyName, "RCBD_2026")
  expect_identical(payload$studies[[1]]$studyDescription, trial$trial$description)
  expect_identical(payload$studies[[1]]$seasons, list("2026"))
  expect_identical(payload$studies[[1]]$experimentalDesign$PUI, "RCBD")
  sheet <- phenotype_template(trial, traits = "fresh root yield")
  expect_equal(nrow(sheet), 6)
  expect_identical(sheet$plot_name, vapply(payload$plots[[1]], `[[`, "", "observationUnitName"))
  expect_equal(sheet$rep, as.integer(as.character(trial$book$block)))
  expect_equal(sum(sheet$is_check), 2)
  expect_true(all(is.na(sheet$`fresh root yield`)))
})

test_that("prepared defaults can be overridden explicitly when saving", {
  trial <- prepare_trial_fixture()
  payload <- suppressMessages(save_breedbase(trial, trial_name = "Another_2027", year = 2027,
    location = "Ubiaja", description = "Another season", trial_type = "Preliminary Yield Trial",
    dry_run = TRUE))
  expect_identical(payload$studies[[1]]$studyName, "Another_2027")
  expect_identical(payload$studies[[1]]$studyDescription, "Another season")
  expect_identical(payload$studies[[1]]$seasons, list("2027"))
  expect_identical(payload$studies[[1]]$studyType, "Preliminary Yield Trial")
  expect_identical(trial$trial$trial_name, "RCBD_2026")
  sheet <- phenotype_template(trial, trial_name = "Another_2027", location = "Ubiaja")
  expect_identical(sheet$plot_name, vapply(payload$plots[[1]], `[[`, "", "observationUnitName"))
})

test_that("a bare field book and an already converted design can be prepared", {
  raw <- trial_fixture()
  from_book <- prepare_trial_fixture(raw$book, accession = "entries", design_type = "RCBD", checks = "C1")
  converted <- as_breedbase_design(raw, checks = "C1")
  from_converted <- prepare_trial_fixture(converted, folder = "Existing folder")
  expect_identical(from_book$design, from_converted$design)
  expect_identical(from_converted$design, converted$design)
  expect_identical(from_converted$trial$folder, "Existing folder")
  expect_error(prepare_trial_fixture(converted, checks = "C1"), "already converted")
  mapped <- raw$book
  names(mapped)[1] <- "plot_id"
  expect_identical(prepare_trial_fixture(mapped, plot = "plot_id", design_type = "RCBD")$design,
                   prepare_trial_fixture()$design)
})

test_that("missing trial information can be requested through the console", {
  labels <- character()
  answers <- c("PROMPTED_2026", "2026", "Ibadan", "IITA", "RCBD with two complete blocks")
  local_mocked_bindings(bb_prompt_trial_field = function(label) {
    labels <<- c(labels, label)
    answers[length(labels)]
  }, .package = "breedtools")
  trial <- prepare_breedbase_trial(trial_fixture(), prompt = TRUE)
  expect_length(labels, 5)
  expect_identical(trial$trial$trial_name, "PROMPTED_2026")
  expect_identical(trial$trial$year, 2026L)
  expect_identical(trial$trial$program, "IITA")
  expect_identical(trial$trial$location, "Ibadan")
  # Supplying the metadata leaves no fields to ask for.
  prepare_trial_fixture(prompt = TRUE)
  expect_length(labels, 5)
})

test_that("required metadata and invalid metadata fail before any database operation", {
  expect_error(prepare_breedbase_trial(trial_fixture(), prompt = FALSE),
               "trial_name, year, location, program, description")
  expect_error(prepare_trial_fixture(location = NULL), "location")
  expect_error(prepare_trial_fixture(program = " "), "program")
  expect_error(prepare_trial_fixture(description = 1), "description")
  expect_error(prepare_trial_fixture(location = c("A", "B")), "location")
  expect_error(prepare_trial_fixture(folder = ""), "folder")
  expect_error(prepare_trial_fixture(trial_type = ""), "trial_type")
  expect_error(prepare_trial_fixture(prompt = NA), "prompt")
  for (year in list(0, -1, 2026.5, Inf, "not a year", TRUE, c(2025, 2026))) {
    expect_error(prepare_trial_fixture(year = year), "year")
  }
  expect_identical(prepare_trial_fixture(year = "2026")$trial$year, 2026L)
  expect_message(trial <- prepare_trial_fixture(trial_name = "RCBD 2026 (Ibadan)"), "trial name")
  expect_identical(trial$trial$trial_name, "RCBD_2026_Ibadan")
})
