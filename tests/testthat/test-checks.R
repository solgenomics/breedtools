test_that("check_result, check_summary and check_reset keep a tally", {
  suppressMessages(check_reset())
  expect_message(check_reset(), "Ready to start: no earlier checks")
  expect_message(expect_true(check_result("passes", TRUE, section = "s")), "\\[OK\\] +passes")
  expect_message(expect_false(check_result("fails", FALSE)), "\\[FAIL\\] +fails")
  expect_message(check_result("missing counts as a failure", NA), "\\[FAIL\\]")
  expect_message(check_result("informational", FALSE, note = TRUE), "\\[NOTE\\]")
  log <- NULL
  expect_output(log <- check_summary(), "1 of 3 checks passed \\(plus 1 notes\\)")
  expect_output(check_summary(), "By section")
  expect_output(check_summary(), "Failed checks")
  expect_output(check_summary(details = TRUE), "Every check")
  expect_equal(log$status, c("OK", "FAIL", "FAIL", "NOTE"))
  expect_equal(log$section[1], "s")
  cleared <- NULL
  expect_message(cleared <- check_reset(), "Earlier checks cleared \\(3 checks: 1 passed, 2 failed\\)")
  expect_equal(nrow(cleared), 4)
  expect_equal(nrow(check_summary(print = FALSE)), 0)
  expect_output(check_summary(), "No checks recorded")
  expect_message(check_result("ok", TRUE), "OK")
  expect_output(check_summary(), "All checks are correct")
  suppressMessages(check_reset())
})

test_that("prediction_accuracy correlates adjusted means with the truth", {
  truth <- c(A = 1, B = 2, C = 3, D = 4)
  means <- data.frame(genotype = c("D", "C", "B", "A"), adjusted_mean = c(8, 6, 4, 2))
  expect_equal(prediction_accuracy(means, truth), 1)
  expect_equal(prediction_accuracy(list(adjusted_means = means), truth, method = "spearman"), 1)
  expect_equal(prediction_accuracy(means, -truth), -1)
  expect_equal(prediction_accuracy(means, truth, genotypes = c("A", "B", "C")), 1)
  expect_error(prediction_accuracy(means, unname(truth)), "named")
  expect_error(prediction_accuracy(means, truth, genotypes = c("A", "B")), "Fewer than 3")
})

test_that("check_design runs the standard checks for every design type", {
  suppressMessages(check_reset())
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                    n_candidates = 5, seed = 1)
  prep <- augmented_prep_design(paste0("G", 1:60), c("A", "B"), block_size = 8, prep_prop = 0.3,
                                resolvable = TRUE, seed = 1)
  rrc <- resolvable_row_column_design(paste0("T", 1:20), 3, 4, 5, latinize = TRUE,
                                      n_starts = 1, n_iter = 300, seed = 1)
  expect_message(expect_true(check_design(rc)), "every field row has a check")
  expect_message(expect_true(check_design(prep)), "every location connected")
  expect_message(expect_true(check_design(rrc)), "latinized")
  log <- check_summary(print = FALSE)
  expect_false(any(log$status == "FAIL"))
  expect_true(all(c("augmented row-column design", "augmented p-rep design",
                    "resolvable row-column design") %in% log$section))
  # A broken design is caught.
  broken <- rc
  broken$design$all_entries[broken$design$is_control == 0][1:2] <- "G1"
  expect_message(expect_false(check_design(broken)), "\\[FAIL\\] +every entry planted exactly once")
  expect_error(check_design(data.frame()), "must come from")
  suppressMessages(check_reset())
})

test_that("check_section groups the following checks", {
  suppressMessages(check_reset())
  expect_message(check_section("part A"), "== part A ==")
  expect_message(check_result("first", TRUE), "OK")
  rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 2, seed = 1)
  expect_message(check_design(rc), "Checking the augmented row-column design")
  expect_message(check_section("part B"), "part B")
  expect_message(check_result("second", TRUE), "OK")
  log <- check_summary(print = FALSE)
  expect_equal(unique(log$section), c("part A", "part B"))
  suppressMessages(check_reset())
  expect_message(check_result("no section", TRUE), "OK")
  expect_equal(check_summary(print = FALSE)$section, "")
  suppressMessages(check_reset())
})

test_that("check_setup checks the setup and starts a fresh tally", {
  suppressMessages(check_reset())
  expect_message(check_result("left over from an earlier run", FALSE), "FAIL")
  out_dir <- file.path(tempdir(), "breedtools-setup-test")
  out <- capture_msgs(ok <- check_setup(out_dir))
  expect_true(ok)
  expect_true(dir.exists(out_dir))
  expect_true(any(grepl("breedtools .* is loaded", out)))
  expect_true(any(grepl("files can be saved", out)))
  expect_true(any(grepl("Everything is ready to start a new design", out)))
  log <- check_summary(print = FALSE)
  expect_false("left over from an earlier run" %in% log$check)
  expect_true(all(log$section == "0. Setup"))
  expect_message(check_result("after setup", TRUE), "OK")
  expect_equal(utils::tail(check_summary(print = FALSE)$section, 1), "")
  expect_false(any(log$status == "FAIL"))
  expect_false(any(list.files(out_dir, all.files = TRUE, no.. = TRUE) == ".breedtools_write_test"))
  unlink(out_dir, recursive = TRUE)
  suppressMessages(check_reset())
})
