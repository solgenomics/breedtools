met_rc <- function(seed = 1) {
  sites <- c("North", "Central", "South")
  designs <- lapply(seq_along(sites), function(i) {
    augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4, n_candidates = 3, seed = i)
  })
  met <- do.call(rbind, lapply(seq_along(sites), function(i) {
    d <- designs[[i]]$design
    d$location <- sites[i]
    d
  }))
  set.seed(seed)
  g <- stats::setNames(stats::rnorm(112, sd = 3), c(paste0("G", 1:108), paste0("C", 1:4)))
  met$yield <- 50 + c(North = 0, Central = 5, South = -3)[met$location] + g[met$all_entries] +
    0.6 * met$row + stats::rnorm(nrow(met), sd = ifelse(met$location == "South", 3, 1))
  list(data = met, designs = designs, g = g)
}

test_that("outliers: IQR fences within locations, values set to NA and reported", {
  d <- data.frame(location = rep(c("A", "B"), each = 20), plots = 1:40, all_entries = paste0("G", 1:40),
                  yield = c(50 + (1:20) / 10, 100 + (1:20) / 10))
  d$yield[c(3, 25)] <- c(80, 10)
  found <- find_outliers(d, "yield")
  expect_setequal(found$plots, c(3, 25))
  q <- stats::quantile(d$yield[d$location == "A"], c(0.25, 0.75), names = FALSE)
  expect_equal(found$upper[found$plots == 3], q[2] + 1.5 * diff(q))
  expect_equal(found$side, c("high", "low")[match(found$plots, c(3, 25))])
  # Without grouping, location B's normal values would not be outliers but the low one in A is hidden.
  expect_false(3 %in% find_outliers(d, "yield", by = NULL)$plots)

  out <- capture.output(clean <- remove_outliers(d, "yield"))
  expect_true(all(is.na(clean$yield[c(3, 25)])))
  expect_equal(sum(is.na(clean$yield)), 2)
  expect_true(any(grepl("Removed 2 of 40 values", out)))
  expect_equal(nrow(attr(clean, "outliers")), 2)
  expect_equal(nrow(find_outliers(d, "yield", multiplier = 100)), 0)
  expect_error(find_outliers(d, "yield", multiplier = -1), "positive")
  expect_true(any(grepl("No outliers", capture.output(remove_outliers(d[-c(3, 25), ], "yield")))))
})

test_that("plot_trait draws and summarises per location", {
  m <- met_rc()$data
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  s <- plot_trait(m, "yield", outliers = find_outliers(m, "yield"))
  expect_equal(s$group, c("North", "Central", "South"))
  expect_equal(s$plots, rep(144L, 3))
  expect_silent(plot_trait(m, "yield", type = "histogram", by = NULL))
  expect_error(plot_trait(m, "nope"), "not found")
})

test_that("augmented row-column analyses return the covariance of the adjusted means", {
  m <- met_rc()$data
  d <- m[m$location == "North", ]
  fit <- analyze_augmented_row_column_design(d, "yield")
  vc <- fit$adjusted_means_vcov
  am <- fit$adjusted_means
  expect_equal(unname(sqrt(diag(vc))[am$genotype]), am$std_error, tolerance = 1e-8)
  skip_if_not_installed("SpATS")
  sp <- analyze_augmented_row_column_design(d, "yield", spatial = TRUE)
  # Residuals come from the spatial model, so adjusted mean + residual is the spatially corrected value.
  expect_equal(sp$plot_results$fitted + sp$plot_results$residual, sp$plot_results$yield)
  expect_false(isTRUE(all.equal(sp$plot_results$residual, fit$plot_results$residual)))
  expect_equal(dim(sp$adjusted_means_vcov), c(112L, 112L))
})

test_that("two-stage analysis for augmented row-column, resolvable row-column and p-rep designs", {
  sim <- met_rc()
  t1 <- analyze_two_stage(sim$data, "yield", stage1 = "design")
  expect_s3_class(t1, "two_stage_analysis")
  expect_equal(nrow(t1$locations), 3)
  expect_equal(nrow(t1$adjusted_means), 112)
  expect_gt(cor(t1$adjusted_means$adjusted_mean, sim$g[t1$adjusted_means$genotype]), 0.85)
  expect_true(any(grepl("Two-stage analysis", capture.output(print(t1)))))
  # The noisy South location gets a larger residual variance.
  expect_equal(t1$locations$location[which.max(t1$locations$residual_variance)], "South")

  r <- resolvable_row_column_design(paste0("T", 1:42), n_reps = 3, rows_per_rep = 6, cols_per_rep = 7, seed = 1)
  rm <- do.call(rbind, lapply(c("A", "B"), function(s) {
    d <- r$design
    d$location <- s
    d
  }))
  set.seed(2)
  rm$y <- 20 + rnorm(nrow(rm))
  t2 <- analyze_two_stage(rm, "y", stage1 = "design")
  expect_equal(t2$locations$stage1_model, rep("resolvable row-column (RC)", 2))
  expect_equal(nrow(t2$adjusted_means), 42)

  prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), field_cols = 10, seed = 1)
  pd <- prep$design
  pd$y <- rnorm(nrow(pd), 25, 3)
  t3 <- analyze_two_stage(pd, "y", stage1 = "design")
  ref <- analyze_augmented_prep_two_stage(pd, "y")
  # Same stage 1 and stage 2 as the p-rep two-stage analysis.
  expect_equal(t3$adjusted_means$adjusted_mean, ref$adjusted_means$adjusted_mean, tolerance = 1e-8)

  skip_if_not_installed("SpATS")
  ts <- analyze_two_stage(sim$data, "yield")
  expect_equal(ts$locations$stage1_model, rep("spatial (SpATS)", 3))
  # With a row trend, the spatial stage 1 is at least as accurate as the design model.
  expect_gte(cor(ts$adjusted_means$adjusted_mean, sim$g[ts$adjusted_means$genotype]),
             cor(t1$adjusted_means$adjusted_mean, sim$g[t1$adjusted_means$genotype]) - 0.01)
  v <- first_stage_values(ts$stage1_fits[["North"]], sim$designs[[1]], "AYT_North", trait = "fresh root yield",
                          location = "North")
  expect_equal(nrow(v), 144)
  expect_false(anyNA(v$`fresh root yield`))
  expect_equal(attr(v, "method"), "spatial")
  vp <- first_stage_values(analyze_two_stage(pd, "y"), prep, "PREP", trait = "y")
  expect_equal(nrow(vp), nrow(pd))
  vp2 <- first_stage_values(analyze_augmented_prep_two_stage(pd, "y", stage1 = "spatial"), prep, "PREP", trait = "y")
  expect_false(anyNA(vp2$y))
})

test_that("compare_first_stage puts raw and adjusted values side by side", {
  sim <- met_rc()
  d <- sim$data[sim$data$location == "North", ]
  fit <- analyze_augmented_row_column_design(d, "yield")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- capture.output(cmp <- compare_first_stage(fit))
  expect_equal(nrow(cmp$plots), 144)
  expect_equal(cmp$plots$raw, d$yield[match(cmp$plots$plots, d$plots)])
  expect_equal(cmp$plots$change, cmp$plots$adjusted - cmp$plots$raw)
  expect_true(any(grepl("Raw versus first-stage", out)))
  # The row trend is in the raw values, not in the adjusted ones.
  expect_lt(abs(coef(lm(adjusted ~ row, cmp$plots))[2]), abs(coef(lm(raw ~ row, cmp$plots))[2]))
  t1 <- analyze_two_stage(sim$data, "yield", stage1 = "design")
  out <- capture.output(cmp2 <- compare_first_stage(t1, plot = FALSE))
  expect_equal(cmp2$summary$location, c("North", "Central", "South"))
})

test_that("first-stage comparisons group Breedbase trials with missing or shared location labels", {
  sim <- met_rc()
  d <- sim$data
  d$trial_name <- paste0("TRIAL_", d$location)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (label in c(NA_character_, "Shared_station")) {
    d$location <- label
    fit <- analyze_two_stage(d, "yield", stage1 = "blocks", location_col = "trial_name")
    invisible(capture.output(cmp <- compare_first_stage(fit)))
    expect_equal(cmp$summary$location, unique(d$trial_name))
    expect_equal(cmp$summary$plots, rep(144L, 3))
    expect_false(anyNA(cmp$plots$location))
    expect_true(all(is.finite(cmp$summary$cor_raw_adjusted)))
  }
})
