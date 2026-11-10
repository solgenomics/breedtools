test_that("augmented_row_column_capacity calculates plot counts", {
  capacity <- augmented_row_column_capacity(
    rows_in_field = 6,
    cols_in_field = 6,
    rows_per_block = 3,
    cols_per_block = 3,
    n_controls = 2
  )

  expect_equal(capacity$n_blocks, 4)
  expect_equal(capacity$total_plots, 36)
  expect_equal(capacity$check_plots, 8)
  expect_equal(capacity$entry_plots, 28)
})

test_that("augmented_row_column_design returns complete and breedbase outputs", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 10,
    seed = 123
  )

  expect_s3_class(design, "augmented_row_column_design")
  expect_equal(nrow(design$design), 144)
  expect_equal(nrow(design$breedbase_design), 144)
  expect_equal(
    names(design$design),
    c(
      "plots",
      "row",
      "col",
      "block",
      "rowgroup",
      "colgroup",
      "all_entries",
      "type",
      "rep",
      "is_control"
    )
  )
  expect_equal(
    names(design$breedbase_design),
    c("plots", "block", "all_entries", "rep", "is_control")
  )
  expect_equal(sum(design$design$is_control), 36)
  expect_equal(sum(design$design$type == "entry"), 108)
  expect_setequal(design$input$controls, paste0("Check", seq_len(4)))
  expect_setequal(design$input$treatments, paste0("G", seq_len(108)))
  expect_setequal(design$design$all_entries[design$design$is_control == 0], paste0("G", seq_len(108)))
  expect_equal(design$summary$error_df, 6)
})

test_that("augmented_row_column_design prevents duplicate check rows and columns within blocks", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 10,
    seed = 321
  )

  checks <- design$design[design$design$is_control == 1, ]

  by_block <- split(checks, checks$block)
  for (block_checks in by_block) {
    expect_equal(anyDuplicated(block_checks$row), 0)
    expect_equal(anyDuplicated(block_checks$col), 0)
    expect_equal(anyDuplicated(block_checks$all_entries), 0)
  }
})

test_that("augmented_row_column_design handles large layouts without recursive search", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(576)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 36,
    cols_in_field = 18,
    rows_per_block = 6,
    cols_per_block = 6,
    n_candidates = 3,
    seed = 123
  )

  checks <- design$design[design$design$is_control == 1, ]

  expect_equal(nrow(design$design), 648)
  expect_equal(nrow(checks), 72)
  expect_setequal(checks$row, seq_len(36))
  expect_setequal(checks$col, seq_len(18))
  expect_false(any(table(checks$all_entries, checks$row) > 1))
  expect_false(any(table(checks$all_entries, checks$col) > 1))

  by_block <- split(checks, checks$block)
  for (block_checks in by_block) {
    expect_equal(anyDuplicated(block_checks$row), 0)
    expect_equal(anyDuplicated(block_checks$col), 0)
    expect_equal(anyDuplicated(block_checks$all_entries), 0)
  }
})


test_that("augmented_row_column_design randomizes check positions across seeds", {
  check_signature <- function(seed) {
    design <- augmented_row_column_design(
      treatments = paste0("G", seq_len(108)),
      controls = paste0("Check", seq_len(4)),
      rows_in_field = 12,
      cols_in_field = 12,
      rows_per_block = 4,
      cols_per_block = 4,
      n_candidates = 1,
      seed = seed
    )
    checks <- design$design[design$design$is_control == 1, ]
    paste(sort(paste(checks$all_entries, checks$row, checks$col, sep = ":")), collapse = "|")
  }

  signatures <- vapply(c(101, 202, 303), check_signature, character(1))

  expect_gt(length(unique(signatures)), 1)
})

test_that("plot_augmented_row_column_design draws design objects and data frames", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 3,
    seed = 111
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  plotted_object <- plot_augmented_row_column_design(
    design,
    label = "none",
    fill = "block",
    legend = FALSE
  )
  plotted_data <- plot_augmented_row_column_design(
    design$design,
    label = "all_entries",
    fill = "type",
    legend = FALSE
  )
  grDevices::dev.off()

  expect_s3_class(plotted_object, "tbl_df")
  expect_s3_class(plotted_data, "tbl_df")
  expect_true(all(c(".plot_row", ".plot_col", ".fill_value") %in% names(plotted_object)))
  expect_true(file.exists(plot_file))
})

test_that("augmented_row_column_design validates treatment count", {
  expect_error(
    augmented_row_column_design(
      treatments = paste0("G", seq_len(27)),
      controls = c("Check1", "Check2"),
      rows_in_field = 6,
      cols_in_field = 6,
      rows_per_block = 3,
      cols_per_block = 3,
      n_candidates = 5
    ),
    "Expected 28 treatments"
  )
})

test_that("analyze_augmented_row_column_design returns model outputs", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 123
  )

  design_data <- design$design
  design_data$yield <- 100 +
    as.numeric(factor(design_data$rowgroup)) * 2 +
    as.numeric(factor(design_data$colgroup)) * 3 +
    ifelse(design_data$is_control == 1, 5, 0) +
    stats::rnorm(nrow(design_data), sd = 1)

  analysis <- analyze_augmented_row_column_design(
    design_data,
    response = "yield"
  )

  expect_s3_class(analysis, "augmented_row_column_analysis")
  expect_s3_class(analysis$model, "lm")
  expect_false(any(is.na(stats::coef(analysis$model))))
  expect_equal(nrow(analysis$plot_results), nrow(design_data))
  expect_equal(nrow(analysis$genotype_summary), dplyr::n_distinct(design_data$all_entries))
  expect_true(all(c("rowgroup:colgroup", "rowgroup:row", "colgroup:col", "genotype") %in% analysis$anova$term))
  expect_equal(
    analysis$anova_formatted$Source,
    c(
      "Row groups",
      "Rows, nested within row groups",
      "Column groups",
      "Columns, nested within column groups",
      "Row groups x column groups (blocks)",
      "Genotypes",
      "Error",
      "Corrected total"
    )
  )
  expect_equal(analysis$anova_formatted$Df, c(2, 9, 2, 9, 4, 111, 6, 143))
  expect_true(all(analysis$anova_formatted$Significance %in% c("", "*", "**")))
  expect_true(analysis$diagnostics$residual_df > 0)
  expect_equal(analysis$diagnostics$anova_error_df, 6)
})

test_that("augmented row-column analysis print method shows formatted ANOVA", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 3,
    seed = 456
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_output(print(analysis), "Augmented row-column design analysis")
  expect_output(print(analysis), "Corrected total")
})

test_that("analyze_augmented_row_column_design validates required columns", {
  expect_error(
    analyze_augmented_row_column_design(
      data.frame(rowgroup = 1, colgroup = 1, row = 1, col = 1, all_entries = "G1"),
      response = "yield"
    ),
    "Missing required column"
  )
})

test_that("analyze_augmented_row_column_design flags a perfect fit", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 123
  )

  design_data <- design$design
  design_data$phenotype <- design_data$row * 2 + design_data$col

  expect_warning(
    analysis <- analyze_augmented_row_column_design(design_data, response = "phenotype"),
    NA
  )

  expect_true(analysis$diagnostics$saturated_or_near_perfect_fit)
  expect_false(any(is.nan(unlist(analysis$anova[vapply(analysis$anova, is.numeric, logical(1))]))))
  expect_false(any(is.nan(unlist(analysis$anova_formatted[vapply(
    analysis$anova_formatted,
    is.numeric,
    logical(1)
  )]))))
  expect_match(analysis$message, "essentially perfect fit")
})

test_that("analyze_augmented_row_column_design reports when no error df remain", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 123
  )

  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  # Keep only the checks of one block, so no genotype is replicated.
  design_data$yield[design_data$is_control == 1 & design_data$block != 1] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_equal(analysis$diagnostics$anova_error_df, 0)
  expect_match(analysis$message, "no positive error degrees of freedom")
})

test_that("augmented ANOVA df come from the observed plots when values are missing", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 3,
    seed = 789
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  entry_plots <- which(design_data$is_control == 0)[1:5]
  design_data$yield[entry_plots] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_equal(unique(analysis$anova_df$k), 12)
  expect_equal(unique(analysis$anova_df$s), 12)
  expect_equal(analysis$diagnostics$design_error_df, 6)
  # Missing entries remove genotype df only; error df is unchanged.
  expect_equal(analysis$anova_formatted$Df, c(2, 9, 2, 9, 4, 106, 6, 138))
})

test_that("augmented ANOVA df follow the design formulas", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(320)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 24,
    cols_in_field = 16,
    rows_per_block = 6,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 987
  )
  design_data <- design$design
  design_data$phenotype <- stats::rnorm(nrow(design_data))

  analysis <- analyze_augmented_row_column_design(design_data, response = "phenotype")

  expect_equal(unique(analysis$anova_df$k), 24)
  expect_equal(unique(analysis$anova_df$s), 16)
  expect_equal(unique(analysis$anova_df$gk), 4)
  expect_equal(unique(analysis$anova_df$gs), 4)
  expect_equal(unique(analysis$anova_df$vc), 4)
  expect_equal(unique(analysis$anova_df$ve), 320)
  # rows within groups = k - gk = 20; columns within groups = s - gs = 12
  expected_df <- c(3, 20, 3, 12, 9, 323, 13, 383)
  expect_equal(analysis$anova_formatted$Df, expected_df)
  expect_equal(analysis$anova_df$display_df, expected_df[1:7])
})

test_that("augmented ANOVA matches a full lm() fit", {
  design <- augmented_row_column_design(
    treatments = paste0("G", seq_len(108)),
    controls = paste0("Check", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 42
  )
  d <- design$design
  d$y <- stats::rnorm(nrow(d)) + d$row * 0.3
  d$y[c(3, 50, 77)] <- NA_real_

  analysis <- analyze_augmented_row_column_design(d, response = "y")

  ref_data <- data.frame(
    y = d$y,
    rg = factor(d$rowgroup),
    rn = interaction(d$rowgroup, d$row, drop = TRUE),
    cg = factor(d$colgroup),
    cn = interaction(d$colgroup, d$col, drop = TRUE),
    bk = interaction(d$rowgroup, d$colgroup, drop = TRUE),
    g = factor(d$all_entries)
  )
  ref <- stats::anova(stats::lm(y ~ rg + rn + cg + cn + bk + g, data = ref_data))

  expect_equal(analysis$anova$df, ref$Df)
  expect_equal(analysis$anova$sum_sq, ref$`Sum Sq`, tolerance = 1e-8)
  expect_equal(analysis$anova$p_value[1:6], ref$`Pr(>F)`[1:6], tolerance = 1e-6)
})

test_that("augmented_row_column_design rejects layouts that cannot be analyzed", {
  expect_error(
    augmented_row_column_design(
      treatments = paste0("G", seq_len(28)),
      controls = c("Check1", "Check2"),
      rows_in_field = 6,
      cols_in_field = 6,
      rows_per_block = 3,
      cols_per_block = 3,
      n_candidates = 5
    ),
    "fewer check plots"
  )

  capacity <- augmented_row_column_capacity(10, 10, 5, 5, 2)
  expect_error(
    augmented_row_column_design(
      treatments = paste0("G", seq_len(capacity$entry_plots)),
      controls = c("Check1", "Check2"),
      rows_in_field = 10,
      cols_in_field = 10,
      rows_per_block = 5,
      cols_per_block = 5,
      n_candidates = 5
    ),
    "Not every field row can contain a check"
  )
})

test_that("augmented_row_column_capacity reports error df", {
  expect_equal(augmented_row_column_capacity(12, 12, 4, 4, 4)$error_df, 6)
  expect_equal(augmented_row_column_capacity(6, 6, 3, 3, 2)$error_df, -5)
})

test_that("augmented row-column adjusted means match a full lm fit", {
  design <- augmented_row_column_design(
    treatments = paste0("E", seq_len(108)),
    controls = paste0("C", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 9
  )
  set.seed(3)
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data), 50, 5) + 0.5 * design_data$row
  design_data$yield[c(3, 40, 77)] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  observed <- design_data[!is.na(design_data$yield), ]
  observed$genotype <- factor(observed$all_entries)
  observed$rowgroup <- factor(observed$rowgroup)
  observed$colgroup <- factor(observed$colgroup)
  observed$row_nested <- interaction(observed$rowgroup, observed$row, drop = TRUE)
  observed$col_nested <- interaction(observed$colgroup, observed$col, drop = TRUE)
  fit <- stats::lm(
    yield ~ rowgroup + row_nested + colgroup + col_nested + rowgroup:colgroup + genotype,
    data = observed
  )
  ref <- reference_adjusted_means(fit, observed)
  got <- analysis$adjusted_means[match(ref$genotype, analysis$adjusted_means$genotype), ]

  expect_equal(got$adjusted_mean, ref$mean, tolerance = 1e-8)
  expect_equal(got$std_error, ref$se, tolerance = 1e-8)
  expect_equal(nrow(analysis$adjusted_means), 112)
  expect_true(all(is.na(analysis$adjusted_means$adjusted_mean[analysis$adjusted_means$n_observations == 0])))
  expect_false(is.unsorted(-stats::na.omit(analysis$adjusted_means$adjusted_mean)))
  expect_true(all(c("adjusted_mean", "std_error") %in% names(analysis$genotype_summary)))
  expect_gt(analysis$diagnostics$mean_sed, 0)
  expect_equal(
    analysis$diagnostics$lsd_5pct,
    stats::qt(0.975, analysis$diagnostics$anova_error_df) * analysis$diagnostics$mean_sed
  )
  expect_output(print(analysis), "Top adjusted means")
})

test_that("augmented row-column adjusted means are NA when not estimable", {
  design <- augmented_row_column_design(
    treatments = paste0("E", seq_len(108)),
    controls = paste0("C", seq_len(4)),
    rows_in_field = 12,
    cols_in_field = 12,
    rows_per_block = 4,
    cols_per_block = 4,
    n_candidates = 5,
    seed = 9
  )
  design_data <- design$design
  design_data$yield <- stats::rnorm(nrow(design_data))
  design_data$yield[design_data$is_control == 1 & design_data$row == 1] <- NA_real_

  analysis <- analyze_augmented_row_column_design(design_data, response = "yield")

  expect_true(analysis$diagnostics$rank_deficient)
  expect_true(all(is.na(analysis$adjusted_means$adjusted_mean)))
  expect_match(analysis$message, "not estimable")
})

# --- Piepho and Williams (2016) options --------------------------------------

pw_example3 <- function(seed = 1) {
  augmented_row_column_design(
    treatments = paste0("E", 1:173),
    controls = c("A", "B", "C", "D"),
    rows_per_rowgroup = c(6, 5, 6, 5),
    cols_per_colgroup = c(3, 4, 3),
    checks_per_block = rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3)),
    n_candidates = 3,
    seed = seed
  )
}

test_that("capacity reproduces Table 2 of Piepho and Williams (2016)", {
  ex1 <- augmented_row_column_capacity(6, 9, 3, 3, n_controls = 3)
  expect_equal(ex1$entry_plots, 36)
  expect_equal(ex1$error_df, 0)
  expect_equal(ex1$min_rowgroup_slack, 2)   # 7 <= 9
  ex2 <- augmented_row_column_capacity(36, 18, 6, 6, n_controls = 4)
  expect_equal(ex2$entry_plots, 576)
  expect_equal(ex2$error_df, 6)
  expect_equal(ex2$min_rowgroup_slack, 1)   # 11 <= 12
  expect_equal(ex2$min_colgroup_slack, 10)  # 14 <= 24
  ex3 <- augmented_row_column_capacity(
    n_controls = 4,
    rows_per_rowgroup = c(6, 5, 6, 5),
    cols_per_colgroup = c(3, 4, 3),
    checks_per_block = rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
  )
  expect_equal(ex3$total_plots, 220)
  expect_equal(ex3$check_plots, 47)
  expect_equal(ex3$entry_plots, 173)
  expect_equal(ex3$error_df, 7)
  expect_true(is.na(ex3$rows_per_block))
})

test_that("Example 1 of the paper (Fig. 2) is built without warnings despite zero error df", {
  # 6 x 9 field, 3 x 3 blocks, each check once per block, 36 entries.
  expect_no_warning(
    fig2 <- augmented_row_column_design(paste0("E", 1:36), LETTERS[1:3], 6, 9, 3, 3,
                                        method = "latin", seed = 1)
  )
  expect_equal(fig2$summary$error_df, 0)
  ck <- fig2$design[fig2$design$is_control == 1, ]
  expect_true(all(table(ck$block, ck$all_entries) == 1))
  expect_true(all(table(ck$row) == 3))        # one check per column group
  expect_true(all(table(ck$row, ck$all_entries) == 1))
  # Analyzable with random rows, columns and blocks.
  skip_if_not_installed("lme4")
  d <- fig2$design
  set.seed(2)
  d$y <- rnorm(nrow(d), 10)
  fit <- suppressWarnings(analyze_augmented_row_column_design(d, "y", random_blocks = TRUE))
  expect_false(anyNA(fit$adjusted_means$std_error))
})

test_that("Example 3 of the paper: unequal groups and blocks with 3 to 5 checks", {
  design <- pw_example3()
  d <- design$design
  checks <- d[d$is_control == 1, ]
  expect_equal(design$summary$method, "optimal")
  expect_equal(nrow(d), 220)
  expect_equal(nrow(checks), 47)
  expect_equal(sum(d$type == "entry"), 173)
  expect_equal(as.integer(table(d$rowgroup) / 10), c(6L, 5L, 6L, 5L))
  expect_equal(as.integer(table(d$colgroup) / 22), c(3L, 4L, 3L))
  pb <- rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
  expect_equal(as.vector(table(factor(checks$block, levels = 1:12))), as.vector(t(pb)))
  expect_setequal(checks$row, 1:22)
  expect_setequal(checks$col, 1:10)
  expect_equal(anyDuplicated(paste(d$row, d$col)), 0)
  expect_gt(design$summary$d_criterion, 0)

  d$y <- stats::rnorm(nrow(d))
  analysis <- analyze_augmented_row_column_design(d, "y")
  expect_equal(analysis$diagnostics$residual_df, 7)
  expect_false(anyNA(analysis$adjusted_means$adjusted_mean))
})

test_that("incomplete blocks and extra check replication are supported", {
  incomplete <- augmented_row_column_design(
    paste0("E", 1:279), paste0("C", 1:6), 18, 18, 6, 6,
    checks_per_block = 5, n_candidates = 2, seed = 1
  )
  ck <- incomplete$design[incomplete$design$is_control == 1, ]
  expect_true(all(table(ck$block) == 5))
  expect_true(all(table(ck$block, ck$all_entries) <= 1))
  expect_equal(incomplete$summary$error_df, 1)

  extra <- augmented_row_column_design(
    paste0("E", 1:90), LETTERS[1:4], 12, 12, 4, 4,
    checks_per_block = 6, n_candidates = 2, seed = 1
  )
  ck <- extra$design[extra$design$is_control == 1, ]
  expect_true(all(table(ck$block) == 6))
  expect_true(all(table(ck$block, ck$all_entries) >= 1))
  expect_equal(extra$summary$error_df, 24)
})

test_that("method optimal can be forced and meets the same placement rules as latin", {
  latin <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                       n_candidates = 50, seed = 2)
  optimal <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                         method = "optimal", n_candidates = 3, seed = 2)
  expect_equal(latin$summary$method, "latin")
  expect_equal(optimal$summary$method, "optimal")
  # Both methods meet the placement rules; the optimal search trades a little
  # D-efficiency for them and is never far behind the Latin construction.
  expect_gte(optimal$summary$d_criterion, 0.9 * latin$summary$d_criterion)
  for (d in list(latin, optimal)) {
    expect_equal(d$summary$adjacent_checks_in_blocks, 0)
    suppressMessages(check_reset())
    expect_message(expect_true(check_design(d)), "even numbers of check plots")
  }
  suppressMessages(check_reset())
  expect_error(
    augmented_row_column_design(paste0("G", 1:117), paste0("C", 1:4), 12, 12, 4, 4,
                                checks_per_block = 3, method = "latin"),
    "latin"
  )
})

test_that("checks_per_block is validated", {
  expect_error(
    augmented_row_column_capacity(12, 12, 4, 4, 4, checks_per_block = matrix(4, 2, 2)),
    "3 x 3 matrix"
  )
  expect_error(augmented_row_column_capacity(12, 12, 4, 4, 4, checks_per_block = 17), "larger than")
  expect_error(
    augmented_row_column_capacity(n_controls = 4, rows_per_rowgroup = c(6, 5), cols_per_colgroup = 3,
                                  rows_in_field = 12),
    "add up"
  )
})

rc_mixed_trial <- function() {
  design <- augmented_row_column_design(paste0("E", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                        n_candidates = 20, seed = 1)
  x <- design$design
  set.seed(3)
  g <- stats::setNames(stats::rnorm(112, sd = 2), unique(x$all_entries))
  rn <- interaction(x$rowgroup, x$row)
  cn <- interaction(x$colgroup, x$col)
  x$y <- 20 + g[x$all_entries] + stats::rnorm(nlevels(rn))[as.integer(rn)] +
    stats::rnorm(nlevels(cn))[as.integer(cn)] + 2 * sin(x$row / 3) + stats::rnorm(nrow(x))
  x$rn <- rn
  x$cn <- cn
  list(data = x, truth = g)
}

test_that("random blocks recover inter-block information and match lme4", {
  skip_if_not_installed("lme4")
  trial <- rc_mixed_trial()
  x <- trial$data
  fit <- analyze_augmented_row_column_design(x, "y", random_blocks = TRUE)
  expect_equal(fit$diagnostics$model_type, "mixed")
  expect_true(inherits(fit$mixed_model, "lmerMod"))
  expect_true("rows within row groups" %in% fit$variance_components$component)

  x$gen <- factor(x$all_entries)
  x$bl <- factor(x$block)
  ref <- lme4::lmer(y ~ factor(rowgroup) + factor(colgroup) + gen + (1 | rn) + (1 | cn) + (1 | bl), data = x)
  fe <- lme4::fixef(ref)
  lv <- levels(x$gen)
  got <- fit$adjusted_means$adjusted_mean[match(lv, fit$adjusted_means$genotype)]
  expect_equal(got - got[1], unname(c(0, fe[paste0("gen", lv[-1])])), tolerance = 1e-4)
  expect_equal(
    sort(fit$variance_components$variance),
    sort(as.data.frame(lme4::VarCorr(ref))$vcov),
    tolerance = 1e-4
  )
  # The ANOVA is still the fixed-effects table.
  fixed <- analyze_augmented_row_column_design(x, "y")
  expect_equal(fit$anova_formatted, fixed$anova_formatted)
})

test_that("random entries give BLUPs that match lme4", {
  skip_if_not_installed("lme4")
  trial <- rc_mixed_trial()
  x <- trial$data
  fit <- analyze_augmented_row_column_design(x, "y", random_blocks = TRUE, entries = "random")
  expect_match(fit$message, "BLUP")
  expect_true("entries" %in% fit$variance_components$component)
  ent <- x$all_entries[x$is_control == 0]

  x$conv <- as.numeric(x$is_control == 0)
  x$chk <- factor(ifelse(x$is_control == 1, x$all_entries, ".e"))
  x$ent <- factor(ifelse(x$is_control == 1, ".c", x$all_entries))
  x$bl <- factor(x$block)
  ref <- lme4::lmer(y ~ factor(rowgroup) + factor(colgroup) + chk + (1 | rn) + (1 | cn) + (1 | bl) +
                      (0 + conv | ent), data = x)
  u <- lme4::ranef(ref)$ent[ent, 1]
  ours <- fit$adjusted_means$adjusted_mean[match(ent, fit$adjusted_means$genotype)]
  expect_equal(ours - mean(ours), u - mean(u), tolerance = 1e-4)
  expect_true(all(fit$adjusted_means$std_error > 0))

  fixed_blocks <- analyze_augmented_row_column_design(x, "y", entries = "random")
  expect_equal(fixed_blocks$diagnostics$model_type, "mixed")
})

test_that("spatial analysis with SpATS runs and improves on a field trend", {
  skip_if_not_installed("SpATS")
  trial <- rc_mixed_trial()
  fixed <- analyze_augmented_row_column_design(trial$data, "y")
  spatial <- analyze_augmented_row_column_design(trial$data, "y", spatial = TRUE)
  expect_equal(spatial$diagnostics$model_type, "spatial")
  expect_gt(spatial$diagnostics$spatial_effective_dim, 0)
  ent <- trial$data$all_entries[trial$data$is_control == 0]
  r <- function(f) stats::cor(f$adjusted_means$adjusted_mean[match(ent, f$adjusted_means$genotype)], trial$truth[ent])
  expect_gt(r(spatial), r(fixed))
  expect_error(
    analyze_augmented_row_column_design(trial$data, "y", spatial = TRUE, entries = "random"),
    "only available"
  )
})

test_that("augmented_row_column_layouts lists the analyzable layouts", {
  layouts <- augmented_row_column_layouts(12, 12, n_controls = 2:5)
  expect_s3_class(layouts, "tbl_df")
  expect_true(all(layouts$error_df >= 1))
  expect_true(all(diff(layouts$entry_plots) <= 0))
  expect_true(all(layouts$rows_per_block >= layouts$n_controls))
  expect_true(all(layouts$cols_per_block >= layouts$n_controls))
  # Every listed layout agrees with augmented_row_column_capacity().
  for (i in seq_len(nrow(layouts))) {
    cap <- augmented_row_column_capacity(12, 12, layouts$rows_per_block[i],
                                         layouts$cols_per_block[i], layouts$n_controls[i])
    expect_equal(layouts$entry_plots[i], cap$entry_plots)
    expect_equal(layouts$error_df[i], cap$error_df)
  }
  row <- layouts[layouts$rows_per_block == 4 & layouts$cols_per_block == 4 & layouts$n_controls == 4, ]
  expect_equal(row$entry_plots, 108)
  expect_equal(row$method, "latin")

  fit100 <- augmented_row_column_layouts(12, 12, n_controls = 2:6, min_error_df = 6, n_entries = 100)
  expect_true(all(fit100$entry_plots >= 100 & fit100$error_df >= 6))
  expect_equal(fit100$spare_plots, fit100$entry_plots - 100)

  # Each listed layout can actually be generated.
  pick <- fit100[1, ]
  design <- augmented_row_column_design(paste0("G", seq_len(pick$entry_plots)), paste0("C", seq_len(pick$n_controls)),
                                        12, 12, pick$rows_per_block, pick$cols_per_block, n_candidates = 5, seed = 1)
  expect_equal(design$summary$error_df, pick$error_df)

  expect_message(none <- augmented_row_column_layouts(6, 6, n_controls = 1, min_error_df = 50), "No analyzable layout")
  expect_equal(nrow(none), 0)
  for (dims in list(c(3, 6), c(6, 3))) {
    expect_message(narrow <- augmented_row_column_layouts(12, 12,
      rows_per_block = dims[1], cols_per_block = dims[2], n_controls = 4), "No analyzable layout")
    expect_equal(nrow(narrow), 0)
  }
})

test_that("checks are never placed next to each other inside a block", {
  neighbours <- function(design) {
    ck <- design$design[design$design$is_control == 1, ]
    key <- paste(ck$row, ck$col, ck$block)
    sum(paste(ck$row + 1, ck$col, ck$block) %in% key) + sum(paste(ck$row, ck$col + 1, ck$block) %in% key)
  }
  latin <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4), 12, 12, 4, 4,
                                       n_candidates = 50, seed = 2026)
  pw <- augmented_row_column_design(paste0("E", 1:173), LETTERS[1:4], rows_per_rowgroup = c(6, 5, 6, 5),
                                    cols_per_colgroup = c(3, 4, 3),
                                    checks_per_block = rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3)),
                                    n_candidates = 2, seed = 1)
  dense <- augmented_row_column_design(paste0("E", 1:90), LETTERS[1:4], 12, 12, 4, 4,
                                       checks_per_block = 6, n_candidates = 2, seed = 1)
  for (d in list(latin, pw, dense)) {
    expect_equal(neighbours(d), 0)
    expect_equal(d$summary$adjacent_checks_in_blocks, 0)
  }
  expect_equal(latin$summary$adjacent_checks_across_blocks, 0)
  expect_message(check_design(pw), "no checks next to each other inside a block \\(0 pairs\\)")
  suppressMessages(check_reset())
})

test_that("blocks with at least as many check plots as checks hold every check", {
  pb <- rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
  for (sd in 1:3) {
    pw <- augmented_row_column_design(paste0("E", 1:173), LETTERS[1:4], rows_per_rowgroup = c(6, 5, 6, 5),
                                      cols_per_colgroup = c(3, 4, 3), checks_per_block = pb,
                                      n_candidates = 2, seed = sd)
    ck <- pw$design[pw$design$is_control == 1, ]
    tab <- table(factor(ck$block, 1:12), ck$all_entries)
    full <- as.vector(t(pb)) >= 4
    expect_true(all(tab[full, ] >= 1))                       # every check in 4- and 5-check blocks
    expect_true(all(tab[full, ] <= 2))                       # at most one repeat
    expect_true(all(tab[!full, ] <= 1))                      # incomplete blocks: distinct checks
  }
  expect_message(expect_true(check_design(pw)), "checks spread evenly in each block")
  broken <- pw
  i <- which(broken$design$is_control == 1 & broken$design$block == 2)
  broken$design$all_entries[i] <- "A"
  expect_message(expect_false(check_design(broken)), "\\[FAIL\\] +checks spread evenly")
  suppressMessages(check_reset())
})


test_that("Figures 3 and 4 of Piepho and Williams (2016) are reproduced", {
  # Figure 3 (Example 2): 6 x 6 blocks are complete blocks of the 4 checks,
  # every column is a complete block with 4 check plots, every row has 2.
  for (m in c("latin", "optimal")) {
    f3 <- augmented_row_column_design(paste0("E", 1:576), LETTERS[1:4], 36, 18, 6, 6,
                                      method = m, n_candidates = if (m == "latin") 100 else 2, seed = 3)
    ck <- f3$design[f3$design$is_control == 1, ]
    expect_true(all(table(factor(ck$block, 1:18), ck$all_entries) == 1))
    expect_true(all(table(factor(ck$col, 1:18), ck$all_entries) == 1))
    expect_true(all(table(factor(ck$row, 1:36)) == 2))
    expect_equal(f3$summary$error_df, 6)
    expect_equal(f3$summary$adjacent_checks_in_blocks, 0)
  }
  # Figure 4 (Example 3): 7 complete blocks and 5 with three check plots, four
  # central blocks with five check plots (one check twice), every check in
  # every column, rows with two or three check plots.
  pb <- rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
  f4 <- augmented_row_column_design(paste0("E", 1:173), LETTERS[1:4], rows_per_rowgroup = c(6, 5, 6, 5),
                                    cols_per_colgroup = c(3, 4, 3), checks_per_block = pb,
                                    n_candidates = 2, seed = 1)
  ck <- f4$design[f4$design$is_control == 1, ]
  tb <- table(factor(ck$block, 1:12), ck$all_entries)
  expect_equal(sum(apply(tb, 1, function(x) all(x >= 1))), 7)
  expect_equal(sum(rowSums(tb) == 3), 5)
  expect_true(all(sapply(c(2, 5, 8, 11), function(b) all(tb[b, ] >= 1) && max(tb[b, ]) == 2)))
  expect_true(all(table(factor(ck$col, 1:10), ck$all_entries) >= 1))
  expect_true(all(table(factor(ck$row, 1:22)) %in% 2:3))
  expect_equal(f4$summary$error_df, 7)
})
