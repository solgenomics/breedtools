# =============================================================================
# breedtools: a step-by-step tour of every design and analysis
#
# Run this script section by section in RStudio (Ctrl+Enter line by line, or
# Ctrl+Alt+T for one section; the outline pane lists the sections). Each
# section designs a trial, looks at it, saves the field book, simulates
# yields, analyses them, and ends with a few checks that print
#   [OK]   when the result is what we expect, or
#   [FAIL] when something is wrong.
# Where a method comes from a paper, the checks compare with the published
# numbers.
#
#   0. Setup
#   1. Augmented row-column design (checks, one location)
#   2. Flexible augmented row-column designs (Piepho and Williams, 2016)
#   3. Augmented p-rep design (no checks, several locations)
#   4. p-rep designs built from alpha-arrays (Williams, Piepho and Whitaker, 2011)
#   5. Two-stage analysis of a multi-location trial
#   6. Resolvable row-column design (Williams, John and Whitaker, 2006)
#   7. Multi-environment trials from Breedbase
#   8. Save a design to Breedbase (BrAPI), including agricolae field books
#   9. Summary of the checks
#
# The yields are simulated, so we know the true genotype values and can see
# how well each analysis recovers them.
#
# For the full automated test battery (thousands of checks, used in CI), run
#   Rscript scripts/run_all_checks.R
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

# Install once (pick the option that applies to you):
# install.packages("remotes")
# remotes::install_github("solgenomics/breedtools")    # from GitHub
devtools::install("~/Documents/codeR/breedtools")  # from a local copy
library(breedtools)

# Folder for files this script writes (field books, plots).
# Change it to a project folder of your own, e.g. "~/trials/2026".
out_dir <- file.path(tempdir(), "breedtools_example")

# Check that everything is ready: breedtools loaded, optional packages for
# mixed (lme4) and spatial (SpATS) models, and a folder to save files in.
check_setup(out_dir)
has_lme4 <- requireNamespace("lme4", quietly = TRUE)
has_spats <- requireNamespace("SpATS", quietly = TRUE)

# From here on, checks print [OK] / [FAIL] and are counted:
#   check_design(design)      standard checks of a design
#   check_result("what", ok)  one check of your own
#   check_summary()           what was checked, and what is correct
# prediction_accuracy(fit, truth) is the correlation between adjusted means
# and the true genotype values (known here because the data are simulated).


# ---- 1. Augmented row-column design -----------------------------------------
check_section("1. Augmented row-column design")
# Unreplicated entries plus a few replicated checks. Rows and columns are
# grouped into blocks; every block gets one plot of each check, and the checks
# are used to estimate row, column and block effects.

# Start with the entry count when field dimensions or checks are undecided:
planning_options <- design_options(100)
print(planning_options[, c("design", "n_checks", "rows_in_field", "cols_in_field", "n_reps",
                           "plots_per_location", "error_df")])
# Use a selected row's generator and arguments; see ?design_options and the
# worked handoff in scripts/validate_designs.R, section 0a.

# 1a. Plan: which layouts of a 12 x 12 field can be analysed? Each suggested
#     block has at least as many rows and columns as checks. More checks mean
#     fewer entries but more error df (more reliable tests); aim for at least ~10.
print(augmented_row_column_layouts(rows_in_field = 12, cols_in_field = 12, n_controls = 2:5))

# Only layouts that fit 100 entries with at least 6 error df:
print(augmented_row_column_layouts(12, 12, n_controls = 2:6, min_error_df = 6, n_entries = 100))

# Details of the layout we pick: 4 x 4 blocks with 4 checks.
print(augmented_row_column_capacity(
  rows_in_field  = 12,
  cols_in_field  = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  n_controls     = 4
))
# 108 entries fit, with 36 check plots and 6 error df.

# 1b. Generate the design.
entries <- sprintf("LINE_%03d", 1:108)
check_names <- c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D")

rc <- augmented_row_column_design(
  treatments     = entries,
  controls       = check_names,
  rows_in_field  = 12,
  cols_in_field  = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  n_candidates   = 200,
  seed           = 2026
)
print(rc$summary)
print(head(rc$design))
check_design(rc)          # standard checks of the design

# 1c. Look at it.
plot_augmented_row_column_design(rc, label = "all_entries", fill = "type")
rc_checks <- rc$design[rc$design$is_control == 1, ]
print(table(rc_checks$all_entries, rc_checks$block))   # one of each check per block

# 1d. Save the field book and the Breedbase upload file.
# write.csv(rc$design, file.path(out_dir, "rc_fieldbook.csv"), row.names = FALSE)
# write.csv(rc$breedbase_design, file.path(out_dir, "rc_breedbase.csv"), row.names = FALSE)

# 1e. Simulate yields: a fertility gradient, true genotype values, noise,
#     and three lost plots.
set.seed(1)
genotype_value <- setNames(rnorm(112, sd = 4), c(entries, check_names))
rc_trial <- rc$design
rc_trial$yield <- 60 + 0.8 * rc_trial$row - 0.3 * rc_trial$col +
  genotype_value[rc_trial$all_entries] + rnorm(nrow(rc_trial), sd = 2)
rc_trial$yield[c(7, 50, 101)] <- NA

#Plotting yield
rc_trial$yield
hist(rc_trial$yield)


# 1f. Analyse.
rc_fit <- analyze_augmented_row_column_design(rc_trial, response = "yield")

# Results: the ANOVA table, the best entries and the precision.
print(rc_fit$anova_formatted)                  # ANOVA (rows, columns, blocks, genotypes)
print(head(rc_fit$adjusted_means[, c("rank", "genotype", "raw_mean", "adjusted_mean", "std_error")], 10))
print(rc_fit$diagnostics[, c("sigma", "mean_sed", "lsd_5pct")])   # precision of comparisons
# The full ranking, best first, for selection decisions:
# write.csv(rc_fit$adjusted_means, file.path(out_dir, "rc_ranking.csv"), row.names = FALSE)

# Optional: random rows, columns and blocks (recover inter-block
# information), random entries (BLUPs), and a SpATS spatial model.
if (has_lme4) {
  rc_random <- analyze_augmented_row_column_design(rc_trial, "yield", random_blocks = TRUE)
  rc_blup <- analyze_augmented_row_column_design(rc_trial, "yield", random_blocks = TRUE,
                                                 entries = "random")
  print(rc_random$variance_components)          # variance of rows, columns, blocks
  print(head(rc_blup$adjusted_means[, c("rank", "genotype", "adjusted_mean", "std_error")], 10))
}
if (has_spats) {
  rc_spatial <- analyze_augmented_row_column_design(rc_trial, "yield", spatial = TRUE)
  print(head(rc_spatial$adjusted_means[, c("rank", "genotype", "adjusted_mean", "std_error")], 10))
}
# Which model compares entries most precisely? (smaller mean SED is better)
print(c(fixed = rc_fit$diagnostics$mean_sed,
        random_blocks = if (has_lme4) rc_random$diagnostics$mean_sed else NA,
        spatial = if (has_spats) rc_spatial$diagnostics$mean_sed else NA))

# 1g. Checks of the analysis.
check_result("each check exactly once in every block",
             all(table(rc_checks$all_entries, rc_checks$block) == 1))
ref <- lm(yield ~ factor(rowgroup) + factor(rowgroup):factor(row) + factor(colgroup) +
            factor(colgroup):factor(col) + factor(rowgroup):factor(colgroup) + all_entries,
          data = rc_trial)
check_result("error df and residual variance equal a plain lm() fit",
             rc_fit$diagnostics$residual_df == df.residual(ref) &&
               isTRUE(all.equal(rc_fit$diagnostics$sigma, summary(ref)$sigma)))
entry_means <- rc_fit$adjusted_means[rc_fit$adjusted_means$genotype %in% entries, ]
# Unreplicated entries have one plot each, so their adjusted means carry the
# full plot error: a correlation of about 0.7-0.9 with the truth is typical.
print(prediction_accuracy(entry_means, genotype_value))
check_result("adjusted means track the true values (r > 0.7)", prediction_accuracy(entry_means, genotype_value) > 0.7)
if (has_lme4) {
  check_result("random blocks: variance components estimated",
               nrow(rc_random$variance_components) == 4)
  check_result("random entries: entry values are BLUPs", grepl("BLUP", rc_blup$message))
}
if (has_spats) {
  check_result("spatial model fitted", rc_spatial$diagnostics$model_type == "spatial")
}


# ---- 2. Flexible augmented row-column designs (Piepho and Williams, 2016) ---
check_section("2. Piepho and Williams (2016) layouts")
# Row and column groups can have different sizes, and blocks can hold fewer
# check plots than there are checks (incomplete blocks) or more (extra
# replication). Checks are then placed by a D-optimal search.

# 2a. Table 2 of the paper: entries and error df of the three examples.
pw_ex1 <- augmented_row_column_capacity(6, 9, 3, 3, n_controls = 3)
pw_ex2 <- augmented_row_column_capacity(36, 18, 6, 6, n_controls = 4)
checks_per_block <- rbind(c(4, 5, 4),     # row group 1: 4, 5 and 4 check plots
                          c(3, 5, 3),     # row group 2
                          c(3, 5, 4),     # row group 3
                          c(3, 5, 3))     # row group 4
pw_ex3 <- augmented_row_column_capacity(
  n_controls        = 4,
  rows_per_rowgroup = c(6, 5, 6, 5),
  cols_per_colgroup = c(3, 4, 3),
  checks_per_block  = checks_per_block
)
print(rbind(pw_ex1, pw_ex2, pw_ex3)[, c("total_plots", "check_plots", "entry_plots", "error_df")])

# 2b. Example 3: 173 entries and 4 checks on 22 rows x 10 columns.
pw <- augmented_row_column_design(
  treatments        = paste0("E", 1:173),
  controls          = c("A", "B", "C", "D"),
  rows_per_rowgroup = c(6, 5, 6, 5),
  cols_per_colgroup = c(3, 4, 3),
  checks_per_block  = checks_per_block,
  n_candidates      = 3,
  seed              = 1
)
pw$breedbase_design$all_entries

print(pw$summary[, c("method", "check_plots", "entry_plots", "error_df", "d_criterion")])
check_design(pw)
plot_augmented_row_column_design(pw, label = "type", fill = "type")

# 2c. Simulate and analyse.
set.seed(2)
pw_value <- setNames(rnorm(177, sd = 3), c(paste0("E", 1:173), LETTERS[1:4]))
pw_trial <- pw$design
pw_trial$yield <- 40 + 0.5 * pw_trial$row + pw_value[pw_trial$all_entries] + rnorm(nrow(pw_trial))

#Plotting yield
library(ggplot2)

ggplot(pw_trial, aes(x = yield)) +
  geom_histogram(
    bins = 20,
    fill = "steelblue",
    color = "white",
    alpha = 0.85
  ) +
  labs(
    title = "Yield Distribution",
    x = "Yield",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)


pw_fit <- analyze_augmented_row_column_design(pw_trial, response = "yield")

# Results: ANOVA (7 error df, as in Table 2 of the paper) and the best entries.
print(pw_fit$anova_formatted)
print(head(pw_fit$adjusted_means[, c("rank", "genotype", "adjusted_mean", "std_error")], 10))
print(pw_fit$diagnostics[, c("sigma", "mean_sed", "lsd_5pct")])
# write.csv(pw_fit$adjusted_means, file.path(out_dir, "pw_ranking.csv"), row.names = FALSE)

# 2d. Figure 3 of the paper (Example 2): 4 checks on 36 rows x 18 columns in
#     6 x 6 blocks. Every block and every column is a complete block of the
#     checks, and every row has two check plots.
fig3 <- augmented_row_column_design(
  treatments     = paste0("E", 1:576),
  controls       = c("A", "B", "C", "D"),
  rows_in_field  = 36,
  cols_in_field  = 18,
  rows_per_block = 6,
  cols_per_block = 6,
  n_candidates   = 200,
  seed           = 3
)
check_design(fig3)
plot_augmented_row_column_design(fig3, label = "type", fill = "type")
fig3_checks <- fig3$design[fig3$design$is_control == 1, ]
print(table(fig3_checks$all_entries, fig3_checks$col))   # each check once in every column

# 2e. Checks against the paper.
check_result("Table 2: Example 1 has 36 entries and 0 error df (cannot be used)",
             pw_ex1$entry_plots == 36 && pw_ex1$error_df == 0)
check_result("Table 2: Example 2 has 576 entries and 6 error df",
             pw_ex2$entry_plots == 576 && pw_ex2$error_df == 6)
check_result("Table 2: Example 3 has 173 entries and 7 error df",
             pw_ex3$entry_plots == 173 && pw_ex3$error_df == 7)
pw_check_plots <- pw$design[pw$design$is_control == 1, ]
check_result("Example 3: 3 to 5 check plots per block, as requested",
             all(as.vector(table(pw_check_plots$block)) == as.vector(t(checks_per_block))))
check_result("Example 3: the analysis has 7 error df", pw_fit$diagnostics$residual_df == 7)
check_result("Figure 3: every block and every column holds each check once",
             all(table(fig3_checks$all_entries, fig3_checks$block) == 1) &&
               all(table(fig3_checks$all_entries, fig3_checks$col) == 1))
check_result("Figure 3: every row has two check plots", all(table(fig3_checks$row) == 2))
check_result("Figure 4 (Example 3): every check in every column",
             all(table(pw_check_plots$all_entries, pw_check_plots$col) >= 1))
check_result("Example 1 is rejected (no error df)",
             inherits(try(augmented_row_column_design(paste0("E", 1:36), LETTERS[1:3], 6, 9, 3, 3),
                   silent = TRUE), "try-error"))


# ---- 3. Augmented p-rep design -----------------------------------------------
check_section("3. Augmented p-rep design")
# A multi-location trial without checks: every entry is planted once at every
# location, and 20% of the entries get a second plot at each location.

# 3a. Generate.
prep_entries <- sprintf("G%03d", 1:120)
sites <- c("North", "Central", "South")

prep <- augmented_prep_design(
  treatments   = prep_entries,
  locations    = sites,
  block_size   = 10,           # at most 10 plots per incomplete block
  prep_prop    = 0.2,          # 20% of the entries repeated at each location
  field_cols   = 10,           # field columns, for row/column coordinates
  n_candidates = 200,
  seed         = 2026
)
print(prep$summary[, c("n_plots", "n_blocks", "max_pairwise_concurrence",
                       "n_disconnected_locations", "efficiency_factor")])
print(table(table(prep$design$all_entries)))   # 3 plots for most entries, 4 for the repeated ones
check_design(prep)

# 3b. Look at one location and save the field book.
# A p-rep design has no checks: the repeated entries (orange) replace them.
plot_augmented_prep_design(prep, locations = "North")
# write.csv(prep$design, file.path(out_dir, "prep_fieldbook.csv"), row.names = FALSE)

# 3c. Simulate yields with location and block effects.
set.seed(3)
prep_value <- setNames(rnorm(120, sd = 4), prep_entries)
block_id <- interaction(prep$design$location, prep$design$block)
block_effect <- setNames(rnorm(nlevels(block_id), sd = 2), levels(block_id))
prep_trial <- prep$design
prep_trial$yield <- 50 + c(North = 0, Central = 8, South = -5)[prep_trial$location] +
  block_effect[as.character(block_id)] + prep_value[prep_trial$all_entries] +
  rnorm(nrow(prep_trial), sd = 3)

# 3d. Analyse all locations together.
prep_fit <- analyze_augmented_prep_design(prep_trial, response = "yield")

# Results: ANOVA (locations, blocks, genotypes) and the best entries.
print(prep_fit$anova_formatted)
print(head(prep_fit$adjusted_means[, c("rank", "genotype", "n_observations", "raw_mean",
                                       "adjusted_mean", "std_error")], 10))
print(prep_fit$diagnostics[, c("sigma", "mean_sed", "lsd_5pct")])
write.csv(prep_fit$adjusted_means, file.path(out_dir, "prep_ranking.csv"), row.names = FALSE)

# One location on its own (possible because every location is connected).
north_fit <- analyze_augmented_prep_design(prep_trial[prep_trial$location == "North", ], "yield")
print(north_fit$anova_formatted)

# Optional: genotype x location interaction as a random effect.
if (has_lme4) {
  prep_gxl <- analyze_augmented_prep_design(prep_trial, "yield", genotype_by_location = "random")
  print(prep_gxl$variance_components)            # genotype x location and residual variance
  print(head(prep_gxl$adjusted_means[, c("rank", "genotype", "adjusted_mean", "std_error")], 10))
}

# 3e. Checks.
check_result("the North location can be analysed on its own",
             !anyNA(north_fit$adjusted_means$adjusted_mean))
check_result("adjusted means beat raw means at finding the best entries",
             prediction_accuracy(prep_fit$adjusted_means, prep_value) >
               cor(prep_fit$adjusted_means$raw_mean, prep_value[prep_fit$adjusted_means$genotype]))


# ---- 4. p-rep designs built from alpha-arrays --------------------------------
check_section("4. p-rep designs from alpha-arrays")
# Williams, Piepho and Whitaker (2011) build p-rep designs from alpha-arrays:
# every entry is repeated at exactly one location, every location has two
# replicates, and every block holds m repeated entries.

# 4a. Which settings are possible? (compare with Table 3 of the paper)
print(augmented_prep_parameters(n_treatments = 42, n_locations = 3, max_block_size = 6))

# 4b. Example 1 of the paper: 20 entries, 2 locations, blocks of 3.
alpha_ex1 <- augmented_prep_design(paste0("G", 1:20), c("L1", "L2"),
                                   method = "alpha", block_size = 3, seed = 1)
# 4c. Example 3 of the paper: 42 entries, 3 locations, blocks of 4 (alpha2).
alpha_ex3 <- augmented_prep_design(paste0("G", 1:42), c("L1", "L2", "L3"),
                                   method = "alpha2", block_size = 4, field_cols = 7, seed = 1)
print(rbind(alpha_ex1$summary, alpha_ex3$summary)[, c("method", "m", "s", "n_plots", "n_blocks",
                                                      "max_pairwise_concurrence", "efficiency_factor")])
# Orange = entries repeated at this location (they replace checks); thick
# black lines = blocks; the blue line separates the two replicates.
plot_augmented_prep_design(alpha_ex3, locations = "L1")
check_design(alpha_ex1)
check_design(alpha_ex3)

# 4d. Simulate and analyse Example 3. Each location has two replicates, so
#     the ANOVA has replicates within locations and blocks within replicates.
set.seed(7)
alpha_value <- setNames(rnorm(42, sd = 3), paste0("G", 1:42))
alpha_trial <- alpha_ex3$design
alpha_trial$yield <- 30 + c(L1 = 0, L2 = 4, L3 = -2)[alpha_trial$location] +
  alpha_value[alpha_trial$all_entries] + rnorm(nrow(alpha_trial))
alpha_fit <- analyze_augmented_prep_design(alpha_trial, response = "yield")
print(alpha_fit$anova_formatted)
print(head(alpha_fit$adjusted_means[, c("rank", "genotype", "adjusted_mean", "std_error")], 10))
print(alpha_fit$diagnostics[, c("sigma", "mean_sed", "lsd_5pct")])

# 4e. Checks against the paper.
check_result("Example 1: m = 2, 5 blocks per replicate, blocks of 3",
             alpha_ex1$summary$m == 2 && alpha_ex1$summary$s == 5 && all(alpha_ex1$design$block_size == 3))
check_result("Example 3: 14 blocks of 4 plots at each location",
             alpha_ex3$summary$n_blocks == 42 && all(alpha_ex3$design$block_size == 4))
check_result("every entry repeated at exactly one location",
             all(table(alpha_ex3$design$all_entries) == 4))
check_result("no pair of entries shares a block more than once",
             alpha_ex1$summary$max_pairwise_concurrence == 1 && alpha_ex3$summary$max_pairwise_concurrence == 1)
check_result("the ANOVA has replicates within locations",
             "Replicates, nested within locations" %in% alpha_fit$anova_formatted$Source)
check_result("adjusted means track the true values (r > 0.8)",
             prediction_accuracy(alpha_fit, alpha_value) > 0.8)

### Second Part: Two Stage Analysis P-Rep
# ---- 5. Two-stage analysis of a multi-location trial -------------------------
check_section("5. Two-stage analysis")
# Locations rarely have the same error variance. The two-stage analysis
# analyses every location on its own and then combines them, giving less
# weight to noisy locations.

# 5a. Make the South location much noisier (error SD 12 instead of 3).
#     Note: each location has only 10 error df here (from its 24 repeated
#     entries), so a test of equal variances only detects large differences.
het_trial <- prep_trial
south <- het_trial$location == "South"
set.seed(6)
het_trial$yield[south] <- het_trial$yield[south] + rnorm(sum(south), sd = 12)

# 5b. Two-stage analysis, compared with treating all locations alike.
two <- analyze_augmented_prep_two_stage(het_trial, response = "yield")
unweighted <- analyze_augmented_prep_two_stage(het_trial, "yield", weighting = "none")

# Results: per-location error variances, the test of equal variances, the
# genotype x location variance and the best entries.
print(two$locations[, c("location", "n_plots", "residual_variance", "residual_df", "mean_sed")])
print(two$heterogeneity)
print(two$variance_components)
print(head(two$adjusted_means, 10))
write.csv(two$adjusted_means, file.path(out_dir, "two_stage_ranking.csv"), row.names = FALSE)
# Accuracy (correlation with the true values): weighting by precision helps.
print(c(weighted = prediction_accuracy(two, prep_value),
        unweighted = prediction_accuracy(unweighted, prep_value)))

# Optional: SpATS spatial model at every location in stage 1.
if (has_spats) {
  two_spatial <- analyze_augmented_prep_two_stage(het_trial, "yield", stage1 = "spatial")
  print(two_spatial$locations[, c("location", "residual_variance", "spatial_effective_dim")])
  print(head(two_spatial$adjusted_means, 10))
}

# 5c. Checks.
check_result("Hartley's test flags the noisy location", isTRUE(two$heterogeneity$heterogeneous))
check_result("the weighted analysis is more accurate than the unweighted one",
             prediction_accuracy(two$adjusted_means, prep_value) > prediction_accuracy(unweighted$adjusted_means, prep_value))


# ---- 6. Resolvable row-column design ------------------------------------------
check_section("6. Resolvable row-column design")
# When every entry can be replicated: each replicate is a grid of rows and
# columns that holds every treatment once (Williams, John and Whitaker, 2006).

# 6a. Generate: 42 treatments, 3 replicates of 6 rows x 7 columns, latinized
#     so that no treatment appears twice in a long column.
rrc <- resolvable_row_column_design(
  treatments   = paste0("T", 1:42),
  n_reps       = 3,
  rows_per_rep = 6,
  cols_per_rep = 7,
  latinize     = TRUE,
  n_starts     = 1,
  seed         = 1
)
print(rrc$summary[, c("efficiency_rc", "efficiency_rows", "efficiency_cols",
                      "efficiency_spatial", "upper_bound_spatial", "latin_ok")])
plot_resolvable_row_column_design(rrc)
check_design(rrc)
write.csv(rrc$design, file.path(out_dir, "rrc_fieldbook.csv"), row.names = FALSE)

# 6b. Evaluate the design of Table 1 of the paper (made with CycDesigN).
#     Each block of 6 rows is one replicate; columns 1 to 7.
table1 <- rbind(
  c(19, 34, 27,  3, 15, 12,  2), c(40,  4,  7, 20, 21, 33, 18), c(23, 35, 17, 29, 22, 28, 38),
  c(10, 14, 31, 39, 32,  5, 42), c(13,  6, 16, 36, 11, 37, 41), c(25,  8,  1, 26,  9, 30, 24),
  c(38, 34, 26, 33,  7, 10, 11), c( 1, 13, 15, 35,  5,  3, 18), c(32, 20, 17, 41,  8, 24, 19),
  c(12, 30,  6, 39, 23, 21, 16), c( 4, 31, 42, 25,  2, 37, 29), c(36, 22, 40, 27,  9, 14, 28),
  c(20, 26, 23, 37, 32, 27, 35), c(10, 22, 16, 18, 42,  8, 12), c(28,  4, 15, 39,  7, 13, 24),
  c( 1, 41, 14, 34, 30, 29, 40), c( 2,  5, 36, 17, 25, 21, 11), c( 6, 19, 33,  9,  3, 38, 31)
)
paper_design <- data.frame(
  replicate   = rep(rep(1:3, each = 6), times = 7),
  row_in_rep  = rep(rep(1:6, 3), times = 7),
  col_in_rep  = rep(1:7, each = 18),
  all_entries = paste0("T", as.vector(table1))
)
paper_eff <- resolvable_row_column_efficiency(paper_design)
print(unlist(paper_eff[c("efficiency_rc", "efficiency_spatial", "upper_bound_spatial")]))

# 6c. Simulate yields with row and column effects inside each replicate.
set.seed(5)
rrc_value <- setNames(rnorm(42, sd = 2), paste0("T", 1:42))
rrc_trial <- rrc$design
rrc_trial$yield <- 50 + rrc_value[rrc_trial$all_entries] +
  rnorm(18, sd = 1.2)[(rrc_trial$replicate - 1) * 6 + rrc_trial$row_in_rep] +
  rnorm(21, sd = 0.8)[(rrc_trial$replicate - 1) * 7 + rrc_trial$col_in_rep] +
  rnorm(nrow(rrc_trial))

# 6d. Analyse: the row-column mixed model, compared with spatial models.
rrc_fit <- analyze_resolvable_row_column(rrc_trial, response = "yield", model = "rc",
                                         long_col = "long_col")

# Results: the intra-block ANOVA, the comparison of variance models (lower
# AIC is better), the variance parameters and the best treatments.
print(rrc_fit$anova)
print(rrc_fit$model_comparison[, c("model", "deviance", "AIC", "lrt_p_value")])
print(rrc_fit$variance_parameters)
print(head(rrc_fit$adjusted_means, 10))
print(rrc_fit$diagnostics[, c("mean_sed", "lsd_5pct")])
write.csv(rrc_fit$adjusted_means, file.path(out_dir, "rrc_ranking.csv"), row.names = FALSE)

# 6e. Checks against the paper.
check_result("Table 1 design: row-column efficiency 0.6941 (as published)",
             round(paper_eff$efficiency_rc, 4) == 0.6941)
check_result("Table 1 design: spatial efficiency 0.9204, upper bound 0.9268 (as published)",
             round(paper_eff$efficiency_spatial, 4) == 0.9204 && round(paper_eff$upper_bound_spatial, 4) == 0.9268)
check_result("generated design within 2% of the paper's row-column efficiency",
             rrc$summary$efficiency_rc >= 0.98 * 0.6941)
print(prediction_accuracy(rrc_fit, rrc_value))
check_result("adjusted means track the true values (r > 0.8)", prediction_accuracy(rrc_fit$adjusted_means, rrc_value) > 0.8)


# ---- 7. Multi-environment trials from Breedbase -------------------------------
check_section("7. Multi-environment trials")
# Before combining trials: are their error variances similar, and which
# environments rank genotypes alike? The package ships an example file with
# four trials.

pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
locations_file <- system.file("extdata", "breedbase_locations.csv", package = "breedtools")
pheno_raw <- read_breedbase_file(pheno_file)   # also reads .xlsx and .txt files
print(head(pheno_raw))

# 7a. Hartley's Fmax test of the error variances, and compatible groups.
variance_check <- run_breedbase_variance_check(pheno_file,
                                               trait_col = "fresh root yield|CO_334:0000013",
                                               n_sim = 20000)
print(variance_check$env_anova[, c("environment", "ms_error", "mean", "cv_percent")])
print(variance_check$overall_hartley)
print(variance_check$maximal_compatible_groups)

# 7b. Group environments by genotype-by-environment interaction (Lin, 1982).
strata <- run_environment_stratification(pheno_file, trait = "fresh root yield")
print(strata$group_summary[, c("group_id", "locations", "n_env", "f_value", "p_value")])

# 7c. Draw the groups on a map (opens in the RStudio viewer).
locations <- read_locations(locations_file, phenotype_data = strata)
print(map_environment_groups(strata, locations))

# 7d. Checks.
check_result("four environments read from the Breedbase file", nrow(variance_check$env_anova) == 4)
check_result("Fmax = largest / smallest error mean square",
             isTRUE(all.equal(variance_check$overall_hartley$fmax_observed,
                       max(variance_check$env_anova$ms_error) / min(variance_check$env_anova$ms_error))))
check_result("two environment groups: Geneva + Ithaca and Ibadan + Ubiaja",
             setequal(strata$group_summary$locations, c("Geneva, Ithaca", "Ibadan, Ubiaja")))


# ---- 8. Save a design to Breedbase -----------------------------------------------
check_section("8. Save to Breedbase")
# Designs are saved through BrAPI: log in, choose the breeding program and the
# folder, and save. Nothing is written until every name has been checked and
# you confirm. Breedbase words: a "folder" holds trials; a "trial" holds plots.

# 8a. Without a server: see exactly what would be sent.
to_send <- save_breedbase(rc, trial_name = "AYT_2026_Ibadan", year = 2026,
                          location = "Ibadan", dry_run = TRUE,
                          description = "Advanced yield trial of 108 lines and 4 checks, Ibadan 2026")
print(to_send$studies[[1]][c("studyName", "studyType", "experimentalDesign")])
str(to_send$plots[[1]][[1]])     # the first plot: name, accession, row/column, rep, block

# Plot names are made automatically as trialName_accession_block_rep:
print(head(vapply(to_send$plots[[1]], function(p) p$observationUnitName, character(1))))
# Or choose another pattern, e.g. with the plot number and field position:
#   save_breedbase(..., plot_name_pattern = "{trial}_{plot}_r{row}c{col}")
# Or give your own names, one per plot (or the name of a column you added):
#   save_breedbase(..., plot_names = my_plot_names)

check_result("one Breedbase plot per field plot", length(to_send$plots[[1]]) == nrow(rc$design))
check_result("every plot name is unique",
             !anyDuplicated(vapply(to_send$plots[[1]], function(p) p$observationUnitName, character(1))))
check_result("checks are marked as checks",
             sum(vapply(to_send$plots[[1]], function(p) p$observationUnitPosition$entryType == "CHECK",
                        logical(1))) == sum(rc$design$is_control))

# 8b. Save to YOUR Breedbase, step by step.
#     1) Fill in your details below.
#     2) Set save_to_breedbase <- TRUE and run this section.
#     Your password is asked for in a pop-up (or read from the
#     BREEDBASE_PASSWORD environment variable): never write it in the script.
#     Nothing is written until every name is checked (step 3 is a dry run).

save_to_breedbase <- FALSE                      # TRUE to run this section

breedbase_url   <- "https://cassavabase.org"    # your Breedbase site
breedbase_user  <- "your_username"              # your Breedbase username
my_program      <- "IITA"                       # breeding program
my_folder       <- "AYT 2026"                   # folder inside the program, or NULL
                                                #   to save the trial directly in the program
my_location     <- "Ibadan"                     # Breedbase location name
my_trial_name   <- "AYT_2026_Ibadan"            # name of the NEW trial
my_description  <- "Advanced yield trial of 108 lines and 4 checks, Ibadan 2026"

if (save_to_breedbase) {
  # Step 1. Log in (the password is asked for in a pop-up).
  con <- breedbase_login(breedbase_url, username = breedbase_user)

  # Step 2. Look up the exact names you can use.
  print(breedbase_programs(con))                    # breeding programs
  print(breedbase_folders(con, program = my_program))  # folders in your program
  print(head(breedbase_locations(con), 20))         # location names

  # Step 2b. Check the names in one go: location, every entry and check
  #          (as accessions), program, folder and trial name. Problems are
  #          listed together with close matches ("ibadan" -> "Ibadan").
  #          Tip: run this BEFORE designing, with your list of names:
  #            breedbase_check(con, accessions = c(my_entries, my_checks), location = my_location)
  ready <- breedbase_check(con, accessions = rc, location = my_location,
                           program = my_program, folder = my_folder, trial_name = my_trial_name)
  if (!attr(ready, "ok")) stop("Fix the names listed above, then run this section again.")

  # Step 3. Dry run: checks the program, folder, location, that every entry and
  #         check is an accession in Breedbase, and that the trial and plot
  #         names are new. Nothing is written.
  save_breedbase(rc, con, program = my_program, folder = my_folder,
                 trial_name = my_trial_name, year = 2026, location = my_location,
                 description = my_description, create_folder = TRUE, dry_run = TRUE)

  # Step 4. Save. Setting save_to_breedbase <- TRUE is your 'yes', so no
  #         question is asked (confirm = FALSE). The plots are read back to verify.
  saved <- save_breedbase(rc, con, program = my_program, folder = my_folder,
                          trial_name = my_trial_name, year = 2026, location = my_location,
                          description = my_description, create_folder = TRUE, confirm = FALSE)
  print(saved)                                      # trial id and plots saved
  saved$trial_page                                  # link: open it in your browser

  # A p-rep design becomes one trial per location; map your location names
  # to Breedbase locations:
  save_breedbase(prep, con, program = my_program, folder = my_folder,
                 trial_name = "PREP_2026", year = 2026,
                 location = c(North = "Ibadan", Central = "Ikenne", South = "Ubiaja"),
                 description = "p-rep trial of 120 entries in three locations, 2026",
                 create_folder = TRUE, dry_run = TRUE)
} else {
  message("Saving to Breedbase skipped: fill in your details above and set ",
          "save_to_breedbase <- TRUE to run it.")
}


# 8c. Phenotypic data to Breedbase.
# After harvest, the values go to the plots saved above, matched by plot name
# and trait name (the trait name as in Breedbase, e.g. "fresh root yield").

# Step 1. A sheet with the plot names of the trial and one empty column per trait.
pheno_sheet <- phenotype_template(rc, trial_name = my_trial_name, location = my_location,
                                  traits = c("fresh root yield", "dry matter content percentage"))
print(head(pheno_sheet))
write.csv(pheno_sheet, file.path(out_dir, "phenotyping_sheet.csv"), row.names = FALSE)

# Step 2. Fill it in. Here the values are simulated; with real data, fill the
#         sheet in a spreadsheet and read it back:
#           pheno_sheet <- read.csv("phenotyping_sheet.csv", check.names = FALSE)
set.seed(2026)
pheno_sheet$`fresh root yield` <- round(rnorm(nrow(pheno_sheet), 25, 4), 1)
pheno_sheet$`dry matter content percentage` <- round(rnorm(nrow(pheno_sheet), 32, 3), 1)

# Step 3. See what would be sent (no server needed): one row per value.
to_send_pheno <- save_breedbase_phenotypes(pheno_sheet, dry_run = TRUE)
print(head(to_send_pheno))
check_result("one value per plot and trait", nrow(to_send_pheno) == 2 * nrow(rc$design))
check_result("phenotype plot names match the saved plot names",
             setequal(to_send_pheno$plot_name,
                      vapply(to_send$plots[[1]], function(p) p$observationUnitName, character(1))))

# Step 4. Save to YOUR Breedbase. The trial must be saved first (above), and
#         your account needs "submitter" rights. Set
#         save_phenotypes_to_breedbase <- TRUE to run it.
save_phenotypes_to_breedbase <- FALSE

if (save_phenotypes_to_breedbase) {
  if (!exists("con")) con <- breedbase_login(breedbase_url, username = breedbase_user)
  # Checks that every plot and trait is in Breedbase; writes nothing.
  save_breedbase_phenotypes(pheno_sheet, con, dry_run = TRUE)
  # Saves (no question: the TRUE above is your 'yes') and reads the values back.
  saved_pheno <- save_breedbase_phenotypes(pheno_sheet, con, confirm = FALSE)
  print(saved_pheno)                                # values sent and found, per trait
} else {
  message("Saving phenotypes to Breedbase skipped: set ",
          "save_phenotypes_to_breedbase <- TRUE to run it.")
}


# 8d. Analysis results to Breedbase (first and second stage).
# First stage: one value per field plot (the genotype's adjusted mean plus the
# plot's residual: blocks, rows, columns and the spatial trend removed). It is
# saved on the field plots as its own TRAIT, so raw and adjusted values stay
# apart: the measured trait composed with a statistic, e.g.
# "fresh root yield|Adjusted Means from Spatial Correction using SpATS R".
# That trait is found in Breedbase, or created when it does not exist yet.
# Second stage: one value per genotype (adjusted means), saved as a new
# Breedbase ANALYSIS, like the Breedbase mixed-models tool does.

field <- rc$design
field$yield <- pheno_sheet$`fresh root yield`[match(field$plots, pheno_sheet$plot_number)]
use_spats <- requireNamespace("SpATS", quietly = TRUE)
fit_yield <- analyze_augmented_row_column_design(field, "yield", spatial = use_spats)

stage1 <- first_stage_values(fit_yield, rc, trial_name = my_trial_name, location = my_location, trait = "fresh root yield")
print(head(stage1))
check_result("first stage: one value per field plot", nrow(stage1) == nrow(rc$design))

stage2 <- second_stage_values(fit_yield, trait = "fresh root yield")
print(head(stage2))
check_result("second stage: one value per genotype", nrow(stage2) == length(unique(rc$design$all_entries)))

# Without a server: see what would be sent (trait = measured trait + statistic).
stage1_obs <- save_breedbase_first_stage(stage1, dry_run = TRUE)
print(head(stage1_obs))
stage2_sent <- save_breedbase_analysis(stage2, analysis_name = paste0(my_trial_name, "_means"),
                                       program = my_program, description = "Adjusted means, row-column model",
                                       dry_run = TRUE)

# Save to YOUR Breedbase (the trial and its phenotypes must be saved first).
# Set save_analysis_to_breedbase <- TRUE to run it.
save_analysis_to_breedbase <- FALSE

if (save_analysis_to_breedbase) {
  if (!exists("con")) con <- breedbase_login(breedbase_url, username = breedbase_user)
  # First stage on the field plots, as a trait. The adjusted trait is looked
  # up and created if needed (see it first with breedbase_adjusted_trait(con,
  # "fresh root yield")); then plots and values are checked and saved.
  save_breedbase_first_stage(stage1, con, confirm = FALSE)
  # Second stage as an analysis (checks accessions, traits and the name first).
  saved_analysis <- save_breedbase_analysis(stage2, con, analysis_name = paste0(my_trial_name, "_means"),
                                            program = my_program,
                                            description = "Adjusted means, row-column model",
                                            confirm = FALSE)
  print(saved_analysis)
  saved_analysis$analysis_page                     # link: open it in your browser
} else {
  message("Saving analysis results to Breedbase skipped: set ",
          "save_analysis_to_breedbase <- TRUE to run it.")
}


# ---- 8e. Agricolae field books to a new Breedbase trial ------------------------
check_section("8e. Agricolae to Breedbase")

# Generate the design with agricolae, then prepare its field book and trial
# information with prepare_breedbase_trial(). Save using save_breedbase().
# Install once with install.packages("agricolae").
# Replace these example names with accessions already in your Breedbase.
# This section reuses the site/program/folder/location settings above.
save_agricolae_to_breedbase <- FALSE
ag_trial_name <- "AGRICOLAE_RCBD_2026_Ibadan"
ag_description <- "RCBD of three entries and one check generated with agricolae"

if (requireNamespace("agricolae", quietly = TRUE)) {
  ag_accessions <- c("G1", "G2", "G3", "CHECK_A")
  ag_randomized <- agricolae::design.rcbd(ag_accessions, r = 3, seed = 2026)
  print(head(ag_randomized$book))   # agricolae columns: plots, block, ag_accessions

  ag_design <- prepare_breedbase_trial(
    ag_randomized, trial_name = ag_trial_name, year = 2026,
    location = my_location, program = my_program, folder = my_folder,
    description = ag_description, checks = "CHECK_A", prompt = FALSE
  )
  print(ag_design$trial)            # trial name, location, breeding program, year
  print(head(ag_design$design))     # plots, all_entries, block, rep, is_control

  # A bare field book (also zigzag() output) needs its Breedbase design type.
  ag_from_book <- as_breedbase_design(
    ag_randomized$book, accession = "ag_accessions",
    checks = "CHECK_A", design_type = "RCBD"
  )

  # Preview the new trial and plots offline.
  ag_payload <- save_breedbase(ag_design, dry_run = TRUE)
  print(ag_payload$studies[[1]])
  str(ag_payload$plots[[1]][[1]])
  ag_sheet <- phenotype_template(
    ag_design, traits = "fresh root yield"
  )
  print(head(ag_sheet))
  check_result("agricolae: one Breedbase plot per field-book row",
               length(ag_payload$plots[[1]]) == nrow(ag_randomized$book))
  check_result("agricolae: original and bare field books give the same design",
               identical(ag_design$design, ag_from_book$design))

  # Set the flag above to TRUE after entering your settings and accession
  # names. The connected dry run checks every name before creating the trial.
  if (save_agricolae_to_breedbase) {
    ag_con <- breedbase_login(breedbase_url, username = breedbase_user)
    save_breedbase(ag_design, ag_con, create_folder = TRUE, dry_run = TRUE)
    ag_saved <- save_breedbase(ag_design, ag_con, create_folder = TRUE, confirm = FALSE)
    print(ag_saved)
    print(ag_saved$trial_page)
  }
} else {
  message("agricolae example skipped: install.packages(\"agricolae\") to run it.")
}

# For multiple treatment factors, select the genotype column explicitly:
#   as_breedbase_design(your_design, accession = "genotype")
# Other treatment factors and main-plot/subplot hierarchy remain in $book;
# saving uploads individual plots. See ?as_breedbase_design for other types.


# ---- 9. Summary of the checks --------------------------------------------------
check_summary()               # what was checked, by section, and what is correct
# check_summary(details = TRUE) lists every single check.
print(list.files(out_dir))    # files written by this script
