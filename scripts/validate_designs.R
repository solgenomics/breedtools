# =============================================================================
# breedtools: a step-by-step tour, one part per trial design
#
# Run it part by part in RStudio (Ctrl+Enter line by line, or Ctrl+Alt+T for
# one section; the outline pane lists the parts). Every design part follows
# the same steps, from the field plan to the results in Breedbase:
#
#   0. suggest trial dimensions with design_options() when field rows,
#      columns or check counts are unknown: start with the entry count and
#      selected design families (worked examples in section 0a)
#   1. plan and generate the design (one randomization per location)
#   2. look at it, and save the field books
#   3. save the trials to Breedbase
#   4. phenotypes: the sheet with the plot names, and the field data
#   5. plot the trait
#   6. find and remove outliers (IQR rule), listing what was removed
#   7. first stage at every location with SpATS, compared with the raw data
#   8. second stage: combine the locations, weighted by precision
#   9. save phenotypes, first stage (plot values) and second stage
#      (genotype means) to Breedbase
#  10. checks: [OK] when a result is what we expect, [FAIL] when not
#  11. analyze a second trait and add it to the existing analysis projects
#      (worked example in A11)
#  12. alternatively, fit all locations in one single-stage analysis
#      (run A1-A6, then jump to A12; no first/second-stage fit is needed)
#
#   0. Setup
#   0a. Suggest trial dimensions and layout options from entry counts
#   A. Augmented row-column design (unreplicated entries + checks)
#   B. Flexible augmented row-column layouts (Piepho and Williams, 2016)
#   C. Augmented p-rep design (no checks, several locations)
#   D. p-rep designs from alpha-arrays (Williams, Piepho and Whitaker, 2011)
#   E. Resolvable row-column design (Williams, John and Whitaker, 2006)
#   F. Multi-environment trials from Breedbase files
#   H. From Breedbase trials to analyses (choose single, two_stage or both;
#      add subsequent traits to the same analysis projects)
#   I. Agricolae: generate RCBD, prepare trial information and save;
#      optional phenotype collection and analysis (also supports alpha)
#   G. Summary of the checks
#
# The field data are simulated, so we know the true genotype values and can
# see how well each analysis recovers them. With your own data, read the
# filled phenotyping sheet instead (shown in each part).
#
# For the full automated test battery (used in CI), run
#   Rscript scripts/run_all_checks.R
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

# Install once (pick the option that applies to you):
# install.packages("remotes")
# remotes::install_github("solgenomics/breedtools")    # from GitHub
# devtools::install("~/Documents/codeR/breedtools")  # from a local copy
library(breedtools)

# Folder for the files this script writes (field books, sheets, rankings).
out_dir <- file.path(tempdir(), "breedtools_example")

# Check that everything is ready: breedtools, the optional packages for mixed
# (lme4) and spatial (SpATS) models, and the output folder.
check_setup(out_dir)
has_lme4 <- requireNamespace("lme4", quietly = TRUE)
has_spats <- requireNamespace("SpATS", quietly = TRUE)
# First stage: SpATS spatial model when installed, else the design's own model.
stage1_model <- if (has_spats) "spatial" else "design"

# Breedbase. Three switches, so each step can be run on its own (e.g. only
# the analyses, once the trials and phenotypes are saved):
#   save_designs     save the trials (designs) of every part
#   save_phenotypes  save the field data (raw phenotypes)
#   save_analyses    save the analyses: first stage (one value per field
#                    plot) and second stage (one value per genotype)
# With all FALSE, every part only shows what would be sent (dry runs).
# Passwords: set my_password / my_site_password below to run the script with
# Rscript, or leave them NULL to read BREEDBASE_PASSWORD /
# BREEDBASE_SITE_PASSWORD (e.g. in ~/.Renviron) or be asked in a pop-up.
save_designs    <- FALSE
save_phenotypes <- FALSE
save_analyses   <- FALSE
breedbase_url   <- "https://cassavabase.org"     # your Breedbase site
breedbase_user  <- "your_username"                    # your Breedbase username
my_program      <- "IITA"                             # breeding program
my_folder       <- "AYT 2026"                   # folder in the program, or NULL
my_locations    <- c("Ibadan", "Ikenne", "Ubiaja")   # Breedbase location names
my_trait        <- "fresh root yield"                # Breedbase trait of the yields
my_second_trait <- "dry matter content percentage"  # second analyzed trait (A11)
my_year         <- 2026
my_trial_prefix <- "AYT"                            # trials are <prefix>_<year>_<location>
# Analyses need a unique name: a run tag (date and time) keeps them unique
# when the script is run again.
my_run_tag      <- format(Sys.time(), "%Y%m%d_%H%M")
# Test sites such as cassava-test.sgn.cornell.edu have an extra web password
# (the browser pop-up "Restricted Test Site"); analyses are saved through the
# web pages, so give its user name here. Its password is read from
# BREEDBASE_SITE_PASSWORD or asked in a pop-up. NULL for sites without it
# (cassavabase.org).
my_site_user    <- NULL
my_password      <- NULL    # Breedbase password; NULL = BREEDBASE_PASSWORD or pop-up
my_site_password <- NULL    # web password of a test site; NULL = BREEDBASE_SITE_PASSWORD or pop-up

# Part H can read and analyze existing trials independently of the save switches.
read_breedbase_trials <- FALSE
breedbase_analysis_route <- "single"  # "single", "two_stage" or "both" (Part H)

use_breedbase <- save_designs || save_phenotypes || save_analyses || read_breedbase_trials
if (use_breedbase) {
  con <- breedbase_login(breedbase_url, username = breedbase_user, password = my_password,
                         site_username = my_site_user, site_password = my_site_password)
  print(con)                                    # "site (web) password: set" when the web pages open
}
# Analyses are saved only for trials that are in Breedbase.
trials_saved <- function(trial_names) {
  ok <- vapply(trial_names, function(t) !inherits(try(breedbase_trial_plots(con, t), silent = TRUE), "try-error"),
               logical(1))
  if (!all(ok)) message("Analyses not saved: these trials are not in Breedbase yet: ",
                        paste(trial_names[!ok], collapse = ", "), " (save them with save_designs <- TRUE).")
  all(ok)
}

# Only for this tour: field data from true genotype values, a fertility
# gradient over rows and columns, and plot noise.
simulate_yield <- function(d, values, mean = 50, row_trend = 0, col_trend = 0, sd = 2) {
  mean + values[as.character(d$all_entries)] + row_trend * d$row + col_trend * d$col +
    stats::rnorm(nrow(d), sd = sd)
}

# Checks print [OK] / [FAIL] and are counted:
#   check_design(design)      standard checks of a design
#   check_result("what", ok)  one check of your own
#   check_summary()           what was checked, and what is correct
# prediction_accuracy(fit, truth): correlation of adjusted means with the
# true genotype values (known here because the data are simulated).


# =============================================================================
# ---- 0a. Suggest trial dimensions from entry counts ---------------------------
check_section("0a. Plan from entry counts")
# Choose the design families and number of test entries (excluding checks).
# design_options() suggests field rows and columns, checks and replication.
# For augmented row-column, every block has at least n_checks rows and
# n_checks columns: four checks require blocks of at least 4 x 4.
# Leave rows, cols or n_checks unspecified to search, or supply one or
# several candidate values for each known detail.
planning_options <- design_options(100, n_locations = length(my_locations))
print(planning_options[, c("option", "design", "n_checks", "rows_in_field", "cols_in_field",
                           "n_reps", "plots_per_location", "error_df")])
planning_fixed <- design_options(108, design = "augmented_row_column", n_checks = 4,
                                 rows = 12, cols = 12, n_candidates = 200)
print(planning_fixed[, c("option", "rows_in_field", "cols_in_field", "rows_per_block",
                         "cols_per_block", "n_candidates", "check_plots", "error_df", "method")])
# Select a row, add your accession names and generate it:
#   selected <- planning_fixed[1, ]
#   trial <- augmented_row_column_design(
#     treatments = your_entries, controls = your_checks,
#     rows_in_field = selected$rows_in_field, cols_in_field = selected$cols_in_field,
#     rows_per_block = selected$rows_per_block, cols_per_block = selected$cols_per_block,
#     checks_per_block = selected$checks_per_block[[1]], method = selected$method,
#     n_candidates = selected$n_candidates, seed = 2026
#   )
#   check_design(trial)
# For other families, selected$arguments[[1]] contains their layout arguments.
# For p-rep replace its locations with your locations. For fully replicated
# designs include both entries and checks in treatments. See ?design_options.
check_result("planner keeps the requested entry count",
             all(planning_fixed$plots_per_location - planning_fixed$check_plots == 108))
check_result("planner respects the known field dimensions",
             all(planning_fixed$rows_in_field == 12 & planning_fixed$cols_in_field == 12))


# =============================================================================
# ---- A. Augmented row-column design ------------------------------------------
# =============================================================================
check_section("A. Augmented row-column design")
# Unreplicated entries plus a few replicated checks. Rows and columns are
# grouped into blocks; every block gets one plot of each check, and the
# checks estimate the row, column and block effects. The trial is grown at
# three locations, each with its own randomization.

# ---- A1. Plan and generate --------------------------------------------------
# Which layouts of a 12 x 12 field can be analysed? Most entry plots first;
# more checks mean fewer entries but more error df (aim for ~10 or more).
# Every suggested block has at least as many rows and columns as checks.
print(augmented_row_column_layouts(rows_in_field = 12, cols_in_field = 12, n_controls = 2:5))
print(augmented_row_column_capacity(12, 12, rows_per_block = 4, cols_per_block = 4, n_controls = 4))
# 108 entries fit, with 36 check plots and 6 error df.

entries     <- sprintf("LINE_%03d", 1:108)
check_names <- c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D")
rc_sites    <- my_locations
rc_trials   <- paste0(my_trial_prefix, "_", my_year, "_", rc_sites)   # one trial per location

# Suggest the field AND block dimensions from your entry and check counts.
# The fixed field below reproduces this example. Omit rows and/or cols to
# let the planner suggest the dimensions you do not know.
rc_options <- design_options(
  n_entries = length(entries), n_checks = length(check_names),
  design = "augmented_row_column", rows = 12, cols = 12, n_candidates = 200
)
print(rc_options[, c("option", "rows_in_field", "cols_in_field", "rows_per_block",
                     "cols_per_block", "n_candidates", "check_plots", "error_df", "method")])
rc_layout <- rc_options[1, ]       # choose a row after reviewing the options

rc_designs <- lapply(seq_along(rc_sites), function(i) {
  augmented_row_column_design(
    treatments       = entries,
    controls         = check_names,
    rows_in_field    = rc_layout$rows_in_field,
    cols_in_field    = rc_layout$cols_in_field,
    rows_per_block   = rc_layout$rows_per_block,
    cols_per_block   = rc_layout$cols_per_block,
    checks_per_block = rc_layout$checks_per_block[[1]],
    method           = rc_layout$method,
    n_candidates     = rc_layout$n_candidates,
    seed             = 2026 + i,          # a different randomization per location
    trial_name       = rc_trials[i],      # plot names made now: see rc$design$plot_name
    location         = rc_sites[i]
  )
})
names(rc_designs) <- rc_sites

# ---- A2. Look at it and save the field books ------------------------------------
# (The write.csv() lines that save files are commented out: remove the "# " to save them.)
rc <- rc_designs[[1]]
print(rc$summary)
check_design(rc)
plot_augmented_row_column_design(rc, label = "all_entries", fill = "type")
rc_checks <- rc$design[rc$design$is_control == 1, ]
print(table(rc_checks$all_entries, rc_checks$block))     # one of each check per block
for (s in rc_sites) {
  # write.csv(rc_designs[[s]]$design, file.path(out_dir, paste0("rc_fieldbook_", s, ".csv")), row.names = FALSE)
}

# ---- A3. Save the trials to Breedbase ---------------------------------------------
# Without a server: what would be sent for the first location.
rc_sent <- save_breedbase(rc, trial_name = rc_trials[1], year = my_year, location = rc_sites[1],
                          dry_run = TRUE, description = "Advanced yield trial, 108 lines and 4 checks")
print(head(vapply(rc_sent$plots[[1]], function(p) p$observationUnitName, character(1))))  # plot names
if (save_designs) {
  # Check the names first: locations, entries and checks, program, folder,
  # trial names and the plot names of every location's design.
  ready <- breedbase_check(con, accessions = rc_designs, location = rc_sites, program = my_program,
                           folder = my_folder, trial_name = rc_trials)
  if (!attr(ready, "ok")) stop("Fix the names listed above, then run this part again.")
  for (i in seq_along(rc_sites)) {
    save_breedbase(rc_designs[[i]], con, program = my_program, folder = my_folder,
                   trial_name = rc_trials[i], year = my_year, location = rc_sites[i],
                   description = "Advanced yield trial, 108 lines and 4 checks",
                   create_folder = TRUE, confirm = FALSE)
  }
}

# ---- A4. Phenotypes ------------------------------------------------------------
# The phenotyping sheets carry the Breedbase plot names of each trial. With
# real data: fill them in and read them back, e.g.
#   sheet <- read.csv("rc_sheet_Ibadan.csv", check.names = FALSE)
rc_sheets <- lapply(seq_along(rc_sites), function(i) {
  sheet <- phenotype_template(rc_designs[[i]], trial_name = rc_trials[i], traits = my_trait)
  # write.csv(sheet, file.path(out_dir, paste0("rc_sheet_", rc_sites[i], ".csv")), row.names = FALSE)
  sheet
})

# Field data for the analysis: the design of every location with its yields.
set.seed(1)
rc_value <- setNames(rnorm(112, sd = 4), c(entries, check_names))        # true genotype values
rc_met <- do.call(rbind, lapply(seq_along(rc_sites), function(i) {
  d <- rc_designs[[i]]$design
  d$location <- rc_sites[i]
  d$yield <- simulate_yield(d, rc_value, mean = c(60, 66, 55)[i], row_trend = 0.8, col_trend = -0.3,
                            sd = c(2, 2, 4)[i])                        # the last location is noisier
  d
}))
rc_met$yield[c(7, 50, 301)] <- NA                                       # three lost plots
rc_met$yield[120] <- rc_met$yield[120] * 10                             # a typing error (outlier)

# The same values in the sheets (what you would upload to Breedbase).
for (i in seq_along(rc_sites)) {
  here <- rc_met[rc_met$location == rc_sites[i], ]
  rc_sheets[[i]][[my_trait]] <- here$yield[match(rc_sheets[[i]]$plot_number, here$plots)]
}

# ---- A5. Plot the trait ------------------------------------------------------------
# Histogram and field map per location: gradients and odd plots stand out.
print(plot_trait(rc_met, "yield"))

# ---- A6. Outliers (IQR rule) -----------------------------------------------------
# Values beyond Q1 - 1.5 IQR or Q3 + 1.5 IQR within each location are set to
# NA (the plots stay in the design). The list shows every plot and value
# removed. Raw values include genetic and field effects, so the rule can also
# flag real extreme lines on a strong gradient: multiplier = 3 keeps only the
# extreme values.
rc_clean <- remove_outliers(rc_met, "yield", multiplier = 1.5)
rc_removed <- attr(rc_clean, "outliers")
plot_trait(rc_met, "yield", type = "field", outliers = rc_removed)      # removed plots circled

# ---- A7. First stage at every location (SpATS), compared with the raw data ----
# The design's own analysis at one location: ANOVA and ranking.
rc_fit <- analyze_augmented_row_column_design(rc_clean[rc_clean$location == rc_sites[1], ], "yield")
print(rc_fit$anova_formatted)
print(head(rc_fit$adjusted_means[, c("rank", "genotype", "raw_mean", "adjusted_mean", "std_error")], 10))

# The first stage of the two-stage analysis: SpATS at every location (genotype
# means with their covariance, and adjusted plot values).
rc_two <- analyze_two_stage(rc_clean, "yield", stage1 = stage1_model)
print(rc_two$locations[, c("location", "stage1_model", "n_plots", "residual_variance", "mean_sed")])
# Raw versus adjusted plot values: field maps before and after, and how much
# the ranking of the genotypes changed.
rc_compare <- compare_first_stage(rc_two)

# ---- A8. Second stage: combine the locations --------------------------------------
print(rc_two)                                # weights, heterogeneity, G x L variance, best genotypes
rc_unweighted <- analyze_two_stage(rc_clean, "yield", stage1 = stage1_model, weighting = "none")
rc_raw_means <- tapply(rc_clean$yield, rc_clean$all_entries, mean, na.rm = TRUE)
print(c(raw_means = cor(rc_raw_means, rc_value[names(rc_raw_means)]),
        two_stage = prediction_accuracy(rc_two, rc_value),
        unweighted = prediction_accuracy(rc_unweighted, rc_value)))
# write.csv(rc_two$adjusted_means, file.path(out_dir, "rc_two_stage_ranking.csv"), row.names = FALSE)

# ---- A9. Results to Breedbase -----------------------------------------------------
# Phenotypes (raw values) on the field plots. First stage (one adjusted value
# per plot), in either or both ways:
#   - as an analysis of the field plots ("Analysis of: plots", like the
#     "Upload analysis" form): a project of type analysis with one plot per
#     field plot;
#   - as its own trait on the field plots ("trait|statistic", created if
#     needed).
# Second stage (one value per genotype) as an analysis of accessions.
# Analyses are created by the web pages of Breedbase. On sites whose web
# pages need their own web password (e.g. cassava-test) and without it, each
# analysis is checked through BrAPI and written as the file of the website's
# "Upload analysis" form (in out_dir), with the values of every form field
# printed: upload it in your browser.
rc_stage1 <- do.call(rbind, lapply(seq_along(rc_sites), function(i) {
  # With Breedbase: check the design against the saved trial (same accession
  # at every plot number) and take the plot names exactly as saved.
  saved_plots <- if (use_breedbase) breedbase_trial_plots(con, rc_trials[i]) else NULL
  first_stage_values(rc_two$stage1_fits[[rc_sites[i]]], rc_designs[[i]], rc_trials[i], trait = my_trait,
                     breedbase_plots = saved_plots)
}))
rc_stage2 <- second_stage_values(rc_two, trait = my_trait)
rc_stage1_name <- paste0(my_trial_prefix, "_", my_year, "_stage1_", my_run_tag)
rc_stage2_name <- paste0(my_trial_prefix, "_", my_year, "_stage2_", my_run_tag)
print(head(rc_stage1))
print(head(rc_stage2))
invisible(save_breedbase_first_stage(rc_stage1, dry_run = TRUE))
invisible(save_breedbase_analysis(rc_stage2, analysis_name = paste0("AYT_", my_year, "_two_stage"),
                                  program = my_program, description = "Two-stage analysis, SpATS stage 1",
                                  dry_run = TRUE))
if (save_phenotypes) {
  for (sheet in rc_sheets) save_breedbase_phenotypes(sheet, con, confirm = FALSE)
}
if (save_analyses && trials_saved(rc_trials)) {
  # First stage: an analysis of the field plots (one value per plot).
  rc_saved1 <- save_breedbase_analysis(rc_stage1, con,
                                       analysis_name = rc_stage1_name,
                                       program = my_program, description = "First stage: SpATS at every location",
                                       protocol = "SpATS(yield ~ genotype + PSANOVA(col, row)) per location",
                                       level = "plot", confirm = FALSE, folder = my_folder, create_folder = TRUE,
                                       upload_file = file.path(out_dir, "rc_stage1_upload.csv"))
  # Second stage: an analysis of the accessions (one value per genotype).
  rc_saved2 <- save_breedbase_analysis(rc_stage2, con,
                                       analysis_name = rc_stage2_name,
                                       program = my_program, description = "Second stage: weighted means across locations",
                                       protocol = "two-stage: SpATS per location; weighted means, G x L random",
                                       confirm = FALSE, folder = my_folder, create_folder = TRUE,
                                       upload_file = file.path(out_dir, "rc_stage2_upload.csv"))
  print(rc_saved1)                                  # analysis id and page
  print(rc_saved2)
  # A11 below analyzes dry matter and adds it to these same analysis projects.
  # Optional: the first-stage values also as their own trait on the field
  # plots ("trait|statistic", created if needed):
  # save_breedbase_first_stage(rc_stage1, con, confirm = FALSE)
}

# ---- A10. Checks ---------------------------------------------------------------
check_result("each check exactly once in every block", all(table(rc_checks$all_entries, rc_checks$block) == 1))
check_result("the typing error was removed as an outlier",
             any(rc_removed$location == rc_met$location[120] & rc_removed$plots == rc_met$plots[120]))
check_result("first stage: one value per field plot", nrow(rc_stage1) == 3 * nrow(rc$design))
check_result("second stage: one value per genotype", nrow(rc_stage2) == 112)
check_result("the noisiest location has the largest residual variance",
             rc_two$locations$location[which.max(rc_two$locations$residual_variance)] == rc_sites[3])
check_result("two-stage means track the true values better than raw means",
             prediction_accuracy(rc_two, rc_value) > cor(rc_raw_means, rc_value[names(rc_raw_means)]))


# ---- A11. Add a second analyzed trait to the existing analyses -----------------
# Run after A1-A10. With real data, read dry matter from the same field books.
# Here it is simulated separately from yield, then analyzed in its own right.
check_section("A11. A second trait in the same analyses")
set.seed(84)
rc_dm_value <- setNames(rnorm(length(rc_value), sd = 1.5), names(rc_value))
rc_dm_met <- rc_met
rc_dm_met$dry_matter <- 32 + rc_dm_value[rc_dm_met$all_entries] +
  c(0, 2, -1)[match(rc_dm_met$location, rc_sites)] + 0.08 * rc_dm_met$row +
  rnorm(nrow(rc_dm_met), sd = c(1, 1.5, 2)[match(rc_dm_met$location, rc_sites)])
rc_dm_met$dry_matter[30] <- 95       # a typing error, reported by the IQR check
plot_trait(rc_dm_met, "dry_matter")
rc_dm_clean <- remove_outliers(rc_dm_met, "dry_matter")
rc_dm_two <- analyze_two_stage(rc_dm_clean, "dry_matter", stage1 = stage1_model)
rc_dm_stage1 <- do.call(rbind, lapply(seq_along(rc_sites), function(i) {
  saved_plots <- if (use_breedbase) breedbase_trial_plots(con, rc_trials[i]) else NULL
  first_stage_values(rc_dm_two$stage1_fits[[rc_sites[i]]], rc_designs[[i]], rc_trials[i],
                     trait = my_second_trait, breedbase_plots = saved_plots)
}))
rc_dm_stage2 <- second_stage_values(rc_dm_two, trait = my_second_trait)
print(head(rc_dm_stage1))
print(head(rc_dm_stage2))

# Store measured dry matter on the same field plots, using their saved names.
rc_dm_sheets <- lapply(seq_along(rc_sites), function(i) {
  sheet <- phenotype_template(rc_designs[[i]], rc_trials[i], traits = my_second_trait)
  here <- rc_dm_met[rc_dm_met$location == rc_sites[i], ]
  sheet[[my_second_trait]] <- here$dry_matter[match(sheet$plot_number, here$plots)]
  sheet
})
if (save_phenotypes) {
  for (sheet in rc_dm_sheets) save_breedbase_phenotypes(sheet, con, confirm = FALSE)
}

# A9 must have saved these analyses, or set these names to existing projects
# with the same plots/accessions. Verify first, then add the dry matter values.
# add_breedbase_analysis_values() requires a connection, even for its dry run.
if (save_analyses && trials_saved(rc_trials)) {
  add_breedbase_analysis_values(rc_dm_stage1, con, analysis_name = rc_stage1_name,
                                level = "plot", dry_run = TRUE)
  rc_dm_added1 <- add_breedbase_analysis_values(rc_dm_stage1, con,
    analysis_name = rc_stage1_name, level = "plot", confirm = FALSE)
  add_breedbase_analysis_values(rc_dm_stage2, con, analysis_name = rc_stage2_name,
                                level = "accession", dry_run = TRUE)
  rc_dm_added2 <- add_breedbase_analysis_values(rc_dm_stage2, con,
    analysis_name = rc_stage2_name, level = "accession", confirm = FALSE)
  print(rc_dm_added1)
  print(rc_dm_added2)
}
check_result("second trait uses exactly the same field plot names",
             setequal(rc_dm_stage1$plot_name, rc_stage1$plot_name))
check_result("second trait uses the same accessions in the means analysis",
             setequal(rc_dm_stage2$accession, rc_stage2$accession))
check_result("the dry matter typing error was excluded",
             is.na(rc_dm_clean$dry_matter[30]))


# ---- A12. Alternative: single-stage analysis of all locations ------------------
# Run A1-A6, then jump here to skip the first and second stages entirely.
# This section needs rc_clean and rc_designs, not rc_two or its stage-1 fits.
check_section("A12. Single-stage analysis")
rc_single <- analyze_single_stage(
  rc_clean, "yield", genotype_col = "all_entries", location_col = "location",
  spatial = has_spats
)
print(rc_single)
print(rc_single$variance_components)

# Despite its name, second_stage_values() also formats a single-stage fit's
# genotype means for upload. It does not perform another analysis.
rc_single_values <- second_stage_values(rc_single, trait = my_trait)
rc_single_name <- paste0(my_trial_prefix, "_", my_year, "_single_stage_", my_run_tag)
rc_single_payload <- save_breedbase_analysis(
  rc_single_values, analysis_name = rc_single_name, program = my_program,
  description = "Single-stage analysis of all locations together", year = my_year,
  estimate = if (has_spats) "spatial" else "blue", dry_run = TRUE
)
if (save_analyses && trials_saved(rc_trials)) {
  rc_single_saved <- save_breedbase_analysis(
    rc_single_values, con, analysis_name = rc_single_name, program = my_program,
    description = "Single-stage analysis of all locations together", year = my_year,
    protocol = "All plot observations in one model, with genotype x location random",
    estimate = if (has_spats) "spatial" else "blue",
    folder = my_folder, create_folder = TRUE, confirm = FALSE,
    upload_file = file.path(out_dir, "rc_single_stage_upload.csv")
  )
  print(rc_single_saved)
}
check_result("single stage: one result per genotype",
             nrow(rc_single_values) == length(unique(rc_clean$all_entries)))
check_result("single stage: all locations in the same model",
             rc_single$diagnostics$n_locations == length(rc_sites))
check_result("single-stage standard errors are finite",
             all(is.finite(rc_single$adjusted_means$std_error)))
# This model shares one residual variance across locations. The weighted
# two-stage route above allows each location to have its own precision.


# =============================================================================
# ---- B. Flexible augmented row-column layouts (Piepho and Williams, 2016) ----
# =============================================================================
check_section("B. Piepho and Williams (2016) layouts")
# Row and column groups can have different sizes, and blocks can hold fewer
# or more check plots than there are checks. Checks are then placed by a
# D-optimal search.

# ---- B1. Plan: Table 2 of the paper ----------------------------------------------
pw_ex1 <- augmented_row_column_capacity(6, 9, 3, 3, n_controls = 3)
pw_ex2 <- augmented_row_column_capacity(36, 18, 6, 6, n_controls = 4)
checks_per_block <- rbind(c(4, 5, 4),     # row group 1: 4, 5 and 4 check plots
                          c(3, 5, 3),     # row group 2
                          c(3, 5, 4),     # row group 3
                          c(3, 5, 3))     # row group 4
pw_ex3 <- augmented_row_column_capacity(n_controls = 4, rows_per_rowgroup = c(6, 5, 6, 5),
                                        cols_per_colgroup = c(3, 4, 3), checks_per_block = checks_per_block)
print(rbind(pw_ex1, pw_ex2, pw_ex3)[, c("total_plots", "check_plots", "entry_plots", "error_df")])

# ---- B2. Generate Example 3 (Figure 4) and Example 2 (Figure 3) ------------------------
pw <- augmented_row_column_design(
  treatments        = paste0("E", 1:173),
  controls          = c("A", "B", "C", "D"),
  rows_per_rowgroup = c(6, 5, 6, 5),
  cols_per_colgroup = c(3, 4, 3),
  checks_per_block  = checks_per_block,
  n_candidates      = 3,
  seed              = 1
)
print(pw$summary[, c("method", "check_plots", "entry_plots", "error_df", "d_criterion")])
check_design(pw)
plot_augmented_row_column_design(pw, label = "type", fill = "type")

fig3 <- augmented_row_column_design(paste0("E", 1:576), c("A", "B", "C", "D"), 36, 18, 6, 6,
                                    n_candidates = 200, seed = 3)
check_design(fig3)
plot_augmented_row_column_design(fig3, label = "type", fill = "type")
fig3_checks <- fig3$design[fig3$design$is_control == 1, ]
print(table(fig3_checks$all_entries, fig3_checks$col))     # each check once in every column

# ---- B3. Field data, plot, outliers ------------------------------------------------
set.seed(2)
pw_value <- setNames(rnorm(177, sd = 3), c(paste0("E", 1:173), LETTERS[1:4]))
pw_trial <- pw$design
pw_trial$yield <- simulate_yield(pw_trial, pw_value, mean = 40, row_trend = 0.5, sd = 1)
plot_trait(pw_trial, "yield")
pw_clean <- remove_outliers(pw_trial, "yield")

# ---- B4. First stage: design model and SpATS, compared with the raw data --------------
pw_fit <- analyze_augmented_row_column_design(pw_clean, response = "yield")
print(pw_fit$anova_formatted)                       # 7 error df, as in Table 2 of the paper
print(head(pw_fit$adjusted_means[, c("rank", "genotype", "adjusted_mean", "std_error")], 10))
if (has_spats) {
  pw_spatial <- analyze_augmented_row_column_design(pw_clean, "yield", spatial = TRUE)
  pw_compare <- compare_first_stage(pw_spatial)        # raw vs spatially adjusted
  print(c(design_model = prediction_accuracy(pw_fit, pw_value),
          spatial = prediction_accuracy(pw_spatial, pw_value)))
}
# One location: the first stage is the final analysis (no second stage).

# ---- B5. Checks against the paper ---------------------------------------------------
check_result("Table 2: Example 1 has 36 entries and 0 error df (cannot be used)",
             pw_ex1$entry_plots == 36 && pw_ex1$error_df == 0)
check_result("Table 2: Example 2 has 576 entries and 6 error df", pw_ex2$entry_plots == 576 && pw_ex2$error_df == 6)
check_result("Table 2: Example 3 has 173 entries and 7 error df", pw_ex3$entry_plots == 173 && pw_ex3$error_df == 7)
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
             inherits(try(augmented_row_column_design(paste0("E", 1:36), LETTERS[1:3], 6, 9, 3, 3), silent = TRUE),
                      "try-error"))


# =============================================================================
# ---- C. Augmented p-rep design ----------------------------------------------
# =============================================================================
check_section("C. Augmented p-rep design")
# A multi-location trial without checks: every entry is planted once at
# every location, and 20% of the entries get a second plot at each location;
# the repeated entries take the place of checks.

# ---- C1. Plan and generate ----------------------------------------------------------
prep_entries <- sprintf("G%03d", 1:120)
prep_sites   <- c("North", "Central", "South")
prep_trial_name <- paste0("PREP_", my_year)            # Breedbase trials PREP_2026_North, ...

prep <- augmented_prep_design(
  treatments   = prep_entries,
  locations    = prep_sites,
  block_size   = 10,           # at most 10 plots per incomplete block
  prep_prop    = 0.2,          # 20% of the entries repeated at each location
  field_cols   = 10,           # field columns: row/column coordinates for SpATS
  n_candidates = 200,
  seed         = 2026
)

# ---- C2. Look at it and save the field book -------------------------------------------
print(prep$summary[, c("n_plots", "n_blocks", "max_pairwise_concurrence", "n_disconnected_locations",
                       "efficiency_factor")])
print(table(table(prep$design$all_entries)))     # 3 plots for most entries, 4 for the repeated ones
check_design(prep)
plot_augmented_prep_design(prep, locations = "North")     # orange: repeated entries
# write.csv(prep$design, file.path(out_dir, "prep_fieldbook.csv"), row.names = FALSE)

# ---- C3. Save the trials to Breedbase (one trial per location) -----------------------------
prep_map <- setNames(my_locations, prep_sites)       # design location -> Breedbase location
invisible(save_breedbase(prep, trial_name = prep_trial_name, year = my_year, location = prep_map,
                         dry_run = TRUE, description = "p-rep trial of 120 entries in three locations"))
if (save_designs) {
  save_breedbase(prep, con, program = my_program, folder = my_folder, trial_name = prep_trial_name,
                 year = my_year, location = prep_map, create_folder = TRUE, confirm = FALSE,
                 description = "p-rep trial of 120 entries in three locations")
}

# ---- C4. Phenotypes -------------------------------------------------------------------
set.seed(3)
prep_value <- setNames(rnorm(120, sd = 4), prep_entries)
block_id <- interaction(prep$design$location, prep$design$block)
block_effect <- setNames(rnorm(nlevels(block_id), sd = 2), levels(block_id))
prep_met <- prep$design
prep_met$yield <- simulate_yield(prep_met, prep_value, mean = 50, row_trend = 0.3, sd = 3) +
  c(North = 0, Central = 8, South = -5)[prep_met$location] + block_effect[as.character(block_id)]
south <- prep_met$location == "South"
prep_met$yield[south] <- prep_met$yield[south] + rnorm(sum(south), sd = 9)     # a noisy location
prep_met$yield[15] <- prep_met$yield[15] / 10                                  # a typing error

# The sheet lists the plots by location and plot number, as in Breedbase.
prep_sheet <- phenotype_template(prep, trial_name = prep_trial_name, traits = my_trait)
in_sheet_order <- prep_met[order(match(prep_met$location, prep_sites), prep_met$plots), ]
prep_sheet[[my_trait]] <- in_sheet_order$yield

# ---- C5. Plot the trait -------------------------------------------------------------------
print(plot_trait(prep_met, "yield"))

# ---- C6. Outliers --------------------------------------------------------------------------
prep_clean <- remove_outliers(prep_met, "yield")
prep_removed <- attr(prep_clean, "outliers")

# ---- C7. First stage at every location (SpATS), compared with the raw data ----------------
# One-stage analysis of all locations (ANOVA and ranking), for reference.
prep_fit <- analyze_augmented_prep_design(prep_clean, response = "yield")
print(prep_fit$anova_formatted)
# Two-stage: SpATS (or the p-rep block model) at every location ...
prep_stage1_model <- if (has_spats) "spatial" else "blocks"
prep_two <- analyze_augmented_prep_two_stage(prep_clean, response = "yield", stage1 = prep_stage1_model)
print(prep_two$locations[, c("location", "n_plots", "residual_variance", "residual_df", "mean_sed")])
prep_compare <- compare_first_stage(prep_two)

# ---- C8. Second stage: combine the locations --------------------------------------------------
print(prep_two)                               # Hartley test, G x L variance, best genotypes
prep_unweighted <- analyze_augmented_prep_two_stage(prep_clean, "yield", stage1 = prep_stage1_model,
                                                    weighting = "none")
print(c(one_stage = prediction_accuracy(prep_fit, prep_value),
        two_stage = prediction_accuracy(prep_two, prep_value),
        unweighted = prediction_accuracy(prep_unweighted, prep_value)))
if (has_lme4) {
  prep_gxl <- analyze_augmented_prep_design(prep_clean, "yield", genotype_by_location = "random")
  print(prep_gxl$variance_components)        # one-stage genotype x location variance
}

# ---- C9. Results to Breedbase -------------------------------------------------------------------
prep_stage1 <- first_stage_values(prep_two, prep, trial_name = prep_trial_name, trait = my_trait)
prep_stage2 <- second_stage_values(prep_two, trait = my_trait)
invisible(save_breedbase_first_stage(prep_stage1, dry_run = TRUE))
if (save_phenotypes) save_breedbase_phenotypes(prep_sheet, con, confirm = FALSE)
if (save_analyses && trials_saved(paste0(prep_trial_name, "_", prep_sites))) {
  save_breedbase_analysis(prep_stage1, con, analysis_name = paste0(prep_trial_name, "_stage1_", my_run_tag),
                          program = my_program, description = "p-rep first stage per location", level = "plot",
                          confirm = FALSE, folder = my_folder, create_folder = TRUE, upload_file = file.path(out_dir, "prep_stage1_upload.csv"))
  save_breedbase_analysis(prep_stage2, con, analysis_name = paste0(prep_trial_name, "_stage2_", my_run_tag),
                          program = my_program, description = "p-rep two-stage analysis", confirm = FALSE, folder = my_folder, create_folder = TRUE,
                          upload_file = file.path(out_dir, "prep_stage2_upload.csv"))
}

# ---- C10. Checks ------------------------------------------------------------------------------
check_result("every location can be analysed on its own (all connected)",
             prep$summary$n_disconnected_locations == 0)
check_result("the typing error was removed as an outlier",
             any(prep_removed$location == prep_met$location[15] & prep_removed$plots == prep_met$plots[15]))
check_result("first stage: one value per field plot", nrow(prep_stage1) == nrow(prep$design))
check_result("the sheet holds the same values as the field data",
             isTRUE(all.equal(sort(prep_sheet[[my_trait]]), sort(prep_met$yield))))
check_result("Hartley's test flags the noisy location", isTRUE(prep_two$heterogeneity$heterogeneous))
check_result("weighting by precision beats the unweighted analysis",
             prediction_accuracy(prep_two, prep_value) > prediction_accuracy(prep_unweighted, prep_value))


# =============================================================================
# ---- D. p-rep designs from alpha-arrays (Williams, Piepho and Whitaker, 2011) -
# =============================================================================
check_section("D. p-rep designs from alpha-arrays")
# Every entry repeated at exactly one location, two replicates per location,
# and every block holds m repeated entries.

# ---- D1. Plan and generate ----------------------------------------------------------
print(augmented_prep_parameters(n_treatments = 42, n_locations = 3, max_block_size = 6))   # Table 3
alpha_ex1 <- augmented_prep_design(paste0("G", 1:20), c("L1", "L2"), method = "alpha", block_size = 3, seed = 1)
alpha_ex3 <- augmented_prep_design(paste0("G", 1:42), c("L1", "L2", "L3"), method = "alpha2",
                                   block_size = 4, field_cols = 7, seed = 1)
print(rbind(alpha_ex1$summary, alpha_ex3$summary)[, c("method", "m", "s", "n_plots", "n_blocks",
                                                      "max_pairwise_concurrence", "efficiency_factor")])

# ---- D2. Look at it ---------------------------------------------------------------------
plot_augmented_prep_design(alpha_ex3, locations = "L1")   # blue line: the two replicates
check_design(alpha_ex1)
check_design(alpha_ex3)

# ---- D3. Field data, plot, outliers --------------------------------------------------------
set.seed(7)
alpha_value <- setNames(rnorm(42, sd = 3), paste0("G", 1:42))
alpha_met <- alpha_ex3$design
alpha_met$yield <- simulate_yield(alpha_met, alpha_value, mean = 30, row_trend = 0.4, sd = 1) +
  c(L1 = 0, L2 = 4, L3 = -2)[alpha_met$location]
plot_trait(alpha_met, "yield", type = "histogram")
alpha_clean <- remove_outliers(alpha_met, "yield")

# ---- D4. First stage per location, compared with the raw data; second stage ----------------
alpha_fit <- analyze_augmented_prep_design(alpha_clean, response = "yield")    # one stage, for reference
print(alpha_fit$anova_formatted)                   # replicates within locations, blocks within replicates
alpha_two <- analyze_two_stage(alpha_clean, "yield", stage1 = stage1_model)
alpha_compare <- compare_first_stage(alpha_two, plot = FALSE)
print(alpha_two)
print(c(one_stage = prediction_accuracy(alpha_fit, alpha_value), two_stage = prediction_accuracy(alpha_two, alpha_value)))

# ---- D5. Checks against the paper -------------------------------------------------------------
check_result("Example 1: m = 2, 5 blocks per replicate, blocks of 3",
             alpha_ex1$summary$m == 2 && alpha_ex1$summary$s == 5 && all(alpha_ex1$design$block_size == 3))
check_result("Example 3: 14 blocks of 4 plots at each location",
             alpha_ex3$summary$n_blocks == 42 && all(alpha_ex3$design$block_size == 4))
check_result("every entry repeated at exactly one location", all(table(alpha_ex3$design$all_entries) == 4))
check_result("no pair of entries shares a block more than once",
             alpha_ex1$summary$max_pairwise_concurrence == 1 && alpha_ex3$summary$max_pairwise_concurrence == 1)
check_result("the ANOVA has replicates within locations",
             "Replicates, nested within locations" %in% alpha_fit$anova_formatted$Source)
check_result("adjusted means track the true values (r > 0.8)", prediction_accuracy(alpha_two, alpha_value) > 0.8)


# =============================================================================
# ---- E. Resolvable row-column design (Williams, John and Whitaker, 2006) -----
# =============================================================================
check_section("E. Resolvable row-column design")
# When every entry can be replicated: each replicate is a grid of rows and
# columns that holds every treatment once. Grown at two locations.

# ---- E1. Plan and generate --------------------------------------------------------------
rrc_entries <- paste0("T", 1:42)
rrc_sites   <- my_locations[1:2]
rrc_trials  <- paste0("RRC_", my_year, "_", rrc_sites)
rrc_designs <- lapply(seq_along(rrc_sites), function(i) {
  resolvable_row_column_design(
    treatments   = rrc_entries,
    n_reps       = 3,
    rows_per_rep = 6,
    cols_per_rep = 7,
    latinize     = TRUE,        # no treatment twice in a long column
    n_starts     = 1,
    seed         = i
  )
})
names(rrc_designs) <- rrc_sites
rrc <- rrc_designs[[1]]

# ---- E2. Look at it, and compare with the published design ----------------------------------
print(rrc$summary[, c("efficiency_rc", "efficiency_rows", "efficiency_cols", "efficiency_spatial",
                      "upper_bound_spatial", "latin_ok")])
plot_resolvable_row_column_design(rrc)
check_design(rrc)
# write.csv(rrc$design, file.path(out_dir, "rrc_fieldbook.csv"), row.names = FALSE)

# Table 1 of the paper (made with CycDesigN): each block of 6 rows is one replicate.
table1 <- rbind(
  c(19, 34, 27,  3, 15, 12,  2), c(40,  4,  7, 20, 21, 33, 18), c(23, 35, 17, 29, 22, 28, 38),
  c(10, 14, 31, 39, 32,  5, 42), c(13,  6, 16, 36, 11, 37, 41), c(25,  8,  1, 26,  9, 30, 24),
  c(38, 34, 26, 33,  7, 10, 11), c( 1, 13, 15, 35,  5,  3, 18), c(32, 20, 17, 41,  8, 24, 19),
  c(12, 30,  6, 39, 23, 21, 16), c( 4, 31, 42, 25,  2, 37, 29), c(36, 22, 40, 27,  9, 14, 28),
  c(20, 26, 23, 37, 32, 27, 35), c(10, 22, 16, 18, 42,  8, 12), c(28,  4, 15, 39,  7, 13, 24),
  c( 1, 41, 14, 34, 30, 29, 40), c( 2,  5, 36, 17, 25, 21, 11), c( 6, 19, 33,  9,  3, 38, 31)
)
paper_design <- data.frame(replicate = rep(rep(1:3, each = 6), times = 7), row_in_rep = rep(rep(1:6, 3), times = 7),
                           col_in_rep = rep(1:7, each = 18), all_entries = paste0("T", as.vector(table1)))
paper_eff <- resolvable_row_column_efficiency(paper_design)
print(unlist(paper_eff[c("efficiency_rc", "efficiency_spatial", "upper_bound_spatial")]))

# ---- E3. Save the trials to Breedbase ----------------------------------------------------------
invisible(save_breedbase(rrc, trial_name = rrc_trials[1], year = my_year, location = rrc_sites[1],
                         dry_run = TRUE, description = "Resolvable row-column trial, 42 lines"))
if (save_designs) {
  for (i in seq_along(rrc_sites)) {
    save_breedbase(rrc_designs[[i]], con, program = my_program, folder = my_folder, trial_name = rrc_trials[i],
                   year = my_year, location = rrc_sites[i], description = "Resolvable row-column trial, 42 lines",
                   create_folder = TRUE, confirm = FALSE)
  }
}

# ---- E4. Phenotypes ---------------------------------------------------------------------------
set.seed(5)
rrc_value <- setNames(rnorm(42, sd = 2), rrc_entries)
rrc_met <- do.call(rbind, lapply(seq_along(rrc_sites), function(i) {
  d <- rrc_designs[[i]]$design
  d$location <- rrc_sites[i]
  d$yield <- simulate_yield(d, rrc_value, mean = c(50, 45)[i], row_trend = 0.15, sd = 1) +
    rnorm(18, sd = 1.2)[(d$replicate - 1) * 6 + d$row_in_rep] +      # row effects within replicates
    rnorm(21, sd = 0.8)[(d$replicate - 1) * 7 + d$col_in_rep]        # column effects within replicates
  d
}))
rrc_sheets <- lapply(seq_along(rrc_sites), function(i) {
  sheet <- phenotype_template(rrc_designs[[i]], trial_name = rrc_trials[i], location = rrc_sites[i],
                              traits = my_trait)
  here <- rrc_met[rrc_met$location == rrc_sites[i], ]
  sheet[[my_trait]] <- here$yield[match(sheet$plot_number, here$plots)]
  sheet
})

# ---- E5. Plot the trait, and outliers ------------------------------------------------------------
print(plot_trait(rrc_met, "yield"))
rrc_clean <- remove_outliers(rrc_met, "yield")

# ---- E6. First stage: variance models at one location; SpATS at every location -----------------
rrc_fit <- analyze_resolvable_row_column(rrc_clean[rrc_clean$location == rrc_sites[1], ], response = "yield",
                                         model = "rc", long_col = "long_col")
print(rrc_fit$anova)
print(rrc_fit$model_comparison[, c("model", "deviance", "AIC", "lrt_p_value")])   # RC, LV, AR models
print(head(rrc_fit$adjusted_means, 10))
rrc_two <- analyze_two_stage(rrc_clean, "yield", stage1 = stage1_model)
rrc_compare <- compare_first_stage(rrc_two)

# ---- E7. Second stage: combine the locations -------------------------------------------------------
print(rrc_two)
rrc_design_two <- analyze_two_stage(rrc_clean, "yield", stage1 = "design")    # RC model in stage 1
print(c(rc_one_location = prediction_accuracy(rrc_fit, rrc_value),
        two_stage_design = prediction_accuracy(rrc_design_two, rrc_value),
        two_stage = prediction_accuracy(rrc_two, rrc_value)))

# ---- E8. Results to Breedbase -----------------------------------------------------------------------
rrc_stage1 <- do.call(rbind, lapply(seq_along(rrc_sites), function(i) {
  first_stage_values(rrc_two$stage1_fits[[rrc_sites[i]]], rrc_designs[[i]], rrc_trials[i], trait = my_trait,
                     location = rrc_sites[i])
}))
rrc_stage2 <- second_stage_values(rrc_two, trait = my_trait)
if (save_phenotypes) for (sheet in rrc_sheets) save_breedbase_phenotypes(sheet, con, confirm = FALSE)
if (save_analyses && trials_saved(rrc_trials)) {
  save_breedbase_analysis(rrc_stage1, con, analysis_name = paste0("RRC_", my_year, "_stage1_", my_run_tag),
                          program = my_program, description = "Resolvable row-column first stage", level = "plot",
                          confirm = FALSE, folder = my_folder, create_folder = TRUE, upload_file = file.path(out_dir, "rrc_stage1_upload.csv"))
  save_breedbase_analysis(rrc_stage2, con, analysis_name = paste0("RRC_", my_year, "_stage2_", my_run_tag),
                          program = my_program, description = "Resolvable row-column, two-stage", confirm = FALSE, folder = my_folder, create_folder = TRUE,
                          upload_file = file.path(out_dir, "rrc_stage2_upload.csv"))
}

# ---- E9. Checks against the paper ----------------------------------------------------------------------
check_result("Table 1 design: row-column efficiency 0.6941 (as published)", round(paper_eff$efficiency_rc, 4) == 0.6941)
check_result("Table 1 design: spatial efficiency 0.9204, upper bound 0.9268 (as published)",
             round(paper_eff$efficiency_spatial, 4) == 0.9204 && round(paper_eff$upper_bound_spatial, 4) == 0.9268)
check_result("generated design within 2% of the paper's row-column efficiency",
             rrc$summary$efficiency_rc >= 0.98 * 0.6941)
check_result("first stage: one value per field plot", nrow(rrc_stage1) == 2 * nrow(rrc$design))
check_result("two-stage means track the true values (r > 0.8)", prediction_accuracy(rrc_two, rrc_value) > 0.8)


# =============================================================================
# ---- F. Multi-environment trials from Breedbase files -------------------------
# =============================================================================
check_section("F. Multi-environment trials")
# Before combining trials: are their error variances similar, and which
# environments rank genotypes alike? Example file with four trials.

pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
locations_file <- system.file("extdata", "breedbase_locations.csv", package = "breedtools")
pheno_raw <- read_breedbase_file(pheno_file)        # also reads .xlsx and .txt files
print(head(pheno_raw))

# F1. Hartley's Fmax test of the error variances, and compatible groups.
variance_check <- run_breedbase_variance_check(pheno_file, trait_col = "fresh root yield|CO_334:0000013",
                                               n_sim = 20000)
print(variance_check$env_anova[, c("environment", "ms_error", "mean", "cv_percent")])
print(variance_check$overall_hartley)
print(variance_check$maximal_compatible_groups)

# F2. Group environments by genotype-by-environment interaction (Lin, 1982).
strata <- run_environment_stratification(pheno_file, trait = "fresh root yield")
print(strata$group_summary[, c("group_id", "locations", "n_env", "f_value", "p_value")])

# F3. The groups on a map (opens in the RStudio viewer).
locations <- read_locations(locations_file, phenotype_data = strata)
print(map_environment_groups(strata, locations))

# F4. Checks.
check_result("four environments read from the Breedbase file", nrow(variance_check$env_anova) == 4)
check_result("Fmax = largest / smallest error mean square",
             isTRUE(all.equal(variance_check$overall_hartley$fmax_observed,
                              max(variance_check$env_anova$ms_error) / min(variance_check$env_anova$ms_error))))
check_result("two environment groups: Geneva + Ithaca and Ibadan + Ubiaja",
             setequal(strata$group_summary$locations, c("Geneva, Ithaca", "Ibadan, Ubiaja")))


# =============================================================================
# ---- H. From Breedbase trials to analyses -------------------------------------
# =============================================================================
check_section("H. From Breedbase trials to analyses")
# Set read_breedbase_trials <- TRUE in Setup to run this part. Choose
# breedbase_analysis_route = "single", "two_stage" or "both". The single
# route calls analyze_single_stage() directly on the cleaned field data.
# save_analyses controls writing independently of reading and analyzing.

my_trials <- paste0(my_trial_prefix, "_", my_year, "_", my_locations)
my_traits <- c(my_trait, my_second_trait)       # read and analyze both traits
h_route <- match.arg(breedbase_analysis_route, c("single", "two_stage", "both"))
h_prefix <- paste0(my_trial_prefix, "_", my_year, "_readback_", my_run_tag)
h_stage1_name <- paste0(h_prefix, "_stage1")
h_stage2_name <- paste0(h_prefix, "_stage2")
h_single_name <- paste0(h_prefix, "_single_stage")

if (read_breedbase_trials) {
  # H1. Read the trials: one row per plot, one column per trait.
  pheno <- breedbase_phenotypes(con, trials = my_trials, traits = my_traits)
  print(head(pheno))
  # Include all saved stocks when creating an analysis, even if the first
  # trait is missing on some plots or accessions. Later traits can use them.
  h_plots <- pheno[, c("plot_name", "accession")]
  h_accessions <- data.frame(accession = unique(pheno$accession))

  for (trait_index in seq_along(my_traits)) {
    tr <- my_traits[trait_index]
    # H2. Look at the trait, and remove outliers within each trial (listed).
    plot_trait(pheno, tr, by = "trial_name")
    clean <- remove_outliers(pheno, tr, by = "trial_name")

    if (h_route %in% c("two_stage", "both")) {
      # H3/H4. First stage per trial, then precision-weighted second stage.
      two <- analyze_two_stage(clean, tr, stage1 = if (has_spats) "spatial" else "blocks",
                               genotype_col = "accession", location_col = "trial_name")
      compare_first_stage(two)
      print(two)
      first <- first_stage_values(two, trait = tr)
      h_first <- h_plots
      h_first[[tr]] <- first[[tr]][match(h_first$plot_name, first$plot_name)]
      means <- second_stage_values(two, trait = tr)
      h_second <- h_accessions
      h_second[[tr]] <- means[[tr]][match(h_second$accession, means$accession)]
      if (save_analyses) {
        if (trait_index == 1L) {
          save_breedbase_analysis(h_first, con, analysis_name = h_stage1_name,
            program = my_program, description = "First-stage values per trial", year = my_year,
            level = "plot", folder = my_folder, create_folder = TRUE, confirm = FALSE)
          save_breedbase_analysis(h_second, con, analysis_name = h_stage2_name,
            program = my_program, description = "Weighted two-stage means", year = my_year,
            folder = my_folder, create_folder = TRUE, confirm = FALSE)
        } else {
          add_breedbase_analysis_values(h_first, con, analysis_name = h_stage1_name,
                                        level = "plot", dry_run = TRUE)
          add_breedbase_analysis_values(h_first, con, analysis_name = h_stage1_name,
                                        level = "plot", confirm = FALSE)
          add_breedbase_analysis_values(h_second, con, analysis_name = h_stage2_name,
                                        level = "accession", dry_run = TRUE)
          add_breedbase_analysis_values(h_second, con, analysis_name = h_stage2_name,
                                        level = "accession", confirm = FALSE)
        }
      }
    }

    if (h_route %in% c("single", "both")) {
      # H5. One fit of all cleaned plots; independent of the two-stage route.
      one <- analyze_single_stage(clean, tr, spatial = has_spats)
      print(one)
      means <- second_stage_values(one, trait = tr)    # format means, no stage-2 fit
      h_single <- h_accessions
      h_single[[tr]] <- means[[tr]][match(h_single$accession, means$accession)]
      check_result(paste("single stage: one mean per accession for", tr),
                   !anyDuplicated(h_single$accession))
      if (save_analyses) {
        if (trait_index == 1L) {
          save_breedbase_analysis(h_single, con, analysis_name = h_single_name,
            program = my_program, description = "Single-stage model of all trials", year = my_year,
            estimate = if (has_spats) "spatial" else "blue",
            folder = my_folder, create_folder = TRUE, confirm = FALSE)
        } else {
          # A second analyzed trait goes into the SAME single-stage analysis.
          add_breedbase_analysis_values(h_single, con, analysis_name = h_single_name,
                                        level = "accession", dry_run = TRUE)
          add_breedbase_analysis_values(h_single, con, analysis_name = h_single_name,
                                        level = "accession", confirm = FALSE)
        }
      }
    }
    if (h_route == "both") {
      both <- merge(two$adjusted_means, one$adjusted_means, by = "genotype")
      cat("Two-stage vs single-stage means: r =",
          round(cor(both$adjusted_mean.x, both$adjusted_mean.y, use = "complete.obs"), 3), "\n")
    }
  }
} else {
  message("Part H skipped: set read_breedbase_trials <- TRUE in Setup to read and analyze saved trials.")
}


# =============================================================================
# ---- I. Agricolae designs: field books, Breedbase and analysis -----------------
# =============================================================================
check_section("I. Agricolae to Breedbase")
# Three steps: (1) generate a design with agricolae, (2) prepare its field
# book AND trial information with prepare_breedbase_trial(), (3) save the
# prepared object with save_breedbase(). RCBD means randomized complete
# block design. The remaining steps demonstrate phenotypes and analysis.
# Install once: install.packages("agricolae").
run_agricolae_analysis <- FALSE   # TRUE to run the optional I4-I10 examples

if (requireNamespace("agricolae", quietly = TRUE)) {
  # I1. Generate your design with agricolae (or use one you already have).
  ag_kind <- "RCBD"              # "RCBD" or "Alpha"
  ag_kind <- match.arg(ag_kind, c("RCBD", "Alpha"))
  ag_accessions <- c(sprintf("LINE_%03d", 1:22), "CHECK_A", "CHECK_B")
  ag_sites <- my_locations
  ag_trials <- paste0("AGRICOLAE_", ag_kind, "_", my_year, "_", ag_sites)
  ag_originals <- lapply(seq_along(ag_sites), function(i) {
    original <- switch(ag_kind,
      RCBD = agricolae::design.rcbd(ag_accessions, r = 3, seed = 2026 + i, kinds = RNGkind()[1L]),
      Alpha = agricolae::design.alpha(ag_accessions, k = 4, r = 3, seed = 2026 + i, kinds = RNGkind()[1L])
    )
    # Plan a field with six columns in field-book row order. These are
    # explicit physical positions, not inferred from agricolae's plot codes
    # or its within-block `cols` column. Adjust them to your planting plan.
    position <- seq_len(nrow(original$book)) - 1L
    original$book$row <- position %/% 6L + 1L
    original$book$col <- position %% 6L + 1L
    original
  })

  # I2. Convert and record the information Breedbase needs for each trial.
  # For one existing RCBD this is simply:
  #   trial <- prepare_breedbase_trial(your_rcbd,
  #     trial_name = "RCBD_2026_Ibadan", year = 2026, location = "Ibadan",
  #     program = "IITA", description = "RCBD with three complete blocks",
  #     checks = your_check_names, prompt = FALSE)
  # In RStudio, prepare_breedbase_trial(your_rcbd) asks for missing trial
  # name, year, location, breeding program and description.
  ag_designs <- lapply(seq_along(ag_sites), function(i) {
    prepare_breedbase_trial(
      design      = ag_originals[[i]],
      trial_name  = ag_trials[i],
      year        = my_year,
      location    = ag_sites[i],
      program     = my_program,
      folder      = my_folder,
      description = paste(ag_kind, "design generated with agricolae, three replicates"),
      checks      = c("CHECK_A", "CHECK_B"),
      prompt      = FALSE
    )
  })
  names(ag_designs) <- ag_sites
  ag_trials <- unname(vapply(ag_designs, function(x) x$trial$trial_name, character(1)))

  ag <- ag_designs[[1]]
  print(ag$trial)                  # stored program, location, trial name, year, description
  print(head(ag$book))             # agricolae format
  print(head(ag$design))           # standardized for the Breedbase uploader
  ag_from_book <- as_breedbase_design(ag$book, accession = "ag_accessions",
    checks = c("CHECK_A", "CHECK_B"), design_type = ag$design_type)
  plot(ag$design$col, ag$design$row, type = "n", asp = 1,
       xlab = "Field column", ylab = "Field row", main = ag_trials[1])
  text(ag$design$col, ag$design$row, labels = ag$design$all_entries, cex = 0.55)
  for (s in ag_sites) {
    # write.csv(ag_designs[[s]]$book, file.path(out_dir, paste0("agricolae_original_", s, ".csv")), row.names = FALSE)
    # write.csv(ag_designs[[s]]$design, file.path(out_dir, paste0("agricolae_converted_", s, ".csv")), row.names = FALSE)
  }

  # I3. Save through the EXISTING function. Trial information is already
  # recorded: there is no need to enter program/location/name/year again.
  ag_payload <- save_breedbase(ag, dry_run = TRUE)   # preview offline
  str(ag_payload$plots[[1]][[1]])
  if (save_designs) {
    for (i in seq_along(ag_sites)) {
      save_breedbase(ag_designs[[i]], con, create_folder = TRUE, dry_run = TRUE)
      save_breedbase(ag_designs[[i]], con, create_folder = TRUE, confirm = FALSE)
    }
  }

  check_result("agricolae: original and bare books give the same converted design",
               identical(ag$design, ag_from_book$design))
  check_result("agricolae: plot numbers preserve the generated randomization",
               identical(as.integer(ag$book$plots), ag$design$plots))
  check_result("agricolae: one uploaded plot per field-book row",
               length(ag_payload$plots[[1]]) == nrow(ag$book))
  check_result("agricolae: trial information is recorded for every location",
               identical(unname(vapply(ag_designs, function(x) x$trial$trial_name, character(1))), ag_trials))
} else {
  message("Part I skipped: install.packages(\"agricolae\") to generate and prepare its trials.")
}


# ---- I4. Optional phenotypes and analysis after preparing the trials ---------
# I1-I3 above complete the conversion and save workflow. Set
# run_agricolae_analysis <- TRUE above to run this separate post-planting
# example using simulated field data. With real data, read your filled sheets.
if (isTRUE(run_agricolae_analysis) && requireNamespace("agricolae", quietly = TRUE)) {
  # I4. Phenotype sheets with exactly the plot names used by the uploader.
  # With real data, fill the sheets in the field and read them back.
  ag_sheets <- lapply(seq_along(ag_sites), function(i) {
    phenotype_template(ag_designs[[i]], traits = my_trait)
  })
  set.seed(85)
  ag_value <- setNames(rnorm(length(ag_accessions), sd = 3), ag_accessions)
  ag_met <- do.call(rbind, lapply(seq_along(ag_sites), function(i) {
    d <- ag_designs[[i]]$design
    d$location <- ag_sites[i]
    d$yield <- simulate_yield(d, ag_value, mean = 45 + 3 * i,
                              row_trend = 0.3, col_trend = -0.1, sd = 2)
    d
  }))
  ag_met$yield[15] <- 500       # a typing error
  for (i in seq_along(ag_sites)) {
    here <- ag_met[ag_met$location == ag_sites[i], ]
    ag_sheets[[i]][[my_trait]] <- here$yield[match(ag_sheets[[i]]$plot_number, here$plots)]
    # write.csv(ag_sheets[[i]], file.path(out_dir, paste0("agricolae_sheet_", ag_sites[i], ".csv")), row.names = FALSE)
  }

  # I5/I6. Look at the trait; list and remove IQR outliers per location.
  plot_trait(ag_met, "yield")
  ag_clean <- remove_outliers(ag_met, "yield")
  ag_removed <- attr(ag_clean, "outliers")
  print(ag_removed)

  # I7/I8. SpATS per location, then precision-weighted means across locations.
  # With no SpATS, use replicate/block models for these imported designs.
  ag_analysis <- ag_clean
  # In RCBD, blocks ARE replicates: include this factor once in the model.
  if (ag_kind == "RCBD") ag_analysis$block <- NULL
  ag_two <- analyze_two_stage(ag_analysis, "yield",
    stage1 = if (has_spats) "spatial" else "blocks", replicate_col = "rep",
    block_col = if (ag_kind == "RCBD") NULL else "block")
  compare_first_stage(ag_two)
  print(ag_two)
  ag_stage1 <- do.call(rbind, lapply(seq_along(ag_sites), function(i) {
    saved_plots <- if (use_breedbase) breedbase_trial_plots(con, ag_trials[i]) else NULL
    first_stage_values(ag_two$stage1_fits[[ag_sites[i]]], ag_designs[[i]],
                       trait = my_trait, breedbase_plots = saved_plots)
  }))
  ag_stage2 <- second_stage_values(ag_two, trait = my_trait)

  # I9. Store phenotypes on field plots and results as analysis projects.
  ag_stage1_name <- paste0("AGRICOLAE_", ag_kind, "_", my_year, "_stage1_", my_run_tag)
  ag_stage2_name <- paste0("AGRICOLAE_", ag_kind, "_", my_year, "_stage2_", my_run_tag)
  invisible(save_breedbase_analysis(ag_stage1, analysis_name = ag_stage1_name,
    program = my_program, description = paste("Agricolae", ag_kind, "first-stage plot values"),
    level = "plot", year = my_year, dry_run = TRUE))
  invisible(save_breedbase_analysis(ag_stage2, analysis_name = ag_stage2_name,
    program = my_program, description = paste("Agricolae", ag_kind, "weighted two-stage means"),
    year = my_year, dry_run = TRUE))
  if (save_phenotypes) {
    for (sheet in ag_sheets) save_breedbase_phenotypes(sheet, con, confirm = FALSE)
  }
  if (save_analyses && trials_saved(ag_trials)) {
    save_breedbase_analysis(ag_stage1, con, analysis_name = ag_stage1_name,
      program = my_program, description = paste("Agricolae", ag_kind, "first-stage plot values"), year = my_year,
      level = "plot", folder = my_folder, create_folder = TRUE, confirm = FALSE,
      upload_file = file.path(out_dir, "agricolae_stage1_upload.csv"))
    save_breedbase_analysis(ag_stage2, con, analysis_name = ag_stage2_name,
      program = my_program, description = paste("Agricolae", ag_kind, "weighted two-stage means"), year = my_year,
      folder = my_folder, create_folder = TRUE, confirm = FALSE,
      upload_file = file.path(out_dir, "agricolae_stage2_upload.csv"))
  }

  # I10. Checks on outliers and result dimensions.
  check_result("agricolae: the typing error was excluded", is.na(ag_clean$yield[15]))
  check_result("agricolae: first-stage values match the phenotype sheet plot names",
               setequal(ag_stage1$plot_name, unlist(lapply(ag_sheets, `[[`, "plot_name"))))
  check_result("agricolae: one combined result per accession",
               setequal(ag_stage2$accession, ag_accessions))
} else if (!isTRUE(run_agricolae_analysis)) {
  message("Optional agricolae phenotype/analysis examples skipped: set run_agricolae_analysis <- TRUE to run I4-I10.")
}


# =============================================================================
# ---- G. Summary of the checks ------------------------------------------------
# =============================================================================
check_summary()               # what was checked, by part, and what is correct
# check_summary(details = TRUE) lists every single check.
print(list.files(out_dir))    # files written by this script
