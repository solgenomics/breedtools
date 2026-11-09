# =============================================================================
# breedtools: example workflow for trial design, analysis and Breedbase
#
# Open this file in RStudio with:
#   file.edit(system.file("examples", "breedtools_workflow.R", package = "breedtools"))
#
# Run it section by section (Ctrl/Cmd + Enter, or Ctrl/Cmd + Alt + T to run
# the current section). Use the outline pane (Ctrl/Cmd + Shift + O) to jump
# between sections.
#
#   Part A - Augmented row-column design (unreplicated entries + checks)
#   Part B - Augmented p-rep design (no checks, several locations)
#   Part C - Resolvable row-column design (every entry replicated)
#   Part D - Multi-environment trials from Breedbase (variances, groups)
#   Part E - Save a design to Breedbase, including agricolae field books (E5)
#
# The phenotypes in this script are SIMULATED so the example runs anywhere.
# Replace those steps with your own field data.
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

# Install once (pick the option that applies to you):
# install.packages("remotes")
# install.packages("breedtools")                    # from CRAN (once released)
# remotes::install_github("solgenomics/breedtools")    # development version
# devtools::install("path/to/breedtools")           # from a local copy
# install.packages(c("lme4", "SpATS", "httr2"))     # optional: mixed and spatial
#                                                   # models, saving to Breedbase

library(breedtools)

# Folder for files this script writes (field books, maps, results).
# Change it to a project folder of your own, e.g. "~/trials/2026".
out_dir <- file.path(tempdir(), "breedtools_example")

# Check that everything is ready (package, optional packages, output folder):
check_setup(out_dir)


# =============================================================================
# PART A - AUGMENTED ROW-COLUMN DESIGN
#
# Use when you have many unreplicated entries (e.g. early generation lines)
# and a few replicated checks. The field is split into rectangular blocks;
# every block gets plots of the checks, and row, column and block effects are
# estimated from the checks.
# =============================================================================


# ---- A1. Plan the layout: will it be analyzable? ----------------------------

# If only your entry count is known, start here. Restrict `design` to the
# families you want and add rows_in_field, cols_in_field, n_checks or n_reps as they become known.
planning_options <- design_options(100)
print(planning_options[, c("design", "n_checks", "rows_in_field", "cols_in_field", "n_reps",
                           "plots_per_location", "error_df")])
# For an augmented row-column design, pass every chosen setting explicitly:
#   selected <- design_options(108, "augmented_row_column", n_checks = 4, rows_in_field = 12, cols_in_field = 12)[1, ]
#   trial <- augmented_row_column_design(
#     treatments = your_entries, controls = your_checks,
#     rows_in_field = selected$rows_in_field, cols_in_field = selected$cols_in_field,
#     rows_per_block = selected$rows_per_block, cols_per_block = selected$cols_per_block,
#     checks_per_block = selected$checks_per_block[[1]], method = selected$method,
#     n_candidates = selected$n_candidates, seed = 2026
#   )
#   check_design(trial)

# Suggested analyzable layouts of a 12 x 12 field, most entry plots first.
# Every block has at least as many rows and columns as checks. There is
# a trade-off: more check plots mean fewer entries but more error df, i.e.
# more reliable F-tests and spatial adjustment. Aim for at least ~10 error df.
print(augmented_row_column_layouts(rows_in_field = 12, cols_in_field = 12, n_controls = 2:5))

# Only the layouts that fit your entries (here 100) with enough error df:
print(augmented_row_column_layouts(12, 12, n_controls = 2:6, min_error_df = 6, n_entries = 100))


# ---- A2. Generate the design ------------------------------------------------

# Your entry and check names (usually read from a seed list / spreadsheet).
entries <- sprintf("LINE_%03d", 1:108)
checks  <- c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D")

rc <- augmented_row_column_design(
  treatments     = entries,
  controls       = checks,
  rows_in_field  = 12,
  cols_in_field  = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  plot_type      = "serpentine",   # plot numbering; or "cartesian"
  n_candidates   = 1000,           # candidate check layouts to compare
  seed           = 2026            # keep the seed to reproduce the design
)

print(rc$summary)          # size, error df, method, D-criterion, neighbouring checks
print(head(rc$design))     # plot-level field book


# ---- A3. Check the design ------------------------------------------------------

# Standard checks: every plot used once, every entry planted once, every row
# and column has a check, check plots per block, error df, estimability, and
# no checks next to each other inside a block.
check_design(rc)


# ---- A4. Plot the field map --------------------------------------------------

plot_augmented_row_column_design(rc, label = "all_entries", fill = "type")
plot_augmented_row_column_design(rc, label = "none", fill = "block",
                                 show_block_labels = TRUE)

# Save the map as a PDF to take to the field:
pdf(file.path(out_dir, "row_column_field_map.pdf"), width = 11, height = 9)
plot_augmented_row_column_design(rc, label = "all_entries", fill = "type")
dev.off()


# ---- A5. Export the design ---------------------------------------------------

# Full field book (for data collection):
write.csv(rc$design, file.path(out_dir, "row_column_fieldbook.csv"), row.names = FALSE)

# Breedbase upload columns (plots, block, all_entries, rep, is_control).
# To save the trial directly into Breedbase instead, see Part E.
write.csv(rc$breedbase_design, file.path(out_dir, "row_column_breedbase.csv"),
          row.names = FALSE)


# ---- A6. Flexible layouts (Piepho and Williams, 2016) -------------------------

# Row and column groups of different sizes, and blocks with fewer check plots
# than checks (incomplete blocks) or more (extra replication). Example 3 of
# the paper: 173 entries, 4 checks, 22 rows x 10 columns, 3-5 checks per block.
checks_per_block <- rbind(c(4, 5, 4),   # one row per row group,
                          c(3, 5, 3),   # one column per column group
                          c(3, 5, 4),
                          c(3, 5, 3))
print(augmented_row_column_capacity(n_controls = 4,
                                    rows_per_rowgroup = c(6, 5, 6, 5),
                                    cols_per_colgroup = c(3, 4, 3),
                                    checks_per_block = checks_per_block))
flex <- augmented_row_column_design(
  treatments        = paste0("E", 1:173),
  controls          = c("A", "B", "C", "D"),
  rows_per_rowgroup = c(6, 5, 6, 5),
  cols_per_colgroup = c(3, 4, 3),
  checks_per_block  = checks_per_block,
  n_candidates      = 3,
  seed              = 1
)
check_design(flex)
plot_augmented_row_column_design(flex, label = "type", fill = "type")


# ---- A7. After harvest: bring in phenotypes ----------------------------------

# In practice you read the file you collected, e.g.
#   pheno <- read_breedbase_file("my_trial_phenotypes.xlsx")
# Here we SIMULATE a yield file with a fertility gradient down the rows and
# a few missing plots, write it, and read it back the same way.
set.seed(1)
sim <- rc$design[, c("plots", "all_entries")]
sim$yield <- round(
  60 + 0.8 * rc$design$row - 0.3 * rc$design$col +       # field trend
    ifelse(rc$design$is_control == 1, 5, 0) +            # checks yield more
    rnorm(nrow(sim), sd = 3), 1
)
sim$yield[c(7, 50, 101)] <- NA                           # lost plots
write.csv(sim, file.path(out_dir, "row_column_phenotypes.csv"), row.names = FALSE)

pheno <- read_breedbase_file(file.path(out_dir, "row_column_phenotypes.csv"))

# Join phenotypes to the design by plot number:
rc_data <- merge(rc$design, pheno[, c("plots", "yield")], by = "plots")
print(head(rc_data))


# ---- A8. Analyze --------------------------------------------------------------

rc_fit <- analyze_augmented_row_column_design(rc_data, response = "yield")

print(rc_fit$anova_formatted)   # ANOVA (rows, columns, blocks, genotypes)
print(rc_fit$diagnostics[, c("anova_error_df", "sigma", "mean_sed", "lsd_5pct")])

# Other models (they need lme4 / SpATS): random rows, columns and blocks
# recover inter-block information; random entries give BLUPs; SpATS fits a
# smooth spatial surface. Smaller mean SED = more precise comparisons.
if (requireNamespace("lme4", quietly = TRUE)) {
  rc_random <- analyze_augmented_row_column_design(rc_data, "yield", random_blocks = TRUE)
  rc_blup <- analyze_augmented_row_column_design(rc_data, "yield", random_blocks = TRUE,
                                                 entries = "random")
  print(rc_random$variance_components)
}
if (requireNamespace("SpATS", quietly = TRUE)) {
  rc_spatial <- analyze_augmented_row_column_design(rc_data, "yield", spatial = TRUE)
}
print(c(fixed = rc_fit$diagnostics$mean_sed,
        random_blocks = if (exists("rc_random")) rc_random$diagnostics$mean_sed else NA,
        spatial = if (exists("rc_spatial")) rc_spatial$diagnostics$mean_sed else NA))


# ---- A9. Adjusted entry means: rank and select -------------------------------

# Adjusted means (BLUEs) remove the row, column, and block effects of the plot
# each entry landed in, estimated from the checks. Rank on these, not on the
# raw plot values: an entry in a poor corner of the field is not penalized.
rc_means <- rc_fit$adjusted_means
print(head(rc_means[, c("rank", "genotype", "raw_mean", "adjusted_mean", "std_error")], 10))

# How much did spatial adjustment reorder the entries?
entries_only <- rc_means[rc_means$is_control == 0 & !is.na(rc_means$adjusted_mean), ]
print(cor(entries_only$raw_mean, entries_only$adjusted_mean, method = "spearman"))

# Entries that beat the best check by more than one LSD (5%):
best_check <- max(rc_means$adjusted_mean[rc_means$is_control == 1], na.rm = TRUE)
lsd <- rc_fit$diagnostics$lsd_5pct
print(entries_only[entries_only$adjusted_mean > best_check + lsd,
                   c("rank", "genotype", "adjusted_mean", "std_error")])

# Residual checks:
op <- par(mfrow = c(1, 2))
plot(rc_fit$plot_results$fitted, rc_fit$plot_results$residual,
     xlab = "Fitted", ylab = "Residual", main = "Residuals vs fitted")
abline(h = 0, lty = 2)
qqnorm(rc_fit$plot_results$residual[rc_fit$plot_results$is_control == 1],
       main = "Q-Q plot (check plots)")
qqline(rc_fit$plot_results$residual[rc_fit$plot_results$is_control == 1])
par(op)

# Save results:
write.csv(rc_fit$anova_formatted, file.path(out_dir, "row_column_anova.csv"), row.names = FALSE)
write.csv(rc_means, file.path(out_dir, "row_column_adjusted_means.csv"), row.names = FALSE)
write.csv(rc_fit$plot_results, file.path(out_dir, "row_column_plot_results.csv"), row.names = FALSE)


# =============================================================================
# PART B - AUGMENTED P-REP DESIGN
#
# Use for multi-location trials without replicated checks: every entry is
# planted once at every location, and a subset of entries gets extra plots
# at one location. These repeated entries replace the checks. Plots are
# grouped into incomplete blocks so that pairs of entries share a block at
# most once.
# =============================================================================


# ---- B1. Generate the design --------------------------------------------------

prep_entries <- sprintf("G%03d", 1:120)
sites        <- c("North", "Central", "South")

prep <- augmented_prep_design(
  treatments      = prep_entries,
  locations       = sites,
  block_size      = 10,          # max plots per incomplete block
  prep_prop       = 0.2,         # 20% of entries repeated at each location
  replications_per_repeated_entry = 2,
  field_cols      = 10,          # field columns, for row/col coordinates
  serpentine      = TRUE,
  n_candidates    = 1000,        # stops early if an optimal design is found
  seed            = 2026
)

print(prep$summary)          # concurrence, connectedness, efficiency factor
print(head(prep$design))
check_design(prep)

# Algebraic constructions of Williams, Piepho and Whitaker (2011): every entry
# duplicated at exactly one location, two replicates per location, m p-rep
# entries in every block. See which settings are feasible first:
print(augmented_prep_parameters(n_treatments = 120, n_locations = 3, max_block_size = 10))

prep_alpha <- augmented_prep_design(
  treatments = prep_entries,
  locations  = sites,
  method     = "alpha",        # or "alpha2" (even m, odd s)
  block_size = 10,
  field_cols = 10,
  seed       = 2026
)
print(prep_alpha$summary[, c("method", "m", "s", "n_plots", "max_pairwise_concurrence",
                             "efficiency_factor")])
check_design(prep_alpha)


# ---- B2. Inspect the design ---------------------------------------------------

# Plots per entry (one per location + extra plots where repeated):
print(table(table(prep$design$all_entries)))

# Pairwise concurrence: how often two entries share a block (ideally <= 1)
conc <- prep$concurrence_matrix
print(table(conc[upper.tri(conc)]))


# ---- B3. Plot and export ------------------------------------------------------

# Orange = entries repeated at this location (they replace the checks);
# thick black lines = blocks; a blue line separates replicates.
plot_augmented_prep_design(prep, locations = "North")
plot_augmented_prep_design(prep_alpha, locations = "North")

pdf(file.path(out_dir, "prep_field_maps.pdf"), width = 11, height = 9)
plot_augmented_prep_design(prep)     # all locations
dev.off()

write.csv(prep$design, file.path(out_dir, "prep_fieldbook.csv"), row.names = FALSE)
write.csv(prep$breedbase_design, file.path(out_dir, "prep_breedbase.csv"), row.names = FALSE)


# ---- B4. After harvest: phenotypes -------------------------------------------

# SIMULATED: location and block effects + genotype differences + noise.
set.seed(2)
genotype_effect <- setNames(rnorm(length(prep_entries), sd = 4), prep_entries)
location_effect <- c(North = 0, Central = 8, South = -5)
block_id        <- interaction(prep$design$location, prep$design$block)
block_effect    <- setNames(rnorm(nlevels(block_id), sd = 2), levels(block_id))

prep_data <- prep$design
prep_data$yield <- 50 +
  location_effect[prep_data$location] +
  block_effect[as.character(block_id)] +
  genotype_effect[prep_data$all_entries] +
  rnorm(nrow(prep_data), sd = 3)
prep_data$yield[sample(nrow(prep_data), 6)] <- NA     # a few lost plots


# ---- B5. Analyze --------------------------------------------------------------

prep_fit <- analyze_augmented_prep_design(prep_data, response = "yield")

print(prep_fit$anova_formatted)
print(prep_fit$diagnostics[, c("anova_error_df", "sigma", "mean_sed", "lsd_5pct")])

# Adjusted means (BLUEs) correct for the location and block of every plot,
# so entries are compared fairly even though only some are repeated.
prep_means <- prep_fit$adjusted_means
print(head(prep_means[, c("rank", "genotype", "n_observations", "adjusted_mean", "std_error")], 10))

# With simulated data we know the truth: adjusted means track it better.
print(c(raw = cor(prep_means$raw_mean, genotype_effect[prep_means$genotype], use = "complete.obs"),
        adjusted = prediction_accuracy(prep_fit, genotype_effect)))

# A single location can be analyzed on its own:
north_fit <- analyze_augmented_prep_design(
  prep_data[prep_data$location == "North", ], response = "yield"
)
print(north_fit$anova_formatted)

# Genotype-by-location interaction added back at the analysis stage
# (random effect by REML; needs the lme4 package):
if (requireNamespace("lme4", quietly = TRUE)) {
  gxl_fit <- analyze_augmented_prep_design(prep_data, response = "yield",
                                           genotype_by_location = "random")
  print(gxl_fit$variance_components)
}

write.csv(prep_fit$anova_formatted, file.path(out_dir, "prep_anova.csv"), row.names = FALSE)
write.csv(prep_means, file.path(out_dir, "prep_adjusted_means.csv"), row.names = FALSE)


# ---- B6. Two-stage analysis: locations with different error variances --------

# Each location is analysed on its own, then the location means are combined
# with weights from their precision, so noisy locations count less.
# Use stage1 = "spatial" for a SpATS spatial model at each location.
two <- analyze_augmented_prep_two_stage(prep_data, response = "yield")
print(two$locations[, c("location", "residual_variance", "residual_df", "mean_sed")])
print(two$heterogeneity)        # Hartley's test of equal error variances
print(head(two$adjusted_means, 10))
write.csv(two$adjusted_means, file.path(out_dir, "prep_two_stage_means.csv"), row.names = FALSE)


# =============================================================================
# PART C - RESOLVABLE ROW-COLUMN DESIGN
#
# Use when every entry can be planted in every replicate. Each replicate is a
# grid of rows and columns holding every entry once; the design is optimized
# for row-column and spatial efficiency (Williams, John and Whitaker, 2006).
# =============================================================================


# ---- C1. Generate, check and plot ----------------------------------------------

rrc <- resolvable_row_column_design(
  treatments   = paste0("T", 1:42),
  n_reps       = 3,
  rows_per_rep = 6,
  cols_per_rep = 7,
  latinize     = TRUE,       # no entry twice in a long column (contiguous reps)
  n_starts     = 1,
  seed         = 1
)
print(rrc$summary[, c("efficiency_rc", "efficiency_rows", "efficiency_cols",
                      "efficiency_spatial", "upper_bound_spatial", "latin_ok")])
check_design(rrc)
plot_resolvable_row_column_design(rrc)
write.csv(rrc$design, file.path(out_dir, "rrc_fieldbook.csv"), row.names = FALSE)


# ---- C2. Simulate phenotypes and analyze ----------------------------------------

set.seed(5)
rrc_effect <- setNames(rnorm(42, sd = 2), paste0("T", 1:42))
rrc_data <- rrc$design
rrc_data$yield <- 50 + rrc_effect[rrc_data$all_entries] +
  rnorm(18, sd = 1.2)[(rrc_data$replicate - 1) * 6 + rrc_data$row_in_rep] +
  rnorm(21, sd = 0.8)[(rrc_data$replicate - 1) * 7 + rrc_data$col_in_rep] +
  rnorm(nrow(rrc_data))

# Row-column mixed model, compared with spatial models (lower AIC is better):
rrc_fit <- analyze_resolvable_row_column(rrc_data, response = "yield", model = "rc",
                                         long_col = "long_col")
print(rrc_fit$model_comparison[, c("model", "deviance", "AIC", "lrt_p_value")])
print(rrc_fit$variance_parameters)
print(head(rrc_fit$adjusted_means, 10))
write.csv(rrc_fit$adjusted_means, file.path(out_dir, "rrc_adjusted_means.csv"), row.names = FALSE)


# =============================================================================
# PART D - MULTI-ENVIRONMENT TRIALS FROM BREEDBASE
#
# Before combining trials: are their error variances similar, and which
# environments rank genotypes alike? The package ships an example file.
# =============================================================================


# ---- D1. Read the phenotypes ------------------------------------------------------

pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
pheno_raw <- read_breedbase_file(pheno_file)     # also .xlsx and .txt files
print(head(pheno_raw))


# ---- D2. Error variances (Hartley's Fmax) and environment groups (Lin) -----------

variance_check <- run_breedbase_variance_check(pheno_file,
                                               trait_col = "fresh root yield|CO_334:0000013",
                                               n_sim = 20000)
print(variance_check$overall_hartley)
print(variance_check$maximal_compatible_groups)

strata <- run_environment_stratification(pheno_file, trait = "fresh root yield")
print(strata$group_summary[, c("group_id", "locations", "n_env", "f_value", "p_value")])

# Map of the environment groups (opens in the RStudio viewer):
locations <- read_locations(system.file("extdata", "breedbase_locations.csv", package = "breedtools"),
                            phenotype_data = strata)
map_environment_groups(strata, locations)


# =============================================================================
# PART E - SAVE A DESIGN TO BREEDBASE
#
# Any design above can be saved as a new trial, with all its plots, into a
# breeding program and folder of your Breedbase, through BrAPI.
# Breedbase words: a breeding program holds folders, a folder holds trials,
# a trial holds plots.
# =============================================================================


# ---- E1. See what would be sent (no server needed) --------------------------------

to_send <- save_breedbase(rc, trial_name = "AYT_2026_Ibadan", year = 2026,
                          location = "Ibadan", dry_run = TRUE,
                          description = "Advanced yield trial of 108 lines and 4 checks, Ibadan 2026")
str(to_send$plots[[1]][[1]])   # first plot: name, accession, row/column, rep, block

# Plot names are made automatically as trialName_accession_block_rep:
print(head(vapply(to_send$plots[[1]], function(p) p$observationUnitName, character(1))))
# Other options:
#   save_breedbase(..., plot_name_pattern = "{trial}_{plot}_r{row}c{col}")
#   save_breedbase(..., plot_names = my_plot_names)   # your own names


# ---- E2. Save to your Breedbase -------------------------------------------------

# Save to YOUR Breedbase, step by step:
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


# ---- E3. Phenotypic data to Breedbase -------------------------------------------

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


# ---- E4. Analysis results to Breedbase (first and second stage) ----------------

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


# ---- E5. Agricolae field books to a new Breedbase trial ------------------------
check_section("E5. Agricolae to Breedbase")

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


# ---- Done ---------------------------------------------------------------------

check_summary()          # all design checks run in this script
print(list.files(out_dir))
# browseURL(out_dir)     # open the output folder
