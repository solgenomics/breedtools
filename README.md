# breedtools

Design and analysis of augmented plant breeding trials.

- **Augmented row-column designs** (Piepho and Williams, 2016):
  unreplicated entries plus replicated checks in rows, columns and
  superimposed blocks. Row and column groups can be unequal, blocks can
  hold any number of check plots, and checks are placed by a D-optimal
  search (or covering Latin rectangles for regular layouts). Analyses can
  use random blocks, random entries (BLUPs), or a 'SpATS' spatial model.
- **Augmented p-rep designs** (Williams, Piepho and Whitaker, 2011):
  multi-location, partially replicated trials in incomplete blocks where
  pairs of entries share a block at most once. Choose the construction with
  `method`: a flexible randomized `"search"`, the paper's partial
  alpha-array construction (`"alpha"`), or its augmented alpha(2)-array
  construction (`"alpha2"`). `augmented_prep_parameters()` lists the
  feasible settings, and every design reports its efficiency factor.
- **Resolvable row-column designs** (Williams, John and Whitaker, 2006):
  fully replicated trials with rows and columns within replicates,
  optimized for row-column and spatial (linear variance) efficiency, with
  latinized options; analyses with RC, LV + RC and AR + RC models by REML.
- **Analysis**: ANOVA with degrees of freedom that are estimable from the
  observed plots, plus **adjusted entry means (BLUEs)** with standard
  errors, the mean SED and the 5% LSD. p-rep analyses fit replicates within
  locations and can add a fixed or random (REML) genotype-by-location
  interaction. A two-stage analysis (`analyze_augmented_prep_two_stage()`)
  handles heterogeneous error variances across locations and optional
  spatial adjustment with 'SpATS'.
- **Field maps and 'Breedbase'**: a plot-level field book and upload files,
  or save the design directly into a Breedbase breeding program and folder
  through BrAPI with `save_breedbase()`, after checking the location and
  accessions with `breedbase_check()`. After harvest, save the phenotypic
  data by plot and trait name with `save_breedbase_phenotypes()` (sheet from
  `phenotype_template()`), and the analysis results: first-stage plot values
  as a trait of the field plots (`save_breedbase_first_stage()`) and
  second-stage genotype means as a Breedbase analysis
  (`save_breedbase_analysis()`).
- **agricolae designs**: generate a design with the CRAN package `agricolae`,
  convert and record trial information with `prepare_breedbase_trial()`,
  then use the same `save_breedbase()` and `phenotype_template()` workflow.
- **Multi-environment tools**: read 'Breedbase' phenotype files, run
  Hartley's Fmax variance check, and group environments by GxE.

## Installation

```r
# install.packages("remotes")
remotes::install_github("solgenomics/breedtools")
```

## Documentation

- **Website**: <https://solgenomics.github.io/breedtools/>. The
  [function reference](https://solgenomics.github.io/breedtools/reference/)
  groups every function by task and shows its arguments with examples you
  can run.
- **Tutorial**: `vignette("breedtools")` walks through every design and
  analysis with simulated data.
- **In R**: `help(package = "breedtools")` lists every help page.
  `?augmented_row_column_design` opens one page, and
  `example(augmented_row_column_design)` runs its examples.
- **All functions** are listed by task [at the end of this
  page](#all-functions).

## Find layout options from your entry count

```r
library(breedtools)
design_options(n_entries = 100)  # field dimensions and checks not yet known
design_options(108, design = "augmented_row_column", n_checks = 4, rows_in_field = 12, cols_in_field = 12)
design_options(42, design = "resolvable", n_reps = 3, cols_in_field = 7)
```

Options report field size, checks, replication, plot counts, unused
positions and nominal error degrees of freedom. Augmented row-column
options expose `rows_in_field`, `cols_in_field`, `rows_per_block`,
`cols_per_block`, `checks_per_block`, `method` and `n_candidates` for an
explicit generator call. See `?design_options` for selecting and generating
one. These are parameter options to check with the generator,
its summary and `check_design()`.

## Example

```r
library(breedtools)

rc <- augmented_row_column_design(
  treatments = sprintf("LINE_%03d", 1:108),
  controls = c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D"),
  rows_in_field = 12, cols_in_field = 12,
  rows_per_block = 4, cols_per_block = 4,
  seed = 2026
)
plot_augmented_row_column_design(rc)

field <- rc$design
field$yield <- rnorm(nrow(field), 60, 3)   # your phenotypes here
fit <- analyze_augmented_row_column_design(field, response = "yield")
fit                      # ANOVA and top adjusted means
fit$adjusted_means       # all entries ranked by adjusted mean
```

To see the full workflow, run `vignette("breedtools")`, or open the RStudio
script:

```r
file.edit(system.file("examples", "breedtools_workflow.R", package = "breedtools"))
```

## Save an agricolae design as a new Breedbase trial

```r
accessions <- c("G1", "G2", "G3", "CHECK_A")  # existing Breedbase accession names
randomized <- agricolae::design.rcbd(accessions, r = 3, seed = 2026)
trial <- prepare_breedbase_trial(
  randomized, trial_name = "RCBD_2026", year = 2026,
  location = "Ibadan", program = "Your breeding program",
  description = "RCBD with three complete blocks generated with agricolae",
  checks = "CHECK_A", prompt = FALSE
)

# Inspect the trial and plot payload without a connection.
payload <- save_breedbase(trial, dry_run = TRUE)

# Check against your database, then save as a new trial.
con <- breedbase_login("https://your-breedbase.org", "your_username")
save_breedbase(trial, con, dry_run = TRUE)
save_breedbase(trial, con)
```

`prepare_breedbase_trial()` records the required trial information once
and converts the field book. In an interactive RStudio session, it can ask
for missing name, year, location, program and description. Use
`phenotype_template(trial)` to make a field sheet with the recorded trial's
plot names. Preparation is offline; saving uses your existing connection.

The converter preserves the randomization. It recognizes CRD, RCBD, alpha,
lattice, augmented block, split-plot and strip-plot types. For a field book
alone, including `agricolae::zigzag()` output, supply `design_type`. For
multiple treatment factors, select the accession column explicitly, e.g.
`accession = "genotype"`. Other factors and the main-plot/subplot hierarchy
remain in `trial$book`; the save function uploads individual plots.
See `?as_breedbase_design` for column mapping and other design types.

The complete tour in `scripts/validate_designs.R` includes a second analyzed
trait added to existing analysis projects (A11), an independent single-stage
analysis (run A1-A6, then A12), and the full agricolae workflow (Part I).
Part H reads existing Breedbase trials and lets you choose `"single"`,
`"two_stage"` or `"both"` through `breedbase_analysis_route`; multiple traits
are stored in the same analysis projects.

## All functions

Each function has a help page with runnable examples: `?function_name`,
then `example(function_name)`. Functions that connect to Breedbase have
examples in `\dontrun{}`. To run them, you need a server and an account.

### 1. Choose a design

| Function | What it does |
|---|---|
| `design_options()` | Suggests layouts for any design from the number of entries and the field size, if you know it |
| `augmented_row_column_capacity()` | Counts the entries and checks that fit in an augmented row-column field |
| `augmented_row_column_layouts()` | Lists the analyzable block layouts of an augmented row-column field |
| `augmented_prep_parameters()` | Lists feasible settings for the alpha and alpha(2) p-rep constructions |

```r
design_options(n_entries = 100)
augmented_row_column_layouts(rows_in_field = 12, cols_in_field = 12)
augmented_prep_parameters(n_treatments = 20, n_locations = 2)
```

### 2. Augmented row-column designs (one plot per entry, replicated checks)

| Function | What it does |
|---|---|
| `augmented_row_column_design()` | Generates the design with D-optimal or Latin check placement |
| `plot_augmented_row_column_design()` | Draws the field map |
| `analyze_augmented_row_column_design()` | ANOVA and adjusted means using fixed or random blocks, random entries (BLUPs), or a SpATS model |

```r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  rows_in_field = 12, cols_in_field = 12,
                                  rows_per_block = 4, cols_per_block = 4, seed = 1)
plot_augmented_row_column_design(rc)
d <- rc$design; d$yield <- rnorm(nrow(d), 50, 3)
fit <- analyze_augmented_row_column_design(d, "yield")
fit$adjusted_means
```

### 3. Augmented p-rep designs (several locations, part of the entries repeated)

| Function | What it does |
|---|---|
| `augmented_prep_design()` | Generates the design (`method = "search"`, `"alpha"` or `"alpha2"`) |
| `augmented_prep_concurrence()` | Counts how often each pair of entries shares a block |
| `plot_augmented_prep_design()` | Draws the field map of each location |
| `analyze_augmented_prep_design()` | Single-stage analysis with optional genotype-by-location interaction |
| `analyze_augmented_prep_two_stage()` | Two-stage analysis with heterogeneous error variances and optional SpATS |

```r
prep <- augmented_prep_design(paste0("G", 1:60), locations = c("North", "South"),
                              block_size = 6, field_cols = 10, seed = 1)
plot_augmented_prep_design(prep)
d <- prep$design; d$yield <- rnorm(nrow(d), 25, 3)
analyze_augmented_prep_design(d, "yield")
analyze_augmented_prep_two_stage(d, "yield")
```

### 4. Resolvable row-column designs (fully replicated)

| Function | What it does |
|---|---|
| `resolvable_row_column_design()` | Generates a spatially optimized design, including latinized options |
| `resolvable_row_column_efficiency()` | Calculates row-column, linear variance and AR efficiencies of a design |
| `plot_resolvable_row_column_design()` | Draws the field map |
| `analyze_resolvable_row_column()` | Fits RC, LV + RC or AR + RC models by REML |

```r
rr <- resolvable_row_column_design(paste0("T", 1:42), n_reps = 3,
                                   rows_per_rep = 6, cols_per_rep = 7, seed = 1)
resolvable_row_column_efficiency(rr)
plot_resolvable_row_column_design(rr)
d <- rr$design; d$yield <- rnorm(nrow(d), 40, 4)
analyze_resolvable_row_column(d, "yield", model = "rc")   # or "lv", "ar"
```

### 5. Analysis for any design

| Function | What it does |
|---|---|
| `plot_trait()` | Plots a trait's distribution and field map |
| `find_outliers()` / `remove_outliers()` | Flags or removes outlying plots with the IQR rule |
| `analyze_two_stage()` | Two-stage analysis of several trials of any design |
| `analyze_single_stage()` | Single-stage analysis of several trials, as an alternative |
| `compare_first_stage()` | Compares raw and first-stage adjusted plot values |
| `first_stage_values()` | Returns the adjusted plot values to store as a trait |
| `second_stage_values()` | Returns the genotype values to store as a Breedbase analysis |
| `prediction_accuracy()` | Correlates adjusted means with known (simulated) genotype values |

```r
plot_trait(d, "yield")
find_outliers(d, "yield")
fit2 <- analyze_two_stage(met, "yield", stage1 = "design")  # met: trials stacked with a location column
second_stage_values(fit2, trait = "fresh root yield")
```

### 6. Design checks

| Function | What it does |
|---|---|
| `check_setup()` | Checks that packages and the output folder are ready |
| `check_design()` | Runs the standard checks on any design, such as plot use, entry counts, check placement and estimable effects |
| `check_section()` / `check_result()` | Labels a section and records your own checks |
| `check_summary()` / `check_reset()` | Prints or clears the recorded checklist |

```r
check_reset(); check_design(rc); check_summary()
```

### 7. Breedbase: connect, look up and save designs

| Function | What it does |
|---|---|
| `breedbase_login()` | Logs in to a Breedbase server through BrAPI |
| `breedbase_programs()` / `breedbase_folders()` / `breedbase_locations()` | Lists breeding programs, trial folders and locations |
| `breedbase_check()` | Checks the program, folder, location, accessions, trial and plot names before you save, and says whether everything is READY |
| `add_plot_names()` | Adds the Breedbase plot names to a design, to see, check or edit them before saving |
| `add_plot_coordinates()` | Adds plot GPS coordinates, from the field corner and plot size, or one point per plot |
| `save_breedbase_coordinates()` | Saves plot GPS coordinates to a trial already in Breedbase |
| `prepare_breedbase_trial()` | Converts an agricolae design and records the trial information |
| `as_breedbase_design()` | Converts an agricolae field book only (column mapping) |
| `save_breedbase()` | Saves a design as a new trial with its plots and their GPS coordinates (`dry_run = TRUE` previews it) |
| `breedbase_trial_plots()` | Reads the plots of a saved trial |
| `phenotype_template()` | Makes a phenotyping sheet with the Breedbase plot names |

```r
con <- breedbase_login("https://your-breedbase.org", "your_username")
breedbase_programs(con)
rc <- add_plot_names(rc, "AYT_2026", location = "Ibadan")       # see rc$design$plot_name
rc <- add_plot_coordinates(rc, origin = c(7.4931, 3.9035), plot_width = 3, plot_length = 5)
breedbase_check(con, accessions = rc, location = "Ibadan", program = "Cassava",
                trial_name = "AYT_2026")
save_breedbase(rc, con, program = "Cassava", folder = "AYT 2026",
               trial_name = "AYT_2026", year = 2026, location = "Ibadan",
               description = "Advanced yield trial", dry_run = TRUE)
phenotype_template(rc, traits = "fresh root yield")      # uses the stored plot names
```

### 8. Breedbase: phenotypes and analysis results

| Function | What it does |
|---|---|
| `breedbase_phenotypes()` | Reads phenotypes of trials by trial and trait name |
| `save_breedbase_phenotypes()` | Saves phenotypes by plot and trait name |
| `save_breedbase_first_stage()` | Saves first-stage plot values as a trait of the field plots |
| `save_breedbase_analysis()` | Saves results as a Breedbase analysis, of plots (first stage) or accessions (second stage) |
| `analysis_plot_names()` | Shows the plot and accession names an analysis will have, before saving |
| `breedbase_analysis_plots()` | Reads the plots, accessions and values of an analysis already in Breedbase |
| `add_breedbase_analysis_values()` | Adds another trait to an existing analysis |
| `breedbase_statistics()` | Lists the statistic terms of the SGNSTAT ontology |
| `breedbase_adjusted_trait()` | Finds or creates the trait used for adjusted values |
| `breedbase_refresh()` | Refreshes the server's phenotype cache |

### 9. Multi-environment trials from Breedbase files

| Function | What it does |
|---|---|
| `read_breedbase_file()` / `read_breedbase_excel()` | Reads Breedbase phenotype or design files (Excel, CSV, TSV) |
| `prepare_breedbase_pheno()` | Prepares phenotypes for multi-environment analysis |
| `run_breedbase_variance_check()` | Analyzes each environment, runs Hartley's Fmax test and finds compatible groups |
| `run_environment_stratification()` | Groups environments by GxE using Lin (1982) |
| `run_breedbase_environment_stratification()` | Runs the same grouping directly from a Breedbase file |
| `read_locations()` / `map_environment_groups()` | Loads location coordinates and draws the groups on a map |

```r
pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
loc_file   <- system.file("extdata", "breedbase_locations.csv", package = "breedtools")
vc <- run_breedbase_variance_check(pheno_file, trait_col = "fresh root yield|CO_334:0000013")
vc$maximal_compatible_groups
st <- run_environment_stratification(pheno_file, trait = "fresh root yield")
map_environment_groups(st, read_locations(loc_file, phenotype_data = st))
```

### 10. Building blocks and helpers

The pipeline functions above call these lower-level functions. You can call
them yourself to run one step at a time.

| Group | Functions |
|---|---|
| Variance check | `analyze_environment_anova()`, `analyze_one_environment()`, `hartley_fmax_test()`, `hartley_critical_sim()`, `make_environment_groups()`, `get_compatible_groups()`, `get_maximal_groups()`, `harmonic_mean()` |
| Environment stratification | `prepare_environment_stratification_data()`, `detect_design()`, `prepare_design_factors()`, `calculate_environment_stratification_anova()`, `lin_group_environments()`, `lin_test_group()`, `lin_ss_ge()`, `lin_error_mse()`, `complete_environment_info()`, `add_environment_metadata()`, `add_pairwise_environment_metadata()`, `environment_summary_by_group()` |
| Columns and empty tables | `find_column()`, `guess_column()`, `normalize_column_name()`, `empty_anova()`, `empty_group_membership()`, `empty_group_summary()`, `empty_pairwise()`, `empty_ungrouped()` |
