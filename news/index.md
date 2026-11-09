# Changelog

## breedtools 0.1.0

- [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  lists, before saving and in a dry run, the composed trait each trait
  will be stored under (“trait\|statistic\|COMP:id”) and whether it
  exists; the dry run returns them as `traits`. Choose another statistic
  from `breedbase_statistics(con)` with `statistic`.

- The bundled statistics ontology (`cxgn_statistics.obo`) matches
  Breedbase (solgenomics/sgn, November 2026);
  `data-raw/update_statistics_ontology.R` refreshes it when Breedbase
  adds statistics.

- [`add_breedbase_analysis_values()`](https://solgenomics.github.io/breedtools/reference/add_breedbase_analysis_values.md)
  never stores analysis values under a plain trait: each trait is
  composed with a statistic (“trait\|statistic\|COMP:id”), by default
  the one the analysis already uses, as Breedbase does when it saves an
  analysis. Before, values added without `statistic` went under the
  plain trait (e.g. “dry matter content percentage\|CO_334:0000092”).

- [`analysis_plot_names()`](https://solgenomics.github.io/breedtools/reference/analysis_plot_names.md)
  previews the plot and accession names of an analysis before saving;
  [`breedbase_analysis_plots()`](https://solgenomics.github.io/breedtools/reference/breedbase_analysis_plots.md)
  reads the plots, accessions and values of a saved analysis.

- Plot GPS coordinates: the new
  [`add_plot_coordinates()`](https://solgenomics.github.io/breedtools/reference/add_plot_coordinates.md)
  adds them to a design, calculated from the GPS of the field corner,
  plot size, alleys and field direction, or given per plot (vectors, or
  a table such as the Breedbase GPS upload template).
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  saves them with the plots (Breedbase stores them only through an
  update after the plots are created, so they are sent just after and
  read back to check), and the new
  [`save_breedbase_coordinates()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_coordinates.md)
  adds them to a trial already saved.

- Tests: the fake Breedbase server reads the login form again (the tests
  of saving to Breedbase had stopped running).

- Plot names can be made when the design is made: give
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md)
  or
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md)
  a `trial_name` (and `location`), or call the new
  [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md).
  They are stored in a `plot_name` column, which can be edited, and are
  used by
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md),
  [`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
  and
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md).
  [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
  now also checks plot names (free in the database, not repeated), and
  accepts a list of designs.

- The default plot name pattern now includes the location:
  `"{trial}_{location}_{accession}_{block}_{rep}"` (new field
  `{location}`). For trials saved with the old default, pass
  `plot_name_pattern = "{trial}_{accession}_{block}_{rep}"` to
  [`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
  and
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md).

- First CRAN release.

- [`design_options()`](https://solgenomics.github.io/breedtools/reference/design_options.md):
  the field dimensions are now `rows_in_field` and `cols_in_field`, the
  same names as in its result and in
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md).
  `rows` and `cols` still work, with a deprecation warning. When
  possible, the first augmented row-column option has square blocks of
  `n_checks` x `n_checks` with each check exactly once (e.g. 3 x 3
  blocks for three checks), flagged by the new column
  `square_check_blocks`; other block shapes are still listed. With one
  known dimension, the square field is also tested.

- [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md)
  now builds layouts with zero error df when rows, columns and blocks
  are fixed, such as Fig. 2 of Piepho and Williams (2016) (36 entries,
  three checks once in each 3 x 3 block). Analyze them with
  `random_blocks = TRUE` or `entries = "random"`; `summary$error_df`
  shows 0 and
  [`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md)
  reports it as a NOTE. By default
  [`design_options()`](https://solgenomics.github.io/breedtools/reference/design_options.md)
  lists the Fig. 2 layout (blocks of `n_checks` x `n_checks` plots, each
  check exactly once) first when it fits, as the one exception to its
  positive error-df requirement; `min_error_df = 0` lists other
  saturated layouts too.

- [`design_options()`](https://solgenomics.github.io/breedtools/reference/design_options.md)
  also suggests augmented row-column layouts with unequal row and column
  groups, as in Fig. 4 of Piepho and Williams (2016) (22 x 10 field, row
  groups 6, 5, 6, 5, column groups 3, 4, 3, 47 check plots), with new
  list-columns `rows_per_rowgroup` and `cols_per_colgroup`. Within a
  field, layouts whose blocks hold about one plot per check come first,
  then squarer blocks: Fig. 3 (36 x 18 field, 6 x 6 blocks, four checks)
  is the first suggestion for its field.

- Documentation website (pkgdown) with every function grouped by task,
  and a grouped list of all functions with examples in the README.

- Names saved in Breedbase contain only letters, digits, “*”, ”-” and
  ”.”: in trial, plot, new folder and analysis names, spaces and
  characters such as % \$ ( ) { } \[ \]  / \< \> \* become ”*”
  (e.g. location “Alcobaca (BA)” gives trial
  “EMBRAPA_2026_Alcobaca_BA”). Accession and location names are kept,
  existing folders are used as they are, and trials can still be looked
  up by the original name.

- [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  has `folder` and `create_folder`: the analysis is put in a folder of
  the program (as the website’s “Move to folder”), instead of directly
  in the breeding program.

- New
  [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
  checks the location(s) and every entry and check (as accessions) in
  Breedbase before a design is built, plus the program, folder and trial
  name when given. All problems are listed together with close matches
  (e.g. “ibadan” -\> “Ibadan”).
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  now reports location and accession problems together too.

- Breedbase searches work on slow servers: accessions are searched 25 at
  a time and a batch that times out (HTTP 504/502/503/429) is split and
  tried again, with a line after each batch (names done, percent, time).
  Accessions found are remembered for the R session, so
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  after
  [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
  does not search them again.

- [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  survives server timeouts when writing: after a timeout it checks
  whether the folder, trial or plots were created before trying again
  (no duplicates), and a re-run continues a save that was cut off (an
  existing folder is used; a trial without plots gets its plots).

- [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
  says plainly when the trial already exists (“choose another
  trial_name”), and that an existing folder is fine: the trial is saved
  inside it.
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  also uses a folder that Breedbase reports as existing when it tries to
  create it.

- The `folder` of
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  is optional: without it the trial is saved directly in the breeding
  program.

- [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  sends the `germplasmDbId` of each accession with the plots, which
  Breedbase requires (“Verify Germplasm Id.” error before).

- New
  [`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md)
  saves phenotypic data to Breedbase by plot name and trait name (wide:
  one column per trait, or long: trait and value columns; optional time
  stamps for repeated measurements). Plots and traits are checked first,
  values are sent in batches (surviving server timeouts) and read back
  to verify. New
  [`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
  makes the phenotyping sheet with the Breedbase plot names of a design.

- Analysis results to Breedbase. First stage:
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)
  gives one value per field plot (adjusted mean of the genotype + plot
  residual, design and spatial effects removed) and
  [`save_breedbase_first_stage()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_first_stage.md)
  saves it on the field plots as its own trait: the measured trait
  composed with a statistic (“fresh root yield\|Adjusted Means from
  Spatial Correction using SpATS R”), found or created by
  [`breedbase_adjusted_trait()`](https://solgenomics.github.io/breedtools/reference/breedbase_adjusted_trait.md).
  Second stage:
  [`second_stage_values()`](https://solgenomics.github.io/breedtools/reference/second_stage_values.md)
  collects the genotype adjusted means and
  [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  saves them as a Breedbase analysis, like the Breedbase mixed-models
  tool.

- [`analyze_resolvable_row_column()`](https://solgenomics.github.io/breedtools/reference/analyze_resolvable_row_column.md)
  now returns `plot_results` with the fitted values and plot residuals.

- New
  [`analyze_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_two_stage.md):
  two-stage analysis of a multi-location trial for every design
  (augmented row-column, augmented p-rep, resolvable row-column), with a
  SpATS spatial stage 1 (default) or the design’s own model, and the
  weighted stage 2 of
  [`analyze_augmented_prep_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md).

- New
  [`compare_first_stage()`](https://solgenomics.github.io/breedtools/reference/compare_first_stage.md):
  raw against first-stage adjusted plot values per location
  (correlation, spread, change in genotype ranking), with raw and
  adjusted field maps.

- New
  [`find_outliers()`](https://solgenomics.github.io/breedtools/reference/find_outliers.md)
  and
  [`remove_outliers()`](https://solgenomics.github.io/breedtools/reference/remove_outliers.md):
  IQR rule (multiplier 1.5 by default) within each location; removed
  values are set to `NA` and listed by plot.

- New
  [`plot_trait()`](https://solgenomics.github.io/breedtools/reference/plot_trait.md):
  histogram and field map of a trait per location, with outliers
  circled.

- [`analyze_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_row_column_design.md)
  returns `adjusted_means_vcov`; with `spatial = TRUE` or random
  effects, `plot_results` residuals now come from that model (they came
  from the fixed model before).

- [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)
  has `site_username` / `site_password` for test sites with an extra web
  password in front of the website (BrAPI works without it, the trait
  composer and analyses need it); a clear message says so. Usually not
  needed: at login, when a site’s web pages ask for a password, the
  Breedbase login is tried there and used if it works.

- [`breedbase_adjusted_trait()`](https://solgenomics.github.io/breedtools/reference/breedbase_adjusted_trait.md)
  finds an existing adjusted trait through BrAPI first, so the first
  stage needs no web page when the trait exists.

- Breedbase’s BrAPI upload does not refresh the phenotype cache that its
  trial pages and search wizard read, so saved values did not show
  there.
  [`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md)
  and
  [`save_breedbase_first_stage()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_first_stage.md)
  now ask Breedbase to refresh it; new
  [`breedbase_refresh()`](https://solgenomics.github.io/breedtools/reference/breedbase_refresh.md)
  does it on demand.

- No fixed Breedbase terms: the statistic of first- and second-stage
  results is chosen from the server’s own statistic terms (new
  [`breedbase_statistics()`](https://solgenomics.github.io/breedtools/reference/breedbase_statistics.md)
  lists them), a given statistic is checked against the server, and
  composed traits are matched by their components. A model type unknown
  to the server gives a clear message.

- [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  uses the same upload as Breedbase’s “Upload analysis” form (name,
  description, protocol, breeding program, statistic, optional time
  terms) with Breedbase’s model type for uploaded results, and has
  `level = "plot"` (“Analysis of: plots”): first-stage values become an
  analysis with one plot per field plot, linked by Breedbase to the
  field trial.

- Statistics are checked against the Breedbase statistics ontology
  (`cxgn_statistics.obo` from the Breedbase source, included in the
  package): a statistic must be a real SGNSTAT term, a given id must
  belong to the given name, and misspelled names get suggestions;
  [`breedbase_statistics()`](https://solgenomics.github.io/breedtools/reference/breedbase_statistics.md)
  lists the terms without a server.

- Sites whose web pages need their own web password (e.g. cassava-test):
  the login asks for it; without it,
  [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  checks everything through BrAPI and writes the file of the website’s
  “Upload analysis” form, printing every form field to fill in.

- New
  [`breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/breedbase_phenotypes.md)
  reads the phenotypes of a list of trials and traits from Breedbase
  (plots with accession, rep, block, row, column).

- New
  [`analyze_single_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_single_stage.md):
  all trials in one model (lme4, or SpATS with `spatial = TRUE`: trials
  side by side, each with its own spatial trend), to compare with the
  two-stage analysis.

- `analyze_two_stage(stage1 = "blocks")`: a first stage without a design
  object, for data read from Breedbase;
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)
  works without a design for such data.

- `add_breedbase_analysis_values(statistic = )` stores values under
  “trait\|statistic”, so several analyses (e.g. two-stage and
  single-stage means) can sit in the same analysis project.

- New
  [`add_breedbase_analysis_values()`](https://solgenomics.github.io/breedtools/reference/add_breedbase_analysis_values.md)
  adds results (e.g. a second trait) to an analysis already in
  Breedbase, through the phenotype upload of the website: Breedbase
  verifies the file first, and values already stored are not replaced
  unless `overwrite = TRUE`.

- The script has three Breedbase switches (`save_designs`,
  `save_phenotypes`, `save_analyses`), a trial-name prefix and a run tag
  for unique analysis names, and saves analyses only for trials in
  Breedbase.

- New
  [`breedbase_trial_plots()`](https://solgenomics.github.io/breedtools/reference/breedbase_trial_plots.md)
  reads the plots of a saved trial;
  `first_stage_values(breedbase_plots = )` checks that the design in R
  is the one saved (same accession at every plot number) and uses the
  saved plot names.

- `scripts/validate_designs.R` is organized in one part per design, each
  from the field plan to the results in Breedbase. See
  `options(breedtools.germplasm_batch)` and
  `options(breedtools.breedbase_timeout)`.

- New
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  saves any design (augmented row-column, p-rep, resolvable row-column)
  as a trial with its plots in a ‘Breedbase’ database through BrAPI v2,
  in a chosen breeding program and folder. It checks the program, folder
  (optionally created), locations, accessions and trial name before
  writing, asks for confirmation, and reads the plots back to verify.
  `dry_run = TRUE` makes all checks without writing. Plot names are made
  automatically as `trialName_accession_block_rep` (`plot_name_pattern`)
  or given by you (`plot_names`), and are checked to be new in the
  database.
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)
  logs in without storing the password;
  [`breedbase_programs()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md),
  [`breedbase_folders()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)
  and
  [`breedbase_locations()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)
  list what exists.

- New resolvable row-column designs (Williams, John and Whitaker, 2006):

  - [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md)
    builds designs for `v <= k * s` treatments in `r` replicates of
    `k x s`, optimizing the row-column, row and column efficiency
    factors and then the spatial (LV + RC) efficiency factor. Latinized
    and t-latinized designs for contiguous replicates, filler plots, and
    improving or laying out an existing design (`start_design`) are
    supported.
  - [`resolvable_row_column_efficiency()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_efficiency.md)
    reports the efficiency factors, the spatial upper bound, the
    AR-model criterion z(X), and treatment neighbour and edge counts.
  - [`analyze_resolvable_row_column()`](https://solgenomics.github.io/breedtools/reference/analyze_resolvable_row_column.md)
    fits the RC, LV + RC and AR + RC models by REML and compares them by
    deviance; long columns of latinized designs can be fitted as fixed.
  - [`plot_resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/plot_resolvable_row_column_design.md)
    draws the field.

- [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md):

  - Check placement now uses covering Latin rectangles, so every field
    row and column gets a check.
  - Layouts whose effects cannot be estimated, or that leave no error
    degrees of freedom, are rejected.
  - The design is much faster to generate.

- New helpers for checking designs and simulations:
  [`check_setup()`](https://solgenomics.github.io/breedtools/reference/check_setup.md)
  checks at the start of a script that breedtools, the optional ‘lme4’
  and ‘SpATS’ packages and the output folder are ready,
  [`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md)
  runs the standard checks of any design and prints `[OK]` / `[FAIL]`
  lines,
  [`check_result()`](https://solgenomics.github.io/breedtools/reference/check_result.md)
  records a check of your own,
  [`check_section()`](https://solgenomics.github.io/breedtools/reference/check_section.md)
  groups the following checks,
  [`check_summary()`](https://solgenomics.github.io/breedtools/reference/check_summary.md)
  reports what was checked by section and what is correct,
  [`check_reset()`](https://solgenomics.github.io/breedtools/reference/check_reset.md)
  clears the tally (saying how many checks it cleared), and
  [`prediction_accuracy()`](https://solgenomics.github.io/breedtools/reference/prediction_accuracy.md)
  gives the correlation between adjusted means and known true values.

- New
  [`augmented_row_column_layouts()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_layouts.md)
  lists every analyzable layout of a field (block sizes and numbers of
  checks), most entry plots first, with optional `n_entries` and
  `min_error_df` filters.

- [`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md)
  reports the error degrees of freedom of a layout, the share of check
  plots, and the row-group and column-group conditions of Piepho and
  Williams (2016).

- Augmented row-column designs now support the full strategy of Piepho
  and Williams (2016): row and column groups of unequal size
  (`rows_per_rowgroup`, `cols_per_colgroup`) and any number of check
  plots per block (`checks_per_block`), giving incomplete blocks or
  extra check replication. The new `method = "optimal"` allocates checks
  by a D-optimal local search that keeps check plots from being
  neighbours (same column in consecutive rows, or same row in
  consecutive columns) and keeps the check allocation balanced: a block
  with at least as many check plots as checks holds every check (extra
  plots repeat checks), an incomplete block holds each check at most
  once, rows and columns of a group get even numbers of check plots, and
  a column (or row) with room for all checks holds every check. Figures
  3 and 4 of Piepho and Williams (2016) are reproduced with all their
  stated properties; `"latin"` keeps the covering Latin rectangle
  construction (now with evenly spaced check rows and columns, and
  polished by swapping rows or columns within groups to remove
  neighbouring checks), and `"auto"` picks it when it applies. `summary`
  reports the D-criterion.

- [`analyze_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_row_column_design.md)
  can take rows, columns and blocks as random to recover inter-block
  information (`random_blocks = TRUE`), model entries as random with
  fixed checks (`entries = "random"`, BLUPs), or add a spatial P-spline
  surface with ‘SpATS’ (`spatial = TRUE`).

- [`analyze_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_row_column_design.md):

  - Degrees of freedom are corrected (rows within row groups `k - gk`,
    columns within column groups `s - gs`), and the ANOVA uses the df
    that are estimable from the observed plots.
  - The model is fitted by absorbing unreplicated entries, so it is fast
    for large fields.

- [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md)
  uses a faster candidate search with incremental concurrence scoring.
  The search stops early once an optimal design is found.

- [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md)
  follows Williams, Piepho and Whitaker (2011) more closely:

  - p-rep plots are spread evenly over the blocks of each location;
  - the blocks of each location are kept connected by the p-rep entries,
    so each location can be analysed on its own (a warning is given when
    there are too few p-rep plots to do so);
  - a swap-based improvement step removes pairs of entries that share
    more than one block, giving (0,1)-designs for cases such as the
    paper’s six-location trial with blocks of eight; repeated
    concurrences are penalized quadratically, so when a (0,1)-design is
    out of reach the excess is spread over pairs instead of one pair
    meeting at every location, and greedy candidates that put an entry
    twice in a block are repaired instead of discarded;
  - new `resolvable = TRUE` option builds a resolvable design at each
    location (every p-rep entry once per replicate), returned in the new
    `replicate` column;
  - `summary` reports the average efficiency factor
    (`efficiency_factor`) and the number of disconnected locations;
    among optimal candidates the most efficient one is returned
    (`n_efficiency_candidates`).
  - new `method` argument chooses the construction: `"search"`
    (default), `"alpha"` (partial alpha-arrays, Sections 2.1-2.2 of the
    paper) or `"alpha2"` (augmented alpha(2)-arrays, Section 2.4), with
    `n_prep_per_block` for the number of p-rep entries per block (`m`).

- New
  [`analyze_augmented_prep_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md):
  two-stage analysis of multi-location p-rep trials. Stage 1 analyses
  each location with its own error variance, either with the
  incomplete-block model or with a spatial P-spline model from the
  suggested ‘SpATS’ package (`stage1 = "spatial"`). Stage 2 combines the
  location means with full or diagonal stage-1 weighting (Möhring and
  Piepho, 2009) and a REML genotype-by-location variance. It reports
  Hartley’s Fmax test of heterogeneous error variances.

- [`analyze_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_design.md)
  returns the covariance matrix of the adjusted means
  (`$adjusted_means_vcov`).

- Methods `"alpha"` and `"alpha2"` stop searching once further random
  starts no longer help, and warn when no (0,1)-design was found.

- New
  [`augmented_prep_parameters()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_parameters.md)
  lists the feasible parameter sets of the algebraic constructions (as
  in Table 3 of the paper).

- [`map_environment_groups()`](https://solgenomics.github.io/breedtools/reference/map_environment_groups.md)
  now uses the same maps as the ‘Breedbase’ Locations page by default
  (`tiles = "breedbase"`): Esri Street, Topographical and Satellite,
  with a switch between them. They need no account or API key and work
  in the RStudio viewer. `tiles` can also be `"street"`, `"topo"`,
  `"satellite"`, `"OpenStreetMap"` or a leaflet provider name.

- [`plot_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/plot_augmented_prep_design.md)
  draws clearer maps: entries repeated at the location in orange with
  bold names (they replace the checks), single plots in light grey,
  thick black lines around blocks, a blue line and labels for
  replicates, a plain-language legend, a caption explaining the map
  (`caption = TRUE`), and locations in the order of the design. It can
  also colour plots by replicate.

- [`analyze_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_design.md):

  - Handles single-location trials.
  - Reports df that are estimable when plots are missing.

- [`analyze_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_design.md)
  fits replicates within locations for resolvable designs
  (`replicate_col`) and can add a genotype-by-location interaction,
  fixed (`genotype_by_location = "fixed"`) or random (`"random"`, REML
  via the suggested ‘lme4’ package, with variance components in
  `$variance_components`).

- Both analyses now return adjusted entry means (BLUEs) with standard
  errors (`$adjusted_means`), and report the mean SED and the 5% LSD in
  `$diagnostics`.

- New tutorial,
  [`vignette("breedtools")`](https://solgenomics.github.io/breedtools/articles/breedtools.md),
  covering every design and analysis and the multi-environment tools,
  and an RStudio workflow script in `inst/examples`. The source package
  includes `scripts/validate_designs.R`, a step-by-step tour of every
  design and analysis with simple checks, and
  `scripts/run_all_checks.R`, the full automated validation against the
  published examples, [`lm()`](https://rdrr.io/r/stats/lm.html), ‘lme4’
  and simulation.
