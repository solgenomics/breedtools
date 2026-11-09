# Package index

## Choose a design

Start from the number of entries and the field space you have.

- [`design_options()`](https://solgenomics.github.io/breedtools/reference/design_options.md)
  : Suggest trial layouts from the number of entries and available field
  space
- [`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md)
  : Calculate augmented row-column design capacity
- [`augmented_row_column_layouts()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_layouts.md)
  : List the analyzable layouts of an augmented row-column field
- [`augmented_prep_parameters()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_parameters.md)
  : Feasible parameters for augmented p-rep constructions

## Augmented row-column designs

Unreplicated entries plus replicated checks in rows, columns and blocks
(Federer, Nair and Raghavarao 1975; Piepho and Williams 2016).

- [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md)
  : Generate an augmented row-column field design
- [`plot_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/plot_augmented_row_column_design.md)
  : Plot an augmented row-column design grid
- [`analyze_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_row_column_design.md)
  : Analyze an augmented row-column design
- [`print(`*`<augmented_row_column_analysis>`*`)`](https://solgenomics.github.io/breedtools/reference/print.augmented_row_column_analysis.md)
  : Print augmented row-column analysis

## Augmented p-rep designs

Multi-location partially replicated designs (Williams, Piepho and
Whitaker 2011).

- [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md)
  : Generate an augmented p-rep design
- [`augmented_prep_concurrence()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_concurrence.md)
  : Calculate pairwise concurrence for a p-rep design
- [`plot_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/plot_augmented_prep_design.md)
  : Plot an augmented p-rep design grid
- [`analyze_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_design.md)
  : Analyze an augmented p-rep design
- [`analyze_augmented_prep_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md)
  : Two-stage analysis of a multi-location augmented p-rep trial
- [`print(`*`<augmented_prep_analysis>`*`)`](https://solgenomics.github.io/breedtools/reference/print.augmented_prep_analysis.md)
  : Print augmented p-rep analysis
- [`print(`*`<augmented_prep_two_stage>`*`)`](https://solgenomics.github.io/breedtools/reference/print.augmented_prep_two_stage.md)
  : Print a two-stage p-rep analysis

## Resolvable row-column designs

Fully replicated spatial designs (Williams, John and Whitaker 2006).

- [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md)
  : Generate a resolvable (spatial) row-column design
- [`resolvable_row_column_efficiency()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_efficiency.md)
  : Efficiency of a resolvable row-column design
- [`plot_resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/plot_resolvable_row_column_design.md)
  : Plot a resolvable row-column design
- [`analyze_resolvable_row_column()`](https://solgenomics.github.io/breedtools/reference/analyze_resolvable_row_column.md)
  : Analyze a resolvable row-column trial with RC, LV + RC or AR + RC
  models
- [`print(`*`<resolvable_row_column_analysis>`*`)`](https://solgenomics.github.io/breedtools/reference/print.resolvable_row_column_analysis.md)
  : Print a resolvable row-column analysis

## Analysis for any design

Trait plots, outliers, two-stage and single-stage analysis of several
trials, and the values to store from them.

- [`plot_trait()`](https://solgenomics.github.io/breedtools/reference/plot_trait.md)
  : Plot a trait: boxplot, distribution and field map
- [`find_outliers()`](https://solgenomics.github.io/breedtools/reference/find_outliers.md)
  : Find outliers in trait data (IQR rule)
- [`remove_outliers()`](https://solgenomics.github.io/breedtools/reference/remove_outliers.md)
  : Remove outliers from trait data (IQR rule)
- [`analyze_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_two_stage.md)
  : Two-stage analysis of a multi-location trial, for any design
- [`analyze_single_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_single_stage.md)
  : Single-stage analysis of several trials
- [`compare_first_stage()`](https://solgenomics.github.io/breedtools/reference/compare_first_stage.md)
  : Compare raw and first-stage adjusted values
- [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)
  : First-stage plot values from an analysis
- [`second_stage_values()`](https://solgenomics.github.io/breedtools/reference/second_stage_values.md)
  : Second-stage genotype values from analyses
- [`prediction_accuracy()`](https://solgenomics.github.io/breedtools/reference/prediction_accuracy.md)
  : Accuracy of adjusted means against known values
- [`print(`*`<two_stage_analysis>`*`)`](https://solgenomics.github.io/breedtools/reference/print.two_stage_analysis.md)
  : Print a two-stage analysis
- [`print(`*`<single_stage_analysis>`*`)`](https://solgenomics.github.io/breedtools/reference/print.single_stage_analysis.md)
  : Print a single-stage analysis

## Design checks

Verify a design and record a checklist of what was tested.

- [`check_setup()`](https://solgenomics.github.io/breedtools/reference/check_setup.md)
  : Check that everything is ready before starting
- [`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md)
  : Standard checks of a design
- [`check_section()`](https://solgenomics.github.io/breedtools/reference/check_section.md)
  : Set the section label for the next checks
- [`check_result()`](https://solgenomics.github.io/breedtools/reference/check_result.md)
  : Record and print a check
- [`check_summary()`](https://solgenomics.github.io/breedtools/reference/check_summary.md)
  : Summary of the recorded checks
- [`check_reset()`](https://solgenomics.github.io/breedtools/reference/check_reset.md)
  : Clear the recorded checks

## Breedbase connection and lookups

Log in through BrAPI, list programs, folders and locations, and check
names before designing.

- [`breedbase_connection`](https://solgenomics.github.io/breedtools/reference/breedbase_connection.md)
  : Connecting to a Breedbase instance
- [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)
  : Log in to a Breedbase server
- [`breedbase_programs()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)
  [`breedbase_folders()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)
  [`breedbase_locations()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)
  : List breeding programs, folders and locations in Breedbase
- [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
  : Check names in Breedbase before building a design
- [`breedbase_statistics()`](https://solgenomics.github.io/breedtools/reference/breedbase_statistics.md)
  : List the statistic terms of Breedbase (SGNSTAT ontology)
- [`breedbase_adjusted_trait()`](https://solgenomics.github.io/breedtools/reference/breedbase_adjusted_trait.md)
  : Find or create the Breedbase trait for adjusted values
- [`breedbase_refresh()`](https://solgenomics.github.io/breedtools/reference/breedbase_refresh.md)
  : Refresh the Breedbase phenotype cache

## Save designs to Breedbase

Convert agricolae designs, add plot names and GPS coordinates, preview
the payload, and save trials with their plots.

- [`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md)
  : Prepare an agricolae design and trial information for Breedbase
- [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md)
  : Convert an agricolae field book for saving to Breedbase
- [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md)
  : Add Breedbase plot names to a design
- [`add_plot_coordinates()`](https://solgenomics.github.io/breedtools/reference/add_plot_coordinates.md)
  : Add GPS coordinates to the plots of a design
- [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  : Save a design to Breedbase
- [`save_breedbase_coordinates()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_coordinates.md)
  : Save plot GPS coordinates of a trial already in Breedbase
- [`breedbase_trial_plots()`](https://solgenomics.github.io/breedtools/reference/breedbase_trial_plots.md)
  : Read the plots of a trial saved in Breedbase
- [`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
  : Phenotyping template with the Breedbase plot names of a design

## Phenotypes and analysis results in Breedbase

Read and save phenotypic data, first-stage plot values and second-stage
analyses.

- [`breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/breedbase_phenotypes.md)
  : Read phenotypic data of trials from Breedbase
- [`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md)
  : Save phenotypic data to Breedbase
- [`save_breedbase_first_stage()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_first_stage.md)
  : Save first-stage results to Breedbase as a trait of the field plots
- [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  : Save analysis results to Breedbase as an analysis
- [`analysis_plot_names()`](https://solgenomics.github.io/breedtools/reference/analysis_plot_names.md)
  : Plot and accession names an analysis will have, before saving
- [`breedbase_analysis_plots()`](https://solgenomics.github.io/breedtools/reference/breedbase_analysis_plots.md)
  : Plots of an analysis saved in Breedbase
- [`add_breedbase_analysis_values()`](https://solgenomics.github.io/breedtools/reference/add_breedbase_analysis_values.md)
  : Add results to an analysis that is already in Breedbase

## Multi-environment trials from Breedbase files

Read Breedbase phenotype files, test homogeneity of error variances
(Hartley’s Fmax) and group environments by GxE (Lin 1982).

- [`read_breedbase_file()`](https://solgenomics.github.io/breedtools/reference/read_breedbase_file.md)
  : Read phenotype or design data from Excel, CSV, or tab-delimited
  files
- [`read_breedbase_excel()`](https://solgenomics.github.io/breedtools/reference/read_breedbase_excel.md)
  : Read a Breedbase Excel phenotype file
- [`prepare_breedbase_pheno()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_pheno.md)
  : Prepare Breedbase phenotype data for multi-environment analysis
- [`run_breedbase_variance_check()`](https://solgenomics.github.io/breedtools/reference/run_breedbase_variance_check.md)
  : Run Breedbase multi-location phenotype variance check
- [`run_environment_stratification()`](https://solgenomics.github.io/breedtools/reference/run_environment_stratification.md)
  : Run Lin-style environment stratification
- [`run_breedbase_environment_stratification()`](https://solgenomics.github.io/breedtools/reference/run_breedbase_environment_stratification.md)
  : Run Lin-style environment stratification from a Breedbase phenotype
  file
- [`read_locations()`](https://solgenomics.github.io/breedtools/reference/read_locations.md)
  : Load location coordinates
- [`map_environment_groups()`](https://solgenomics.github.io/breedtools/reference/map_environment_groups.md)
  : Create a map of environment groups

## Variance check building blocks

Lower-level steps used by
[`run_breedbase_variance_check()`](https://solgenomics.github.io/breedtools/reference/run_breedbase_variance_check.md).

- [`analyze_environment_anova()`](https://solgenomics.github.io/breedtools/reference/analyze_environment_anova.md)
  : Run individual ANOVA by environment
- [`analyze_one_environment()`](https://solgenomics.github.io/breedtools/reference/analyze_one_environment.md)
  : Analyze one environment using ANOVA
- [`hartley_fmax_test()`](https://solgenomics.github.io/breedtools/reference/hartley_fmax_test.md)
  : Hartley Fmax test for residual variance homogeneity
- [`hartley_critical_sim()`](https://solgenomics.github.io/breedtools/reference/hartley_critical_sim.md)
  : Simulate Hartley Fmax critical value
- [`make_environment_groups()`](https://solgenomics.github.io/breedtools/reference/make_environment_groups.md)
  : Test all compatible environment groups using Hartley Fmax
- [`get_compatible_groups()`](https://solgenomics.github.io/breedtools/reference/get_compatible_groups.md)
  : Get compatible environment groups
- [`get_maximal_groups()`](https://solgenomics.github.io/breedtools/reference/get_maximal_groups.md)
  : Select maximal compatible environment groups
- [`harmonic_mean()`](https://solgenomics.github.io/breedtools/reference/harmonic_mean.md)
  : Calculate harmonic mean

## Environment stratification building blocks

Lower-level steps used by
[`run_environment_stratification()`](https://solgenomics.github.io/breedtools/reference/run_environment_stratification.md).

- [`prepare_environment_stratification_data()`](https://solgenomics.github.io/breedtools/reference/prepare_environment_stratification_data.md)
  : Prepare phenotype data for Lin environment stratification
- [`detect_design()`](https://solgenomics.github.io/breedtools/reference/detect_design.md)
  : Detect field trial design for environment stratification
- [`prepare_design_factors()`](https://solgenomics.github.io/breedtools/reference/prepare_design_factors.md)
  : Prepare design factors
- [`calculate_environment_stratification_anova()`](https://solgenomics.github.io/breedtools/reference/calculate_environment_stratification_anova.md)
  : Calculate design-aware ANOVA
- [`lin_group_environments()`](https://solgenomics.github.io/breedtools/reference/lin_group_environments.md)
  : Group environments using Lin-style GE compatibility testing
- [`lin_test_group()`](https://solgenomics.github.io/breedtools/reference/lin_test_group.md)
  : Test one Lin environment group
- [`lin_ss_ge()`](https://solgenomics.github.io/breedtools/reference/lin_ss_ge.md)
  : Calculate Lin genotype-by-environment sum of squares
- [`lin_error_mse()`](https://solgenomics.github.io/breedtools/reference/lin_error_mse.md)
  : Calculate residual error MSE for a Lin environment group
- [`complete_environment_info()`](https://solgenomics.github.io/breedtools/reference/complete_environment_info.md)
  : Complete environment metadata
- [`add_environment_metadata()`](https://solgenomics.github.io/breedtools/reference/add_environment_metadata.md)
  : Add environment metadata to result table
- [`add_pairwise_environment_metadata()`](https://solgenomics.github.io/breedtools/reference/add_pairwise_environment_metadata.md)
  : Add environment metadata to pairwise Lin results
- [`environment_summary_by_group()`](https://solgenomics.github.io/breedtools/reference/environment_summary_by_group.md)
  : Summarize environment names by group

## Column and table helpers

- [`find_column()`](https://solgenomics.github.io/breedtools/reference/find_column.md)
  : Find the first matching column
- [`guess_column()`](https://solgenomics.github.io/breedtools/reference/guess_column.md)
  : Guess a column name from possible Breedbase column names
- [`normalize_column_name()`](https://solgenomics.github.io/breedtools/reference/normalize_column_name.md)
  : Normalize Breedbase-style column names

## Package

- [`breedtools`](https://solgenomics.github.io/breedtools/reference/breedtools-package.md)
  [`breedtools-package`](https://solgenomics.github.io/breedtools/reference/breedtools-package.md)
  : breedtools: design, analysis and 'Breedbase' storage of plant
  breeding trials
