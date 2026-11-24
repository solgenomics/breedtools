#' breedtools: design, analysis and 'Breedbase' storage of plant breeding trials
#'
#' breedtools designs and analyses plant breeding field trials, from early
#' generations with a single plot per entry to fully replicated trials, and
#' saves the designs into a 'Breedbase' database. The tutorial,
#' `vignette("breedtools")`, walks through every design and analysis.
#' For one complete workflow with offline Field Book collection, start with
#' `vignette("getting-started", package = "breedtools")`.
#'
#' @section Which design?:
#' Start with [design_options()] when the number of entries is known but
#' the field dimensions, checks or replication are still to be decided.
#' \itemize{
#'   \item One plot per entry, one location: augmented row-column design
#'     with replicated checks.
#'   \item One plot per entry per location, several locations: augmented
#'     partially replicated (p-rep) design, without checks.
#'   \item Every entry in every replicate: resolvable (spatial) row-column
#'     design.
#' }
#'
#' @section Augmented row-column designs:
#' Unreplicated entries plus replicated checks in rows, columns and blocks
#' (Federer, Nair and Raghavarao, 1975; Piepho and Williams, 2016), with
#' equal or unequal row and column groups, complete or incomplete blocks of
#' checks, and Latin or D-optimal check placement that keeps checks apart.
#' \itemize{
#'   \item Plan: [augmented_row_column_layouts()], [augmented_row_column_capacity()]
#'   \item Design: [augmented_row_column_design()]
#'   \item Map: [plot_augmented_row_column_design()]
#'   \item Analysis (fixed, random blocks, random entries, SpATS spatial):
#'     [analyze_augmented_row_column_design()]
#' }
#'
#' @section Augmented p-rep designs:
#' Multi-location trials in which a share of the entries is repeated at each
#' location (Williams, Piepho and Whitaker, 2011), built by search or from
#' alpha and alpha(2) arrays.
#' \itemize{
#'   \item Plan: [augmented_prep_parameters()]
#'   \item Design: [augmented_prep_design()], [augmented_prep_concurrence()]
#'   \item Map: [plot_augmented_prep_design()]
#'   \item Analysis (with replicates and fixed or random genotype x location):
#'     [analyze_augmented_prep_design()]
#'   \item Two-stage weighted analysis with block or SpATS spatial stage 1:
#'     [analyze_augmented_prep_two_stage()]
#' }
#'
#' @section Resolvable row-column designs:
#' Fully replicated designs optimized for row-column and spatial (linear
#' variance) efficiency, with latinized options (Williams, John and
#' Whitaker, 2006).
#' \itemize{
#'   \item Design: [resolvable_row_column_design()]
#'   \item Efficiency, neighbours, AR criterion: [resolvable_row_column_efficiency()]
#'   \item Map: [plot_resolvable_row_column_design()]
#'   \item Analysis (RC, LV + RC and AR + RC by REML): [analyze_resolvable_row_column()]
#' }
#'
#' @section Saving to Breedbase:
#' Designs are saved as trials with their plots into a breeding program and
#' folder through BrAPI, after checking every name; after harvest, the
#' phenotypic data are saved to the same plots.
#' \itemize{
#'   \item Log in: [breedbase_login()]
#'   \item Look up names: [breedbase_programs()], [breedbase_folders()], [breedbase_locations()]
#'   \item Check location and accessions before designing: [breedbase_check()]
#'   \item Convert designs generated with agricolae: [as_breedbase_design()]
#'   \item Prepare an agricolae trial with its name, location and program:
#'     [prepare_breedbase_trial()]
#'   \item Save (with a dry run first): [save_breedbase()]
#'   \item Phenotypes: [phenotype_template()] (sheet with the plot names),
#'     [save_breedbase_phenotypes()] (values by plot and trait name)
#'   \item Offline Field Book collection: [phenotype_template()] (field CSV
#'     with no trait columns), [trait_template()] (trait definitions as `.trt`)
#'   \item Analysis results: [first_stage_values()] and
#'     [save_breedbase_first_stage()] (plot values as a trait of the field
#'     plots, trait made by [breedbase_adjusted_trait()]);
#'     [second_stage_values()] and [save_breedbase_analysis()]
#'     (genotype values as a Breedbase analysis)
#' }
#'
#' @section Multi-environment trials from Breedbase:
#' \itemize{
#'   \item Read files: [read_breedbase_file()], [prepare_breedbase_pheno()]
#'   \item Error variances (Hartley's Fmax): [run_breedbase_variance_check()], [hartley_fmax_test()]
#'   \item Environment groups (Lin, 1982): [run_environment_stratification()],
#'     [lin_group_environments()]
#'   \item Maps: [read_locations()], [map_environment_groups()]
#' }
#'
#' @section Trait data and two-stage analysis:
#' \itemize{
#'   \item Look: [plot_trait()]; outliers: [find_outliers()], [remove_outliers()]
#'   \item Two-stage analysis for any design: [analyze_two_stage()];
#'     single stage: [analyze_single_stage()]
#'   \item Phenotypes of trials from Breedbase: [breedbase_phenotypes()]
#'   \item Raw versus first-stage values: [compare_first_stage()]
#' }
#'
#' @section Checks:
#' [check_setup()], [check_design()], [check_result()], [check_summary()],
#' [prediction_accuracy()].
#'
"_PACKAGE"
