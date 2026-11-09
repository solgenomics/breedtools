# breedtools: design, analysis and 'Breedbase' storage of plant breeding trials

breedtools designs and analyses plant breeding field trials, from early
generations with a single plot per entry to fully replicated trials, and
saves the designs into a 'Breedbase' database. The tutorial,
[`vignette("breedtools")`](https://solgenomics.github.io/breedtools/articles/breedtools.md),
walks through every design and analysis.

## Which design?

Start with
[`design_options()`](https://solgenomics.github.io/breedtools/reference/design_options.md)
when the number of entries is known but the field dimensions, checks or
replication are still to be decided.

- One plot per entry, one location: augmented row-column design with
  replicated checks.

- One plot per entry per location, several locations: augmented
  partially replicated (p-rep) design, without checks.

- Every entry in every replicate: resolvable (spatial) row-column
  design.

## Augmented row-column designs

Unreplicated entries plus replicated checks in rows, columns and blocks
(Federer, Nair and Raghavarao, 1975; Piepho and Williams, 2016), with
equal or unequal row and column groups, complete or incomplete blocks of
checks, and Latin or D-optimal check placement that keeps checks apart.

- Plan:
  [`augmented_row_column_layouts()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_layouts.md),
  [`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md)

- Design:
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md)

- Map:
  [`plot_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/plot_augmented_row_column_design.md)

- Analysis (fixed, random blocks, random entries, SpATS spatial):
  [`analyze_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_row_column_design.md)

## Augmented p-rep designs

Multi-location trials in which a share of the entries is repeated at
each location (Williams, Piepho and Whitaker, 2011), built by search or
from alpha and alpha(2) arrays.

- Plan:
  [`augmented_prep_parameters()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_parameters.md)

- Design:
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
  [`augmented_prep_concurrence()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_concurrence.md)

- Map:
  [`plot_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/plot_augmented_prep_design.md)

- Analysis (with replicates and fixed or random genotype x location):
  [`analyze_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_design.md)

- Two-stage weighted analysis with block or SpATS spatial stage 1:
  [`analyze_augmented_prep_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md)

## Resolvable row-column designs

Fully replicated designs optimized for row-column and spatial (linear
variance) efficiency, with latinized options (Williams, John and
Whitaker, 2006).

- Design:
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md)

- Efficiency, neighbours, AR criterion:
  [`resolvable_row_column_efficiency()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_efficiency.md)

- Map:
  [`plot_resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/plot_resolvable_row_column_design.md)

- Analysis (RC, LV + RC and AR + RC by REML):
  [`analyze_resolvable_row_column()`](https://solgenomics.github.io/breedtools/reference/analyze_resolvable_row_column.md)

## Saving to Breedbase

Designs are saved as trials with their plots into a breeding program and
folder through BrAPI, after checking every name; after harvest, the
phenotypic data are saved to the same plots.

- Log in:
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)

- Look up names:
  [`breedbase_programs()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md),
  [`breedbase_folders()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md),
  [`breedbase_locations()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)

- Check location and accessions before designing:
  [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)

- Convert designs generated with agricolae:
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md)

- Prepare an agricolae trial with its name, location and program:
  [`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md)

- Save (with a dry run first):
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)

- Phenotypes:
  [`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
  (sheet with the plot names),
  [`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md)
  (values by plot and trait name)

- Analysis results:
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)
  and
  [`save_breedbase_first_stage()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_first_stage.md)
  (plot values as a trait of the field plots, trait made by
  [`breedbase_adjusted_trait()`](https://solgenomics.github.io/breedtools/reference/breedbase_adjusted_trait.md));
  [`second_stage_values()`](https://solgenomics.github.io/breedtools/reference/second_stage_values.md)
  and
  [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  (genotype values as a Breedbase analysis)

## Multi-environment trials from Breedbase

- Read files:
  [`read_breedbase_file()`](https://solgenomics.github.io/breedtools/reference/read_breedbase_file.md),
  [`prepare_breedbase_pheno()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_pheno.md)

- Error variances (Hartley's Fmax):
  [`run_breedbase_variance_check()`](https://solgenomics.github.io/breedtools/reference/run_breedbase_variance_check.md),
  [`hartley_fmax_test()`](https://solgenomics.github.io/breedtools/reference/hartley_fmax_test.md)

- Environment groups (Lin, 1982):
  [`run_environment_stratification()`](https://solgenomics.github.io/breedtools/reference/run_environment_stratification.md),
  [`lin_group_environments()`](https://solgenomics.github.io/breedtools/reference/lin_group_environments.md)

- Maps:
  [`read_locations()`](https://solgenomics.github.io/breedtools/reference/read_locations.md),
  [`map_environment_groups()`](https://solgenomics.github.io/breedtools/reference/map_environment_groups.md)

## Trait data and two-stage analysis

- Look:
  [`plot_trait()`](https://solgenomics.github.io/breedtools/reference/plot_trait.md);
  outliers:
  [`find_outliers()`](https://solgenomics.github.io/breedtools/reference/find_outliers.md),
  [`remove_outliers()`](https://solgenomics.github.io/breedtools/reference/remove_outliers.md)

- Two-stage analysis for any design:
  [`analyze_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_two_stage.md);
  single stage:
  [`analyze_single_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_single_stage.md)

- Phenotypes of trials from Breedbase:
  [`breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/breedbase_phenotypes.md)

- Raw versus first-stage values:
  [`compare_first_stage()`](https://solgenomics.github.io/breedtools/reference/compare_first_stage.md)

## Checks

[`check_setup()`](https://solgenomics.github.io/breedtools/reference/check_setup.md),
[`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md),
[`check_result()`](https://solgenomics.github.io/breedtools/reference/check_result.md),
[`check_summary()`](https://solgenomics.github.io/breedtools/reference/check_summary.md),
[`prediction_accuracy()`](https://solgenomics.github.io/breedtools/reference/prediction_accuracy.md).

## See also

Useful links:

- <https://github.com/solgenomics/breedtools>

- <https://solgenomics.github.io/breedtools/>

- Report bugs at <https://github.com/solgenomics/breedtools/issues>

## Author

**Maintainer**: Christiano Costa Simoes <ccs263@cornell.edu>

Authors:

- Christiano Costa Simoes <ccs263@cornell.edu>
