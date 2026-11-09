# First-stage plot values from an analysis

Turns analyses of a trial into one value per field plot, ready for
[`save_breedbase_first_stage()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_first_stage.md):
the genotype's adjusted mean plus the plot's residual. This removes the
design effects (blocks, rows, columns, replicates) and the spatial
trend, and keeps the plot-to-plot variation, so the result has as many
plots as the field trial. With SpATS (`spatial = TRUE` or
`stage1 = "spatial"`) these are spatially corrected plot values.

## Usage

``` r
first_stage_values(
  fits,
  design = NULL,
  trial_name = NULL,
  trait = NULL,
  plot_name_pattern = NULL,
  plot_names = NULL,
  breedbase_plots = NULL,
  location = NULL
)
```

## Arguments

- fits:

  One analysis, or a named list of analyses, one per trait (the names
  are the Breedbase trait names). Analyses from
  [`analyze_augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_row_column_design.md),
  [`analyze_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_design.md),
  [`analyze_augmented_prep_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md)
  (stage-1 fits per location) or
  [`analyze_resolvable_row_column()`](https://solgenomics.github.io/breedtools/reference/analyze_resolvable_row_column.md).

- design:

  The design that was saved to Breedbase. Leave it `NULL` for data read
  from Breedbase
  ([`breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/breedbase_phenotypes.md)):
  the plot names are then taken from the data.

- trial_name:

  The trial name used when saving the design (not needed without a
  design).

- trait:

  Trait name when `fits` is a single analysis.

- plot_name_pattern, plot_names, location:

  The plot naming used in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md);
  `location` is the Breedbase location of the trial, needed for the
  `{location}` field when the design has no stored plot names.

- breedbase_plots:

  Optional plots of the saved trial, from
  [`breedbase_trial_plots()`](https://solgenomics.github.io/breedtools/reference/breedbase_trial_plots.md).
  The design is then checked against the saved trial (the same accession
  at every plot number), and the plot names are taken exactly as saved.

## Value

A tibble with `plot_name`, `trial_name`, `accession`, `plot_number`,
`block`, `rep`, `row`, `col`, `is_check` and one column per trait. Its
attribute `"method"` is `"spatial"` when any analysis used SpATS,
otherwise `"mixed"`.

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
d <- rc$design
set.seed(1)
d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
fit <- analyze_augmented_row_column_design(d, "yield")
values <- first_stage_values(fit, rc, "AYT_2026_Ibadan", trait = "fresh root yield",
                            location = "Ibadan")
head(values)
#> # A tibble: 6 × 10
#>   plot_name    trial_name accession plot_number block   rep   row   col is_check
#>   <chr>        <chr>      <chr>           <int> <int> <int> <int> <int> <lgl>   
#> 1 AYT_2026_Ib… AYT_2026_… G51                 1     1     1     1     1 FALSE   
#> 2 AYT_2026_Ib… AYT_2026_… G64                 2     1     1     1     2 FALSE   
#> 3 AYT_2026_Ib… AYT_2026_… G97                 3     1     1     1     3 FALSE   
#> 4 AYT_2026_Ib… AYT_2026_… C3                  4     1     1     1     4 TRUE    
#> 5 AYT_2026_Ib… AYT_2026_… G8                  5     2     1     1     5 FALSE   
#> 6 AYT_2026_Ib… AYT_2026_… C2                  6     2     2     1     6 TRUE    
#> # ℹ 1 more variable: `fresh root yield` <dbl>
```
