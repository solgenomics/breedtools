# Add Breedbase plot names to a design

Adds a `plot_name` column to `design$design` with the names the plots
will have in Breedbase, so they can be seen, checked with
[`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
and changed before
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md).
The design functions call this when given a `trial_name`. Run it again
to make the names for another trial name or pattern, or edit the column
by hand: any unique names are accepted.
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md),
[`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
and
[`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)
use the stored names.

## Usage

``` r
add_plot_names(
  design,
  trial_name = NULL,
  location = NULL,
  plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}",
  plot_names = NULL
)
```

## Arguments

- design:

  A design from
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md),
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md)
  or
  [`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md).

- trial_name:

  The Breedbase trial name. For a multi-location p-rep design, the
  location is added as in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  (one trial per location). Can be omitted for a prepared Breedbase
  trial.

- location:

  The Breedbase location of the trial, used by the `{location}` field of
  the pattern. Not needed for p-rep designs (their own locations are
  used, or give a named vector mapping them to Breedbase locations, as
  in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md))
  or for a prepared Breedbase trial (its recorded location is used).

- plot_name_pattern:

  Pattern for the names, as in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md):
  fields `{trial}`, `{location}`, `{accession}`, `{block}`, `{rep}`,
  `{plot}`, `{row}` and `{col}`.

- plot_names:

  Optional names given by you instead of a pattern: one per row of
  `design$design`, or the name of a column holding them.

## Value

The design with a `plot_name` column (`NA` for filler plots, which are
not saved) and the trial name in `plot_name_trial`.

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
rc <- add_plot_names(rc, "AYT_2026", location = "Ibadan")
head(rc$design$plot_name)
#> [1] "AYT_2026_Ibadan_G31_1_1"  "AYT_2026_Ibadan_G70_1_1" 
#> [3] "AYT_2026_Ibadan_G104_1_1" "AYT_2026_Ibadan_C4_1_1"  
#> [5] "AYT_2026_Ibadan_G66_2_1"  "AYT_2026_Ibadan_G5_2_1"  
# Other names, e.g. with the field position:
rc <- add_plot_names(rc, "AYT_2026_Ibadan", plot_name_pattern = "{trial}_{plot}_r{row}c{col}")
head(rc$design$plot_name)
#> [1] "AYT_2026_Ibadan_1_r1c1" "AYT_2026_Ibadan_2_r1c2" "AYT_2026_Ibadan_3_r1c3"
#> [4] "AYT_2026_Ibadan_4_r1c4" "AYT_2026_Ibadan_5_r1c5" "AYT_2026_Ibadan_6_r1c6"
```
