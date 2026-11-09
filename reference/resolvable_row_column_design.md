# Generate a resolvable (spatial) row-column design

Builds a resolvable row-column design for `v <= k * s` treatments in `r`
replicates, each a `k x s` array of rows and columns that contains every
treatment once (Williams, John and Whitaker, 2006). The design is found
by a search over swaps of treatments within replicates:

1.  The average efficiency factor of the resolvable row-column design
    (rows and columns within replicates), together with the efficiency
    factors of the component row design and column design, is maximized.

2.  With `spatial = TRUE`, the spatial average efficiency factor of the
    two-dimensional linear variance (LV + RC) model is then maximized
    while keeping the three efficiency factors of step 1 within
    `efficiency_tolerance` of their optimized values. This spreads
    treatments over edge and interior plots and over row and column
    neighbours, without needing prior values for spatial parameters.

When replicates are contiguous, `latinize = TRUE` builds a latinized
design: no treatment appears twice in the same long column (replicates
stacked vertically) or long row (replicates side by side), within every
group of `latin_t` adjacent replicates (t-latinized designs; John and
Williams, 1998).

## Usage

``` r
resolvable_row_column_design(
  treatments,
  n_reps,
  rows_per_rep,
  cols_per_rep,
  spatial = TRUE,
  latinize = FALSE,
  replicate_layout = c("vertical", "horizontal"),
  latin_t = NULL,
  n_starts = 3,
  n_iter = NULL,
  efficiency_tolerance = 0.002,
  plot_type = c("serpentine", "cartesian"),
  seed = NULL,
  output_file = NULL,
  start_design = NULL,
  trial_name = NULL,
  location = NULL,
  plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}"
)
```

## Arguments

- treatments:

  Character vector of treatment names (`v <= k * s`).

- n_reps:

  Number of replicates (`r`).

- rows_per_rep:

  Number of rows in each replicate (`k`).

- cols_per_rep:

  Number of columns in each replicate (`s`).

- spatial:

  Logical. If `TRUE` (default), optimize the spatial efficiency factor
  in a second stage.

- latinize:

  Logical. If `TRUE`, build a latinized design for contiguous
  replicates.

- replicate_layout:

  How replicates are placed in the field: `"vertical"` (stacked, giving
  long columns) or `"horizontal"` (side by side, giving long rows).

- latin_t:

  Optional size of the groups of adjacent replicates within which
  latinization holds (t-latinized design). Defaults to all replicates.

- n_starts:

  Number of random starting designs.

- n_iter:

  Number of swap proposals per stage and start. Defaults to
  `min(20000, max(3000, 60 * v * r))`.

- efficiency_tolerance:

  Largest loss in the row-column, row and column efficiency factors
  accepted in the spatial stage.

- plot_type:

  Plot numbering order: `"serpentine"` or `"cartesian"`.

- seed:

  Optional random seed.

- output_file:

  Optional tab-delimited output file for the design.

- start_design:

  Optional existing design to start from (for example one produced by
  CycDesigN), as a data frame with `replicate`, `row_in_rep`,
  `col_in_rep` and `all_entries` columns and complete replicates. The
  search then improves this design; use `n_iter = 1` and
  `spatial = FALSE` to keep it essentially as it is and only lay it out.

- trial_name:

  Optional Breedbase trial name. When given, a `plot_name` column with
  the plot names used by
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  is added to the design (see
  [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md)),
  so the names can be checked with
  [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
  and changed before saving.

- location:

  Breedbase location of the trial, for the `{location}` field of the
  plot names. Used only with `trial_name`.

- plot_name_pattern:

  Pattern for the plot names, as in
  [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md).
  Used only with `trial_name`.

## Value

An object of class `resolvable_row_column_design`: a list with `input`,
`summary` (efficiency factors, spatial upper bound, latinization check),
`design` (one row per plot with `replicate`, `row_in_rep`, `col_in_rep`,
field `row` and `col`, `long_row` / `long_col`, `all_entries`,
`is_filler`), `breedbase_design`, and `efficiency` (see
[`resolvable_row_column_efficiency()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_efficiency.md)).

## Details

When there are fewer treatments than plots per replicate, the remaining
plots are fillers (`all_entries = "FILLER"`, `is_filler = 1`), placed by
the search like an extra treatment.

## References

Williams, E. R., John, J. A. and Whitaker, D. (2006). Construction of
resolvable spatial row-column designs. Biometrics 62, 103-108.
[doi:10.1111/j.1541-0420.2005.00393.x](https://doi.org/10.1111/j.1541-0420.2005.00393.x)

John, J. A. and Williams, E. R. (1998). t-Latinized designs. Australian
and New Zealand Journal of Statistics 40, 111-118.

## See also

[`resolvable_row_column_efficiency()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_efficiency.md),
[`analyze_resolvable_row_column()`](https://solgenomics.github.io/breedtools/reference/analyze_resolvable_row_column.md),
[`plot_resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/plot_resolvable_row_column_design.md)

## Examples

``` r
# Example 2 of Williams, John and Whitaker (2006): 42 treatments in three
# replicates of 6 rows x 7 columns.
design <- resolvable_row_column_design(
  treatments = paste0("T", 1:42),
  n_reps = 3,
  rows_per_rep = 6,
  cols_per_rep = 7,
  n_starts = 1,
  n_iter = 2000,
  seed = 1
)
design$summary
#> # A tibble: 1 × 15
#>   n_treatments n_fillers n_reps rows_per_rep cols_per_rep n_plots efficiency_rc
#>          <int>     <int>  <int>        <int>        <int>   <int>         <dbl>
#> 1           42         0      3            6            7     126         0.684
#> # ℹ 8 more variables: efficiency_rows <dbl>, efficiency_cols <dbl>,
#> #   efficiency_spatial <dbl>, upper_bound_spatial <dbl>, spatial <lgl>,
#> #   latinized <lgl>, latin_ok <lgl>, replicate_layout <chr>
```
