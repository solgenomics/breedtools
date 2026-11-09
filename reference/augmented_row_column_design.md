# Generate an augmented row-column field design

Randomizes unreplicated treatments and replicated controls (checks) in
an augmented row-column design with superimposed blocks (Piepho and
Williams, 2016). Adjacent rows form row groups and adjacent columns form
column groups; each intersection is a block that receives a given number
of check plots. Row, column and block effects are all estimable from the
checks, and the layout must leave positive error degrees of freedom (see
[`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md)).

## Usage

``` r
augmented_row_column_design(
  treatments,
  controls,
  rows_in_field = NULL,
  cols_in_field = NULL,
  rows_per_block = NULL,
  cols_per_block = NULL,
  plot_type = c("serpentine", "cartesian"),
  n_candidates = 1000,
  seed = NULL,
  output_file = NULL,
  rows_per_rowgroup = NULL,
  cols_per_colgroup = NULL,
  checks_per_block = NULL,
  method = c("auto", "latin", "optimal"),
  trial_name = NULL,
  location = NULL,
  plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}"
)
```

## Arguments

- treatments:

  Character vector of unreplicated treatment/entry names. Its length
  must equal the number of entry plots.

- controls:

  Character vector of replicated check/control names.

- rows_in_field:

  Number of rows in the whole field.

- cols_in_field:

  Number of columns in the whole field.

- rows_per_block:

  Number of rows in each row group (equal groups).

- cols_per_block:

  Number of columns in each column group (equal groups).

- plot_type:

  Plot numbering order. Use `"serpentine"` or `"cartesian"`.

- n_candidates:

  Number of candidate layouts (`"latin"`) or random restarts of the
  local search (`"optimal"`, at most 50; the search stops after 3
  restarts without improvement).

- seed:

  Optional random seed.

- output_file:

  Optional path for writing the complete design as a tab-delimited file.

- rows_per_rowgroup, cols_per_colgroup:

  Optional vectors of unequal group sizes, as in
  [`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md).

- checks_per_block:

  Optional number of check plots per block, as in
  [`augmented_row_column_capacity()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_capacity.md).

- method:

  Check allocation method: `"auto"`, `"latin"` or `"optimal"`.

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

A list with:

- input:

  Validated design inputs.

- summary:

  Design size, error df, method, and the D-criterion `d_criterion`
  (`det(X'X)^(1/p)` of the check-plot model; higher is better), and the
  number of pairs of neighbouring check plots (same column in
  consecutive rows, or same row in consecutive columns) inside a block
  and across block boundaries.

- design:

  Complete plot-level design with row, column, block, entry, control
  flag, and plot order.

- breedbase_design:

  Breedbase-compatible columns: `plots`, `block`, `all_entries`, `rep`,
  and `is_control`.

## Details

Two check-allocation methods are available:

- `"latin"`:

  For equal groups where every block is a complete block of checks (one
  copy of each). Checks are placed with covering Latin rectangles, so
  controls never share a row or column within a block, the same control
  never appears twice in a field row or column, and every field row and
  column has a check. Candidates are scored for spatial spread and
  row/column balance.

- `"optimal"`:

  The general strategy of Piepho and Williams (2016): unequal group
  sizes and any number of check plots per block (incomplete blocks or
  extra replication of checks). Check labels are first allocated to
  blocks balancing row groups, column groups and blocks (Step 2); check
  plots are then placed within blocks and the layout improved by a local
  search that maximizes the D-criterion of the check-plot model
  `row + column + rowgroup:colgroup + check` (Steps 3-5). Check plots
  next to each other (same column in consecutive rows, or same row in
  consecutive columns) are avoided, inside blocks above all, and spatial
  spread breaks ties. Entries fill the remaining plots at random (Step
  6).

`method = "auto"` uses `"latin"` when it applies and `"optimal"`
otherwise.

## References

Piepho, H.-P. and Williams, E. R. (2016). Augmented row-column designs
for a small number of checks. Agronomy Journal 108, 2256-2262.
[doi:10.2134/agronj2016.06.0325](https://doi.org/10.2134/agronj2016.06.0325)

## Examples

``` r
treatments <- paste0("G", seq_len(108))
controls <- paste0("Check", seq_len(4))

design <- augmented_row_column_design(
  treatments = treatments,
  controls = controls,
  rows_in_field = 12,
  cols_in_field = 12,
  rows_per_block = 4,
  cols_per_block = 4,
  n_candidates = 50,
  seed = 123
)

head(design$design)
#>   plots row col block rowgroup colgroup all_entries  type rep is_control
#> 1     1   1   1     1        1        1         G16 entry   1          0
#> 2     2   1   2     1        1        1         G74 entry   1          0
#> 3     3   1   3     1        1        1         G63 entry   1          0
#> 4     4   1   4     1        1        1      Check3 check   1          1
#> 5     5   1   5     2        1        2         G88 entry   1          0
#> 6     6   1   6     2        1        2         G44 entry   1          0
design$summary
#> # A tibble: 1 × 25
#>   rows_in_field cols_in_field rows_per_block cols_per_block n_block_rows
#>           <int>         <int>          <int>          <int>        <int>
#> 1            12            12              4              4            3
#> # ℹ 20 more variables: n_block_cols <int>, n_blocks <int>, n_controls <int>,
#> #   total_plots <int>, check_plots <int>, entry_plots <int>, error_df <dbl>,
#> #   check_fraction <dbl>, min_rowgroup_slack <dbl>, min_colgroup_slack <dbl>,
#> #   every_row_has_check <lgl>, every_col_has_check <lgl>, method <chr>,
#> #   plot_type <chr>, n_candidates <int>, n_valid_candidates <int>,
#> #   design_score <dbl>, d_criterion <dbl>, adjacent_checks_in_blocks <int>,
#> #   adjacent_checks_across_blocks <int>

# Example 3 of Piepho and Williams (2016): 173 entries, 4 checks, unequal
# row and column groups, and 3 to 5 check plots per block.
pb <- rbind(c(4, 5, 4), c(3, 5, 3), c(3, 5, 4), c(3, 5, 3))
ex3 <- augmented_row_column_design(
  treatments = paste0("E", 1:173),
  controls = c("A", "B", "C", "D"),
  rows_per_rowgroup = c(6, 5, 6, 5),
  cols_per_colgroup = c(3, 4, 3),
  checks_per_block = pb,
  n_candidates = 3,
  seed = 1
)
ex3$summary[, c("method", "check_plots", "entry_plots", "error_df", "d_criterion")]
#> # A tibble: 1 × 5
#>   method  check_plots entry_plots error_df d_criterion
#>   <chr>         <int>       <int>    <dbl>       <dbl>
#> 1 optimal          47         173        7        1.94
```
