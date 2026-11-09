# Generate an augmented p-rep design

Creates a multi-location augmented partially replicated (p-rep) design
(Williams, Piepho and Whitaker, 2011). Instead of replicated checks, a
subset of the test entries is replicated at each location: every entry
is planted once at every location, and the entries chosen for p-rep at a
location get extra plots there. All plots of a location are arranged in
incomplete blocks, chosen so that pairs of entries share a block at most
once across the whole trial ((0,1)-designs).

## Usage

``` r
augmented_prep_design(
  treatments,
  locations,
  block_size = 8,
  prep_prop = NULL,
  n_duplicate_per_location = NULL,
  replications_per_repeated_entry = 2,
  n_blocks_per_location = NULL,
  resolvable = FALSE,
  n_candidates = 1000,
  seed = NULL,
  field_cols = NULL,
  rows_per_block = NULL,
  blocks_per_row = NULL,
  serpentine = TRUE,
  output_file = NULL,
  method = c("search", "alpha", "alpha2"),
  n_prep_per_block = NULL,
  n_efficiency_candidates = 5,
  trial_name = NULL,
  plot_name_pattern = "{trial}_{location}_{accession}_{block}_{rep}"
)
```

## Arguments

- treatments:

  Character vector of entry names.

- locations:

  Character vector of location names.

- block_size:

  Maximum final incomplete block size. Can be a single integer or a
  named vector indexed by location (the algebraic methods use its
  smallest value). Actual block sizes may differ by one plot when totals
  do not divide evenly.

- prep_prop:

  Optional proportion of entries to duplicate per location
  (`method = "search"` only). If `NULL`, entries are split across
  locations so each entry is duplicated at approximately one location.

- n_duplicate_per_location:

  Optional number of entries duplicated at each location
  (`method = "search"` only). Can be a single integer, a vector matching
  `locations`, or a named vector indexed by location.

- replications_per_repeated_entry:

  Total number of plots for entries selected for p-rep replication at a
  location. The default `2` gives one base plot plus one extra plot. Can
  be a single integer or a named vector indexed by location. Must be 2
  for the algebraic methods.

- n_blocks_per_location:

  Optional number of incomplete blocks at each location
  (`method = "search"` only). If `NULL`, the number of blocks is the
  minimum needed so no block is larger than `block_size`. Can be a
  single integer or a named vector indexed by location.

- resolvable:

  Logical (`method = "search"`). If `TRUE`, each location is laid out as
  a resolvable design: the blocks are grouped into
  `replications_per_repeated_entry` replicates, each replicate holds
  every p-rep entry of that location exactly once plus an even share of
  the unreplicated entries, and all replicates have the same number of
  blocks. If `n_blocks_per_location` is given it must be a multiple of
  the number of replicates. Designs from the algebraic methods are
  always resolvable.

- n_candidates:

  Maximum number of candidate designs (`"search"`) or of random starting
  arrays (`"alpha"`, `"alpha2"`) to generate.
  `summary$n_candidates_evaluated` reports how many were used.

- seed:

  Optional random seed.

- field_cols:

  Optional number of field columns for row/column layout. Can be a
  single integer or a named vector indexed by location.

- rows_per_block:

  Optional number of physical field rows each incomplete block should
  occupy. Can be a single integer or a named vector for selected
  locations. Must be used with `blocks_per_row` for those locations.

- blocks_per_row:

  Optional number of incomplete blocks to place side by side in each
  field row band. Can be a single integer or a named vector for selected
  locations. Must be used with `rows_per_block` for those locations.

- serpentine:

  Logical. If `TRUE`, row/column coordinates use serpentine ordering
  within each location.

- output_file:

  Optional tab-delimited output file for the complete design.

- method:

  Construction method: `"search"`, `"alpha"` or `"alpha2"` (see
  Details).

- n_prep_per_block:

  For the algebraic methods, the number `m` of p-rep entries in every
  block. If `NULL`, the largest feasible `m` whose blocks fit in
  `block_size` is used.

- n_efficiency_candidates:

  Number of optimal candidates (no entry repeated within a block,
  connected locations, and no pair of entries sharing more than one
  block) to compare on the average efficiency factor; the one with the
  highest efficiency is returned. Use `1` to return the first optimal
  candidate. Efficiency is compared for up to 2000 entries.

- trial_name:

  Optional Breedbase trial name. When given, a `plot_name` column with
  the plot names used by
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  is added to the design (see
  [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md)),
  so the names can be checked with
  [`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
  and changed before saving. For a multi-location design the location is
  added, one trial per location.

- plot_name_pattern:

  Pattern for the plot names, as in
  [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md).
  Used only with `trial_name`; `{location}` is the design location.

## Value

A list with design inputs, summary diagnostics, complete plot-level
design, Breedbase-compatible design columns, block summary, replication
summary, repeated-entry groups, and pairwise concurrence matrix. In the
complete design, `repeated_at` is the location where the entry receives
its extra p-rep plot(s), or `"not_repeated"` when the entry is not
selected for extra replication. `replicate` is the replicate
(superblock) of the plot for resolvable designs, and 1 otherwise; the
Breedbase `rep` column is `replicate` for resolvable designs and the
plot's occurrence number otherwise.

## Details

Three construction methods are available through `method`:

- `"search"` (default):

  A flexible randomized search. Any number of entries, locations, block
  sizes, p-rep proportions (`prep_prop`, `n_duplicate_per_location`),
  replication levels (`replications_per_repeated_entry`) and block
  counts can be used, with an optional resolvable layout (`resolvable`).

- `"alpha"`:

  The partial alpha-array construction of Section 2.1 of the paper. An
  alpha-array for `l * m` rows and `2l` replicates is searched, the
  `m x 2` rectangles on its diagonal give the two-replicate design for
  the p-rep entries of each location, and half of the remaining elements
  of each column are deleted (Section 2.2) so that the design is a
  (0,1)-design. Every entry is duplicated at exactly one location, each
  location is resolvable with two replicates of `s` blocks, and each
  block holds `m` p-rep entries; the final block size is `(l + 1) m / 2`
  (blocks of the two replicates differ by one plot when this is not an
  integer).

- `"alpha2"`:

  The augmented alpha(2)-array construction of Section 2.4, for even `m`
  and an odd number of blocks per replicate `s`. Every block has
  `(l + 1) m / 2` plots.

For the algebraic methods, `m` is set with `n_prep_per_block` (or chosen
as the largest value whose blocks fit in `block_size`), the number of
blocks and of duplicated entries follow from the construction, and
`replications_per_repeated_entry` must be 2. When the number of entries
is not a multiple of the array size, the surplus array positions are
dropped as dummy entries, so some blocks are one plot smaller.
[`augmented_prep_parameters()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_parameters.md)
lists the feasible parameter sets.

The incomplete block structure is independent from the field geometry.
By default, `field_cols` lays the ordered plots into rows and columns.
For a rectangular block layout, use `rows_per_block` and
`blocks_per_row` to place multiple blocks side by side within each
location.

With `method = "search"`, the p-rep plots of every location are spread
evenly over its blocks, and the p-rep entries are placed so that they
link all blocks of the location (the location is connected), so entry
effects can also be estimated from each location on its own. A warning
is given when there are too few p-rep plots to connect the blocks. After
the random greedy search, the best candidates are improved by swapping
plots between blocks to remove pairs of entries that share more than one
block.

For all methods, `summary` reports the maximum pairwise concurrence, the
number of locations whose blocks are not connected, and the average
efficiency factor `efficiency_factor` of the block design (Section 2.3
of the paper; higher is better), which can be used to compare designs.
Entry labels, block order within replicates and plot order within blocks
are randomized.

## References

Williams, E., Piepho, H.-P. and Whitaker, D. (2011). Augmented p-rep
designs. Biometrical Journal 53, 19-27.
[doi:10.1002/bimj.201000102](https://doi.org/10.1002/bimj.201000102)

## See also

[`augmented_prep_parameters()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_parameters.md),
[`analyze_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_design.md),
[`plot_augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/plot_augmented_prep_design.md)

## Examples

``` r
design <- augmented_prep_design(
  treatments = paste0("G", seq_len(12)),
  locations = c("Loc1", "Loc2", "Loc3"),
  block_size = c(Loc1 = 5, Loc2 = 5, Loc3 = 6),
  n_blocks_per_location = c(Loc1 = 4, Loc2 = 4, Loc3 = 4),
  n_duplicate_per_location = c(Loc1 = 4, Loc2 = 2, Loc3 = 4),
  replications_per_repeated_entry = c(Loc1 = 2, Loc2 = 3, Loc3 = 2),
  n_candidates = 5,
  seed = 123,
  field_cols = 8
)

design$summary
#> # A tibble: 1 × 20
#>       score max_pairwise_concurr…¹ n_pairs_with_concurr…² total_excess_concurr…³
#>       <dbl>                  <int>                  <int>                  <dbl>
#> 1 130002004                      2                     13                     13
#> # ℹ abbreviated names: ¹​max_pairwise_concurrence,
#> #   ²​n_pairs_with_concurrence_gt_1, ³​total_excess_concurrence
#> # ℹ 16 more variables: same_entry_same_block <dbl>,
#> #   n_disconnected_locations <int>, efficiency_factor <dbl>, resolvable <lgl>,
#> #   min_replication <int>, max_replication <int>, mean_replication <dbl>,
#> #   n_locations <int>, n_entries <int>, n_plots <int>, n_blocks <int>,
#> #   n_candidates_evaluated <int>, method <chr>, m <int>, s <int>, …
head(design$design)
#> # A tibble: 6 × 16
#>   plots location   row   col replicate block block_uid block_size plot_in_block
#>   <int> <chr>    <dbl> <int>     <int> <int> <chr>          <int>         <int>
#> 1     1 Loc1         1     1         1     1 Loc1_B001          4             1
#> 2     2 Loc1         1     2         1     1 Loc1_B001          4             2
#> 3     3 Loc1         1     3         1     1 Loc1_B001          4             3
#> 4     4 Loc1         1     4         1     1 Loc1_B001          4             4
#> 5     5 Loc1         1     5         1     2 Loc1_B002          4             1
#> 6     6 Loc1         1     6         1     2 Loc1_B002          4             2
#> # ℹ 7 more variables: all_entries <chr>, occurrence <int>, entry_type <chr>,
#> #   is_extra_repeat <int>, repeated_at <chr>, rep <int>, is_control <int>

# Example 1 of Williams, Piepho and Whitaker (2011): 20 entries, two
# locations, blocks of 3 built from a partial alpha-array.
alpha <- augmented_prep_design(
  treatments = paste0("G", 1:20),
  locations = c("L1", "L2"),
  method = "alpha",
  block_size = 3,
  seed = 1
)
alpha$summary[, c("method", "m", "s", "max_pairwise_concurrence", "efficiency_factor")]
#> # A tibble: 1 × 5
#>   method     m     s max_pairwise_concurrence efficiency_factor
#>   <chr>  <int> <int>                    <int>             <dbl>
#> 1 alpha      2     5                        1             0.609
```
