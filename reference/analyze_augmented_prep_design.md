# Analyze an augmented p-rep design

Fits the fixed-effects incomplete-block model:
`y = mean + location + block(location) + genotype + error`. For
resolvable designs (a `replicate` column with more than one replicate
per location) the model is
`y = mean + location + replicate(location) + block(replicate, location) + genotype + error`.
With `genotype_by_location = "fixed"` (or `TRUE`) a fixed
genotype-by-location interaction is added after genotypes; the error is
then estimated only from the p-rep entries, as in a per-location
analysis with pooled error (Williams, Piepho and Whitaker, 2011, Section
3). With `genotype_by_location = "random"` the interaction is a random
effect, as recommended by the paper when locations are a sample of
environments, and the model is fitted by REML with 'lme4'. Sums of
squares are sequential, in the order above. Degrees of freedom in the
formatted ANOVA are those estimable from the observed (non-missing)
plots, so F-tests stay valid when plots are missing. For a complete,
connected design they equal the design formulas reported in `anova_df`,
which use total plots (`n`), locations (`l`), incomplete blocks (`b`),
and genotypes (`v`). Terms with a single level in the data (for example
the location term in a single-location trial) are left out of the model
and shown with 0 df.

## Usage

``` r
analyze_augmented_prep_design(
  data,
  response,
  genotype_col = "all_entries",
  location_col = "location",
  block_col = "block",
  entry_type_col = "entry_type",
  extra_repeat_col = "is_extra_repeat",
  repeated_at_col = "repeated_at",
  control_col = "is_control",
  replicate_col = "replicate",
  genotype_by_location = FALSE
)
```

## Arguments

- data:

  A complete augmented p-rep design data frame, usually `design$design`
  from
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
  or an object returned by
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
  with a response column added.

- response:

  Name of the numeric response column.

- genotype_col:

  Name of the genotype/treatment column.

- location_col:

  Name of the location column.

- block_col:

  Name of the incomplete block column. Block numbers may repeat across
  locations; blocks are treated as nested within locations.

- entry_type_col:

  Optional entry type column. Used only for genotype summaries.

- extra_repeat_col:

  Optional extra-repeat indicator column. Used only for genotype
  summaries.

- repeated_at_col:

  Optional column identifying where each entry is repeated. Used only
  for genotype summaries.

- control_col:

  Optional control flag column. Used only for genotype summaries.

- replicate_col:

  Optional replicate (superblock) column. When present and a location
  has more than one replicate, replicates nested within locations are
  fitted, and blocks are nested within replicates. Set to `NULL` to
  ignore it.

- genotype_by_location:

  Genotype-by-location interaction in a multi-location analysis:
  `"none"` (default, or `FALSE`), `"fixed"` (or `TRUE`), or `"random"`
  (needs the 'lme4' package). With `"fixed"`, adjusted means are
  averaged over locations. With `"random"`, adjusted means are the
  genotype BLUEs from the mixed model, `variance_components` gives the
  REML estimates, and the ANOVA table is the fixed-interaction table,
  whose genotype-by-location and error mean squares show the size of the
  interaction.

## Value

An object of class `augmented_prep_analysis`: a list with

- input, formula, model:

  The column mapping, model formula, and the `lm` fit.

- anova, anova_df, anova_formatted:

  Sequential ANOVA (estimable df), design df formulas, and a formatted
  ANOVA table.

- plot_results:

  The input data with fitted values and residuals.

- genotype_summary:

  One row per genotype with raw, fitted, and adjusted means and the
  standard error of the adjusted mean.

- adjusted_means:

  Genotypes ranked by adjusted mean (best first).

- adjusted_means_vcov:

  Covariance matrix of the estimable adjusted means, used by
  [`analyze_augmented_prep_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_augmented_prep_two_stage.md)
  for weighting.

- diagnostics:

  Error df, residual standard deviation, R-squared, `mean_sed` (mean
  standard error of a difference between two adjusted means), `lsd_5pct`
  (the matching least significant difference), and flags for rank
  deficiency or a saturated fit.

- variance_components, mixed_model:

  For `genotype_by_location = "random"`: REML variance components and
  the `lmerMod` fit; otherwise `NULL`.

- message:

  A note when F-tests or adjusted means are not estimable.

## Adjusted means

Genotype effects are fixed, so the adjusted means are best linear
unbiased estimates (BLUEs): the genotype effect plus the average
location and block effect over all observed plots, i.e. the expected
value of the genotype in an average plot of the trial. They remove the
advantage or penalty of the blocks and locations where each entry
happened to be planted, so use `adjusted_mean` (not `raw_mean`) to rank
and select entries. Differences between adjusted means are the same as
those from `emmeans` on the fitted model. With
`genotype_by_location = "fixed"` the genotype-by-location effects are
averaged over locations, weighting each location by its share of the
observed plots. With `"random"`, standard errors come from the REML
covariance of the fixed effects, which includes the genotype-by-location
variance; `lsd_5pct` uses the error df of the fixed-interaction model as
an approximation.

## Per-location analysis

Williams, Piepho and Whitaker (2011) suggest analysing each location
with an incomplete-block model. Pass the rows of one location, e.g.
`analyze_augmented_prep_design(data[data$location == "A", ], "yield")`.
[`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md)
keeps the blocks of every location connected, so this is estimable
whenever the location has enough p-rep plots.

## Examples

``` r
design <- augmented_prep_design(
  treatments = paste0("G", seq_len(8)),
  locations = c("Loc1", "Loc2"),
  block_size = 4,
  n_duplicate_per_location = c(Loc1 = 2, Loc2 = 2),
  n_candidates = 3,
  seed = 123
)

design_data <- design$design
design_data$yield <- rnorm(nrow(design_data), mean = 100, sd = 10)
analysis <- analyze_augmented_prep_design(design_data, response = "yield")
analysis
#> Augmented p-rep design analysis
#> 
#>                            Source Df     Sum Sq   Mean Sq  F value     Pr(>F)
#> 1                       Locations  1   77.50845  77.50845 2.059080 0.19443287
#> 2 Blocks, nested within locations  4  378.47777  94.61944 2.513649 0.13570851
#> 3                       Genotypes  7 1632.20594 233.17228 6.194427 0.01402908
#> 4                           Error  7  263.49586  37.64227       NA         NA
#> 5                 Corrected total 19 2351.68801        NA       NA         NA
#>   Significance
#> 1             
#> 2             
#> 3            *
#> 4             
#> 5             
#> 
#> Top adjusted means (8 genotypes; see $adjusted_means):
#>  rank genotype n_observations adjusted_mean std_error
#>     1       G1              2        120.18     4.814
#>     2       G5              2        108.08     4.964
#>     3       G2              2        105.41     5.303
#>     4       G7              3         98.24     3.911
#>     5       G8              3         95.13     4.037
#>     6       G4              2         93.29     5.147
#>     7       G6              3         93.24     3.924
#>     8       G3              3         85.17     3.888
#> Mean SED = 6.52; LSD (5%) = 15.42
head(analysis$adjusted_means)
#> # A tibble: 6 × 11
#>    rank genotype entry_type    repeated_at is_control n_observations n_locations
#>   <int> <chr>    <chr>         <chr>            <int>          <int>       <int>
#> 1     1 G1       single_entry… not_repeat…          0              2           2
#> 2     2 G5       single_entry… not_repeat…          0              2           2
#> 3     3 G2       single_entry… not_repeat…          0              2           2
#> 4     4 G7       p_rep_entry_… Loc1                 0              3           2
#> 5     5 G8       p_rep_entry_… Loc1                 0              3           2
#> 6     6 G4       single_entry… not_repeat…          0              2           2
#> # ℹ 4 more variables: n_extra_repeats <int>, raw_mean <dbl>,
#> #   adjusted_mean <dbl>, std_error <dbl>
```
