# Plot and accession names an analysis will have, before saving

Shows, without contacting Breedbase, the names under which
[`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
(or
[`add_breedbase_analysis_values()`](https://solgenomics.github.io/breedtools/reference/add_breedbase_analysis_values.md))
stores each row of `data`: Breedbase gives every analysis its own plots,
named `"<analysis_name>_<field plot name>"` ("Analysis of: plots", with
the accession of the field plot) or `"<analysis_name>_<accession>"`
("Analysis of: accessions"). Use it to see what will be saved, and to
compare with the plots of a saved analysis
([`breedbase_analysis_plots()`](https://solgenomics.github.io/breedtools/reference/breedbase_analysis_plots.md)).

## Usage

``` r
analysis_plot_names(
  data,
  analysis_name,
  level = c("plot", "accession"),
  plot = "plot_name",
  accession = "accession",
  traits = NULL
)
```

## Arguments

- data:

  A data frame with one row per field plot (`level = "plot"`, e.g. from
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md))
  or per accession (`level = "accession"`, e.g. from
  [`second_stage_values()`](https://solgenomics.github.io/breedtools/reference/second_stage_values.md))
  and one column per trait.

- analysis_name:

  The analysis name, as given to
  [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  (special characters are replaced with `"_"` the same way).

- level:

  `"plot"` or `"accession"`.

- plot, accession:

  Names of the plot-name and accession columns.

- traits:

  Trait columns (default: all except the name and design columns).

## Value

A tibble with `analysis_plot_name` (the name in the analysis),
`plot_name` (the field plot; for accession analyses, the accession),
`accession`, and the trait values.

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
d <- rc$design
set.seed(1)
d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
fit <- analyze_augmented_row_column_design(d, "yield")
values <- first_stage_values(fit, rc, "AYT_2026", trait = "fresh root yield", location = "Ibadan")
head(analysis_plot_names(values, "AYT_2026_stage1", level = "plot"))
#> # A tibble: 6 × 4
#>   analysis_plot_name                      plot_name accession `fresh root yield`
#>   <chr>                                   <chr>     <chr>                  <dbl>
#> 1 AYT_2026_stage1_AYT_2026_Ibadan_G31_1_1 AYT_2026… G31                     25.0
#> 2 AYT_2026_stage1_AYT_2026_Ibadan_G70_1_1 AYT_2026… G70                     23.9
#> 3 AYT_2026_stage1_AYT_2026_Ibadan_G104_1… AYT_2026… G104                    24.1
#> 4 AYT_2026_stage1_AYT_2026_Ibadan_C4_1_1  AYT_2026… C4                      29.5
#> 5 AYT_2026_stage1_AYT_2026_Ibadan_G66_2_1 AYT_2026… G66                     23.4
#> 6 AYT_2026_stage1_AYT_2026_Ibadan_G5_2_1  AYT_2026… G5                      22.9
```
