# Plots of an analysis saved in Breedbase

Reads back the plots of an existing analysis (from
[`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
or made on the website): the analysis plot names, the field plot (or
accession) each one stands for, its accession and design, and the values
of every trait stored in the analysis. Use it to see what an analysis
holds before adding values with
[`add_breedbase_analysis_values()`](https://solgenomics.github.io/breedtools/reference/add_breedbase_analysis_values.md),
e.g. a second trait.

## Usage

``` r
breedbase_analysis_plots(con, analysis_name)
```

## Arguments

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- analysis_name:

  Name of the analysis.

## Value

A tibble with `analysis_plot_name`, `plot_name` (the field plot, or the
accession for accession analyses: the name to use in the `plot` or
`accession` column of
[`add_breedbase_analysis_values()`](https://solgenomics.github.io/breedtools/reference/add_breedbase_analysis_values.md)),
`accession`, `plot_number`, `rep`, `block`, `row`, `col` and one column
per trait, ordered by plot number. The analysis id is in the attribute
`"analysis_db_id"`.

## Details

Analysis plots are not part of BrAPI, so they are read through the web
pages of Breedbase (the site's web password is needed when it has one;
see
[`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)).

Breedbase lists the plots of an analysis through their values: a plot
without any value (e.g. a lost plot, or an outlier removed before the
analysis) is not listed, although it is part of the analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassava-test.sgn.cornell.edu", "jdoe", site_username = "sgn")
plots <- breedbase_analysis_plots(con, "AYT_2026_stage1")
head(plots)
} # }
```
