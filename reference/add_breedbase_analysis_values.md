# Add results to an analysis that is already in Breedbase

Adds values (for example a second trait, or the same trait with another
statistic) to the plots of an analysis saved before, e.g. with
[`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md).
Breedbase names the plots of an analysis
`"<analysis_name>_<field plot name>"` ("Analysis of: plots") or
`"<analysis_name>_<accession>"` ("Analysis of: accessions"); the values
are uploaded to those plots with the phenotype upload of the website
(simple spreadsheet), which first verifies the file and then stores it.

## Usage

``` r
add_breedbase_analysis_values(
  data,
  con,
  analysis_name,
  level = c("plot", "accession"),
  traits = NULL,
  plot = "plot_name",
  accession = "accession",
  statistic = NULL,
  estimate = NULL,
  overwrite = FALSE,
  dry_run = FALSE,
  confirm = interactive()
)
```

## Arguments

- data:

  A data frame with the field plot names (`level = "plot"`) or the
  accession names (`level = "accession"`) and one numeric column per
  trait, e.g. from
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)
  or
  [`second_stage_values()`](https://solgenomics.github.io/breedtools/reference/second_stage_values.md).

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)
  (the site's web pages must be open; see
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)).

- analysis_name:

  Name of the existing analysis.

- level:

  `"plot"` or `"accession"`, as the analysis was saved.

- traits:

  Trait columns (default: all except the name and design columns). Their
  names are Breedbase trait names, including composed traits such as
  `"dry matter content percentage|Adjusted Means from Spatial Correction using SpATS R"`
  (see
  [`breedbase_adjusted_trait()`](https://solgenomics.github.io/breedtools/reference/breedbase_adjusted_trait.md)).

- plot, accession:

  Names of the plot-name and accession columns.

- statistic:

  Optional statistic (see
  [`breedbase_statistics()`](https://solgenomics.github.io/breedtools/reference/breedbase_statistics.md)).
  Values are always stored under "trait\|statistic" (a composed trait,
  created if needed), never under the plain trait, as Breedbase does for
  analyses. By default the statistic the analysis already uses (from its
  stored traits); give another to keep, e.g., single-stage means next to
  two-stage means in the same analysis.

- estimate:

  When the analysis has no statistic yet: `"blue"`, `"blup"` or
  `"spatial"`, as in
  [`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)
  (by default from the data).

- overwrite:

  Logical. Replace values already stored.

- dry_run:

  Logical. Only verify (Breedbase checks the file), store nothing.

- confirm:

  Logical. Ask before storing (default in interactive sessions).

## Value

Invisibly, a list with the analysis name, the number of plots and values
sent, and Breedbase's messages.

## Details

Values already present for a plot and trait are not replaced unless
`overwrite = TRUE`: the verification lists them and the upload stops.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassava-test.sgn.cornell.edu", "jdoe", site_username = "web_user_name")
# dm1: first-stage values of dry matter, one row per field plot
add_breedbase_analysis_values(dm1, con, analysis_name = "FALL_2026_stage1", level = "plot")
} # }
```
