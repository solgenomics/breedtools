# Read phenotypic data of trials from Breedbase

Reads, through BrAPI, the plots of the given trials (plot name and
number, accession, rep, block, row, column, check) and their values of
the given traits, and returns one table with one row per plot and one
column per trait, ready for
[`plot_trait()`](https://solgenomics.github.io/breedtools/reference/plot_trait.md),
[`remove_outliers()`](https://solgenomics.github.io/breedtools/reference/remove_outliers.md),
[`analyze_two_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_two_stage.md)
and
[`analyze_single_stage()`](https://solgenomics.github.io/breedtools/reference/analyze_single_stage.md).
Repeated measurements of a trait on a plot are averaged (by Breedbase).

## Usage

``` r
breedbase_phenotypes(con, trials, traits)
```

## Arguments

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- trials:

  Trial names, as in Breedbase.

- traits:

  Trait names, as in Breedbase, with or without the ontology id (e.g.
  `"fresh root yield"` or `"fresh root yield|CO_334:0000013"`).

## Value

A tibble: `trial_name`, `location`, `year`, `plot_name`, `plots` (plot
number), `accession`, `rep`, `block`, `row`, `col`, `is_check` and one
numeric column per trait (named as given in `traits`). The trait names
in Breedbase are in `attr(result, "trait_names")`.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassava-test.sgn.cornell.edu", "jdoe", site_username = "web_user_name")
pheno <- breedbase_phenotypes(con, trials = c("FALL_2026_Ibadan", "FALL_2026_Ikenne"),
                              traits = c("fresh root yield", "dry matter content percentage"))
plot_trait(pheno, "fresh root yield", by = "trial_name")
} # }
```
