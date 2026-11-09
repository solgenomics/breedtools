# List the statistic terms of Breedbase (SGNSTAT ontology)

Breedbase saves analysis results as a trait combined with a statistic
term of its statistics ontology (SGNSTAT), e.g. "Adjusted Means from
Spatial Correction using SpATS R". Breedbase loads these terms from the
file `ontology/cxgn_statistics.obo` of its source code, which breedtools
includes
(`system.file("extdata", "cxgn_statistics.obo", package = "breedtools")`),
so the terms can be listed and checked without a server. With a
connection, the terms of that Breedbase are listed and compared with the
ontology file. When Breedbase adds statistics, the bundled copy is
refreshed with the script `data-raw/update_statistics_ontology.R` of the
package source.

## Usage

``` r
breedbase_statistics(con = NULL, pattern = "")
```

## Arguments

- con:

  Optional connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).
  Without it, the terms of the ontology file are listed.

- pattern:

  Optional text the term must contain (any case).

## Value

A tibble with `statistic` (the name), `accession` (e.g.
`"SGNSTAT:0000054"`) and `definition` (from the ontology file); with a
connection also `db_id` (the server's id) and `in_ontology` (`TRUE` when
name and accession match the ontology file).

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
breedbase_statistics(con)
breedbase_statistics(con, "BLUE")
} # }
breedbase_statistics(pattern = "spatial")      # from the ontology file
#> # A tibble: 5 × 3
#>   statistic                                                 accession definition
#>   <chr>                                                     <chr>     <chr>     
#> 1 Adjusted Means from Spatial Correction using SpATS R      SGNSTAT:… Values of…
#> 2 Multivariate linear mixed model 2D spline spatial BLUPs … SGNSTAT:… NA        
#> 3 Multivariate linear mixed model genetic BLUPs using gene… SGNSTAT:… NA        
#> 4 Multivariate linear mixed model residual using genetic r… SGNSTAT:… NA        
#> 5 Univariate linear mixed model spatial BLUPs using geneti… SGNSTAT:… NA        
```
