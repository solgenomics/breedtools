# Run individual ANOVA by environment

Run individual ANOVA by environment

## Usage

``` r
analyze_environment_anova(pheno)
```

## Arguments

- pheno:

  Standardized phenotype data from
  [`prepare_breedbase_pheno()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_pheno.md).

## Value

A tibble with one ANOVA summary row per environment.

## Examples

``` r
df <- data.frame(
  environment = factor(rep(c("Loc1", "Loc2"), each = 8)),
  genotype = factor(rep(rep(c("G1", "G2"), each = 4), times = 2)),
  replicate = factor(rep(rep(c(1, 2), each = 2), times = 4)),
  block = factor(1),
  value = c(10, 11, 12, 13, 14, 15, 16, 17,
            15, 16, 17, 18, 19, 20, 21, 22)
)

analyze_environment_anova(df)
#> # A tibble: 2 × 12
#>   environment n_genotypes r_eff df_genotype ms_genotype df_error ms_error  mean
#>   <chr>             <int> <dbl>       <int>       <dbl>    <int>    <dbl> <dbl>
#> 1 Loc1                  2     4           1        32.0        5    0.400  13.5
#> 2 Loc2                  2     4           1        32          5    0.400  18.5
#> # ℹ 4 more variables: cv_percent <dbl>, sigma2_g <dbl>, f_value <dbl>,
#> #   p_value <dbl>
```
