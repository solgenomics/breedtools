# Detect field trial design for environment stratification

Detects whether the data look like CRD, RCBD, row-column, alpha-lattice,
or augmented/block-only design.

## Usage

``` r
detect_design(data)
```

## Arguments

- data:

  Standardized phenotype data.

## Value

A list with `design`, `formula`, and `message`.

## Examples

``` r
pheno <- read_breedbase_file(
  system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
)
prepared <- prepare_environment_stratification_data(pheno, trait = "fresh root yield")

design <- detect_design(prepared$data)
design$design
#> [1] "RCBD"
design$formula
#> phenotype ~ environment + environment:block_number + accession_name + 
#>     environment:accession_name
#> <environment: 0x623d45e7ebc8>
```
