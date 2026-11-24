# Create a Field Book trait file

Makes trait definitions for offline import into the Field Book Android
app. Use
[`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
with no trait columns for the companion field CSV. No connection to
Breedbase is needed to create either file.

## Usage

``` r
trait_template(
  traits,
  format = "numeric",
  default_value = "",
  minimum = NA_real_,
  maximum = NA_real_,
  details = "",
  categories = "",
  is_visible = TRUE,
  output_file = NULL
)
```

## Arguments

- traits:

  Character vector of unique trait names. Use names that Breedbase can
  resolve when uploading observations later.

- format:

  Field Book format: `"numeric"`, `"categorical"`, `"percent"`,
  `"date"`, `"boolean"`, `"text"`, `"photo"`, `"counter"`, `"location"`,
  or `"angle"`. A single value is repeated for all traits.

- default_value, minimum, maximum, details:

  Optional defaults, numeric bounds and descriptions, respectively. Each
  must have length one or the number of traits. `NA` becomes an empty
  cell. Units can be included in `details`.

- categories:

  For categorical traits, a character vector of slash-separated choices,
  or a list of character vectors, one per trait (e.g.
  `list(character(), c("low", "medium", "high"))`). A single value or
  list element is repeated. Other formats must leave this empty.

- is_visible:

  Logical; show the trait in Field Book's Collect screen. Length one or
  the number of traits.

- output_file:

  Optional path to a `.trt` file. Written as UTF-8 CSV with empty
  missing cells and no row names.

## Value

A tibble with `trait`, `format`, `defaultValue`, `minimum`, `maximum`,
`details`, `categories`, `isVisible`, and `realPosition`. The table is
also returned when `output_file` is supplied.

## Details

The columns follow Field Book's bundled [sample trait
file](https://github.com/PhenoApps/Field-Book/blob/main/app/src/main/assets/trait/trait_sample.trt).
Import the resulting `.trt` through Traits \> Import/Export in the app.
See the [Field Book trait
documentation](https://field-book.readthedocs.io/en/latest/traits.html).
These local definitions do not create ontology traits in Breedbase;
those must already exist before observations can be saved there.

## See also

[`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md),
[`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md)

## Examples

``` r
traits <- trait_template(
  c("fresh root yield", "notes"), format = c("numeric", "text"),
  minimum = c(0, NA), details = c("Yield in t/ha", "Field observations")
)
trait_template("severity", format = "categorical",
               categories = list(c("low", "medium", "high")))
#> # A tibble: 1 × 9
#>   trait    format      defaultValue minimum maximum details categories isVisible
#>   <chr>    <chr>       <chr>        <chr>   <chr>   <chr>   <chr>      <chr>    
#> 1 severity categorical ""           ""      ""      ""      low/mediu… true     
#> # ℹ 1 more variable: realPosition <int>
path <- tempfile(fileext = ".trt")
trait_template("fresh root yield", minimum = 0, output_file = path)
#> # A tibble: 1 × 9
#>   trait         format defaultValue minimum maximum details categories isVisible
#>   <chr>         <chr>  <chr>        <chr>   <chr>   <chr>   <chr>      <chr>    
#> 1 fresh root y… numer… ""           0       ""      ""      ""         true     
#> # ℹ 1 more variable: realPosition <int>
```
