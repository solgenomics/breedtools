# Create a map of environment groups

Creates an interactive `leaflet` map where locations belonging to the
same environment group use the same marker color. Ungrouped environments
are shown with `ungrouped_color`.

## Usage

``` r
map_environment_groups(
  stratification_result,
  locations,
  ignore_case = FALSE,
  group_palette = NULL,
  ungrouped_color = "#7A7A7A",
  match_by = "location",
  strict = TRUE,
  tiles = "open"
)
```

## Arguments

- stratification_result:

  Result from
  [`run_environment_stratification()`](https://solgenomics.github.io/breedtools/reference/run_environment_stratification.md)
  or
  [`run_breedbase_environment_stratification()`](https://solgenomics.github.io/breedtools/reference/run_breedbase_environment_stratification.md).

- locations:

  Location coordinates from
  [`read_locations()`](https://solgenomics.github.io/breedtools/reference/read_locations.md),
  a location file path, or a data frame with location, latitude, and
  longitude columns.

- ignore_case:

  Logical. If `TRUE`, location matching ignores case.

- group_palette:

  Optional vector of colors used for environment groups.

- ungrouped_color:

  Color for ungrouped environments.

- match_by:

  Name from the stratification result used to match the `location`
  column in `locations`. Use `"location"` when coordinates are named by
  location, `"environment"` when coordinates are named by the
  constructed environment identifier, or `"environment_label"` when
  coordinates are named by display labels.

- strict:

  Logical. If `TRUE`, stop when any environment location has no matching
  coordinates.

- tiles:

  Background map. All maps are open: free to use, with no account, API
  key or commercial service. The default, `"open"`, offers three maps
  with a switch in the top-right corner: Street (OpenStreetMap, shown
  first), Topographical (OpenTopoMap) and Satellite (Sentinel-2
  cloudless 2016 by EOX, from Copernicus Sentinel data). `"street"`,
  `"topo"` or `"satellite"` shows only that map. They work in the
  RStudio viewer and in saved HTML files; an internet connection is
  needed to load them. Any other value is passed to
  [`leaflet::addProviderTiles()`](https://rstudio.github.io/leaflet/reference/addProviderTiles.html)
  as a provider name; check that provider's terms, as some need an API
  key.

## Value

A `leaflet` htmlwidget.

## Examples

``` r
pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
loc_file <- system.file("extdata", "breedbase_locations.csv", package = "breedtools")

result <- run_environment_stratification(pheno_file, trait = "fresh root yield")
locations <- read_locations(loc_file, phenotype_data = result)

group_map <- map_environment_groups(result, locations)
if (interactive()) group_map
```
