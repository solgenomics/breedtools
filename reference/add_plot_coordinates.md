# Add GPS coordinates to the plots of a design

Adds the GPS position of every plot to `design$design`, so that
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
saves it with the plots and Breedbase shows the trial on its field map.
Coordinates are WGS84 decimal degrees (as from a phone or GPS). There
are three ways to give them:

## Usage

``` r
add_plot_coordinates(
  design,
  coordinates = NULL,
  latitude = NULL,
  longitude = NULL,
  origin = NULL,
  plot_width = NULL,
  plot_length = NULL,
  col_gap = 0,
  row_gap = 0,
  col_bearing = 90,
  row_bearing = 0
)
```

## Arguments

- design:

  A design from
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md),
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md)
  or
  [`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md).

- coordinates:

  Optional table of coordinates, as described above.

- latitude, longitude:

  Optional point per plot, as described above.

- origin:

  Optional `c(latitude, longitude)` of the outer corner of the plot in
  row 1, column 1. For a p-rep design with several locations, a list
  with one such pair per location, named by location.

- plot_width, plot_length:

  Plot size in metres: `plot_width` across columns (from column 1
  towards column 2), `plot_length` across rows (from row 1 towards row
  2). Used with `origin`.

- col_gap, row_gap:

  Space in metres between plots of neighbouring columns and rows
  (alleys). Used with `origin`.

- col_bearing, row_bearing:

  Compass direction (degrees clockwise from north) in which column
  numbers and row numbers increase; they must be at right angles. The
  default puts column 1 at the west, columns running east, and row 1 at
  the south, rows running north. Used with `origin`.

## Value

The design with the coordinate columns added.

## Details

- **From the field corner** (`origin`): give the GPS of the outer corner
  of the plot in row 1, column 1, the plot size and the directions in
  which rows and columns run. The four corners of every plot are
  calculated from its row and column.

- **One point per plot** (`latitude`, `longitude`): vectors with one
  value per row of `design$design`, or the names of columns of the
  design that hold them.

- **A table** (`coordinates`), e.g. read with
  [`read_breedbase_file()`](https://solgenomics.github.io/breedtools/reference/read_breedbase_file.md):
  one row per plot, matched by `plot_name` (when both have it) or by
  plot number (`plots`, `plot_number` or `plot`; for p-rep designs also
  `location`). It holds either `latitude` and `longitude` (or `WGS84_x`,
  `WGS84_y`) for a point, or the eight corner columns of the Breedbase
  GPS upload template (`WGS84_bottom_left_x`, `WGS84_bottom_left_y`,
  ..., `WGS84_top_left_y`, where `x` is the longitude and `y` the
  latitude) for a polygon.

With `origin`, positions are calculated on a local flat approximation of
the Earth, accurate to centimetres over a field. The design needs `row`
and `col` columns (designs converted from agricolae need them in the
field book).

The design gets `latitude` and `longitude` columns (the centre of the
plot), plus the eight corner columns for polygons. They can be edited.
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
sends a polygon when all eight corner columns of a plot are filled, and
otherwise a point from `latitude` and `longitude`. Plots without
coordinates are saved without them. For a trial that is already in
Breedbase, use
[`save_breedbase_coordinates()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_coordinates.md).

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
# Plots of 3 m x 5 m, 1 m alleys between rows, columns running east:
rc <- add_plot_coordinates(rc, origin = c(7.4931, 3.9035),
                           plot_width = 3, plot_length = 5, row_gap = 1)
head(rc$design[, c("plots", "row", "col", "latitude", "longitude")])
#>   plots row col latitude longitude
#> 1     1   1   1 7.493123  3.903514
#> 2     2   1   2 7.493123  3.903541
#> 3     3   1   3 7.493123  3.903568
#> 4     4   1   4 7.493123  3.903595
#> 5     5   1   5 7.493123  3.903622
#> 6     6   1   6 7.493123  3.903649

# A field laid out at an angle: columns running north-east.
rc <- add_plot_coordinates(rc, origin = c(7.4931, 3.9035), plot_width = 3,
                           plot_length = 5, col_bearing = 45, row_bearing = 315)

# One GPS point per plot, measured in the field:
pts <- data.frame(plots = rc$design$plots,
                  latitude = 7.4931 + rc$design$row * 5e-5,
                  longitude = 3.9035 + rc$design$col * 3e-5)
rc <- add_plot_coordinates(rc, coordinates = pts)
```
