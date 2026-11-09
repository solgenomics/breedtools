# Save plot GPS coordinates of a trial already in Breedbase

Sends the plot coordinates of a design (see
[`add_plot_coordinates()`](https://solgenomics.github.io/breedtools/reference/add_plot_coordinates.md))
to a trial that is already saved in Breedbase, e.g. when the plots were
measured after the trial was saved. Plots are matched by plot number,
after checking that every plot has the same accession in Breedbase as in
the design. Coordinates already stored for these plots are replaced.
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
does this by itself for a design that has coordinates.

## Usage

``` r
save_breedbase_coordinates(
  design,
  con,
  trial_name = NULL,
  dry_run = FALSE,
  confirm = interactive()
)
```

## Arguments

- design:

  A design with coordinates, from
  [`add_plot_coordinates()`](https://solgenomics.github.io/breedtools/reference/add_plot_coordinates.md).

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- trial_name:

  The trial name in Breedbase. For a multi-location p-rep design, the
  location is added as in
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md).
  Can be omitted for a prepared Breedbase trial or a design with plot
  names.

- dry_run:

  Logical. Check the trial and plots but write nothing.

- confirm:

  Logical. Ask before writing (default in interactive sessions).

## Value

Invisibly, a tibble with one row per trial: `trial_name`,
`plots_with_gps` (sent) and `gps_found` (read back from Breedbase).

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
rc <- add_plot_coordinates(rc, origin = c(7.4931, 3.9035), plot_width = 3, plot_length = 5)
save_breedbase_coordinates(rc, con, trial_name = "AYT_2026_Ibadan")
} # }
```
