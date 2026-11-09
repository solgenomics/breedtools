# Save first-stage results to Breedbase as a trait of the field plots

Saves plot-level first-stage results (from
[`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md))
on the plots of the field trial, as values of a trait made for adjusted
values, so raw and adjusted values stay apart. By default that trait is
the measured trait composed with a statistic, e.g. "fresh root
yield\|Adjusted Means from Spatial Correction using SpATS R": it is
found, or created (after you confirm) when it does not exist yet; see
[`breedbase_adjusted_trait()`](https://solgenomics.github.io/breedtools/reference/breedbase_adjusted_trait.md).
Or name your own traits with `as_trait`. The values are saved through
BrAPI, like
[`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md).

## Usage

``` r
save_breedbase_first_stage(
  values,
  con = NULL,
  as_trait = NULL,
  statistic = NULL,
  create_trait = TRUE,
  dry_run = FALSE,
  confirm = interactive(),
  overwrite = FALSE,
  remove_empty = FALSE
)
```

## Arguments

- values:

  A table from
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)
  (or any table with a `plot_name` column and one numeric column per
  trait).

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md);
  `NULL` with `dry_run = TRUE` returns what would be sent.

- as_trait:

  Optional named vector: for each trait column of `values`, an existing
  Breedbase trait to save it under, e.g.
  `c("fresh root yield" = "fresh root yield adjusted|COMP:0000123")`.

- statistic:

  The statistic composed with each trait when `as_trait` is not given.
  By default chosen from this Breedbase's own statistics (the spatial
  correction term); see
  [`breedbase_statistics()`](https://solgenomics.github.io/breedtools/reference/breedbase_statistics.md).

- create_trait:

  Logical. Create the composed traits that do not exist yet (you are
  asked first when `confirm = TRUE`).

- dry_run, confirm, overwrite, remove_empty:

  As in
  [`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md):
  run it again after a new analysis with `overwrite = TRUE` to replace
  the stored values.

## Value

As
[`save_breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_phenotypes.md):
one row per trait with the values sent and found.

## Slow servers

Accessions are searched 25 names at a time, because Breedbase looks up
the pedigree of each one. If the server times out (HTTP 504, 502, 503 or
429), the batch is split in half and tried again, so a busy server is
slower but still works. A line is printed after each batch with the
names done and the time. Accessions found are remembered until R is
restarted, so
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
after
[`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
does not search them again.

When creating the folder, the trial or its plots times out, the save
does not stop: the server may still finish the job, so it checks for the
result every 20 seconds for 3 minutes, and tries once more if it is not
there (Breedbase refuses names that already exist, so this cannot make
duplicates). If a save is cut off anyway, run
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
again: an existing folder is used, and a trial left without plots gets
its plots. Settings:
`options(breedtools.write_wait = 30, breedtools.write_checks = 10)`.
Settings: `options(breedtools.germplasm_batch = 10)` for smaller
batches, `options(breedtools.breedbase_timeout = 600)` to wait longer
(seconds) for each answer.

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
d <- rc$design
set.seed(1)
d$yield <- 25 + 0.4 * d$row + rnorm(nrow(d), sd = 2)
fit <- analyze_augmented_row_column_design(d, "yield")
values <- first_stage_values(fit, rc, "AYT_2026_Ibadan", trait = "fresh root yield",
                            location = "Ibadan")
obs <- save_breedbase_first_stage(values, dry_run = TRUE)
#> Dry run without a connection: returning the observations that would be sent.
head(obs)
#> # A tibble: 6 × 4
#>   plot_name                      trait                           value timestamp
#>   <chr>                          <chr>                           <chr> <chr>    
#> 1 AYT_2026_Ibadan_Ibadan_G51_1_1 fresh root yield|Adjusted Mean… 28.7… NA       
#> 2 AYT_2026_Ibadan_Ibadan_G64_1_1 fresh root yield|Adjusted Mean… 29.9… NA       
#> 3 AYT_2026_Ibadan_Ibadan_G97_1_1 fresh root yield|Adjusted Mean… 29.8… NA       
#> 4 AYT_2026_Ibadan_Ibadan_C3_1_1  fresh root yield|Adjusted Mean… 27.3… NA       
#> 5 AYT_2026_Ibadan_Ibadan_G8_2_1  fresh root yield|Adjusted Mean… 28.2… NA       
#> 6 AYT_2026_Ibadan_Ibadan_C2_2_2  fresh root yield|Adjusted Mean… 28.0… NA       

if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
save_breedbase_first_stage(values, con, dry_run = TRUE)   # check; shows the trait to create
save_breedbase_first_stage(values, con)                   # creates the trait if needed, saves
} # }
```
