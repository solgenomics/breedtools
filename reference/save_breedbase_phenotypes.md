# Save phenotypic data to Breedbase

Uploads trait values of plots that are already in Breedbase (for example
saved with
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md))
through BrAPI. Give the plot name and the trait name with each value;
the plots and traits are looked up in Breedbase first, and every problem
is listed before anything is written.

## Usage

``` r
save_breedbase_phenotypes(
  data,
  con = NULL,
  plot = "plot_name",
  traits = NULL,
  trait = NULL,
  value = NULL,
  timestamp = NULL,
  overwrite = FALSE,
  remove_empty = FALSE,
  dry_run = FALSE,
  confirm = interactive(),
  refresh = TRUE
)
```

## Arguments

- data:

  A data frame with the plot names and the values.

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md);
  `NULL` with `dry_run = TRUE` returns the observations without checking
  them.

- plot:

  Name of the column with the plot names.

- traits:

  For wide data: the trait columns. By default every column except
  `plot`, `timestamp` and the design columns of
  [`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)
  (`trial_name`, `accession`, `plot_number`, `block`, `rep`, `row`,
  `col`, `is_check`).

- trait, value:

  For long data: the names of the trait and value columns.

- timestamp:

  Optional name of a column with the date and time of each measurement
  (e.g. `"2026-08-15 10:30:00"`). Breedbase uses it to tell repeated
  measurements apart.

- overwrite:

  Logical. Replace values that Breedbase already has with different ones
  (see "Saving again").

- remove_empty:

  Logical. Delete stored values whose cell is empty (`NA` or `""`) in
  `data`; needs `overwrite = TRUE` (see "Saving again").

- dry_run:

  Logical. Check everything but write nothing.

- confirm:

  Logical. Ask before writing (default in interactive sessions).

- refresh:

  Logical. After saving, refresh the phenotype cache of Breedbase so the
  values show on the trial pages (see
  [`breedbase_refresh()`](https://solgenomics.github.io/breedtools/reference/breedbase_refresh.md)).

## Value

Invisibly, a tibble with one row per trait: `trait`,
`observation_variable_db_id`, `values_sent`, `values_found` (read back
from Breedbase after saving), and how the values were handled: `new`,
`replaced`, `unchanged` and `removed`. Returns invisible `NULL` when
saving is declined, including changed values with `overwrite = FALSE`.

## Details

The data can be wide, one row per plot and one column per trait (as made
by
[`phenotype_template()`](https://solgenomics.github.io/breedtools/reference/phenotype_template.md)):

|                         |                  |                               |
|-------------------------|------------------|-------------------------------|
| plot_name               | fresh root yield | dry matter content percentage |
| AYT_2026_Ibadan_G31_1_1 | 25.4             | 31.2                          |

or long, one row per value, with `trait` and `value` columns
(`trait = "trait", value = "value"`).

Trait names are the Breedbase names, with or without the ontology id:
`"fresh root yield"` or `"fresh root yield|CO_334:0000013"` (the column
names of a Breedbase phenotype download). Values are checked by
Breedbase against the trait (e.g. numbers within the allowed range).

## Saving again (values already in Breedbase)

Before writing, the values already stored for these plots and traits are
read, and every value is sorted into new, unchanged (same value, not
sent again) or changed:

- `overwrite = FALSE` (default):

  If changed values exist, a message is printed and the function returns
  without an error; nothing is written.

- `overwrite = TRUE`:

  Changed values replace the stored ones.

Empty cells (`NA` or `""`) are never sent. When Breedbase has a value
for an empty cell (e.g. an outlier set to `NA` with
[`remove_outliers()`](https://solgenomics.github.io/breedtools/reference/remove_outliers.md)):

- `remove_empty = FALSE` (default):

  The stored value is kept (and reported).

- `remove_empty = TRUE`:

  The stored value is deleted. Like the "remove values" option of
  Breedbase's upload page, this needs `overwrite = TRUE`, a curator
  account and access to the web pages (the site password, see
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)).

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
sheet <- phenotype_template(rc, "AYT_2026", location = "Ibadan", traits = "fresh root yield")
sheet$`fresh root yield` <- round(rnorm(nrow(sheet), 25, 4), 1)

# What would be sent (no connection needed):
obs <- save_breedbase_phenotypes(sheet, dry_run = TRUE)
#> Dry run without a connection: returning the observations that would be sent.
head(obs)
#> # A tibble: 6 × 4
#>   plot_name               trait            value timestamp
#>   <chr>                   <chr>            <chr> <chr>    
#> 1 AYT_2026_Ibadan_G51_1_1 fresh root yield 21.9  NA       
#> 2 AYT_2026_Ibadan_G64_1_1 fresh root yield 17.4  NA       
#> 3 AYT_2026_Ibadan_G97_1_1 fresh root yield 33.6  NA       
#> 4 AYT_2026_Ibadan_C3_1_1  fresh root yield 24.9  NA       
#> 5 AYT_2026_Ibadan_G8_2_1  fresh root yield 23.8  NA       
#> 6 AYT_2026_Ibadan_C2_2_2  fresh root yield 22.4  NA       

if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
save_breedbase_phenotypes(sheet, con, dry_run = TRUE)   # check first
save_breedbase_phenotypes(sheet, con)

# Saving a corrected sheet again: replace changed values, and delete the
# stored values of cells now empty (e.g. outliers set to NA).
save_breedbase_phenotypes(sheet, con, overwrite = TRUE, remove_empty = TRUE)

# Long data: one row per value
long <- data.frame(plot_name = sheet$plot_name, trait = "fresh root yield",
                   value = sheet$`fresh root yield`)
save_breedbase_phenotypes(long, con, trait = "trait", value = "value")
} # }
```
