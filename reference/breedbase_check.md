# Check names in Breedbase before building a design

Run this first, before designing the trial: it checks that the
location(s) and every entry and check are already in Breedbase (and,
when given, the breeding program, the folder and that the trial name is
free). All problems are reported together, with close matches as
suggestions (for example `"ibadan"` -\> `"Ibadan"`), so they can be
fixed in one go.
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
repeats the same checks before writing.

## Usage

``` r
breedbase_check(
  con,
  accessions = NULL,
  location = NULL,
  program = NULL,
  folder = NULL,
  trial_name = NULL,
  plot_names = NULL
)
```

## Arguments

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- accessions:

  Entry and check names: a character vector, a design from
  [`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
  [`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
  [`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md)
  or
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md),
  or a list of such designs (e.g. one per location).

- location:

  Breedbase location name(s).

- program, folder:

  Optional breeding program and folder (the folder is optional for
  saving too: without it the trial goes directly in the program).

- trial_name:

  Optional trial name(s), e.g. one per location. A name already used by
  a trial (with plots) is reported as `"TAKEN"`: change it before
  saving.

- plot_names:

  Optional plot names to check: a character vector, a design, or a list
  of designs. By default the plot names stored in the design(s) given in
  `accessions` are checked (see
  [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md)).
  A design without stored plot names gets the default names of
  [`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
  for the check when it can be matched to one trial name and location
  (one design and one `trial_name`, or a list of designs and one trial
  name and location per design). Plot names must be unique in the whole
  database: a name already used is reported as `"TAKEN"`, a name given
  to two plots as `"REPEATED"`. To save with other names, change the
  `plot_name` column of the design or use
  [`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md).

## Value

Invisibly, a tibble with one row per name checked: `what`, `name`,
`status` (`"OK"`, `"MISSING"`, `"NEW"`, `"TAKEN"`, `"EMPTY"` or
`"REPEATED"`) and `suggestion`. Its attribute `"ok"` is `TRUE` when
nothing needs fixing.

## Slow servers

Accessions are searched 25 names at a time, because Breedbase looks up
the pedigree of each one. If the server times out (HTTP 504, 502, 503 or
429), the batch is split in half and tried again, so a busy server is
slower but still works. A line is printed after each batch with the
names done and the time. Accessions found are remembered until R is
restarted, so
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
after `breedbase_check()` does not search them again.

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
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
entries <- c(paste0("IITA-TMS-IBA", 1:108))
checks <- c("TMEB419", "IITA-TMS-IBA30572", "TMEB693", "IITA-TMS-IBA000070")
breedbase_check(con, accessions = c(entries, checks), location = "Ibadan",
                program = "IITA", folder = "AYT 2026", trial_name = "AYT_2026_Ibadan")
} # }
```
