# Save a design to Breedbase

Creates the trial and its plots in a Breedbase database through BrAPI,
in a breeding program and, optionally, a folder of your choice. Works
for designs from
[`augmented_row_column_design()`](https://solgenomics.github.io/breedtools/reference/augmented_row_column_design.md),
[`augmented_prep_design()`](https://solgenomics.github.io/breedtools/reference/augmented_prep_design.md),
[`resolvable_row_column_design()`](https://solgenomics.github.io/breedtools/reference/resolvable_row_column_design.md),
and agricolae designs converted with
[`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md)
or prepared with
[`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md).

## Usage

``` r
save_breedbase(
  design,
  con = NULL,
  program = NULL,
  folder = NULL,
  trial_name,
  year,
  location = NULL,
  description,
  trial_type = "phenotyping_trial",
  design_type = NULL,
  plot_names = NULL,
  plot_name_pattern = NULL,
  create_folder = FALSE,
  dry_run = FALSE,
  confirm = interactive()
)
```

## Arguments

- design:

  A design object, including the result of
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md)
  or
  [`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md).
  For a prepared trial, recorded trial name, year, location, program,
  description, folder and trial type are used when their arguments are
  omitted. Explicit arguments override them.

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md).

- program:

  Breeding program name.

- folder:

  Optional folder inside the program. Without a folder (the default,
  `NULL`) the trial is saved directly in the breeding program. An
  existing folder is used as it is; a new one is created when
  `create_folder = TRUE`.

- trial_name:

  Name of the new trial (must not exist yet). In the names saved in
  Breedbase (trial, plot and new folder names, and the location part of
  p-rep trial names) spaces and special characters such as
  `% $ ( ) { } [ ] \\ / < > *` are replaced with `"_"` (only letters,
  digits, `_`, `-` and `.` are kept), e.g. `"Alcobaca (BA)"` gives
  `"Alcobaca_BA"`. Accession and location names are kept as they are, to
  match Breedbase.

  \[ \]: R:%20

- year:

  Year of the trial.

- location:

  Breedbase location name. For p-rep designs, a named vector mapping the
  design's locations to Breedbase locations, e.g.
  `c(North = "Ibadan", South = "Ubiaja")`; by default the design's
  location names are used.

- description:

  Description of the trial (required), e.g. its purpose, material and
  season:
  `"Advanced yield trial of 108 lines and 4 checks, Ibadan 2026"`. It is
  shown on the trial page in Breedbase.

- trial_type:

  Breedbase trial type, e.g. `"Preliminary Yield Trial"`. Types the
  database does not know are saved as a miscellaneous trial.

- design_type:

  Breedbase design type. Defaults to `"Augmented"` for augmented
  row-column designs, `"p-rep"` for p-rep designs and `"Alpha"` (rows as
  incomplete blocks within replicates) for resolvable row-column
  designs. Converted agricolae designs use the type selected by
  [`as_breedbase_design()`](https://solgenomics.github.io/breedtools/reference/as_breedbase_design.md).

- plot_names:

  Optional plot names given by you: either a character vector with one
  name per plot, in the row order of `design$design`, or the name of a
  column of `design$design` that holds them. Names must be unique.

- plot_name_pattern:

  Pattern for automatic plot names, used when `plot_names` is `NULL`.
  With `NULL` (the default) the `plot_name` column of the design is used
  when it has one, and otherwise the pattern
  `"{trial}_{location}_{accession}_{block}_{rep}"`, which gives names
  like `"AYT_2026_Ibadan_LINE_007_3_1"`. Available fields: `{trial}`
  (the Breedbase trial name), `{location}` (the Breedbase location, with
  spaces and special characters replaced by `"_"`), `{accession}`,
  `{block}`, `{rep}`, `{plot}` (plot number), `{row}` and `{col}`. If a
  pattern gives the same name to two plots (the same accession twice in
  one block and rep), the plot number is added to those names.

- create_folder:

  Logical. Create `folder` if it does not exist.

- dry_run:

  Logical. Make all checks but write nothing.

- confirm:

  Logical. Ask before writing (default in interactive sessions).

## Value

Invisibly, a tibble with one row per Breedbase trial: location, trial
name, trial id, number of plots sent and found after saving, number of
plots with GPS coordinates sent and found, and the address of the trial
page. With `dry_run = TRUE` and no connection, a list with the trial and
plot data that would be sent.

## Details

Before anything is written, every name is checked: the breeding program,
the folder if given (created when `create_folder = TRUE`), the
location(s), that every entry and check already exists as an accession,
and that the trial name is free. In an interactive session you are then
asked to confirm. After saving, the plots are read back from Breedbase
to verify the upload.

Each plot is saved with its name, plot number, row and column, rep,
block and whether it is a check, and its GPS coordinates when the design
has them (see
[`add_plot_coordinates()`](https://solgenomics.github.io/breedtools/reference/add_plot_coordinates.md);
they are sent just after the plots and read back to check). Plot names
stored in the design (a `plot_name` column, made by the design functions
when given a `trial_name` or by
[`add_plot_names()`](https://solgenomics.github.io/breedtools/reference/add_plot_names.md))
are used as they are, so they can be reviewed and changed before saving.
Otherwise they are made automatically from the trial name, accession,
block and rep (see `plot_name_pattern`), or you can give your own with
`plot_names`. Plot names must be unique in the whole database, which is
checked before saving. A multi-location p-rep design becomes one
Breedbase trial per location, named `"<trial_name>_<location>"`. Filler
plots of resolvable row-column designs are not saved.

Run first with `dry_run = TRUE`: all the checks are made, but nothing is
written. Without a connection, `dry_run = TRUE` returns the data that
would be sent, without contacting any server.

## Slow servers

Accessions are searched 25 names at a time, because Breedbase looks up
the pedigree of each one. If the server times out (HTTP 504, 502, 503 or
429), the batch is split in half and tried again, so a busy server is
slower but still works. A line is printed after each batch with the
names done and the time. Accessions found are remembered until R is
restarted, so `save_breedbase()` after
[`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md)
does not search them again.

When creating the folder, the trial or its plots times out, the save
does not stop: the server may still finish the job, so it checks for the
result every 20 seconds for 3 minutes, and tries once more if it is not
there (Breedbase refuses names that already exist, so this cannot make
duplicates). If a save is cut off anyway, run `save_breedbase()` again:
an existing folder is used, and a trial left without plots gets its
plots. Settings:
`options(breedtools.write_wait = 30, breedtools.write_checks = 10)`.
Settings: `options(breedtools.germplasm_batch = 10)` for smaller
batches, `options(breedtools.breedbase_timeout = 600)` to wait longer
(seconds) for each answer.

## See also

[`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md),
[`breedbase_programs()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)

## Examples

``` r
rc <- augmented_row_column_design(paste0("G", 1:108), paste0("C", 1:4),
                                  12, 12, 4, 4, n_candidates = 5, seed = 1)
# What would be sent (no server needed):
payload <- save_breedbase(rc, trial_name = "AYT_2026_Ithaca", year = 2026,
                          location = "Ithaca",
                          description = "Advanced yield trial of 108 lines and 4 checks",
                          dry_run = TRUE)
#> Dry run without a connection: returning the data that would be sent.
str(payload$studies[[1]], max.level = 1)
#> List of 6
#>  $ studyName         : chr "AYT_2026_Ithaca"
#>  $ studyDescription  : chr "Advanced yield trial of 108 lines and 4 checks"
#>  $ studyType         : chr "phenotyping_trial"
#>  $ seasons           :List of 1
#>  $ experimentalDesign:List of 2
#>  $ additionalInfo    :List of 1
payload$plots[[1]][[1]]
#> $observationUnitName
#> [1] "AYT_2026_Ithaca_Ithaca_G51_1_1"
#> 
#> $germplasmName
#> [1] "G51"
#> 
#> $observationUnitPosition
#> $observationUnitPosition$observationLevel
#> $observationUnitPosition$observationLevel$levelName
#> [1] "plot"
#> 
#> $observationUnitPosition$observationLevel$levelCode
#> [1] "1"
#> 
#> 
#> $observationUnitPosition$observationLevelRelationships
#> $observationUnitPosition$observationLevelRelationships[[1]]
#> $observationUnitPosition$observationLevelRelationships[[1]]$levelName
#> [1] "rep"
#> 
#> $observationUnitPosition$observationLevelRelationships[[1]]$levelCode
#> [1] "1"
#> 
#> 
#> $observationUnitPosition$observationLevelRelationships[[2]]
#> $observationUnitPosition$observationLevelRelationships[[2]]$levelName
#> [1] "block"
#> 
#> $observationUnitPosition$observationLevelRelationships[[2]]$levelCode
#> [1] "1"
#> 
#> 
#> 
#> $observationUnitPosition$entryType
#> [1] "TEST"
#> 
#> $observationUnitPosition$positionCoordinateX
#> [1] "1"
#> 
#> $observationUnitPosition$positionCoordinateXType
#> [1] "GRID_COL"
#> 
#> $observationUnitPosition$positionCoordinateY
#> [1] "1"
#> 
#> $observationUnitPosition$positionCoordinateYType
#> [1] "GRID_ROW"
#> 
#> 

if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
               trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
               description = "Advanced yield trial of 108 lines and 4 checks, Ibadan 2026",
               create_folder = TRUE, dry_run = TRUE)   # check first
save_breedbase(rc, con, program = "IITA", folder = "AYT 2026",
               trial_name = "AYT_2026_Ibadan", year = 2026, location = "Ibadan",
               description = "Advanced yield trial of 108 lines and 4 checks, Ibadan 2026")
} # }
```
