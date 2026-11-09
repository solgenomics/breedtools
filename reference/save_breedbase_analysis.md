# Save analysis results to Breedbase as an analysis

Saves analysis results as a Breedbase analysis (a project of type
analysis), through the same upload as the "Upload analysis" form of
Breedbase:

- `level = "accession"`:

  "Analysis of: accessions": one value per genotype, e.g. second-stage
  adjusted means from
  [`second_stage_values()`](https://solgenomics.github.io/breedtools/reference/second_stage_values.md).

- `level = "plot"`:

  "Analysis of: plots": one value per field plot, e.g. first-stage
  values from
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md).
  Breedbase finds the field trial(s) of the plots and gives the analysis
  one plot per field plot (named `"<analysis_name>_<field plot name>"`),
  with the same accession, rep, block, row and column.

Each trait is saved together with a statistic term, as Breedbase does
for every analysis. The analysis then appears on the Analyses page and
can be used in later analyses. Every name (analysis, accessions or
plots, traits, statistic) is checked first, and all problems are listed
before anything is written. The field trials do not change.

## Usage

``` r
save_breedbase_analysis(
  data,
  con = NULL,
  analysis_name,
  program,
  description,
  year = format(Sys.Date(), "%Y"),
  level = c("accession", "plot"),
  traits = NULL,
  accession = "accession",
  plot = "plot_name",
  estimate = NULL,
  statistic = NULL,
  protocol = NULL,
  time_terms = NULL,
  upload_file = NULL,
  dry_run = FALSE,
  confirm = interactive(),
  folder = NULL,
  create_folder = FALSE
)
```

## Arguments

- data:

  A data frame with one row per accession (`level = "accession"`) or per
  plot (`level = "plot"`) and one numeric column per trait.

- con:

  A connection from
  [`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md);
  `NULL` with `dry_run = TRUE` returns the file that would be uploaded.

- analysis_name:

  Name of the new analysis (must not exist yet). Spaces and special
  characters are replaced with `"_"` (only letters, digits, `_`, `-` and
  `.` are kept).

- program:

  Breeding program of the analysis.

- description:

  Description of the analysis (required).

- year:

  Year of the analysis (default: this year).

- level:

  `"accession"` or `"plot"` (Breedbase's "Analysis of").

- traits:

  Trait columns (default: all except the name and design columns). Their
  names are the Breedbase trait names.

- accession, plot:

  Names of the accession and plot-name columns.

- estimate:

  For the default statistic: `"blue"` or `"blup"` (adjusted means of
  that kind), or `"spatial"` (spatial correction, the default for
  `level = "plot"`). By default taken from the data
  ([`second_stage_values()`](https://solgenomics.github.io/breedtools/reference/second_stage_values.md),
  [`first_stage_values()`](https://solgenomics.github.io/breedtools/reference/first_stage_values.md)).

- statistic:

  Optional statistic term: its name, `"name|SGNSTAT:id"` or the id, as
  listed by
  [`breedbase_statistics()`](https://solgenomics.github.io/breedtools/reference/breedbase_statistics.md).
  It must be a term of the Breedbase statistics ontology (checked
  against the ontology file, with suggestions for misspelled names, and
  then on the server). By default chosen by `estimate`.

- protocol:

  The model, e.g. a formula: "SpATS(yield ~ genotype + PSANOVA(col,
  row))". Default: the description.

- time_terms:

  Optional time terms per trait, as in the form's "Compose Trait Info
  With Time Term", e.g.
  `c("fresh root yield|CO_334:0000013" = "day 300|TIME:0000300")`.

- upload_file:

  Where to write the upload file when the site's web pages are closed by
  a web password (see Details). Default: the analysis name +
  `"_upload.csv"` in the working folder.

- dry_run:

  Logical. Check everything but write nothing.

- confirm:

  Logical. Ask before writing (default in interactive sessions).

- folder:

  Optional folder of the program to put the analysis in, e.g. the folder
  of its trials. Without a folder (the default, `NULL`) the analysis is
  in the breeding program itself. An existing folder is used as it is; a
  new one (special characters replaced with `"_"`) is created when
  `create_folder = TRUE`.

- create_folder:

  Logical. Create `folder` when it does not exist.

## Value

Invisibly, a tibble with the analysis name and id, the number of rows
(accessions or plots) sent and read back, and the analysis page; or,
when the web pages are closed, a list with the upload `file` and the
`form` values.

## Details

The form's fields map to the arguments: Analysis Name (`analysis_name`,
unique), Analysis Description (`description`), Analysis Protocol
(`protocol`), Analysis Of (`level`), Breeding Program (`program`),
Statistical Ontology Term (`statistic`) and the optional time terms
(`time_terms`). The model type is Breedbase's own one for uploaded
results (`uploaded_generic_analysis_model`).

## Sites with a web password

Breedbase creates analyses through its web pages; BrAPI cannot create an
analysis project. Test sites such as cassava-test.sgn.cornell.edu close
their web pages with an extra web password. With it (see
[`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)),
the analysis is uploaded directly. Without it, everything is checked
through BrAPI and the file of the "Upload analysis" form is written,
with the value of every form field printed: upload it in the browser and
Breedbase creates the same analysis project.

## Examples

``` r
prep <- augmented_prep_design(paste0("E", 1:60), locations = c("North", "South"), seed = 1)
d <- prep$design
set.seed(1)
d$yield <- rnorm(nrow(d), 25, 3)
fit <- analyze_augmented_prep_two_stage(d, "yield")
means <- second_stage_values(fit, trait = "fresh root yield")
sent <- save_breedbase_analysis(means, analysis_name = "PREP_2026_stage2", program = "IITA",
                                description = "Two-stage analysis, weighted", dry_run = TRUE)
#> Dry run without a connection: returning the upload that would be sent (trait names and the statistic are completed from Breedbase when connected).
head(sent$file)
#>   stock_names fresh root yield
#> 1         E50 30.0238650149583
#> 2         E45 29.8904456685267
#> 3         E39 28.5504753824034
#> 4         E11 27.8069635554438
#> 5          E2 27.5851896298825
#> 6         E28 27.4592783525859

if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", "jdoe")
# Second stage: one value per accession.
save_breedbase_analysis(means, con, analysis_name = "PREP_2026_stage2", program = "IITA",
                        description = "Two-stage analysis, weighted")
# First stage: one value per field plot ("Analysis of: plots").
stage1 <- first_stage_values(fit, prep, "PREP_2026", trait = "fresh root yield")
save_breedbase_analysis(stage1, con, analysis_name = "PREP_2026_stage1", program = "IITA",
                        description = "SpATS per location", level = "plot")
} # }
```
