# Check that everything is ready before starting

Run this at the start of a script. It starts a fresh list of checks (so
the final
[`check_summary()`](https://solgenomics.github.io/breedtools/reference/check_summary.md)
only counts this run) and checks that you are ready to design and
analyse trials:

- breedtools is loaded, and which version;

- the optional packages 'lme4' (mixed models with random blocks, random
  entries or random genotype-by-location) and 'SpATS' (spatial models)
  are installed; a missing one is a note, not a failure, and the message
  says what it is needed for;

- the folder for your files exists (it is created if needed) and files
  can be written to it.

## Usage

``` r
check_setup(out_dir = NULL)
```

## Arguments

- out_dir:

  Optional folder where the script will write field books and other
  files.

## Value

Invisibly, `TRUE` if everything needed is in place.

## See also

[`check_summary()`](https://solgenomics.github.io/breedtools/reference/check_summary.md),
[`check_design()`](https://solgenomics.github.io/breedtools/reference/check_design.md),
[`check_result()`](https://solgenomics.github.io/breedtools/reference/check_result.md)

## Examples

``` r
check_setup(out_dir = file.path(tempdir(), "my_trials"))
#> 
#> == 0. Setup ==
#>   [OK]   breedtools 0.1.0 is loaded
#>   [OK]   lme4 is installed (needed for mixed models: random blocks, random entries, random genotype x location)
#>   [OK]   SpATS is installed (needed for spatial models)
#>   [OK]   files can be saved in /tmp/RtmpRMxJvq/my_trials
#> Everything is ready to start a new design.
```
