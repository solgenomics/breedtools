# Log in to a Breedbase server

Logs in through the BrAPI token call and returns a connection used by
the other `breedbase_*()` functions and by
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md).
Only the access token is kept; the password is not stored.

## Usage

``` r
breedbase_login(
  url,
  username,
  password = NULL,
  site_username = NULL,
  site_password = NULL
)
```

## Arguments

- url:

  Address of the Breedbase site, e.g. `"https://cassavabase.org"`.

- username:

  Your Breedbase username.

- password:

  Optional password; see Details.

- site_username, site_password:

  Only for sites with an extra web password in front of the website
  (test servers such as cassava-test.sgn.cornell.edu ask for it in a
  browser pop-up, "Restricted Test Site"): that user name and password.
  They are used for the non-BrAPI calls (trait composer, analyses);
  BrAPI does not need them. Usually not needed: when a site asks for a
  web password, the login first tries your Breedbase user name and
  password there, and uses them if they work. Leave `site_password`
  empty to read it from `BREEDBASE_SITE_PASSWORD` or to be asked in a
  pop-up.

## Value

A `breedbase_connection` object.

## Details

Never write your password in a script. Leave `password` empty: it is
read from the environment variable `BREEDBASE_PASSWORD` if set (for
example in your `.Renviron` file), and otherwise you are asked for it in
a pop-up.

## See also

[breedbase_connection](https://solgenomics.github.io/breedtools/reference/breedbase_connection.md)
for a step-by-step guide,
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md),
[`breedbase_programs()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md),
[`breedbase_folders()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md),
[`breedbase_locations()`](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md)

## Examples

``` r
if (FALSE) { # \dontrun{
con <- breedbase_login("https://cassavabase.org", username = "jdoe")
con
} # }
```
