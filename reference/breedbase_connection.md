# Connecting to a Breedbase instance

How to connect breedtools to a 'Breedbase' site, such as the production
site <https://cassavabase.org> or the test site
<https://cassava-test.sgn.cornell.edu>. Every function that reads from
or writes to Breedbase takes the connection returned by
[`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md)
as its `con` argument.

## What you need

- The address of the site, e.g. `"https://cassavabase.org"`. Any
  Breedbase site with BrAPI v2 works, for example other crop databases
  on the same software.

- A user account on that site, created on the site itself ("Login" then
  "New user"). Accounts are separate for each site: an account on
  cassavabase.org does not exist on cassava-test, and the reverse.

- To save trials, phenotypes or analyses, your account must be allowed
  to write to the breeding program (the site's curators grant this).
  Reading works with any account.

- The 'httr2' package, and optionally 'askpass' for password pop-ups.

## Production and test sites

A production site (e.g. cassavabase.org) holds the real data of the
breeding programs: everything you save there is shared with them. A test
site (e.g. cassava-test.sgn.cornell.edu) is a copy for trying things
out; use it to learn the workflow or to try a new script before saving
to production. Test sites also protect their web pages with an extra web
password (the browser pop-up "Restricted Test Site"). BrAPI calls, used
for trials, plots and phenotypes, do not need it; the trait composer and
the upload of analyses do. Give it with `site_username` (and
`site_password` or `BREEDBASE_SITE_PASSWORD`). Without it, the login
first tries your Breedbase user name and password there and, in an
interactive session, asks for the web password.

## Passwords

Never write a password in a script. Either leave it empty and type it in
the pop-up, or store it once in your personal `~/.Renviron` file (open
it with `usethis::edit_r_environ()`), then restart R:


    BREEDBASE_PASSWORD=your_breedbase_password
    BREEDBASE_SITE_PASSWORD=the_test_site_web_password

Scripts run without a console (Rscript, scheduled jobs) need these
variables, since there is nobody to answer the pop-up. Only the access
token is kept in the connection; the password is not stored.

## Check the connection

Printing the connection shows the site, the user and whether a web
password is set. List the programs, folders and locations to confirm
that the account sees the data you expect, and use `dry_run = TRUE` in
the `save_*()` functions to check everything without writing. When the
token expires (a long session), log in again.

## Working without a connection

Designs, analyses,
[`prepare_breedbase_trial()`](https://solgenomics.github.io/breedtools/reference/prepare_breedbase_trial.md)
and the dry runs of
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md)
work offline, so a whole design can be prepared before logging in.

## See also

[`breedbase_login()`](https://solgenomics.github.io/breedtools/reference/breedbase_login.md),
[breedbase_lists](https://solgenomics.github.io/breedtools/reference/breedbase_lists.md),
[`breedbase_check()`](https://solgenomics.github.io/breedtools/reference/breedbase_check.md),
[`save_breedbase()`](https://solgenomics.github.io/breedtools/reference/save_breedbase.md),
[`breedbase_phenotypes()`](https://solgenomics.github.io/breedtools/reference/breedbase_phenotypes.md),
[`save_breedbase_analysis()`](https://solgenomics.github.io/breedtools/reference/save_breedbase_analysis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Production site
con <- breedbase_login("https://cassavabase.org", username = "your_username")
con
breedbase_programs(con)

# Test site, with its extra web password
test <- breedbase_login("https://cassava-test.sgn.cornell.edu",
                        username = "your_username",
                        site_username = "web_user_name")
breedbase_locations(test)
} # }
```
