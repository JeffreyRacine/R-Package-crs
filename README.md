# crs

![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/crs)

This is the R package `crs' (Categorical Regression Splines) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca) with the invaluable assistance of Zhenghua Nie.

## Installation

You can install the stable version on [CRAN](https://cran.r-project.org/package=crs):

```r
install.packages('crs', dependencies = TRUE)
```

Or download the [zip ball](https://github.com/JeffreyRacine/R-Package-crs/zipball/master) or [tar ball](https://github.com/JeffreyRacine/R-Package-crs/tarball/master), decompress and run `R CMD INSTALL` on it, or install then use the **devtools** package to install the development version:

```r
library(devtools); install_github('R-Package-crs', 'JeffreyRacine')
```

Note Windows users have to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/), while OS X users have to first install [Xcode](https://apps.apple.com/us/app/xcode/id497799835?mt=12) and the command line tools (in OS X 10.9 or higher, once you have Xcode installed, open a terminal and run xcode-select --install).

For more information on this project please visit the maintainer's website (https://www.socialsciences.mcmaster.ca/people/racinej).

## Release Validation

CRAN-ready release claims should use the shared release gate:

```sh
cd /Users/jracine/Development
./release_protocol/run_crs_release_gate.sh
```

The gate builds the source tarball, installs it into a private library, runs
installed smoke tests, runs tarball-first `R CMD check --as-cran`, records
reverse-dependency inventory, and runs the containerized `rchk` native-code
protection check when feasible (`RUN_RCHK=auto` by default).

For changes touching `src/`, NOMAD interface code, registered native
interfaces, or `.C`/`.Call` payload lifetimes, require local `rchk` proof when
infrastructure is available:

```sh
cd /Users/jracine/Development
RUN_RCHK=1 ./release_protocol/run_crs_release_gate.sh
```

Use `RUN_RCHK=auto` for ordinary full release rehearsal; it records a precise
`SKIP` when Docker/rchk infrastructure is unavailable.
