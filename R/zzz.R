.onLoad <- function(libname, pkgname) {
  initScClass()
  packageStartupMessage("Using sciClone version 1.1")
}

# Internal replacement for NORMT3::erf() using base R identity
# erf(x) = 2 * pnorm(x * sqrt(2)) - 1  (exact for real x)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
