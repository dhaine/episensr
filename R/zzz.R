#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Thank you for using episensr!")
    packageStartupMessage("This is version ", packageVersion(pkgname), " of ", pkgname)
    packageStartupMessage("Type 'citation(\"episensr\")' for citing this R package in publications.\n")
}
