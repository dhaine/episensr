.onAttach <- function(libname, pkgname) {
    message("Thank you for using episensr!")
    packageStartupMessage("This is version ", packageVersion(pkgname), " of ", pkgname)
    packageStartupMessage("Type 'citation(\"episensr\")' for citing this R package in publications.\n")
}
