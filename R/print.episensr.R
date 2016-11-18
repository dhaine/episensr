#' Print associations for episensr class
#'
#' Print associations for episensr class.
#'
#' @param x an object of class 'episensr'.
#' @param digits minimal number of _significant_ digits, see 'print.default'. 
#' @param ... other unused arguments.
#'
#' @return print the observed and adjusted measures of association.
#'
#' @export
print.episensr <- function(x, digits = getOption("digits"), ...) {
    cat("--Observed data--",
        "\n         Outcome:", rownames(x$obs.data)[1],
        "\n       Comparing:", colnames(x$obs.data)[1], "vs.", colnames(x$obs.data)[2], "\n\n")
    print.table(x$obs.data, digits = digits, ...)
    cat("\n")
    print.table(x$obs.measures, digits = digits, ...)
    cat("---\n")
    print.table(x$adj.measures, digits = digits, ...)
    invisible(NULL)
}
