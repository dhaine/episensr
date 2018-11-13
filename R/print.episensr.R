#' Print associations for episensr class
#'
#' Print associations for episensr objects.
#'
#' @param x An object of class 'episensr'.
#' @param digits Minimal number of _significant_ digits, see 'print.default'. 
#' @param ... Other unused arguments.
#'
#' @return Print the observed and adjusted measures of association.
#'
#' @export
print.episensr <- function(x, digits = getOption("digits"), ...) {
    if(class(x)[1] == "episensr") {
        cat("--Observed data--",
            "\n         Outcome:", rownames(x$obs.data)[1],
            "\n       Comparing:", colnames(x$obs.data)[1], "vs.",
            colnames(x$obs.data)[2], "\n\n")
        print.table(x$obs.data, digits = digits, ...)
        cat("\n")
        print.table(x$obs.measures, digits = digits, ...)
        cat("---\n")
        print.table(x$adj.measures, digits = digits, ...)
        invisible(NULL)
    } else if(class(x)[1] == "episensr.multiple") {
        cat("\nMultiple bias analysis\n")
        cat("---\n")
        print.table(x$adj.measures, digits = digits, ...)
        invisible(NULL)
    }
}
