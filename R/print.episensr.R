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
    if (class(x)[1] == "episensr") {
        cli::cli_h1("Observed data")
        cli::cli_par()
        cli::cli_ul(c("Outcome: {rownames(x$obs.data)[1]}",
                      "Comparing: {colnames(x$obs.data)[1]} vs. {colnames(x$obs.data)[2]}"))
        cli::cli_end()
        cli::cli_par()
        print.table(x$obs.data, digits = digits, ...)
        cli::cli_end()
        cli::cli_par()
        print.table(x$obs.measures, digits = digits, ...)
        cli::cli_end()
        cli::cli_h2("Bias-adjusted measures")
        print.table(x$adj.measures, digits = digits, ...)
        invisible(NULL)
    } else if (class(x)[1] == "episensr.multiple") {
        cat("\nMultiple bias analysis\n")
        cat("---\n")
        print.table(x$adj.measures, digits = digits, ...)
        invisible(NULL)
    } else if (class(x)[1] == "episensr.evalue") {
        cat("\n--E-value--\n")
        print.table(x, digits = digits, ...)
        invisible(NULL)
    } else if (class(x)[1] == "episensr.confounder") {
        cat("--Input bias parameters--\n")
        print.table(x$bias.parms, digits = digits, ...)
        cat("---\n\n")
        print.table(x$adj.measures, digits = digits, ...)
        invisible(NULL)
    }
}
