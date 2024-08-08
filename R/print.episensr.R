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
        cli::cli_ul(c("Outcome: {rownames(x$obs_data)[1]}",
                      "Comparing: {colnames(x$obs_data)[1]} vs. {colnames(x$obs_data)[2]}"))
        cli::cli_end()
        cli::cli_par()
        print.table(x$obs_data, digits = digits, ...)
        cli::cli_end()
        cli::cli_par()
        print.table(x$obs_measures, digits = digits, ...)
        cli::cli_end()
        cli::cli_h2("Bias-adjusted measures")
        print.table(x$adj_measures, digits = digits, ...)
        invisible(NULL)
    } else if (class(x)[1] == "episensr.multiple") {
        cat("\nMultiple bias analysis\n")
        cat("---\n")
        print.table(x$adj_measures, digits = digits, ...)
        invisible(NULL)
    } else if (class(x)[1] == "episensr.evalue") {
        cli::cli_h1("E-value")
        print.table(x, digits = digits, ...)
        invisible(NULL)
    } else if (class(x)[1] == "episensr.confounder") {
        cli::cli_h2("Input bias parameters")
        print.table(x$bias_parms, digits = digits, ...)
        cli::cli_rule()
        print.table(x$adj_measures, digits = digits, ...)
        invisible(NULL)
    }
}
