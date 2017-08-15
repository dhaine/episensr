#' Print bootstrapped confidence intervals
#'
#' Print bootstrap-ed confidence intervals for selection and misclassification bias functions.
#'
#' @param x An object of class 'episensr.booted'.
#' @param digits Minimal number of _significant_ digits, see 'print.default'. 
#' @param ... Other unused arguments.
#'
#' @return Print the confidence interval of the adjusted measures of association.
#'
#' @export
print.episensr.booted <- function(x, digits = getOption("digits"), ...) {
    cat(100*x$conf,"% confidence interval of the bias adjusted measures:",
        "\n   RR:", x$ci[1, ],
        "\n   OR:", x$ci[2, ],
        "\n---\n",
        "Based on", x$nrep, "bootstrap replicates\n")
    invisible(NULL)
}
