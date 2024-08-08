#' Plot of bootstrap simulation output for selection and misclassification bias
#'
#' This takes an episensr bootstrap object and produces the plot of bootstrap
#' replicates for selection or misclassification bias of the variable of interest,
#' either relative risk or odds ratio. It also draws the confidence interval.
#'
#' @param x An object of class "episensr.booted" returned from the episensr bootstrap generation function.
#' @param association Choice between bias adjusted relative risk (rr) and odds ratio (or).
#' @param ... Other unused arguments.
#'
#' @family visualization
#'
#' @seealso \code{\link{boot.bias}, \link{boot}, \link{selection}, \link{misclass}}
#'
#' @examples
#' misclass_eval <- misclass(matrix(c(215, 1449, 668, 4296),
#' dimnames = list(c("Breast cancer+", "Breast cancer-"),
#' c("Smoker+", "Smoker-")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure",
#' bias_parms = c(.78, .78, .99, .99))
#'
#' set.seed(123)
#' misclass_boot <- boot.bias(misclass_eval)
#' plot(misclass_boot, association = "rr")
#'
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes after_stat geom_density geom_vline ggtitle xlab
plot.episensr.booted <- function(x,
                                 association = c("rr", "or"),
                                 ...) {
    association <- match.arg(association)
    .data <- as.data.frame(x[[2]]$t)
    colnames(.data) <- c("rr", "or")
    bins <- ceiling(length(.data$rr)/25)
    bins <- ifelse(bins < 10, 10,
            ifelse(bins > 100, 100, bins))

    if (association == "rr") {
        ggplot(.data, aes(x = .data$rr)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$ci[1, 1], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$ci[1, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted relative risk") +
            xlab("Relative risk")
    } else if (association == "or") {
        ggplot(.data, aes(x = .data$or)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$ci[2, 1], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$ci[2, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted odds ratio") +
            xlab("Odds ratio")
    }
}
