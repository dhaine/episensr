#' Plot of bootstrap simulation output for selection and misclassification bias
#'
#' This takes an episensr bootstrap object and produces the pot of bootstrap replicates for selection or misclassification bias of the variable of interest, either relative risk or odds ratio.
#'
#' @param x An object of class "episensr.booted" returned from the episensr bootstrap generation function.
#' @param association Choice between bias adjusted relative risk and odds ratio.
#' @param ... Other unused arguments.
#' 
#' @seealso \code{\link{boot.bias}, \link{boot}, \link{selection}, \link{misclassification}}
#'
#' @examples
#' misclass_eval <- misclassification(matrix(c(215, 1449, 668, 4296),
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
#' @importFrom ggplot2 ggplot geom_histogram aes geom_density geom_vline ggtitle xlab
plot.episensr.booted <- function(x,
                                 association = c("rr", "or"),
                                 ...) {
    rr <- NULL
    or <- NULL
    ..density.. <- NULL
    association <- match.arg(association)
    replicates <- as.data.frame(x[[2]]$t)
    colnames(replicates) <- c("rr", "or")
    bins <- ceiling(length(replicates$r)/25)
    bins <- ifelse(bins < 10, 10,
            ifelse(bins > 100, 100, bins))

    if (association == "rr") {
        bias.plot <- ggplot(replicates, aes(x = rr)) +
            geom_histogram(aes(y = ..density..), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() + 
            geom_vline(xintercept = x$ci[1, 1], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$ci[1, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted relative risk") +
            xlab("Relative risk")
    } else if (association == "or") {
        bias.plot <- ggplot(replicates, aes(x = or)) +
            geom_histogram(aes(y = ..density..), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$ci[1, 1], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$ci[1, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted odds ratio") +
            xlab("Odds ratio")
    }
    
    bias.plot
}
