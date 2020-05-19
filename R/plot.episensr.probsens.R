#' Plot(s) of probabilistic bias analyses
#'
#' This takes a \code{probsens}-family object and produces the distribution plot of
#' chosen bias parameters, as well as distribution of adjusted measures.
#'
#' @param x An object of class "episensr.probsens" returned from the
#' \code{episensr probsens}, \code{probsens.conf} functions.
#' @param dist Choice between adjusted relative risk (\code{rr}) and odds ratio (\code{or}),
#' total error relative risk and odds ratio (\code{rr_tot} and \code{or_tot}), \code{seca},
#' \code{seexp}, \code{spca} and \code{spexp}, \code{prev.exp}, \code{prev.nexp} and
#' \code{risk}.
#' @param ... Other unused arguments.
#' 
#' @seealso \code{\link{probsens}, \link{probsens.conf}}
#'
#' @examples
#' set.seed(123)
#' plot(probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 20000,
#' seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
#' spca.parms = list("trapezoidal", c(.75, .85, .95, 1))))
#'
#' set.seed(123)
#' plot(probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 20000,
#' seca.parms = list("beta", c(908, 16)),
#' seexp.parms = list("beta", c(156, 56)),
#' spca.parms = list("beta", c(153, 6)),
#' spexp.parms = list("beta", c(205, 18)),
#' corr.se = .8,
#' corr.sp = .8))
#'
#' set.seed(123)
#' plot(probsens.conf(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")), nrow = 2, byrow = TRUE),
#' reps = 20000,
#' prev.exp = list("triangular", c(.7, .9, .8)),
#' prev.nexp = list("trapezoidal", c(.03, .04, .05, .06)),
#' risk = list("triangular", c(.6, .7, .63)),
#' corr.p = .8))
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes geom_density geom_vline ggtitle xlab
plot.episensr.probsens <- function(x,
                                   dist = c("rr", "or", "rr_tot", "or_tot",
                                            "seca", "seexp", "spca", "spexp",
                                            "prev.exp", "prev.nexp", "risk"),
                                   ...) {
    obj <- x$sim.df
    rr <- NULL
    or <- NULL
    rr_tot <- NULL
    or_tot <- NULL
    seca <- NULL
    seexp <- NULL
    spca <- NULL
    spexp <- NULL
    prev.exp <- NULL
    prev.nexp <- NULL
    risk <- NULL
    ..density.. <- NULL
    dist <- match.arg(dist)

    bins <- ceiling(x$reps/25)
    bins <- ifelse(bins < 10, 10,
            ifelse(bins > 100, 100, bins))

    if (dist == "rr" & dim(obj)[2] == 12) {
        ggplot(obj, aes(x = corr.RR)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted relative risk") +
            xlab("Relative risk")
    } else if (dist == "or" & dim(obj)[2] == 12) {
        ggplot(obj, aes(x = corr.OR)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted odds ratio") +
            xlab("Odds ratio")
    } else if (dist == "rr_tot" & dim(obj)[2] == 12) {
        ggplot(obj, aes(x = tot.RR)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted relative risk (total error)") +
            xlab("Relative risk")
    } else if (dist == "or_tot" & dim(obj)[2] == 12) {
        ggplot(obj, aes(x = tot.OR)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted odds ratio (total error)") +
            xlab("Odds ratio")
    } else if (dist == "rr" & dim(obj)[2] == 13) {
        ggplot(obj, aes(x = RR.SMR.rr)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted relative risk") +
            xlab("Relative risk (SMR)")
    } else if (dist == "or" & dim(obj)[2] == 13) {
        ggplot(obj, aes(x = OR.SMR.rr)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted odds ratio") +
            xlab("Odds ratio (SMR)")
    } else if (dist == "rr_tot" & dim(obj)[2] == 13) {
        ggplot(obj, aes(x = tot.RRadj.SMR)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted relative risk (total error)") +
            xlab("Relative risk (SMR)")
    } else if (dist == "or_tot" & dim(obj)[2] == 13) {
        ggplot(obj, aes(x = tot.ORadj.SMR)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            ggtitle("Bias adjusted odds ratio (total error)") +
            xlab("Odds ratio (SMR)")
    } else if (dist == "seca") {
        ggplot(obj, aes(x = seca)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            xlab("seca")
    } else if (dist == "seexp") {
        ggplot(obj, aes(x = seexp)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            xlab("seexp")
    } else if (dist == "spca") {
        ggplot(obj, aes(x = spca)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            xlab("spca")
    } else if (dist == "spexp") {
        ggplot(obj, aes(x = spexp)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            xlab("spexp")
    } else if (dist == "prev.exp") {
        ggplot(obj, aes(x = p1)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            xlab("prev.exp")
    } else if (dist == "prev.nexp") {
        ggplot(obj, aes(x = p0)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            xlab("prev.nexp")
    } else if (dist == "risk") {
        ggplot(obj, aes(x = RR.cd)) +
            geom_histogram(aes(y = ..density..), bins = bins) +
            geom_density() + 
            xlab("risk")
    }
}
