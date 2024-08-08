#' Plot(s) of probabilistic bias analyses
#'
#' This takes a \code{probsens}-family object and produces the distribution plot of
#' chosen bias parameters, as well as distribution of adjusted measures (with confidence
#' interval). It can also produce a forest plot of relative risks or odds ratios (with
#' `probsens()`, `probsens.conf()`, or `probsens.sel()`)
#'
#' @param x An object of class "episensr.probsens" returned from the
#' \code{episensr probsens}, \code{probsens.sel}, \code{probsens.conf}, \code{probsens.irr},
#' \code{probsens.irr.conf} functions.
#' @param parms Choice between adjusted relative risk (\code{rr}) and odds ratio (\code{or}),
#' total error relative risk and odds ratio (\code{rr_tot} and \code{or_tot}), forest
#' plots (\code{forest_rr} and \code{forest_or}), \code{seca}, \code{seexp},
#' \code{spca}, and \code{spexp}, \code{prev.exp}, \code{prev.nexp} and \code{risk},
#' \code{irr} and \code{irr_tot}.
#' @param ... Other unused arguments.
#'
#' @family visualization
#'
#' @seealso \code{\link{probsens}, \link{probsens.sel}, \link{probsens.conf},
#' \link{probsens.irr}, \link{probsens.irr.conf}}
#'
#' @examples
#' set.seed(123)
#' risk <- probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure", reps = 20000,
#' seca = list("trapezoidal", c(.75, .85, .95, 1)),
#' spca = list("trapezoidal", c(.75, .85, .95, 1)))
#' plot(risk, "rr")
#'
#' set.seed(123)
#' odds <- probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure", reps = 20000,
#' seca = list("beta", c(908, 16)),
#' seexp = list("beta", c(156, 56)),
#' spca = list("beta", c(153, 6)),
#' spexp = list("beta", c(205, 18)),
#' corr_se = .8,
#' corr_sp = .8)
#' plot(odds, "seca")
#'
#' set.seed(123)
#' smoke <- probsens(matrix(c(215, 1449, 668, 4296),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure", reps = 20000,
#' seca = list("uniform", c(.7, .95)),
#' spca = list("uniform", c(.9, .99)))
#' plot(smoke, "forest_or")
#'
#' set.seed(123)
#' conf <- probsens_conf(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")), nrow = 2, byrow = TRUE),
#' reps = 20000,
#' prev_exp = list("triangular", c(.7, .9, .8)),
#' prev_nexp = list("trapezoidal", c(.03, .04, .05, .06)),
#' risk = list("triangular", c(.6, .7, .63)),
#' corr_p = .8)
#' plot(conf, "prev_exp")
#'
#' set.seed(123)
#' inc1 <- probsens.irr(matrix(c(2, 67232, 58, 10539000),
#' dimnames = list(c("GBS+", "Person-time"), c("HPV+", "HPV-")), ncol = 2),
#' reps = 20000,
#' seca = list("trapezoidal", c(.4, .45, .55, .6)),
#' spca = list("constant", 1))
#' plot(inc1, "irr")
#'
#' set.seed(123)
#' inc2 <- probsens.irr.conf(matrix(c(77, 10000, 87, 10000),
#' dimnames = list(c("D+", "Person-time"), c("E+", "E-")), ncol = 2),
#' reps = 20000,
#' prev_exp = list("trapezoidal", c(.01, .2, .3, .51)),
#' prev_nexp = list("trapezoidal", c(.09, .27, .35, .59)),
#' risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
#' corr_p = .8)
#' plot(inc2, "risk")
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes after_stat geom_density geom_point geom_linerange geom_vline ggtitle xlab ylab
#' @importFrom forcats fct_rev
plot.episensr.probsens <- function(x,
                                   parms = c("rr", "or", "rr_tot", "or_tot",
                                             "forest_rr", "forest_or",
                                             "irr", "irr_tot",
                                             "seca", "seexp", "spca", "spexp",
                                             "prev_exp", "prev_nexp", "risk"),
                                   ...) {
    .data <- x$sim_df
    prob.fun <- x$fun

    if (length(parms) != 1)
        stop(cli::format_error(c("x" = "Please provide one parameter to plot.")))

    if (prob.fun == "probsens" & !(parms %in% c("rr", "or", "rr_tot", "or_tot",
                                                "forest_rr", "forest_or",
                                                "seca", "seexp", "spca", "spexp")))
        stop(cli::format_error(c("i" = "Please provide parameters to plot: rr, or,
rr_tot, or_tot, seca, seexp, spca, spexp.")))
    if (prob.fun == "probsens.sel" & !(parms %in% c("rr", "or", "rr_tot", "or_tot",
                                                    "forest_rr", "forest_or",
                                                    "case_exp", "case_nexp", "ncase_exp",
                                                    "ncase_nexp")))
        stop(cli::format_error(c("i" = "Please provide parameters to plot: rr, rr_tot,
or, or_tot, case_exp, case_nexp, ncase_exp, ncase_nexp")))
    if (prob.fun == "probsens_conf" & !(parms %in% c("rr", "or", "rr_tot", "or_tot",
                                                     "forest_rr", "forest_or",
                                                     "prev_exp", "prev_nexp", "risk")))
        stop(cli::format_error(c("i" = "Please provide parameters to plot: rr, or,
rr_tot, or_tot, prev_exp, prev_nexp, risk.")))
    if (prob.fun == "probsens.irr" & !(parms %in% c("irr", "irr_tot",
                                                    "seca", "seexp", "spca", "spexp")))
        stop(cli::format_error(c("i" = "Please provide parameters to plot: irr,
irr_tot, seca, seexp, spca, spexp")))
    if (prob.fun == "probsens.irr.conf" & !(parms %in% c("irr", "irr_tot",
                                                         "prev_exp", "prev_nexp", "risk")))
        stop(cli::format_error(c("x" = "Please provide parameters to plot: irr,
irr_tot, prev_exp, prev_nexp, risk.")))

    parms <- match.arg(parms)

    bins <- ceiling(x$reps/25)
    bins <- ifelse(bins < 10, 10,
            ifelse(bins > 100, 100, bins))

    if (prob.fun %in% c("probsens", "probsens.conf", "probsens.sel")) {
        .rr_estimates <- rbind(x$obs_measures[1, ], x$adj_measures[1:2, ])
        .rr_estimates <- as.data.frame(.rr_estimates)
        .rr_estimates <- cbind(.rr_estimates, c("Random error", "Systematic error", "Total error"))
        colnames(.rr_estimates) <- c("est", "lo", "hi", "txt")
        .or_estimates <- rbind(x$obs_measures[2, ], x$adj_measures[3:4, ])
        .or_estimates <- as.data.frame(.or_estimates)
        .or_estimates <- cbind(.or_estimates, c("Random error", "Systematic error", "Total error"))
        colnames(.or_estimates) <- c("est", "lo", "hi", "txt")
    }

    if (parms == "rr" & (prob.fun %in% c("probsens", "probsens_conf", "probsens.sel"))) {
        ggplot(.data, aes(x = .data$corr_RR)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$adj_measures[1, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$adj_measures[1, 3], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted relative risk") +
            xlab("Relative risk")
    } else if (parms == "or" &
               (prob.fun %in% c("probsens", "probsens_conf", "probsens.sel"))) {
        ggplot(.data, aes(x = .data$corr_OR)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$adj_measures[2, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$adj_measures[2, 3], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted odds ratio") +
            xlab("Odds ratio")
    } else if (parms == "rr_tot" &
               (prob.fun %in% c("probsens", "probsens_conf", "probsens.sel"))) {
        ggplot(.data, aes(x = .data$tot_RR)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$adj_measures[3, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$adj_measures[3, 3], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted relative risk (total error)") +
            xlab("Relative risk")
    } else if (parms == "or_tot" &
               (prob.fun %in% c("probsens", "probsens_conf", "probsens.sel"))) {
        ggplot(.data, aes(x = .data$tot_OR)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$adj_measures[4, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$adj_measures[4, 3], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted odds ratio (total error)") +
            xlab("Odds ratio")
    } else if (parms == "forest_rr" &
               (prob.fun %in% c("probsens", "probsens_conf", "probsens.sel"))) {
        ggplot(.rr_estimates, aes(y = forcats::fct_rev(.rr_estimates$txt))) +
            geom_point(aes(x = .rr_estimates$est), shape = 15, size = 3) +
            geom_linerange(aes(xmin = .rr_estimates$lo, xmax = .rr_estimates$hi)) +
            geom_vline(xintercept = 1, linetype = "dashed") +
            xlab("Risk ratio (95% CI)") +
            ylab("")
    } else if (parms == "forest_or" &
               (prob.fun %in% c("probsens", "probsens_conf", "probsens.sel"))) {
        ggplot(.or_estimates, aes(y = forcats::fct_rev(.or_estimates$txt))) +
            geom_point(aes(x = .or_estimates$est), shape = 15, size = 3) +
            geom_linerange(aes(xmin = .or_estimates$lo, xmax = .or_estimates$hi)) +
            geom_vline(xintercept = 1, linetype = "dashed") +
            xlab("Odds ratio (95% CI)") +
            ylab("")
    } else if (parms == "irr") {
        ggplot(.data, aes(x = .data$corr_IRR)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$adj_measures[1, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$adj_measures[1, 3], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted IRR") +
            xlab("IRR")
    } else if (parms == "irr_tot") {
        ggplot(.data, aes(x = .data$tot_IRR)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            geom_vline(xintercept = x$adj_measures[2, 2], size = .75,
                       colour = "black", linetype = "dashed") +
            geom_vline(xintercept = x$adj_measures[2, 3], size = .75,
                       colour = "black", linetype = "dashed") +
            ggtitle("Bias adjusted IRR (total error)") +
            xlab("IRR")
    } else if (parms == "seca") {
        ggplot(.data, aes(x = .data$seca)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("seca")
    } else if (parms == "seexp") {
        ggplot(.data, aes(x = .data$seexp)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("seexp")
    } else if (parms == "spca") {
        ggplot(.data, aes(x = .data$spca)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("spca")
    } else if (parms == "spexp") {
        ggplot(.data, aes(x = .data$spexp)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("spexp")
    } else if (parms == "case_exp") {
        ggplot(.data, aes(x = .data$S1_1)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("case_exp")
    } else if (parms == "case_nexp") {
        ggplot(.data, aes(x = .data$S0_1)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("case_nexp")
    } else if (parms == "ncase_exp") {
        ggplot(.data, aes(x = .data$S1_0)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("ncase_exp")
    } else if (parms == "ncase_nexp") {
        ggplot(.data, aes(x = .data$S0_0)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("ncase_nexp")
    } else if (parms == "prev_exp") {
        ggplot(.data, aes(x = .data$p1)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("prev_exp")
    } else if (parms == "prev_nexp") {
        ggplot(.data, aes(x = .data$p0)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("prev_nexp")
    } else if (parms == "risk") {
        ggplot(.data, aes(x = .data$RR_cd)) +
            geom_histogram(aes(y = after_stat(.data$density)), bins = bins,
                           colour = "grey", fill = "dimgrey") +
            geom_density() +
            xlab("risk")
    }
}
