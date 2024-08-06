#' Probabilistic sensitivity analysis for unmeasured confounding of person-time data and random error.
#'
#' Probabilistic sensitivity analysis to correct for unmeasured confounding when
#' person-time data has been collected.
#'
#' Correlations between prevalences of exposure classification among cases and
#' controls can be specified and use the NORmal To Anything (NORTA)
#' transformation (Li & Hammond, 1975).
#'
#' @section Updated calculations:
#' episensr 2.0.0 introduced updated calculations of probabilistic bias analyses
#' by (1) using the NORTA transformation to define a correlation between
#' distributions, and (2) sampling true prevalences and then sampling the
#' adjusted cell counts rather than just using the expected cell counts from a
#' simple quantitative bias analysis. This updated version should be preferred
#' but if you need to run an old analysis, you can easily revert to the
#' computation using [probsens.irr.conf_legacy()] as follows:
#'
#' ```
#' library(episensr)
#' probsens.irr.conf <- probsens.irr.conf_legacy
#' ```
#'
#' @param counts A table or matrix where first row contains disease counts and
#' second row contains person-time at risk, and first and second columns are
#' exposed and unexposed observations, as:
#' \tabular{lll}{
#' \tab Exposed \tab Unexposed \cr
#' Cases \tab a \tab b \cr
#' Person-time \tab N1 \tab N0
#' }
#' @param pt A numeric vector of person-time at risk. If provided, \code{counts}
#' must be a numeric vector of disease counts.
#' @param reps Number of replications to run.
#' @param prev_exp List defining the prevalence of exposure among the exposed.
#' The first argument provides the probability distribution function (constant,
#' uniform, triangular, trapezoidal, truncated normal, or beta) and the second
#' its parameters as a vector. Lower and upper bounds for truncated normal distribution
#' cannot be les than zero.
#' \enumerate{
#' \item constant; value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item normal: lower bound, upper bound, mean, sd,
#' \item beta: alpha, beta.
#' }
#' @param prev_nexp List defining the prevalence of exposure among the unexposed.
#' @param risk List defining the confounder-disease relative risk or the
#' confounder-exposure odds ratio. The first argument provides the probability
#' distribution function (constant,uniform, triangular, trapezoidal, log-logistic,
#' or log-normal) and the second its parameters as a vector:
#' \enumerate{
#' \item constant: value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item log-logistic: shape, rate. Must be strictly positive,
#' \item log-normal: meanlog, sdlog. This is the mean and standard deviation on the log scale.
#' }
#' @param corr_p Correlation between the exposure-specific confounder prevalences.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs_measures}{A table of observed incidence rate ratio with exact confidence interval.}
#' \item{adj_measures}{A table of corrected incidence rate ratios.}
#' \item{sim_df}{Data frame of random parameters and computed values.}
#'
#' @family confounding
#'
#' @references
#' Li, S.T., Hammond, J.L., 1975. \emph{Generation of Pseudorandom Numbers
#' with Specified Univariate Distributions and Correlation Coefficients}.
#' IEEE Trans Syst Man Cybern 5:557-561.
#' @examples
#' set.seed(123)
#' # Unmeasured confounding
#' probsens.irr.conf(matrix(c(77, 10000, 87, 10000),
#' dimnames = list(c("D+", "Person-time"), c("E+", "E-")), ncol = 2),
#' reps = 20000,
#' prev_exp = list("trapezoidal", c(.01, .2, .3, .51)),
#' prev_nexp = list("trapezoidal", c(.09, .27, .35, .59)),
#' risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
#' corr_p = .8)
#' @export
#' @importFrom stats binom.test median quantile runif rbeta qbeta
probsens.irr.conf <- function(counts,
                         pt = NULL,
                         reps = 1000,
                         prev_exp = list(dist = c("constant", "uniform",
                                                  "triangular", "trapezoidal",
                                                  "normal", "beta"),
                                         parms = NULL),
                         prev_nexp = list(dist = c("constant", "uniform",
                                                   "triangular", "trapezoidal",
                                                   "normal", "beta"),
                                          parms = NULL),
                         risk = list(dist = c("constant", "uniform", "triangular",
                                              "trapezoidal", "log-logistic", "log-normal"),
                                     parms = NULL),
                         corr_p = NULL,
                         alpha = 0.05) {
    if (reps < 1)
        stop(cli::format_error(c("x" = "Wrong number of replications: reps = {reps}",
                                 "i" = "reps must be >= 1")))

    if (is.null(prev_exp) | is.null(prev_nexp))
        stop(cli::format_error(c("x" = "Missing argument(s) for prev_exp or prev_nexp",
                                 "i" = "Please provide prevalences among the exposed and unexposed.")))
    if (is.null(risk))
        stop(cli::format_error(c("x" = "Missing argument(s) for risk",
                                 "i" = "Please provide risk of acquiring outcome.")))

    if (!is.list(prev_exp))
        stop(cli::format_error(c("i" = "Prevalence of exposure among the exposed should be a list.")))
    else prev_exp <- prev_exp
    if ((length(prev_exp) != 2) | (length(prev_nexp) != 2) | (length(risk) != 2))
        stop(cli::format_error(c("i" = "Check distribution parameters")))
    if ((length(prev_exp[[1]]) != 1) | (length(prev_nexp[[1]]) != 1) | (length(risk[[1]]) != 1))
        stop(cli::format_error(c("x" = "Which distribution?")))
    if (!is.null(corr_p) && (prev_exp[[1]] == "constant" | prev_nexp[[1]] == "constant"))
        stop(cli::format_error(c("x" = "No correlated distributions with constant values.")))
    if (prev_exp[[1]] == "constant" & length(prev_exp[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (prev_exp[[1]] == "uniform" & length(prev_exp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (prev_exp[[1]] == "uniform" & prev_exp[[2]][1] >= prev_exp[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is
greater than upper limit.")))
    if (prev_exp[[1]] == "triangular" & length(prev_exp[[2]]) != 3)
        stop(cli::format_error(c("x" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (prev_exp[[1]] == "triangular" & ((prev_exp[[2]][1] > prev_exp[[2]][3]) |
                                         (prev_exp[[2]][2] < prev_exp[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (prev_exp[[1]] == "trapezoidal" & length(prev_exp[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide
vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (prev_exp[[1]] == "trapezoidal" & ((prev_exp[[2]][1] > prev_exp[[2]][2]) |
                                          (prev_exp[[2]][2] > prev_exp[[2]][3]) |
                                          (prev_exp[[2]][3] > prev_exp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (prev_exp[[1]] == "normal" & (length(prev_exp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide vector of lower and upper bound limits, mean and SD")))
    if (prev_exp[[1]] == "normal" & length(prev_exp[[2]]) == 4 &
        ((prev_exp[[2]][1] >= prev_exp[[2]][2]) | (prev_exp[[2]][1] < 0)))
        stop(cli::format_error(c("x" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if ((prev_exp[[1]] == "constant" | prev_exp[[1]] == "uniform" |
         prev_exp[[1]] == "triangular" | prev_exp[[1]] == "trapezoidal") &
        !all(prev_exp[[2]] >= 0 & prev_exp[[2]] <= 1))
        stop(cli::format_error(c("x" = "Prevalence should be between 0 and 1.")))
    if (!is.null(prev_exp) && prev_exp[[1]] == "beta" && length(prev_exp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(prev_exp) && prev_exp[[1]] == "beta" && (prev_exp[[2]][1] < 0 | prev_exp[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.list(prev_nexp))
        stop(cli::format_error(c("i" = "Prevalence of exposure among the non-exposed should be a list.")))
    else prev_nexp <- prev_nexp
    if (prev_nexp[[1]] == "constant" & length(prev_nexp[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (prev_nexp[[1]] == "uniform" & length(prev_nexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (prev_nexp[[1]] == "uniform" & prev_nexp[[2]][1] >= prev_nexp[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution
is greater than upper limit.")))
    if (prev_nexp[[1]] == "triangular" & length(prev_nexp[[2]]) != 3)
        stop(cli::format_error(c("x" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (prev_nexp[[1]] == "triangular" & ((prev_nexp[[2]][1] > prev_nexp[[2]][3]) |
                                          (prev_nexp[[2]][2] < prev_nexp[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (prev_nexp[[1]] == "trapezoidal" & length(prev_nexp[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide
vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (prev_nexp[[1]] == "trapezoidal" & ((prev_nexp[[2]][1] > prev_nexp[[2]][2]) |
                                           (prev_nexp[[2]][2] > prev_nexp[[2]][3]) |
                                           (prev_nexp[[2]][3] > prev_nexp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (prev_nexp[[1]] == "normal" & (length(prev_nexp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide vector of lower and upper bound limits, mean and SD.")))
    if (prev_nexp[[1]] == "normal" & length(prev_nexp[[2]]) == 4 &
        ((prev_nexp[[2]][1] >= prev_nexp[[2]][2]) | (prev_nexp[[2]][1] < 0)))
        stop(cli::format_error(c("x" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if ((prev_nexp[[1]] == "constant" | prev_nexp[[1]] == "uniform" |
         prev_nexp[[1]] == "triangular" | prev_nexp[[1]] == "trapezoidal") &
        !all(prev_nexp[[2]] >= 0 & prev_nexp[[2]] <= 1))
        stop(cli::format_error(c("x" = "Prevalence should be between 0 and 1.")))
    if (!is.null(prev_nexp) && prev_nexp[[1]] == "beta" && length(prev_nexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(prev_nexp) && prev_nexp[[1]] == "beta" && (prev_nexp[[2]][1] < 0 | prev_nexp[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.list(risk))
        stop(cli::format_error(c("i" = "Risk should be a list.")))
    else risk <- risk
    if (risk[[1]] == "constant" & length(risk[[2]]) != 1)
        stop(cli::format_error(c("x" = "For constant value, please provide a single value.")))
    if (risk[[1]] == "uniform" & length(risk[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (risk[[1]] == "uniform" & risk[[2]][1] >= risk[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution
is greater than upper limit.")))
    if (risk[[1]] == "triangular" & length(risk[[2]]) != 3)
        stop(cli::format_error(c("i" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (risk[[1]] == "triangular" & ((risk[[2]][1] > risk[[2]][3]) | (risk[[2]][2] < risk[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (risk[[1]] == "trapezoidal" & length(risk[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please
provide vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (risk[[1]] == "trapezoidal" & ((risk[[2]][1] > risk[[2]][2]) |
                                      (risk[[2]][2] > risk[[2]][3]) |
                                      (risk[[2]][3] > risk[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (risk[[1]] == "log-logistic" & length(risk[[2]]) != 2)
        stop(cli::format_error(c("x" = "For log-logistic distribution, please provide vector of location and scale.")))
    if (risk[[1]] == "log-normal" & length(risk[[2]]) != 2)
        stop(cli::format_error(c("x" = "For log-logistic distribution, please provide vector of meanlog and sdlog.")))

    if (!is.null(corr_p) && (corr_p == 0 | corr_p == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))

    if (!is.null(pt) && inherits(counts, c("table", "matrix")))
        stop(cli::format_error(c("x" = "pt argument should be NULL.")))
    if (!inherits(counts, c("vector", "table", "matrix")))
        stop(cli::format_error(c("x" = "counts argument should be a vector, a table, or a matrix.")))
    if (is.null(pt) && inherits(counts, c("table", "matrix")))
        tab <- counts
    else tab <- rbind(counts, pt)
    a <- as.numeric(tab[1, 1])
    b <- as.numeric(tab[1, 2])
    c <- as.numeric(tab[2, 1])
    d <- as.numeric(tab[2, 2])

    draws <- matrix(NA, nrow = reps, ncol = 15)
    colnames(draws) <- c("p1", "p0", "RR_cd",
                         "M1", "N1", "A1", "B1",
                         "M0", "N0", "A0", "B0",
                         "corr_IRR", "reps", "tot_IRR", "flag")
    corr_draws <- matrix(NA, nrow = reps, ncol = 2)

    p1 <- c(reps, prev_exp[[2]])
    p0 <- c(reps, prev_nexp[[2]])
    rr_cd <- c(reps, risk[[2]])

    cli::cli_alert_info("Calculating observed measures")
    obs_irr <- (a / c) / (b / d)
    lci_obs_irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1]) * c)
    uci_obs_irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2]) * c)

    cli::cli_progress_step("Assign probability distributions", spinner = TRUE)
    if (is.null(corr_p)) {
        if (prev_exp[[1]] == "constant") {
            draws[, 1] <- prev_exp[[2]]
        }
        if (prev_exp[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(p1))
        }
        if (prev_exp[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(p1))
        }
        if (prev_exp[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(p1))
        }
        if (prev_exp[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::rtruncnorm, as.list(p1))
        }
        if (prev_exp[[1]] == "beta") {
            draws[, 1] <- do.call(rbeta, as.list(p1))
        }
        if (prev_nexp[[1]] == "constant") {
            draws[, 2] <- prev_nexp[[2]]
        }
        if (prev_nexp[[1]] == "uniform") {
            draws[, 2] <- do.call(runif, as.list(p0))
        }
        if (prev_nexp[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::rtriangle, as.list(p0))
        }
        if (prev_nexp[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::rtrapezoid, as.list(p0))
        }
        if (prev_nexp[[1]] == "normal") {
            draws[, 2] <- do.call(truncnorm::rtruncnorm, as.list(p0))
        }
        if (prev_nexp[[1]] == "beta") {
            draws[, 2] <- do.call(rbeta, as.list(p0))
        }
    } else {
        if (is.null(corr_p)) corr_p <- .8
        norta_prev <- matrix(c(1, corr_p, corr_p, 1), ncol = 2)
        corr_draws <- MASS::mvrnorm(reps, c(0, 0), norta_prev)
        corr_draws <- pnorm(corr_draws)

        if (prev_exp[[1]] == "uniform") {
            draws[, 1] <- do.call(qunif, c(list(corr_draws[, 1]), as.list(p1[-1])))
        }
        if (prev_exp[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::qtriangle, c(list(corr_draws[, 1]), as.list(p1[-1])))
        }
        if (prev_exp[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::qtrapezoid, c(list(corr_draws[, 1]), as.list(p1[-1])))
        }
        if (prev_exp[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::qtruncnorm, c(list(corr_draws[, 1]), as.list(p1[-1])))
        }
        if (prev_exp[[1]] == "beta") {
            draws[, 1] <- do.call(qbeta, c(list(corr_draws[, 1]), as.list(p1[-1])))
        }

        if (prev_nexp[[1]] == "uniform") {
            draws[, 2] <- do.call(qunif, c(list(corr_draws[, 2]), as.list(p0[-1])))
        }
        if (prev_nexp[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::qtriangle, c(list(corr_draws[, 2]), as.list(p0[-1])))
        }
        if (prev_nexp[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::qtrapezoid, c(list(corr_draws[, 2]), as.list(p0[-1])))
        }
        if (prev_nexp[[1]] == "normal") {
            draws[, 2] <- do.call(truncnorm::qtruncnorm, c(list(corr_draws[, 2]), as.list(p0[-1])))
        }
        if (prev_nexp[[1]] == "beta") {
            draws[, 2] <- do.call(qbeta, c(list(corr_draws[, 2]), as.list(p0[-1])))
        }
    }

    if (risk[[1]] == "constant") {
        draws[, 3] <- risk[[2]]
    }
    if (risk[[1]] == "uniform") {
        draws[, 3] <- do.call(runif, as.list(rr_cd))
    }
    if (risk[[1]] == "triangular") {
        draws[, 3] <- do.call(triangle::rtriangle, as.list(rr_cd))
    }
    if (risk[[1]] == "trapezoidal") {
        draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(rr_cd))
    }
    if (risk[[1]] == "log-logistic") {
        draws[, 3] <- do.call(actuar::rllogis, as.list(rr_cd))
    }
    if (risk[[1]] == "log-normal") {
        draws[, 3] <- do.call(rlnorm, as.list(rr_cd))
    }

    cli::cli_progress_step("Bias analysis", spinner = TRUE)
    draws[, 13] <- runif(reps)

    draws[, 4] <- c * draws[, 1]
    draws[, 6] <- (draws[, 3] * draws[, 4] * a) / ((draws[, 3] * draws[, 6]) + c - draws[, 6])
    draws[, 5] <- d * draws[, 2]
    draws[, 7] <- (draws[, 3] * draws[, 5] * b) / ((draws[, 3] * draws[, 5]) + d - draws[, 5])
    draws[, 10] <- a - draws[, 6]
    draws[, 8] <- c - draws[, 4]
    draws[, 11] <- b - draws[, 7]
    draws[, 9] <- d - draws[, 5]

    draws[, 12] <- a /
        ((draws[, 4] * draws[, 7] / draws[, 5]) +
             (draws[, 8] * draws[, 11] / draws[, 9]))

    ## Clean up
    draws[, 15] <- apply(draws[, c(4:5, 7, 8, 11)], MARGIN = 1, function(x) sum(x > 0))
    draws[, 15] <- ifelse(draws[, 15] != 5 | is.na(draws[, 15]), NA, 1)
    discard <- sum(is.na(draws[, 15]))
    if (sum(is.na(draws[, 15])) > 0) {
        cli::cli_alert_warning("Chosen prevalence distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
        neg_warn <- paste("Prior prevalence distributions lead to",  discard, "impossible value(s).")
    } else neg_warn <- NULL

    draws <- draws[draws[, 15] == 1 & !is.na(draws[, 15]), ]

    draws[, 14] <- exp(log(draws[, 12]) -
                           qnorm(draws[, 13]) *
                               ((log(uci_obs_irr) - log(lci_obs_irr)) /
                                    (qnorm(.975) * 2)))

    corr_irr <- c(median(draws[, 12], na.rm = TRUE),
                  quantile(draws[, 12], probs = .025, na.rm = TRUE),
                  quantile(draws[, 12], probs = .975, na.rm = TRUE))
    tot_irr <- c(median(draws[, 14], na.rm = TRUE),
                 quantile(draws[, 14], probs = .025, na.rm = TRUE),
                 quantile(draws[, 14], probs = .975, na.rm = TRUE))

    if (is.null(rownames(tab)))
        rownames(tab) <- c("Cases", "Person-time")
    if (is.null(colnames(tab)))
        colnames(tab) <- c("Exposed", "Unexposed")
    rmat <- matrix(c(obs_irr, lci_obs_irr, uci_obs_irr), nrow = 1)
    rownames(rmat) <- " Observed Incidence Rate ratio:"
    colnames(rmat) <- c(" ",
                        paste(100 * (alpha / 2), "%", sep = ""),
                        paste(100 * (1 - alpha / 2), "%", sep = ""))
    rmatc <- rbind(corr_irr, tot_irr)
    rownames(rmatc) <- c("           Incidence Rate Ratio -- systematic error:",
                         "Incidence Rate Ratio -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    res <- list(obs_data = tab,
                obs_measures = rmat,
                adj_measures = rmatc,
                sim_df = as.data.frame(draws[, -c(13, 15)]),
                reps = reps,
                fun = "probsens.irr.conf",
                warnings = neg_warn)
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
