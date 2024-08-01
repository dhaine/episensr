#' Probabilistic sensitivity analysis for exposure misclassification of person-time data and random error.
#'
#' Probabilistic sensitivity analysis to correct for exposure misclassification
#' when person-time data has been collected.
#' Non-differential misclassification is assumed when only the two bias parameters
#' \code{seca} and \code{spca} are provided. Adding the 2 parameters
#' \code{seexp} and \code{spexp} (i.e. providing the 4 bias parameters)
#' evaluates a differential misclassification.
#'
#' Correlations between sensitivity (or specificity) of exposure classification
#' among cases and controls can be specified and use the NORmal To Anything
#' (NORTA) transformation (Li & Hammond, 1975).
#'
#' @section Updated calculations:
#' episensr 2.0.0 introduced updated calculations of probabilistic bias analyses
#' by (1) using the NORTA transformation to define a correlation between
#' distributions, and (2) sampling true prevalences and then sampling the
#' adjusted cell counts rather than just using the expected cell counts from a
#' simple quantitative bias analysis. This updated version should be preferred
#' but if you need to run an old analysis, you can easily revert to the
#' computation using [probsens.irr_legacy()] as follows:
#'
#' ```
#' library(episensr)
#' probsens.irr <- probsens.irr_legacy
#' ```
#'
#' @param counts A table or matrix where first row contains disease counts and
#' second row contains person-time at risk, and first and second columns are exposed
#' and unexposed observations, as:
#' \tabular{lll}{
#' \tab Exposed \tab Unexposed \cr
#' Cases \tab a \tab b \cr
#' Person-time \tab N1 \tab N0
#' }
#' @param pt A numeric vector of person-time at risk. If provided, \code{counts}
#' must be a numeric vector of disease counts.
#' @param reps Number of replications to run.
#' @param seca List defining the sensitivity of exposure classification among
#' those with the outcome. The first argument provides the probability distribution
#' function (uniform, triangular, trapezoidal, truncated normal, or beta) and the second
#' its parameters as a vector. Lower and upper bounds of the truncated normal have
#' to be between 0 and 1.
#' \enumerate{
#' \item constant: constant value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max,
#' \item normal: lower bound, upper bound, mean, sd,
#' \item beta: alpha, beta.
#' }
#' @param seexp List defining the sensitivity of exposure classification among those without the outcome.
#' @param spca List defining the specificity of exposure classification among those with the outcome.
#' @param spexp List defining the specificity of exposure classification among those without the outcome.
#' @param corr_se Correlation between case and non-case sensitivities.
#' @param corr_sp Correlation between case and non-case specificities.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs_measures}{A table of observed incidence rate ratio with exact confidence interval.}
#' \item{adj_measures}{A table of corrected incidence rate ratios.}
#' \item{sim_df}{Data frame of random parameters and computed values.}
#'
#' @references Li, S.T., Hammond, J.L., 1975. \emph{Generation of Pseudorandom Numbers
#' with Specified Univariate Distributions and Correlation Coefficients}.
#' IEEE Trans Syst Man Cybern 5:557-561.
#' @examples
#' set.seed(123)
#' # Exposure misclassification, non-differential
#' probsens.irr(matrix(c(2, 67232, 58, 10539000),
#' dimnames = list(c("GBS+", "Person-time"), c("HPV+", "HPV-")), ncol = 2),
#' reps = 20000,
#' seca = list("trapezoidal", c(.4, .45, .55, .6)),
#' spca = list("constant", 1))
#' @export
#' @importFrom stats binom.test median quantile qnorm runif qbeta rbeta
probsens.irr <- function(counts,
                         pt = NULL,
                         reps = 1000,
                         seca = list(dist = c("constant", "uniform",
                                              "triangular", "trapezoidal",
                                              "normal", "beta"),
                                     parms = NULL),
                         seexp = NULL,
                         spca = list(dist = c("constant", "uniform",
                                              "triangular", "trapezoidal",
                                              "normal", "beta"),
                                     parms = NULL),
                         spexp = NULL,
                         corr_se = NULL,
                         corr_sp = NULL,
                         alpha = 0.05) {
    if (reps < 1)
        stop(cli::format_error(c("x" = "Wrong number of replications: reps = {reps}",
                                 "i" = "reps must be >= 1")))

    if (is.null(seca) | is.null(spca))
        stop(cli::format_error(c("x" = "Missing argument(s) for seca or spca",
                                 "i" = "At least one Se and one Sp should be
provided through outcome parameters.")))
    if (!is.list(seca))
        stop(cli::format_error(c("i" = "Sensitivity of exposure classification among
those with the outcome should be a list.")))
    else seca <- seca
    if ((length(seca) != 2) | (length(spca) != 2))
        stop(cli::format_error(c("i" = "Check distribution parameters")))
    if ((!is.null(seexp) & length(seexp) != 2) |
        (!is.null(spexp) & length(spexp) != 2))
        stop(cli::format_error(c("i" = "Check distribution parameters")))
    if ((length(seca[[1]]) != 1) | (length(spca[[1]]) != 1))
        stop(cli::format_error(c("x" = "Which distribution?")))
    if ((!is.null(seexp[[1]]) & length(seexp[[1]]) != 1) |
        (!is.null(spexp[[1]]) & length(spexp[[1]]) != 1))
        stop(cli::format_error(c("x" = "Which distribution?")))
    if (!is.null(corr_se) && (seca[[1]] == "constant" | seexp[[1]] == "constant"))
        stop(cli::format_error(c("x" = "No correlated distributions with constant values.")))
    if (!is.null(corr_sp) && (spca[[1]] == "constant" | spexp[[1]] == "constant"))
        stop(cli::format_error(c("x" = "No correlated distributions with constant values.")))
    if (seca[[1]] == "constant" & length(seca[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (seca[[1]] == "uniform" & length(seca[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (seca[[1]] == "uniform" & seca[[2]][1] >= seca[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is
greater than upper limit.")))
    if (seca[[1]] == "triangular" & length(seca[[2]]) != 3)
        stop(cli::format_error(c("x" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (seca[[1]] == "triangular" & ((seca[[2]][1] > seca[[2]][3]) |
                                           (seca[[2]][2] < seca[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (seca[[1]] == "trapezoidal" & length(seca[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide
vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (seca[[1]] == "trapezoidal" & ((seca[[2]][1] > seca[[2]][2]) |
                                            (seca[[2]][2] > seca[[2]][3]) |
                                            (seca[[2]][3] > seca[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (seca[[1]] == "normal" & (length(seca[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide vector of lower and upper bound limits, mean and SD")))
    if (seca[[1]] == "normal" & ((seca[[2]][1] >= seca[[2]][2]) |
                                       (!all(seca[[2]][1:2] >= 0 & seca[[2]][1:2] <= 1))))
        stop(cli::format_error(c("x" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if ((seca[[1]] == "constant" | seca[[1]] == "uniform" |
         seca[[1]] == "triangular" | seca[[1]] == "trapezoidal") &
        !all(seca[[2]] >= 0 & seca[[2]] <= 1))
        stop(cli::format_error(c("x" = "Sensitivity of exposure classification
among those with the outcome should be between 0 and 1.")))
    if (seca[[1]] == "beta" & length(seca[[2]]) != 2)
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (seca[[1]] == "beta" & (seca[[2]][1] < 0 | seca[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(seexp) & !is.list(seexp))
        stop(cli::format_error(c("i" = "Sensitivity of exposure classification
among those without the outcome should be a list.")))
    else seexp <- seexp
    if (!is.null(seexp) && seexp[[1]] == "constant" & length(seexp[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (!is.null(seexp) && seexp[[1]] == "uniform" & length(seexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (!is.null(seexp) && seexp[[1]] == "uniform" && seexp[[2]][1] >= seexp[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution
is greater than upper limit.")))
    if (!is.null(seexp) && seexp[[1]] == "triangular" & length(seexp[[2]]) != 3)
        stop(cli::format_error(c("x" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (!is.null(seexp) && seexp[[1]] == "triangular" &&
        ((seexp[[2]][1] > seexp[[2]][3]) | (seexp[[2]][2] < seexp[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (!is.null(seexp) && seexp[[1]] == "trapezoidal" & length(seexp[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide
vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (!is.null(seexp) && seexp[[1]] == "trapezoidal" &&
        ((seexp[[2]][1] > seexp[[2]][2]) |
         (seexp[[2]][2] > seexp[[2]][3]) |
         (seexp[[2]][3] > seexp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (!is.null(seexp) && seexp[[1]] == "normal" & (length(seexp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide vector of lower and upper bound limits, mean and SD.")))
    if (!is.null(seexp) && seexp[[1]] == "normal" &&
        length(seexp[[2]]) == 4 && ((seexp[[2]][1] >= seexp[[2]][2]) |
                                          (!all(seexp[[2]][1:2] >= 0 & seexp[[2]][1:2] <= 1))))
        stop(cli::format_error(c("x" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if (!is.null(seexp) && (seexp[[1]] == "constant" |
                                  seexp[[1]] == "uniform" |
                                  seexp[[1]] == "triangular" |
                                  seexp[[1]] == "trapezoidal") &
        !all(seexp[[2]] >= 0 & seexp[[2]] <= 1))
        stop(cli::format_error(c("x" = "Sensitivity of exposure classification among
those without the outcome should be between 0 and 1.")))
    if (!is.null(seexp) && seexp[[1]] == "beta" && length(seexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(seexp) && seexp[[1]] == "beta" && (seexp[[2]][1] < 0 |
                                                                seexp[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.list(spca))
        stop(cli::format_error(c("i" = "Specificity of exposure classification
among those with the outcome should be a list.")))
    else spca <- spca
    if (spca[[1]] == "constant" & length(spca[[2]]) != 1)
        stop(cli::format_error(c("x" = "For constant value, please provide a single value.")))
    if (spca[[1]] == "uniform" & length(spca[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (spca[[1]] == "uniform" & spca[[2]][1] >= spca[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution
is greater than upper limit.")))
    if (spca[[1]] == "triangular" & length(spca[[2]]) != 3)
        stop(cli::format_error(c("i" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (spca[[1]] == "triangular" & ((spca[[2]][1] > spca[[2]][3]) |
                                           (spca[[2]][2] < spca[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (spca[[1]] == "trapezoidal" & length(spca[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please
provide vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (spca[[1]] == "trapezoidal" & ((spca[[2]][1] > spca[[2]][2]) |
                                            (spca[[2]][2] > spca[[2]][3]) |
                                            (spca[[2]][3] > spca[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (spca[[1]] == "normal" & (length(spca[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide vector of lower and upper bound limits, mean and SD.")))
    if (spca[[1]] == "normal" & ((spca[[2]][1] >= spca[[2]][2]) |
                                       (!all(spca[[2]][1:2] >= 0 & spca[[2]][1:2] <= 1))))
        stop(cli::format_error(c("x" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0
and 1; lower limit < upper limit).")))
    if ((spca[[1]] == "constant" | spca[[1]] == "uniform" |
         spca[[1]] == "triangular" | spca[[1]] == "trapezoidal") &
        !all(spca[[2]] >= 0 & spca[[2]] <= 1))
        stop(cli::format_error(c("x" = "Specificity of exposure classification
among those with the outcome should be between 0 and 1.")))
    if (spca[[1]] == "beta" & length(spca[[2]]) != 2)
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (spca[[1]] == "beta" & (spca[[2]][1] < 0 | spca[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(spexp) & !is.list(spexp))
        stop(cli::format_error(c("i" = "Specificity of exposure classification
among those without the outcome should be a list.")))
    else spexp <- spexp
    if (!is.null(spexp) && spexp[[1]] == "constant" & length(spexp[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (!is.null(spexp) && spexp[[1]] == "uniform" & length(spexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (!is.null(spexp) && spexp[[1]] == "uniform" && spexp[[2]][1] >= spexp[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is
greater than upper limit.")))
    if (!is.null(spexp) && spexp[[1]] == "triangular" & length(spexp[[2]]) != 3)
        stop(cli::format_error(c("i" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (!is.null(spexp) && spexp[[1]] == "triangular" &&
        ((spexp[[2]][1] > spexp[[2]][3]) | (spexp[[2]][2] < spexp[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (!is.null(spexp) && spexp[[1]] == "trapezoidal" & length(spexp[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide
vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (!is.null(spexp) && spexp[[1]] == "trapezoidal" &&
        ((spexp[[2]][1] > spexp[[2]][2]) |
         (spexp[[2]][2] > spexp[[2]][3]) |
         (spexp[[2]][3] > spexp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" & (length(spexp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution, please
provide vector of lower and upper bound limits, meand and SD.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" &&
        ((spexp[[2]][1] >= spexp[[2]][2]) | (!all(spexp[[2]][1:2] >= 0 &
                                                              spexp[[2]][1:2] <= 1))))
        stop(cli::format_error(c("x" = "For truncated normal distribution, please
provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if (!is.null(spexp) && (spexp[[1]] == "constant" | spexp[[1]] == "uniform" |
                                  spexp[[1]] == "triangular" | spexp[[1]] == "trapezoidal") &
        !all(spexp[[2]] >= 0 & spexp[[2]] <= 1))
        stop(cli::format_error(c("x" = "Specificity of exposure classification
among those without the outcome should be between 0 and 1.")))
    if (!is.null(spexp) && spexp[[1]] == "beta" && length(spexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(spexp) && spexp[[1]] == "beta" &&
        (spexp[[2]][1] < 0 | spexp[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(seexp) & (is.null(spca) | is.null(spexp) | is.null(corr_se) | is.null(corr_sp)))
        stop(cli::format_error(c("i" = "For differential misclassification type,
have to provide Se and Sp for among those with and without the outcome as well as
Se and Sp correlations.")))

    if (!is.null(corr_se) && (corr_se == 0 | corr_se == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))
    if (!is.null(corr_sp) && (corr_sp == 0 | corr_sp == 1))
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

    draws <- matrix(NA, nrow = reps, ncol = 12)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A1", "B1", "C1", "D1",
                         "corr_IRR", "tot_IRR", "reps", "flag")
    corr_draws <- matrix(NA, nrow = reps, ncol = 4)

    se1 <- c(reps, seca[[2]])
    se0 <- c(reps, seexp[[2]])
    sp1 <- c(reps, spca[[2]])
    sp0 <- c(reps, spexp[[2]])

    cli::cli_alert_info("Calculating observed measures")
    obs_irr <- (a / c) / (b / d)
    lci_obs_irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1]) * c)
    uci_obs_irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2]) * c)

    cli::cli_progress_step("Assign probability distributions", spinner = TRUE)
    if (is.null(seexp) & !is.null(spca) & is.null(spexp) &
        is.null(corr_se) & is.null(corr_sp)) {
        if (seca[[1]] == "constant") {
            draws[, 1] <- seca[[2]]
        }
        if (seca[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(se1))
        }
        if (seca[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(se1))
        }
        if (seca[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(se1))
        }
        if (seca[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::rtruncnorm, as.list(se1))
        }
        if (seca[[1]] == "beta") {
            draws[, 1] <- do.call(rbeta, as.list(se1))
        }
        draws[, 2] <- draws[, 1]
        if (spca[[1]] == "constant") {
            draws[, 3] <- spca[[2]]
        }
        if (spca[[1]] == "uniform") {
            draws[, 3] <- do.call(runif, as.list(sp1))
        }
        if (spca[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::rtriangle, as.list(sp1))
        }
        if (spca[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(sp1))
        }
        if (spca[[1]] == "normal") {
            draws[, 3] <- do.call(truncnorm::rtruncnorm, as.list(sp1))
        }
        if (spca[[1]] == "beta") {
            draws[, 3] <- do.call(rbeta, as.list(sp1))
        }
        draws[, 4] <- draws[, 3]
    } else {
        norta_se <- matrix(c(1, corr_se, corr_se, 1), ncol = 2)
        norta_sp <- matrix(c(1, corr_sp, corr_sp, 1), ncol = 2)
        corr_draws[, 1:2] <- MASS::mvrnorm(reps, c(0, 0), norta_se)
        corr_draws[, 3:4] <- MASS::mvrnorm(reps, c(0, 0), norta_sp)
        corr_draws <- pnorm(corr_draws)

        if (seca[[1]] == "uniform") {
            draws[, 1] <- do.call(qunif, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::qtriangle, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::qtrapezoid, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::qtruncnorm, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "beta") {
            draws[, 1] <- do.call(qbeta, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seexp[[1]] == "uniform") {
            draws[, 2] <- do.call(qunif, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::qtriangle, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::qtrapezoid, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "normal") {
            draws[, 2] <- do.call(truncnorm::qtruncnorm, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "beta") {
            draws[, 2] <- do.call(qbeta, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (spca[[1]] == "uniform") {
            draws[, 3] <- do.call(qunif, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::qtriangle, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::qtrapezoid, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "normal") {
            draws[, 3] <- do.call(truncnorm::qtruncnorm, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "beta") {
            draws[, 3] <- do.call(qbeta, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spexp[[1]] == "uniform") {
            draws[, 4] <- do.call(qunif, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "triangular") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "trapezoidal") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "normal") {
            draws[, 4] <- do.call(truncnorm::qtruncnorm, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "beta") {
            draws[, 4] <- do.call(qbeta, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
    }

    cli::cli_progress_step("Bias analysis", spinner = TRUE)
    draws[, 11] <- runif(reps)

    draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) / (draws[, 1] - (1 - draws[, 3]))
    draws[, 6] <- (a + b) - draws[, 5]
    draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) / (draws[, 2] - (1 - draws[, 4]))
    draws[, 8] <- (c + d) - draws[, 7]

    draws[, 9] <- (draws[, 5] / (draws[, 5] + draws[, 7])) / (draws[, 6] / (draws[, 6] + draws[, 8]))

    ## Clean up
    draws[, 12] <- apply(draws[, 5:8], MARGIN = 1, function(x) sum(x > 0))
    draws[, 12] <- ifelse(draws[, 12] != 4 | is.na(draws[, 12]), NA, 1)
    discard <- sum(is.na(draws[, 12]))
    if (sum(is.na(draws[, 12])) > 0) {
        cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
        neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
    } else neg_warn <- NULL

    draws <- draws[draws[, 12] == 1 & !is.na(draws[, 12]), ]

    draws[, 10] <- exp(log(draws[, 9]) -
                           qnorm(draws[, 11]) *
                               ((log(uci_obs_irr) - log(lci_obs_irr)) /
                                    (qnorm(.975) * 2)))

    corr_irr <- c(median(draws[, 9], na.rm = TRUE),
                  quantile(draws[, 9], probs = .025, na.rm = TRUE),
                  quantile(draws[, 9], probs = .975, na.rm = TRUE))
    tot_irr <- c(median(draws[, 10], na.rm = TRUE),
                 quantile(draws[, 10], probs = .025, na.rm = TRUE),
                 quantile(draws[, 10], probs = .975, na.rm = TRUE))

    if (is.null(rownames(tab)))
        rownames(tab) <- c("Cases", "Person-time")
    if (is.null(colnames(tab)))
        colnames(tab) <- c("Exposed", "Unexposed")
    rmat <- matrix(c(obs_irr, lci_obs_irr, uci_obs_irr), nrow = 1)
    rownames(rmat) <- " Observed Incidence Rate Ratio:"
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
                sim_df = as.data.frame(draws[, -c(11, 12)]),
                reps = reps,
                fun = "probsens.irr",
                warnings = neg_warn)
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
