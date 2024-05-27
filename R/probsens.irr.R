#' Probabilistic sensitivity analysis for exposure misclassification of person-time data and random error.
#'
#' Probabilistic sensitivity analysis to correct for exposure misclassification
#' when person-time data has been collected.
#' Non-differential misclassification is assumed when only the two bias parameters
#' \code{seca.parms} and \code{spca.parms} are provided. Adding the 2 parameters
#' \code{seexp.parms} and \code{spexp.parms} (i.e. providing the 4 bias parameters)
#' evaluates a differential misclassification.
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
#' @param seca.parms List defining the sensitivity of exposure classification among
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
#' @param seexp.parms List defining the sensitivity of exposure classification among those without the outcome.
#' @param spca.parms List defining the specificity of exposure classification among those with the outcome.
#' @param spexp.parms List defining the specificity of exposure classification among those without the outcome.
#' @param corr.se Correlation between case and non-case sensitivities.
#' @param corr.sp Correlation between case and non-case specificities.
#' @param discard A logical scalar. In case of negative adjusted count, should
#' the draws be discarded? If set to FALSE, negative counts are set to zero.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs.measures}{A table of observed incidence rate ratio with exact confidence interval.}
#' \item{adj.measures}{A table of corrected incidence rate ratios.}
#' \item{sim.df}{Data frame of random parameters and computed values.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#' @examples
#' set.seed(123)
#' # Exposure misclassification, non-differential
#' probsens.irr(matrix(c(2, 67232, 58, 10539000),
#' dimnames = list(c("GBS+", "Person-time"), c("HPV+", "HPV-")), ncol = 2),
#' reps = 20000,
#' seca.parms = list("trapezoidal", c(.4, .45, .55, .6)),
#' spca.parms = list("constant", 1))
#' @export
#' @importFrom stats binom.test median quantile qnorm runif qbeta rbeta
probsens.irr <- function(counts,
                         pt = NULL,
                         reps = 1000,
                         seca.parms = list(dist = c("constant", "uniform",
                                                    "triangular", "trapezoidal",
                                                    "normal", "beta"),
                                           parms = NULL),
                         seexp.parms = NULL,
                         spca.parms = list(dist = c("constant", "uniform",
                                                    "triangular", "trapezoidal",
                                                    "normal", "beta"),
                                           parms = NULL),
                         spexp.parms = NULL,
                         corr.se = NULL,
                         corr.sp = NULL,
                         discard = TRUE,
                         alpha = 0.05) {
    if (reps < 1)
        stop(paste("Invalid argument: reps = ", reps))

    if (is.null(seca.parms) | is.null(spca.parms))
        stop("At least one Se and one Sp should be provided through outcome parameters.")
    if (!is.list(seca.parms))
        stop("Sensitivity of exposure classification among those with the outcome should be a list.")
    else seca.parms <- seca.parms
    if ((length(seca.parms) != 2) | (length(spca.parms) != 2))
        stop("Check distribution parameters")
    if ((!is.null(seexp.parms) & length(seexp.parms) != 2) |
        (!is.null(spexp.parms) & length(spexp.parms) != 2))
        stop("Check distribution parameters")
    if ((length(seca.parms[[1]]) != 1) | (length(spca.parms[[1]]) != 1))
        stop("Which distribution?")
    if ((!is.null(seexp.parms[[1]]) & length(seexp.parms[[1]]) != 1) |
        (!is.null(spexp.parms[[1]]) & length(spexp.parms[[1]]) != 1))
        stop("Which distribution?")
    if (!is.null(corr.se) && (seca.parms[[1]] == "constant" | seexp.parms[[1]] == "constant"))
        stop("No correlated distributions with constant values.")
    if (!is.null(corr.sp) && (spca.parms[[1]] == "constant" | spexp.parms[[1]] == "constant"))
        stop("No correlated distributions with constant values.")
    if (seca.parms[[1]] == "constant" & length(seca.parms[[2]]) != 1)
        stop("For constant value, please provide a single value.")
    if (seca.parms[[1]] == "uniform" & length(seca.parms[[2]]) != 2)
        stop("For uniform distribution, please provide vector of lower and upper limits.")
    if (seca.parms[[1]] == "uniform" & seca.parms[[2]][1] >= seca.parms[[2]][2])
        stop("Lower limit of your uniform distribution is greater than upper limit.")
    if (seca.parms[[1]] == "triangular" & length(seca.parms[[2]]) != 3)
        stop("For triangular distribution, please provide vector of lower, upper limits, and mode.")
    if (seca.parms[[1]] == "triangular" & ((seca.parms[[2]][1] > seca.parms[[2]][3]) |
                                           (seca.parms[[2]][2] < seca.parms[[2]][3])))
        stop("Wrong arguments for your triangular distribution.")
    if (seca.parms[[1]] == "trapezoidal" & length(seca.parms[[2]]) != 4)
        stop("For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")
    if (seca.parms[[1]] == "trapezoidal" & ((seca.parms[[2]][1] > seca.parms[[2]][2]) |
                                            (seca.parms[[2]][2] > seca.parms[[2]][3]) |
                                            (seca.parms[[2]][3] > seca.parms[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (seca.parms[[1]] == "normal" & (length(seca.parms[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")
    if (seca.parms[[1]] == "normal" & ((seca.parms[[2]][1] >= seca.parms[[2]][2]) |
                                       (!all(seca.parms[[2]][1:2] >= 0 & seca.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if ((seca.parms[[1]] == "constant" | seca.parms[[1]] == "uniform" |
         seca.parms[[1]] == "triangular" | seca.parms[[1]] == "trapezoidal") &
        !all(seca.parms[[2]] >= 0 & seca.parms[[2]] <= 1))
        stop("Sensitivity of exposure classification among those with the outcome should be between 0 and 1.")
    if (seca.parms[[1]] == "beta" & length(seca.parms[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (seca.parms[[1]] == "beta" & (seca.parms[[2]][1] < 0 | seca.parms[[2]][2] < 0))
        stop("Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")

    if (!is.null(seexp.parms) & !is.list(seexp.parms))
        stop("Sensitivity of exposure classification among those without the outcome should be a list.")
    else seexp.parms <- seexp.parms
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "constant" & length(seexp.parms[[2]]) != 1)
        stop("For constant value, please provide a single value.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "uniform" & length(seexp.parms[[2]]) != 2)
        stop("For uniform distribution, please provide vector of lower and upper limits.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "uniform" && seexp.parms[[2]][1] >= seexp.parms[[2]][2])
        stop("Lower limit of your uniform distribution is greater than upper limit.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "triangular" & length(seexp.parms[[2]]) != 3)
        stop("For triangular distribution, please provide vector of lower, upper limits, and mode.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "triangular" &&
        ((seexp.parms[[2]][1] > seexp.parms[[2]][3]) | (seexp.parms[[2]][2] < seexp.parms[[2]][3])))
        stop("Wrong arguments for your triangular distribution.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "trapezoidal" & length(seexp.parms[[2]]) != 4)
        stop("For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "trapezoidal" &&
        ((seexp.parms[[2]][1] > seexp.parms[[2]][2]) |
         (seexp.parms[[2]][2] > seexp.parms[[2]][3]) |
         (seexp.parms[[2]][3] > seexp.parms[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "normal" & (length(seexp.parms[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "normal" &&
        length(seexp.parms[[2]]) == 4 && ((seexp.parms[[2]][1] >= seexp.parms[[2]][2]) |
                                          (!all(seexp.parms[[2]][1:2] >= 0 & seexp.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if (!is.null(seexp.parms) && (seexp.parms[[1]] == "constant" |
                                  seexp.parms[[1]] == "uniform" |
                                  seexp.parms[[1]] == "triangular" |
                                  seexp.parms[[1]] == "trapezoidal") &
        !all(seexp.parms[[2]] >= 0 & seexp.parms[[2]] <= 1))
        stop("Sensitivity of exposure classification among those without the outcome should be between 0 and 1.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "beta" && length(seexp.parms[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "beta" && (seexp.parms[[2]][1] < 0 |
                                                                seexp.parms[[2]][2] < 0))
        stop("Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")

    if (!is.list(spca.parms))
        stop("Specificity of exposure classification among those with the outcome should be a list.")
    else spca.parms <- spca.parms
    if (spca.parms[[1]] == "constant" & length(spca.parms[[2]]) != 1)
        stop("For constant value, please provide a single value.")
    if (spca.parms[[1]] == "uniform" & length(spca.parms[[2]]) != 2)
        stop("For uniform distribution, please provide vector of lower and upper limits.")
    if (spca.parms[[1]] == "uniform" & spca.parms[[2]][1] >= spca.parms[[2]][2])
        stop("Lower limit of your uniform distribution is greater than upper limit.")
    if (spca.parms[[1]] == "triangular" & length(spca.parms[[2]]) != 3)
        stop("For triangular distribution, please provide vector of lower, upper limits, and mode.")
    if (spca.parms[[1]] == "triangular" & ((spca.parms[[2]][1] > spca.parms[[2]][3]) |
                                           (spca.parms[[2]][2] < spca.parms[[2]][3])))
        stop("Wrong arguments for your triangular distribution.")
    if (spca.parms[[1]] == "trapezoidal" & length(spca.parms[[2]]) != 4)
        stop("For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")
    if (spca.parms[[1]] == "trapezoidal" & ((spca.parms[[2]][1] > spca.parms[[2]][2]) |
                                            (spca.parms[[2]][2] > spca.parms[[2]][3]) |
                                            (spca.parms[[2]][3] > spca.parms[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (spca.parms[[1]] == "normal" & (length(spca.parms[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")
    if (spca.parms[[1]] == "normal" & ((spca.parms[[2]][1] >= spca.parms[[2]][2]) |
                                       (!all(spca.parms[[2]][1:2] >= 0 & spca.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if ((spca.parms[[1]] == "constant" | spca.parms[[1]] == "uniform" |
         spca.parms[[1]] == "triangular" | spca.parms[[1]] == "trapezoidal") &
        !all(spca.parms[[2]] >= 0 & spca.parms[[2]] <= 1))
        stop("Specificity of exposure classification among those with the outcome should be between 0 and 1.")
    if (spca.parms[[1]] == "beta" & length(spca.parms[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (spca.parms[[1]] == "beta" & (spca.parms[[2]][1] < 0 | spca.parms[[2]][2] < 0))
        stop("Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")

    if (!is.null(spexp.parms) & !is.list(spexp.parms))
        stop("Specificity of exposure classification among those without the outcome should be a list.")
    else spexp.parms <- spexp.parms
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "constant" & length(spexp.parms[[2]]) != 1)
        stop("For constant value, please provide a single value.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "uniform" & length(spexp.parms[[2]]) != 2)
        stop("For uniform distribution, please provide vector of lower and upper limits.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "uniform" && spexp.parms[[2]][1] >= spexp.parms[[2]][2])
        stop("Lower limit of your uniform distribution is greater than upper limit.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "triangular" & length(spexp.parms[[2]]) != 3)
        stop("For triangular distribution, please provide vector of lower, upper limits, and mode.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "triangular" &&
        ((spexp.parms[[2]][1] > spexp.parms[[2]][3]) | (spexp.parms[[2]][2] < spexp.parms[[2]][3])))
        stop("Wrong arguments for your triangular distribution.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "trapezoidal" & length(spexp.parms[[2]]) != 4)
        stop("For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "trapezoidal" &&
        ((spexp.parms[[2]][1] > spexp.parms[[2]][2]) |
         (spexp.parms[[2]][2] > spexp.parms[[2]][3]) |
         (spexp.parms[[2]][3] > spexp.parms[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "normal" & (length(spexp.parms[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "normal" &&
        ((spexp.parms[[2]][1] >= spexp.parms[[2]][2]) | (!all(spexp.parms[[2]][1:2] >= 0 &
                                                              spexp.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if (!is.null(spexp.parms) && (spexp.parms[[1]] == "constant" | spexp.parms[[1]] == "uniform" |
                                  spexp.parms[[1]] == "triangular" | spexp.parms[[1]] == "trapezoidal") &
        !all(spexp.parms[[2]] >= 0 & spexp.parms[[2]] <= 1))
        stop("Specificity of exposure classification among those without the outcome should be between 0 and 1.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "beta" && length(spexp.parms[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "beta" &&
        (spexp.parms[[2]][1] < 0 | spexp.parms[[2]][2] < 0))
        stop("Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")

    if (!is.null(seexp.parms) & (is.null(spca.parms) | is.null(spexp.parms) | is.null(corr.se) | is.null(corr.sp)))
        stop("For differential misclassification type, have to provide Se and Sp for among those with and without the outcome as well as Se and Sp correlations.")

    if (!is.null(corr.se) && (corr.se == 0 | corr.se == 1))
        stop("Correlations should be > 0 and < 1.")
    if (!is.null(corr.sp) && (corr.sp == 0 | corr.sp == 1))
        stop("Correlations should be > 0 and < 1.")

    if (!is.null(pt) && inherits(counts, c("table", "matrix")))
        stop("pt argument should be NULL.")
    if (!inherits(counts, c("vector", "table", "matrix")))
        stop("counts argument should be a vector, a table, or a matrix.")
    if (is.null(pt) && inherits(counts, c("table", "matrix")))
        tab <- counts
    else tab <- rbind(counts, pt)
    a <- as.numeric(tab[1, 1])
    b <- as.numeric(tab[1, 2])
    c <- as.numeric(tab[2, 1])
    d <- as.numeric(tab[2, 2])

    draws <- matrix(NA, nrow = reps, ncol = 11)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A1", "B1", "C1", "D1",
                         "corr.IRR", "tot.IRR", "reps")
    corr.draws <- matrix(NA, nrow = reps, ncol = 4)

    seca <- c(reps, seca.parms[[2]])
    seexp <- c(reps, seexp.parms[[2]])
    spca <- c(reps, spca.parms[[2]])
    spexp <- c(reps, spexp.parms[[2]])

    obs.irr <- (a / c) / (b / d)
    lci.obs.irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1]) * c)
    uci.obs.irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2]) * c)

    if (is.null(seexp.parms) & !is.null(spca.parms) & is.null(spexp.parms) &
        is.null(corr.se) & is.null(corr.sp)) {
        if (seca.parms[[1]] == "constant") {
            draws[, 1] <- seca.parms[[2]]
        }
        if (seca.parms[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(seca))
        }
        if (seca.parms[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(seca))
        }
        if (seca.parms[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(seca))
        }
        if (seca.parms[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::rtruncnorm, as.list(seca))
        }
        if (seca.parms[[1]] == "beta") {
            draws[, 1] <- do.call(rbeta, as.list(seca))
        }
        draws[, 2] <- draws[, 1]
        if (spca.parms[[1]] == "constant") {
            draws[, 3] <- spca.parms[[2]]
        }
        if (spca.parms[[1]] == "uniform") {
            draws[, 3] <- do.call(runif, as.list(spca))
        }
        if (spca.parms[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::rtriangle, as.list(spca))
        }
        if (spca.parms[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(spca))
        }
        if (spca.parms[[1]] == "normal") {
            draws[, 3] <- do.call(truncnorm::rtruncnorm, as.list(spca))
        }
        if (spca.parms[[1]] == "beta") {
            draws[, 3] <- do.call(rbeta, as.list(spca))
        }
        draws[, 4] <- draws[, 3]
    } else {
        norta_se <- matrix(c(1, corr.se, corr.se, 1), ncol = 2)
        norta_sp <- matrix(c(1, corr.sp, corr.sp, 1), ncol = 2)
        corr.draws[, 1:2] <- MASS::mvrnorm(reps, c(0, 0), norta_se)
        corr.draws[, 3:4] <- MASS::mvrnorm(reps, c(0, 0), norta_sp)
        corr.draws <- pnorm(corr.draws)

        if (seca.parms[[1]] == "uniform") {
            draws[, 1] <- do.call(qunif, c(list(corr.draws[, 1]), as.list(seca[-1])))
        }
        if (seca.parms[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::qtriangle, c(list(corr.draws[, 1]), as.list(seca[-1])))
        }
        if (seca.parms[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::qtrapezoid, c(list(corr.draws[, 1]), as.list(seca[-1])))
        }
        if (seca.parms[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::qtruncnorm, c(list(corr.draws[, 1]), as.list(seca[-1])))
        }
        if (seca.parms[[1]] == "beta") {
            draws[, 1] <- do.call(qbeta, c(list(corr.draws[, 1]), as.list(seca[-1])))
        }
        if (seexp.parms[[1]] == "uniform") {
            draws[, 2] <- do.call(qunif, c(list(corr.draws[, 2]), as.list(seexp[-1])))
        }
        if (seexp.parms[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::qtriangle, c(list(corr.draws[, 2]), as.list(seexp[-1])))
        }
        if (seexp.parms[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::qtrapezoid, c(list(corr.draws[, 2]), as.list(seexp[-1])))
        }
        if (seexp.parms[[1]] == "normal") {
            draws[, 2] <- do.call(truncnorm::qtruncnorm, c(list(corr.draws[, 2]), as.list(seexp[-1])))
        }
        if (seexp.parms[[1]] == "beta") {
            draws[, 2] <- do.call(qbeta, c(list(corr.draws[, 2]), as.list(seexp[-1])))
        }
        if (spca.parms[[1]] == "uniform") {
            draws[, 3] <- do.call(qunif, c(list(corr.draws[, 3]), as.list(spca[-1])))
        }
        if (spca.parms[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::qtriangle, c(list(corr.draws[, 3]), as.list(spca[-1])))
        }
        if (spca.parms[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::qtrapezoid, c(list(corr.draws[, 3]), as.list(spca[-1])))
        }
        if (spca.parms[[1]] == "normal") {
            draws[, 3] <- do.call(truncnorm::qtruncnorm, c(list(corr.draws[, 3]), as.list(spca[-1])))
        }
        if (spca.parms[[1]] == "beta") {
            draws[, 3] <- do.call(qbeta, c(list(corr.draws[, 3]), as.list(spca[-1])))
        }
        if (spexp.parms[[1]] == "uniform") {
            draws[, 4] <- do.call(qunif, c(list(corr.draws[, 4]), as.list(spexp[-1])))
        }
        if (spexp.parms[[1]] == "triangular") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid, c(list(corr.draws[, 4]), as.list(spexp[-1])))
        }
        if (spexp.parms[[1]] == "trapezoidal") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid, c(list(corr.draws[, 4]), as.list(spexp[-1])))
        }
        if (spexp.parms[[1]] == "normal") {
            draws[, 4] <- do.call(truncnorm::qtruncnorm, c(list(corr.draws[, 4]), as.list(spexp[-1])))
        }
        if (spexp.parms[[1]] == "beta") {
            draws[, 4] <- do.call(qbeta, c(list(corr.draws[, 4]), as.list(spexp[-1])))
        }
    }

    draws[, 11] <- runif(reps)

    draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) / (draws[, 1] - (1 - draws[, 3]))
    draws[, 6] <- (a + b) - draws[, 5]
    draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) / (draws[, 2] - (1 - draws[, 4]))
    draws[, 8] <- (c + d) - draws[, 7]

    draws[, 9] <- (draws[, 5] / (draws[, 5] + draws[, 7])) / (draws[, 6] / (draws[, 6] + draws[, 8]))

    draws[, 9] <- ifelse(draws[, 5] < 0 |
                             draws[, 6] < 0 |
                                 draws[, 7] < 0 |
                                     draws[, 8] < 0, NA, draws[, 9])
    if (all(is.na(draws[, 9]))) {
        warning("Prior Se/Sp distributions lead to all negative adjusted counts.")
        neg_warn <- "Prior Se/Sp distributions lead to all negative adjusted counts."
    } else neg_warn <- NULL
    if (discard) {
        if (sum(is.na(draws[, 9])) > 0) {
            message("Chosen prior Se/Sp distributions lead to ",
                    sum(is.na(draws[, 9])),
                    " negative adjusted counts which were discarded.")
            discard_mess <- c(paste("Chosen prior Se/Sp distributions lead to ",
                                    sum(is.na(draws[, 9])),
                                    " negative adjusted counts which were discarded."))
        } else discard_mess <- NULL
    } else {
        if (sum(is.na(draws[, 9])) > 0) {
            message("Chosen prior Se/Sp distributions lead to ",
                    sum(is.na(draws[, 9])),
                    " negative adjusted counts which were set to zero.")
            discard_mess <- c(paste("Chosen prior Se/Sp distributions lead to ",
                                    sum(is.na(draws[, 9])),
                                    " negative adjusted counts which were set to zero."))
            draws[, 9] <- ifelse(is.na(draws[, 9]), 0, draws[, 9])
        } else discard_mess <- NULL
    }

    draws[, 10] <- exp(log(draws[, 9]) -
                           qnorm(draws[, 11]) *
                               ((log(uci.obs.irr) - log(lci.obs.irr)) /
                                    (qnorm(.975) * 2)))

    corr.irr <- c(median(draws[, 9], na.rm = TRUE),
                  quantile(draws[, 9], probs = .025, na.rm = TRUE),
                  quantile(draws[, 9], probs = .975, na.rm = TRUE))
    tot.irr <- c(median(draws[, 10], na.rm = TRUE),
                 quantile(draws[, 10], probs = .025, na.rm = TRUE),
                 quantile(draws[, 10], probs = .975, na.rm = TRUE))

    if (is.null(rownames(tab)))
        rownames(tab) <- c("Cases", "Person-time")
    if (is.null(colnames(tab)))
        colnames(tab) <- c("Exposed", "Unexposed")
    rmat <- matrix(c(obs.irr, lci.obs.irr, uci.obs.irr), nrow = 1)
    rownames(rmat) <- " Observed Incidence Rate Ratio:"
    colnames(rmat) <- c(" ",
                        paste(100 * (alpha / 2), "%", sep = ""),
                        paste(100 * (1 - alpha / 2), "%", sep = ""))
    rmatc <- rbind(corr.irr, tot.irr)
    rownames(rmatc) <- c("           Incidence Rate Ratio -- systematic error:",
                         "Incidence Rate Ratio -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    res <- list(obs.data = tab,
                obs.measures = rmat,
                adj.measures = rmatc,
                sim.df = as.data.frame(draws[, -11]),
                reps = reps,
                fun = "probsens.irr",
                warnings = neg_warn,
                message = discard_mess)
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
