#' Probabilistic sensitivity analysis.
#'
#' Probabilistic sensitivity analysis to correct for exposure misclassification or
#' outcome misclassification and random error.
#' Non-differential misclassification is assumed when only the two bias parameters
#' \code{seca.parms} and \code{spca.parms} are provided. Adding the 2 parameters
#' \code{seexp.parms} and \code{spexp.parms} (i.e. providing the 4 bias parameters)
#' evaluates a differential misclassification.
#'
#' Correlations between sensitivity (or specificity) of exposure classification
#' among cases and controls can be specified and use the NORmal To Anything
#' (NORTA) transformation (Li & Hammond, 1975).
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param type Choice of correction for exposure or outcome misclassification.
#' @param reps Number of replications to run.
#' @param seca.parms List defining:
#' \enumerate{
#' \item The sensitivity of exposure classification among those with the outcome
#' (when \code{type = "exposure"}), or
#' \item The sensitivity of outcome classification among those with the exposure
#' (when \code{type = "outcome"}).
#' }
#' The first argument provides the probability distribution function (constant,
#' uniform, triangular, trapezoidal, truncated normal, or beta) and the second
#' its parameters as a vector. Lower and upper bounds of the truncated normal
#' have to be between 0 and 1.
#' \enumerate{
#' \item constant: constant value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max,
#' \item normal: lower bound, upper bound, mean, sd.
#' \item beta: alpha, beta.
#' }
#' @param seexp.parms List defining:
#' \enumerate{
#' \item The sensitivity of exposure classification among those without the outcome (when \code{type = "exposure"}), or
#' \item The sensitivity of outcome classification among those without the exposure (when \code{type = "outcome"}).
#' }
#' @param spca.parms List as above for \code{seca.parms} but for specificity.
#' @param spexp.parms List as above for \code{seexp.parms} but for specificity.
#' @param corr.se Correlation between case and non-case sensitivities.
#' @param corr.sp Correlation between case and non-case specificities.
#' @param discard A logical scalar. In case of negative adjusted count, should the draws be discarded? If set to FALSE, negative counts are set to zero.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs.measures}{A table of observed relative risk and odds ratio with confidence intervals.}
#' \item{adj.measures}{A table of corrected relative risks and odds ratios.}
#' \item{sim.df}{Data frame of random parameters and computed values.}
#' \item{reps}{Number of replications.}
#'
#' @references
#' Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#'
#' Li, S.T., Hammond, J.L., 1975. \emph{Generation of Pseudorandom Numbers
#' with Specified Univariate Distributions and Correlation Coefficients}.
#' IEEE Trans Syst Man Cybern 5:557-561.
#' @examples
#' # The data for this example come from:
#' # Greenland S., Salvan A., Wegman D.H., Hallock M.F., Smith T.J.
#' # A case-control study of cancer mortality at a transformer-assembly facility.
#' # Int Arch Occup Environ Health 1994; 66(1):49-54.
#' set.seed(123)
#' # Exposure misclassification, non-differential
#' probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 20000,
#' seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
#' spca.parms = list("trapezoidal", c(.75, .85, .95, 1)))
#'
#' # Exposure misclassification, differential
#' probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 20000,
#' seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
#' seexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
#' spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
#' spexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
#' corr.se = .8,
#' corr.sp = .8)
#'
#' probsens(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 20000,
#' seca.parms = list("beta", c(908, 16)),
#' seexp.parms = list("beta", c(156, 56)),
#' spca.parms = list("beta", c(153, 6)),
#' spexp.parms = list("beta", c(205, 18)),
#' corr.se = .8,
#' corr.sp = .8)
#'
#' probsens(matrix(c(338, 490, 17984, 32024),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 1000,
#' seca.parms = list("trapezoidal", c(.8, .9, .9, 1)),
#' spca.parms = list("trapezoidal", c(.8, .9, .9, 1)))
#'
#' # Disease misclassification
#' probsens(matrix(c(173, 602, 134, 663),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "outcome",
#' reps = 20000,
#' seca.parms = list("uniform", c(.8, 1)),
#' spca.parms = list("uniform", c(.8, 1)))
#'
#' probsens(matrix(c(338, 490, 17984, 32024),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "outcome",
#' reps = 20000,
#' seca.parms = list("uniform", c(.2, .6)),
#' seexp.parms = list("uniform", c(.1, .5)),
#' spca.parms = list("uniform", c(.99, 1)),
#' spexp.parms = list("uniform", c(.99, 1)),
#' corr.se = .8,
#' corr.sp = .8)
#'
#' probsens(matrix(c(173, 602, 134, 663),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "outcome",
#' reps = 20000,
#' seca.parms = list("beta", c(100, 5)),
#' seexp.parms = list("beta", c(110, 10)),
#' spca.parms = list("beta", c(120, 15)),
#' spexp.parms = list("beta", c(130, 30)),
#' corr.se = .8,
#' corr.sp = .8)
#' @export
#' @importFrom stats median pnorm qnorm quantile qunif runif qbeta rbeta
probsens <- function(case,
                     exposed,
                     type = c("exposure", "outcome"),
                     reps = 1000,
                     seca.parms = list(dist = c("constant", "uniform", "triangular",
                                                "trapezoidal", "normal", "beta"),
                                       parms = NULL),
                     seexp.parms = NULL,
                     spca.parms = list(dist = c("constant", "uniform", "triangular",
                                                "trapezoidal", "normal", "beta"),
                                       parms = NULL),
                     spexp.parms = NULL,
                     corr.se = NULL,
                     corr.sp = NULL,
                     discard = TRUE,
                     alpha = 0.05) {
    if (reps < 1)
        stop(paste("Invalid argument: reps =", reps))

    if (is.null(seca.parms) | is.null(spca.parms))
        stop("At least one Se and one Sp should be provided through outcome parameters.")
    if (!is.list(seca.parms))
        stop("Sensitivity of exposure classification among those with the outcome should be a list.")
    else seca.parms <- seca.parms
    if ((length(seca.parms) != 2) | (length(spca.parms) != 2))
        stop("Check distribution parameters")
    if ((!is.null(seexp.parms) & length(seexp.parms) != 2) | (!is.null(spexp.parms) & length(spexp.parms) != 2))
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
        stop("For truncated normal distribution, please provide vector of lower and upper bound limits, mean and SD")
    if (seca.parms[[1]] == "normal" & ((seca.parms[[2]][1] >= seca.parms[[2]][2]) |
                                       (!all(seca.parms[[2]][1:2] >= 0 & seca.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if ((seca.parms[[1]] == "constant" | seca.parms[[1]] == "uniform" | seca.parms[[1]] == "triangular" |
         seca.parms[[1]] == "trapezoidal") & !all(seca.parms[[2]] >= 0 & seca.parms[[2]] <= 1))
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
        ((seexp.parms[[2]][1] > seexp.parms[[2]][2]) | (seexp.parms[[2]][2] > seexp.parms[[2]][3]) |
         (seexp.parms[[2]][3] > seexp.parms[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "normal" & (length(seexp.parms[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bound limits, mean and SD.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "normal" &&
        ((seexp.parms[[2]][1] >= seexp.parms[[2]][2]) | (!all(seexp.parms[[2]][1:2] >= 0 &
                                                              seexp.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if (!is.null(seexp.parms) && (seexp.parms[[1]] == "constant" | seexp.parms[[1]] == "uniform" |
                                  seexp.parms[[1]] == "triangular" | seexp.parms[[1]] == "trapezoidal") &
        !all(seexp.parms[[2]] >= 0 & seexp.parms[[2]] <= 1))
        stop("Sensitivity of exposure classification among those without the outcome should be between 0 and 1.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "beta" && length(seexp.parms[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (!is.null(seexp.parms) && seexp.parms[[1]] == "beta" && (seexp.parms[[2]][1] < 0 | seexp.parms[[2]][2] < 0))
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
        stop("For truncated normal distribution, please provide vector of lower and upper bound limits, mean and SD.")
    if (spca.parms[[1]] == "normal" & ((spca.parms[[2]][1] >= spca.parms[[2]][2]) |
         (!all(spca.parms[[2]][1:2] >= 0 & spca.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if ((spca.parms[[1]] == "constant" | spca.parms[[1]] == "uniform" |
         spca.parms[[1]] == "triangular" |
         spca.parms[[1]] == "trapezoidal") & !all(spca.parms[[2]] >= 0 & spca.parms[[2]] <= 1))
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
        ((spexp.parms[[2]][1] > spexp.parms[[2]][2]) | (spexp.parms[[2]][2] > spexp.parms[[2]][3]) |
         (spexp.parms[[2]][3] > spexp.parms[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "normal" & (length(spexp.parms[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bound limits, meand and SD.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "normal" &&
        ((spexp.parms[[2]][1] >= spexp.parms[[2]][2]) |
         (!all(spexp.parms[[2]][1:2] >= 0 & spexp.parms[[2]][1:2] <= 1))))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).")
    if (!is.null(spexp.parms) && (spexp.parms[[1]] == "constant" | spexp.parms[[1]] == "uniform" |
                                  spexp.parms[[1]] == "triangular" | spexp.parms[[1]] == "trapezoidal") &
        !all(spexp.parms[[2]] >= 0 & spexp.parms[[2]] <= 1))
        stop("Specificity of exposure classification among those without the outcome should be between 0 and 1.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "beta" && length(spexp.parms[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (!is.null(spexp.parms) && spexp.parms[[1]] == "beta" && (spexp.parms[[2]][1] < 0 | spexp.parms[[2]][2] < 0))
        stop("Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")

    if (!is.null(seexp.parms) & (is.null(spca.parms) | is.null(spexp.parms) | is.null(corr.se) | is.null(corr.sp)))
        stop("For differential misclassification type, have to provide Se and Sp for among those with and without the outcome as well as Se and Sp correlations.")

    if (!is.null(corr.se) && (corr.se == 0 | corr.se == 1))
        stop("Correlations should be > 0 and < 1.")
    if (!is.null(corr.sp) && (corr.sp == 0 | corr.sp == 1))
        stop("Correlations should be > 0 and < 1.")

    if (!inherits(case, "episensr.probsens")) {
        if (inherits(case, c("table", "matrix")))
            tab <- case
        else {
            tab.df <- table(case, exposed)
            tab <- tab.df[2:1, 2:1]
        }

        a <- as.numeric(tab[1, 1])
        b <- as.numeric(tab[1, 2])
        c <- as.numeric(tab[2, 1])
        d <- as.numeric(tab[2, 2])
    } else {
        a <- as.numeric(case[[3]][, 1])
        b <- as.numeric(case[[3]][, 2])
        c <- as.numeric(case[[3]][, 3])
        d <- as.numeric(case[[3]][, 4])

        reps <- case[[4]]
    }

    obs.rr <- (a / (a + c)) / (b / (b + d))
    se.log.obs.rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
    lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha / 2) * se.log.obs.rr)
    uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha / 2) * se.log.obs.rr)

    obs.or <- (a / b) / (c / d)
    se.log.obs.or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha / 2) * se.log.obs.or)
    uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha / 2) * se.log.obs.or)

    draws <- matrix(NA, nrow = reps, ncol = 13)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A1", "B1", "C1", "D1",
                         "corr.RR", "corr.OR",
                         "reps",
                         "tot.RR", "tot.OR")
    corr.draws <- matrix(NA, nrow = reps, ncol = 4)

    seca <- c(reps, seca.parms[[2]])
    seexp <- c(reps, seexp.parms[[2]])
    spca <- c(reps, spca.parms[[2]])
    spexp <- c(reps, spexp.parms[[2]])

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

    type <- match.arg(type)
    if (type == "exposure") {
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) / (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (a + b) - draws[, 5]
        draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) / (draws[, 2] - (1 - draws[, 4]))
        draws[, 8] <- (c + d) - draws[, 7]

        draws[, 9] <- (draws[, 5] / (draws[, 5] + draws[, 7])) / (draws[, 6] / (draws[, 6] + draws[, 8]))
        draws[, 10] <- (draws[, 5] / draws[, 7]) / (draws[, 6] / draws[, 8])

        draws[, 9] <- ifelse(draws[, 5] < 0 |
                             draws[, 6] < 0 |
                             draws[, 7] < 0 |
                             draws[, 8] < 0, NA, draws[, 9])
        draws[, 10] <- ifelse(draws[, 5] < 0 |
                              draws[, 6] < 0 |
                              draws[, 7] < 0 |
                              draws[, 8] < 0, NA, draws[, 10])
        if (all(is.na(draws[, 9])) | all(is.na(draws[, 10]))) {
            warning("Prior Se/Sp distributions lead to all negative adjusted counts.")
            neg_warn <- "Prior Se/Sp distributions lead to all negative adjusted counts."
        } else {
            neg_warn <- NULL
        }
        if (discard) {
            if (sum(is.na(draws[, 9])) > 0) {
                message("Chosen prior Se/Sp distributions lead to ",
                        sum(is.na(draws[, 9])),
                        " negative adjusted counts which were discarded.")
                discard_mess <- c(paste("Chosen prior Se/Sp distributions lead to ",
                                        sum(is.na(draws[, 9])),
                                        " negative adjusted counts which were discarded."))
            } else discard_mess <- NULL
        }
        else {
            if (sum(is.na(draws[, 9])) > 0) {
                message("Chosen prior Se/Sp distributions lead to ",
                        sum(is.na(draws[, 9])),
                        " negative adjusted counts which were set to zero.")
                discard_mess <- c(paste("Chosen prior Se/Sp distributions lead to ",
                                        sum(is.na(draws[, 9])),
                                        " negative adjusted counts which were set to zero."))
                draws[, 9] <- ifelse(is.na(draws[, 9]), 0, draws[, 9])
                draws[, 10] <- ifelse(is.na(draws[, 10]), 0, draws[, 10])
            } else discard_mess <- NULL
        }

        draws[, 12] <- exp(log(draws[, 9]) - qnorm(draws[, 11]) * ((log(uci.obs.rr) - log(lci.obs.rr)) /
                                                                   (qnorm(.975) * 2)))
        draws[, 13] <- exp(log(draws[, 10]) - qnorm(draws[, 11]) * ((log(uci.obs.or) - log(lci.obs.or)) /
                                                                    (qnorm(.975) * 2)))

        corr.rr <- c(median(draws[, 9], na.rm = TRUE),
                     quantile(draws[, 9], probs = .025, na.rm = TRUE),
                     quantile(draws[, 9], probs = .975, na.rm = TRUE))
        corr.or <- c(median(draws[, 10], na.rm = TRUE),
                     quantile(draws[, 10], probs = .025, na.rm = TRUE),
                     quantile(draws[, 10], probs = .975, na.rm = TRUE))
        tot.rr <- c(median(draws[, 12], na.rm = TRUE),
                    quantile(draws[, 12], probs = .025, na.rm = TRUE),
                    quantile(draws[, 12], probs = .975, na.rm = TRUE))
        tot.or <- c(median(draws[, 13], na.rm = TRUE),
                    quantile(draws[, 13], probs = .025, na.rm = TRUE),
                    quantile(draws[, 13], probs = .975, na.rm = TRUE))

        if (!inherits(case, "episensr.probsens")) {
            tab <- tab
            rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr),
                          c(obs.or, lci.obs.or, uci.obs.or))
            rownames(rmat) <- c(" Observed Relative Risk:",
                                "    Observed Odds Ratio:")
            colnames(rmat) <- c(" ",
                                paste(100 * (alpha / 2), "%", sep = ""),
                                paste(100 * (1 - alpha / 2), "%", sep = ""))
        } else {
            tab <- case[[1]]
            rmat <- case[[2]]
        }
        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        rmatc <- rbind(corr.rr, corr.or, tot.rr, tot.or)
        rownames(rmatc) <- c("           Relative Risk -- systematic error:",
                             "              Odds Ratio -- systematic error:",
                             "Relative Risk -- systematic and random error:",
                             "   Odds Ratio -- systematic and random error:")
        colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    }

    if (type == "outcome") {
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + c)) / (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (b - (1 - draws[, 4]) * (b + d)) / (draws[, 2] - (1 - draws[, 4]))
        draws[, 7] <- (a + c) - draws[, 5]
        draws[, 8] <- (b + d) - draws[, 6]

        draws[, 9] <- (draws[, 5] / (draws[, 5] + draws[, 7])) / (draws[, 6] / (draws[, 6] + draws[, 8]))
        draws[, 10] <- (draws[, 5] / draws[, 7]) / (draws[, 6] / draws[, 8])

        draws[, 9] <- ifelse(draws[, 5] < 0 |
                             draws[, 6] < 0 |
                             draws[, 7] < 0 |
                             draws[, 8] < 0, NA, draws[, 9])
        draws[, 10] <- ifelse(draws[, 5] < 0 |
                              draws[, 6] < 0 |
                              draws[, 7] < 0 |
                              draws[, 8] < 0, NA, draws[, 10])

        if (all(is.na(draws[, 9])) | all(is.na(draws[, 10]))) {
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
                draws[, 10] <- ifelse(is.na(draws[, 10]), 0, draws[, 10])
            } else discard_mess <- NULL
        }

        draws[, 12] <- exp(log(draws[, 9]) -
                               qnorm(draws[, 11]) *
                                         ((log(uci.obs.rr) - log(lci.obs.rr)) /
                                              (qnorm(.975) * 2)))
        draws[, 13] <- exp(log(draws[, 10]) -
                               qnorm(draws[, 11]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))

        corr.rr <- c(median(draws[, 9], na.rm = TRUE),
                     quantile(draws[, 9], probs = .025, na.rm = TRUE),
                     quantile(draws[, 9], probs = .975, na.rm = TRUE))
        corr.or <- c(median(draws[, 10], na.rm = TRUE),
                     quantile(draws[, 10], probs = .025, na.rm = TRUE),
                     quantile(draws[, 10], probs = .975, na.rm = TRUE))
        tot.rr <- c(median(draws[, 12], na.rm = TRUE),
                    quantile(draws[, 12], probs = .025, na.rm = TRUE),
                    quantile(draws[, 12], probs = .975, na.rm = TRUE))
        tot.or <- c(median(draws[, 13], na.rm = TRUE),
                    quantile(draws[, 13], probs = .025, na.rm = TRUE),
                    quantile(draws[, 13], probs = .975, na.rm = TRUE))

        if(!inherits(case, "episensr.probsens")) {
            tab <- tab
            rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr),
                          c(obs.or, lci.obs.or, uci.obs.or))
            rownames(rmat) <- c(" Observed Relative Risk:",
                                "    Observed Odds Ratio:")
            colnames(rmat) <- c(" ",
                                paste(100 * (alpha / 2), "%", sep = ""),
                                paste(100 * (1 - alpha / 2), "%", sep = ""))
        } else {
            tab <- case[[1]]
            rmat <- case[[2]]
        }
        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        rmatc <- rbind(corr.rr, corr.or, tot.rr, tot.or)
        rownames(rmatc) <- c("           Relative Risk -- systematic error:",
                             "              Odds Ratio -- systematic error:",
                             "Relative Risk -- systematic and random error:",
                             "   Odds Ratio -- systematic and random error:")
        colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    }
    res <- list(obs.data = tab,
                obs.measures = rmat,
                adj.measures = rmatc,
                sim.df = as.data.frame(draws[, -11]),
                reps = reps,
                fun = "probsens",
                warnings = neg_warn,
                message = discard_mess)
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
