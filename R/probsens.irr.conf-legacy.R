#' Legacy version of `probsens.irr.conf()`.
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' episensr 2.0.0 introduced breaking changes in probabilistic bias analyses by
#' (1) using the NORTA transformation to define a correlation between
#' distributions, and (2) sampling true prevalences and then sampling the
#' adjusted cell counts rather than just using the expected cell counts from a
#' simple quantitative bias analysis. This updated version should be preferred
#' and this legacy version will be deprecated in future versions. However, if
#' you need to quickly roll back to the previous calculations, this function
#' provides the previous interface. To make old code work as is, add the
#' following code to the top of your script:
#'
#' ```
#' library(episensr)
#' probsens.irr.conf <- probsens.irr.conf_legacy
#' ```
#'
#' @param counts A table or matrix where first row contains disease counts and second row contains person-time at risk, and first and second columns are exposed and unexposed observations, as:
#' \tabular{lll}{
#' \tab Exposed \tab Unexposed \cr
#' Cases \tab a \tab b \cr
#' Person-time \tab N1 \tab N0
#' }
#' @param pt A numeric vector of person-time at risk. If provided, \code{counts} must be a numeric vector of disease counts.
#' @param reps Number of replications to run.
#' @param prev.exp List defining the prevalence of exposure among the exposed. The first argument provides the probability distribution function (constant,uniform, triangular, trapezoidal, logit-logistic, logit-normal, or beta) and the second its parameters as a vector. Logit-logistic and logit-normal distributions can be shifted by providing lower and upper bounds. Avoid providing these values if a non-shifted distribution is desired.
#' \enumerate{
#' \item constant; value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item logit-logistic: location, scale, lower bound shift, upper bound shift,
#' \item logit-normal: location, scale, lower bound shift, upper bound shift,
#' \item beta: alpha, beta.
#' }
#' @param prev.nexp List defining the prevalence of exposure among the unexposed.
#' @param risk List defining the confounder-disease relative risk or the confounder-exposure odds ratio. The first argument provides the probability distribution function (constant,uniform, triangular, trapezoidal, log-logistic, or log-normal) and the second its parameters as a vector:
#' \enumerate{
#' \item constant: value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item log-logistic: shape, rate. Must be strictly positive,
#' \item log-normal: meanlog, sdlog. This is the mean and standard deviation on the log scale.
#' }
#' @param corr.p Correlation between the exposure-specific confounder prevalences.
#' @param discard A logical scalar. In case of negative adjusted count, should the draws be discarded? If set to FALSE, negative counts are set to zero.
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
#' \dontrun{
#' set.seed(123)
#' # Unmeasured confounding
#' probsens.irr.conf(matrix(c(77, 10000, 87, 10000),
#' dimnames = list(c("D+", "Person-time"), c("E+", "E-")), ncol = 2),
#' reps = 20000,
#' prev.exp = list("trapezoidal", c(.01, .2, .3, .51)),
#' prev.nexp = list("trapezoidal", c(.09, .27, .35, .59)),
#' risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
#' corr.p = .8)
#' }
#' @export
#' @importFrom stats binom.test median quantile runif rbeta qbeta
#' @rdname probsens.irr.conf_legacy
probsens.irr.conf_legacy <- function(counts,
                                     pt = NULL,
                                     reps = 1000,
                                     prev.exp = list(dist = c("constant", "uniform",
                                                              "triangular", "trapezoidal",
                                                              "logit-logistic", "logit-normal", "beta"),
                                                     parms = NULL),
                                     prev.nexp = list(dist = c("constant", "uniform",
                                                               "triangular", "trapezoidal",
                                                               "logit-logistic", "logit-normal", "beta"),
                                                      parms = NULL),
                                     risk = list(dist = c("constant", "uniform", "triangular",
                                                          "trapezoidal", "log-logistic",
                                                          "log-normal"),
                                                 parms = NULL),
                                     corr.p = NULL,
                                     discard = TRUE,
                                     alpha = 0.05){
    if(reps < 1)
        stop(paste("Invalid argument: reps = ", reps))

    if(is.null(prev.exp) | is.null(prev.nexp))
        stop('Please provide prevalences among the exposed and unexposed.')
    if(is.null(risk))
        stop('Please provide risk of acquiring outcome.')

    if(!is.list(prev.exp))
        stop('Prevalence of exposure among the exposed should be a list.')
    else prev.exp <- prev.exp
    if((length(prev.exp) != 2) | (length(prev.nexp) != 2) | (length(risk) != 2))
        stop('Check distribution parameters.')
    if((length(prev.exp[[1]]) != 1) | (length(prev.nexp[[1]]) != 1) |
       (length(risk[[1]]) != 1))
        stop('Which distribution?')
    if(!is.null(corr.p) && (prev.exp[[1]] == "constant" | prev.nexp[[1]] == "constant"))
        stop('No correlated distributions with constant values.')
    if(prev.exp[[1]] == "constant" & length(prev.exp[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(prev.exp[[1]] == "uniform" & length(prev.exp[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(prev.exp[[1]] == "uniform" & prev.exp[[2]][1] >= prev.exp[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(prev.exp[[1]] == "triangular" & length(prev.exp[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(prev.exp[[1]] == "triangular" & ((prev.exp[[2]][1] > prev.exp[[2]][3]) |
                                        (prev.exp[[2]][2] < prev.exp[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(prev.exp[[1]] == "trapezoidal" & length(prev.exp[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(prev.exp[[1]] == "trapezoidal" & ((prev.exp[[2]][1] > prev.exp[[2]][2]) |
                                         (prev.exp[[2]][2] > prev.exp[[2]][3]) |
                                         (prev.exp[[2]][3] > prev.exp[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')
    if(prev.exp[[1]] == "logit-logistic" & (length(prev.exp[[2]]) < 2 | length(prev.exp[[2]]) == 3 | length(prev.exp[[2]]) > 4))
        stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(prev.exp[[1]] == "logit-logistic" & length(prev.exp[[2]]) == 4 &
       ((prev.exp[[2]][3] >= prev.exp[[2]][4]) | (!all(prev.exp[[2]][3:4] >= 0 & prev.exp[[2]][3:4] <= 1))))
        stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(prev.exp[[1]] == "logit-logistic" & length(prev.exp[[2]]) == 2)
        prev.exp <- list(prev.exp[[1]], c(prev.exp[[2]], c(0, 1)))
    if(prev.exp[[1]] == "logit-normal" & (length(prev.exp[[2]]) < 2 | length(prev.exp[[2]]) == 3 | length(prev.exp[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(prev.exp[[1]] == "logit-normal" & length(prev.exp[[2]]) == 4 &
       ((prev.exp[[2]][3] >= prev.exp[[2]][4]) | (!all(prev.exp[[2]][3:4] >= 0 & prev.exp[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(prev.exp[[1]] == "logit-normal" & length(prev.exp[[2]]) == 2)
        prev.exp <- list(prev.exp[[1]], c(prev.exp[[2]], c(0, 1)))
    if((prev.exp[[1]] == "constant" | prev.exp[[1]] == "uniform" | prev.exp[[1]] == "triangular" | prev.exp[[1]] == "trapezoidal") & !all(prev.exp[[2]] >= 0 & prev.exp[[2]] <= 1))
        stop('Prevalence should be between 0 and 1.')
    if(!is.null(prev.exp) && prev.exp[[1]] == "beta" && length(prev.exp[[2]]) != 2)
        stop('For beta distribution, please provide alpha and beta.')
    if(!is.null(prev.exp) && prev.exp[[1]] == "beta" &&
       (prev.exp[[2]][1] < 0 | prev.exp[[2]][2] < 0))
        stop('Wrong arguments for your beta distribution. Alpha and Beta should be > 0.')

    if(!is.list(prev.nexp))
        stop('Prevalence of exposure among the non-exposed should be a list.')
    else prev.nexp <- prev.nexp
    if(prev.nexp[[1]] == "constant" & length(prev.nexp[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(prev.nexp[[1]] == "uniform" & length(prev.nexp[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(prev.nexp[[1]] == "uniform" & prev.nexp[[2]][1] >= prev.nexp[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(prev.nexp[[1]] == "triangular" & length(prev.nexp[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(prev.nexp[[1]] == "triangular" & ((prev.nexp[[2]][1] > prev.nexp[[2]][3]) |
                                        (prev.nexp[[2]][2] < prev.nexp[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(prev.nexp[[1]] == "trapezoidal" & length(prev.nexp[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(prev.nexp[[1]] == "trapezoidal" & ((prev.nexp[[2]][1] > prev.nexp[[2]][2]) |
                                         (prev.nexp[[2]][2] > prev.nexp[[2]][3]) |
                                         (prev.nexp[[2]][3] > prev.nexp[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')
    if(prev.nexp[[1]] == "logit-logistic" & (length(prev.nexp[[2]]) < 2 | length(prev.nexp[[2]]) == 3 | length(prev.nexp[[2]]) > 4))
        stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(prev.nexp[[1]] == "logit-logistic" & length(prev.nexp[[2]]) == 4 &
       ((prev.nexp[[2]][3] >= prev.nexp[[2]][4]) | (!all(prev.nexp[[2]][3:4] >= 0 & prev.nexp[[2]][3:4] <= 1))))
        stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(prev.nexp[[1]] == "logit-logistic" & length(prev.nexp[[2]]) == 2)
        prev.nexp <- list(prev.nexp[[1]], c(prev.nexp[[2]], c(0, 1)))
    if(prev.nexp[[1]] == "logit-normal" & (length(prev.nexp[[2]]) < 2 | length(prev.nexp[[2]]) == 3 | length(prev.nexp[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(prev.nexp[[1]] == "logit-normal" & length(prev.nexp[[2]]) == 4 &
       ((prev.nexp[[2]][3] >= prev.nexp[[2]][4]) | (!all(prev.nexp[[2]][3:4] >= 0 & prev.nexp[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(prev.nexp[[1]] == "logit-normal" & length(prev.nexp[[2]]) == 2)
        prev.nexp <- list(prev.nexp[[1]], c(prev.nexp[[2]], c(0, 1)))
    if((prev.nexp[[1]] == "constant" | prev.nexp[[1]] == "uniform" | prev.nexp[[1]] == "triangular" | prev.nexp[[1]] == "trapezoidal") & !all(prev.nexp[[2]] >= 0 & prev.nexp[[2]] <= 1))
        stop('Prevalence should be between 0 and 1.')
    if(!is.null(prev.nexp) && prev.nexp[[1]] == "beta" && length(prev.nexp[[2]]) != 2)
        stop('For beta distribution, please provide alpha and beta.')
    if(!is.null(prev.nexp) && prev.nexp[[1]] == "beta" &&
       (prev.nexp[[2]][1] < 0 | prev.nexp[[2]][2] < 0))
        stop('Wrong arguments for your beta distribution. Alpha and Beta should be > 0.')

    if(!is.list(risk))
        stop('Risk should be a list.')
    else risk <- risk
    if(risk[[1]] == "constant" & length(risk[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(risk[[1]] == "uniform" & length(risk[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(risk[[1]] == "uniform" & risk[[2]][1] >= risk[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(risk[[1]] == "triangular" & length(risk[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(risk[[1]] == "triangular" & ((risk[[2]][1] > risk[[2]][3]) |
                                        (risk[[2]][2] < risk[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(risk[[1]] == "trapezoidal" & length(risk[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(risk[[1]] == "trapezoidal" & ((risk[[2]][1] > risk[[2]][2]) |
                                         (risk[[2]][2] > risk[[2]][3]) |
                                         (risk[[2]][3] > risk[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')
    if(risk[[1]] == "log-logistic" & length(risk[[2]]) != 2)
        stop('For log-logistic distribution, please provide vector of location and scale.')
    if(risk[[1]] == "log-normal" & length(risk[[2]]) != 2)
        stop('For log-logistic distribution, please provide vector of meanlog and sdlog.')

    if(!is.null(corr.p) && (corr.p == 0 | corr.p == 1))
        stop('Correlations should be > 0 and < 1.')

    if(!is.null(pt) && inherits(counts, c("table", "matrix")))
        stop("pt argument should be NULL.")
    if(!inherits(counts, c("vector", "table", "matrix")))
        stop("counts argument should be a vector, a table, or a matrix.")
    if(is.null(pt) && inherits(counts, c("table", "matrix")))
        tab <- counts
    else tab <- rbind(counts, pt)
    a <- as.numeric(tab[1, 1])
    b <- as.numeric(tab[1, 2])
    c <- as.numeric(tab[2, 1])
    d <- as.numeric(tab[2, 2])

    draws <- matrix(NA, nrow = reps, ncol = 14)
    colnames(draws) <- c("p1", "p0", "RR.cd",
                         "M1", "N1", "A1", "B1",
                         "M0", "N0", "A0", "B0",
                         "corr.IRR", "reps", "tot.IRR")
    corr.draws <- matrix(NA, nrow = reps, ncol = 5)

    p1 <- c(reps, prev.exp[[2]])
    p0 <- c(reps, prev.nexp[[2]])
    rr.cd <- c(reps, risk[[2]])

    obs.irr <- (a / c) / (b / d)
    lci.obs.irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[1]) * c)
    uci.obs.irr <- (binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2] * d) /
        ((1 - binom.test(a, a + b, conf.level = 1 - alpha)$conf.int[2]) * c)

    if (is.null(corr.p)) {
        if (prev.exp[[1]] == "constant") {
            draws[, 1] <- prev.exp[[2]]
        }
        if (prev.exp[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(p1))
            }
        if (prev.exp[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(p1))
            }
        if (prev.exp[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(p1))
            }
        if (prev.exp[[1]] == "logit-logistic") {
            draws[, 1] <- logitlog.dstr(p1)
            }
        if (prev.exp[[1]] == "logit-normal") {
            draws[, 1] <- logitnorm.dstr(p1)
            }
        if (prev.exp[[1]] == "beta") {
            draws[, 1] <- do.call(rbeta, as.list(p1))
            }
        if (prev.nexp[[1]] == "constant") {
            draws[, 2] <- prev.nexp[[2]]
        }
        if (prev.nexp[[1]] == "uniform") {
            draws[, 2] <- do.call(runif, as.list(p0))
            }
        if (prev.nexp[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::rtriangle, as.list(p0))
            }
        if (prev.nexp[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::rtrapezoid, as.list(p0))
            }
        if (prev.nexp[[1]] == "logit-logistic") {
            draws[, 2] <- logitlog.dstr(p0)
            }
        if (prev.nexp[[1]] == "logit-normal") {
            draws[, 2] <- logitnorm.dstr(p0)
            }
        if (prev.nexp[[1]] == "beta") {
            draws[, 2] <- do.call(rbeta, as.list(p0))
           }
    } else {
        corr.draws[, 1:3] <- apply(corr.draws[, 1:3],
                                   2,
                                   function(x) x = runif(reps))
        corr.draws[, 1:3] <- apply(corr.draws[, 1:3],
                                   2,
                                   function(x) log(x / (1 - x)))
        corr.draws[, 4] <- exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 2]) /
            (1 + (exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 2])))
        corr.draws[, 5] <- exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 3]) /
            (1 + (exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 3])))

    if (prev.exp[[1]] == "uniform") {
        draws[, 1] <- prev.exp[[2]][2] -
            (prev.exp[[2]][2] - prev.exp[[2]][1]) * corr.draws[, 4]
    }
    if (prev.exp[[1]] == "triangular") {
        draws[, 1] <- (corr.draws[, 4] *
            (prev.exp[[2]][2] - prev.exp[[2]][1]) + (prev.exp[[2]][1] + prev.exp[[2]][3])) / 2
        draws[, 1] <- ifelse(draws[, 1] < prev.exp[[2]][3],
                             prev.exp[[2]][1] + sqrt(abs((prev.exp[[2]][3] - prev.exp[[2]][1]) * (2 * draws[, 1] - prev.exp[[2]][1] - prev.exp[[2]][3]))),
                             draws[, 1])
        draws[, 1] <- ifelse(draws[, 1] > prev.exp[[2]][3],
                             prev.exp[[2]][2] - sqrt(abs(2 * (prev.exp[[2]][2] - prev.exp[[2]][3]) * (draws[, 1] - prev.exp[[2]][3]))),
                             draws[, 1])
    }
    if (prev.exp[[1]] == "trapezoidal") {
        draws[, 1] <- (corr.draws[, 4] *
            (prev.exp[[2]][4] + prev.exp[[2]][3] - prev.exp[[2]][1] - prev.exp[[2]][2]) + (prev.exp[[2]][1] + prev.exp[[2]][2])) / 2
        draws[, 1] <- ifelse(draws[, 1] < prev.exp[[2]][2],
                             prev.exp[[2]][1] + sqrt(abs((prev.exp[[2]][2] - prev.exp[[2]][1]) * (2 * draws[, 1] - prev.exp[[2]][1] - prev.exp[[2]][2]))),
                             draws[, 1])
        draws[, 1] <- ifelse(draws[, 1] > prev.exp[[2]][3],
                             prev.exp[[2]][4] - sqrt(abs(2 * (prev.exp[[2]][4] - prev.exp[[2]][3]) * (draws[, 1] - prev.exp[[2]][3]))),
                             draws[, 1])
    }
    if (prev.exp[[1]] == "logit-logistic") {
        pexp.w <- prev.exp[[2]][1] + (prev.exp[[2]][2] * log(corr.draws[, 4] / (1 - corr.draws[, 4])))
        draws[, 1] <- prev.exp[[2]][3] + (prev.exp[[2]][4] - prev.exp[[2]][3]) * exp(pexp.w) / (1 + exp(pexp.w))
    }
    if (prev.exp[[1]] == "logit-normal") {
        pexp.w <- prev.exp[[2]][1] + (prev.exp[[2]][2] * qnorm(corr.draws[, 4]))
        draws[, 1] <- prev.exp[[2]][3] + (prev.exp[[2]][4] - prev.exp[[2]][3]) * exp(pexp.w) / (1 + exp(pexp.w))
    }
    if (prev.exp[[1]] == "beta") {
        draws[, 1] <- qbeta(corr.draws[, 4]/(1 + corr.draws[, 4]),
                            prev.exp[[2]][1],
                            prev.exp[[2]][2])
    }

    if (prev.nexp[[1]] == "uniform") {
        draws[, 2] <- prev.nexp[[2]][2] -
            (prev.nexp[[2]][2] - prev.nexp[[2]][1]) * corr.draws[, 5]
    }
    if (prev.nexp[[1]] == "triangular") {
        draws[, 2] <- (corr.draws[, 5] *
                           (prev.nexp[[2]][2] - prev.nexp[[2]][1]) + (prev.nexp[[2]][1] + prev.nexp[[2]][3])) / 2
        draws[, 2] <- ifelse(draws[, 2] < prev.nexp[[2]][3],
                             prev.nexp[[2]][1] + sqrt(abs((prev.nexp[[2]][3] - prev.nexp[[2]][1]) * (2 * draws[, 2] - prev.nexp[[2]][1] - prev.nexp[[2]][3]))),
                             draws[, 2])
        draws[, 2] <- ifelse(draws[, 2] > prev.nexp[[2]][3],
                             prev.nexp[[2]][2] - sqrt(abs(2 * (prev.nexp[[2]][2] - prev.nexp[[2]][3]) * (draws[, 2] - prev.nexp[[2]][3]))),
                             draws[, 2])
    }
    if (prev.nexp[[1]] == "trapezoidal") {
        draws[, 2] <- (corr.draws[, 5] *
                           (prev.nexp[[2]][4] + prev.nexp[[2]][3] - prev.nexp[[2]][1] - prev.nexp[[2]][2]) + (prev.nexp[[2]][1] + prev.nexp[[2]][2])) / 2
        draws[, 2] <- ifelse(draws[, 2] < prev.nexp[[2]][2],
                             prev.nexp[[2]][1] + sqrt(abs((prev.nexp[[2]][2] - prev.nexp[[2]][1]) * (2 * draws[, 2] - prev.nexp[[2]][1] - prev.nexp[[2]][2]))),
                             draws[, 2])
        draws[, 2] <- ifelse(draws[, 2] > prev.nexp[[2]][3],
                             prev.nexp[[2]][4] - sqrt(abs(2 * (prev.nexp[[2]][4] - prev.nexp[[2]][3]) * (draws[, 2] - prev.nexp[[2]][3]))),
                             draws[, 2])
    }
    if (prev.nexp[[1]] == "logit-logistic") {
        punexp.w <- prev.nexp[[2]][1] + (prev.nexp[[2]][2] * log(corr.draws[, 5] / (1 - corr.draws[, 5])))
        draws[, 2] <- prev.nexp[[2]][3] + (prev.nexp[[2]][4] - prev.nexp[[2]][3]) * exp(punexp.w) / (1 + exp(punexp.w))
    }
    if (prev.nexp[[1]] == "logit-normal") {
        punexp.w <- prev.nexp[[2]][1] + (prev.nexp[[2]][2] * qnorm(corr.draws[, 5]))
        draws[, 2] <- prev.nexp[[2]][3] + (prev.nexp[[2]][4] - prev.nexp[[2]][3]) * exp(punexp.w) / (1 + exp(punexp.w))
    }
    if (prev.nexp[[1]] == "beta") {
        draws[, 2] <- qbeta(corr.draws[, 5]/(1 + corr.draws[, 5]),
                            prev.nexp[[2]][1],
                            prev.nexp[[2]][2])
    }
    }

    if(risk[[1]] == "constant") {
        draws[, 3] <- risk[[2]]
    }
    if (risk[[1]] == "uniform") {
        draws[, 3] <- do.call(runif, as.list(rr.cd))
    }
    if (risk[[1]] == "triangular") {
        draws[, 3] <- do.call(triangle::rtriangle, as.list(rr.cd))
    }
    if (risk[[1]] == "trapezoidal") {
        draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(rr.cd))
    }
    if (risk[[1]] == "log-logistic") {
        draws[, 3] <- do.call(actuar::rllogis, as.list(rr.cd))
    }
    if (risk[[1]] == "log-normal") {
        draws[, 3] <- do.call(rlnorm, as.list(rr.cd))
    }

    draws[, 13] <- runif(reps)

    draws[, 4] <- c * draws[, 1]
    draws[, 6] <- (draws[, 3] * draws[, 4] * a) /
        ((draws[, 3] * draws[, 6]) + c - draws[, 6])
    draws[, 5] <- d * draws[, 2]
    draws[, 7] <- (draws[, 3] * draws[, 5] * b) /
        ((draws[, 3] * draws[, 5]) + d - draws[, 5])
    draws[, 10] <- a - draws[, 6]
    draws[, 8] <- c - draws[, 4]
    draws[, 11] <- b - draws[, 7]
    draws[, 9] <- d - draws[, 5]

    draws[, 12] <- a /
        ((draws[, 4] * draws[, 7] / draws[, 5]) +
             (draws[, 8] * draws[, 11] / draws[, 9]))

    draws[, 12] <- ifelse(draws[, 4] < 0 |
                             draws[, 5] < 0 |
                                 draws[, 7] < 0 |
                                     draws[, 8] < 0 |
                                         draws[, 11] < 0, NA, draws[, 12])

    if(all(is.na(draws[, 12]))) {
        warning('Prior prevalence distributions lead to all negative adjusted values.')
        neg_warn <- "Prior Se/Sp distributions lead to all negative adjusted counts."
    } else neg_warn <- NULL
    if (discard) {
        if(sum(is.na(draws[, 12])) > 0) {
            message('Chosen prior prevalence distributions lead to ',
                    sum(is.na(draws[, 12])),
                    ' negative adjusted values which were discarded.')
            discard_mess <- c(paste('Chosen prior Se/Sp distributions lead to ',
                                    sum(is.na(draws[, 12])),
                                    ' negative adjusted counts which were discarded.'))
        } else discard_mess <- NULL
    }
    else {
        if(sum(is.na(draws[, 12])) > 0) {
            message('Chosen prior Se/Sp distributions lead to ',
                    sum(is.na(draws[, 12])),
                    ' negative adjusted counts which were set to zero.')
                discard_mess <- c(paste('Chosen prior Se/Sp distributions lead to ',
                                        sum(is.na(draws[, 12])),
                                        ' negative adjusted counts which were set to zero.'))
            draws[, 12] <- ifelse(is.na(draws[, 12]), 0, draws[, 12])
        } else discard_mess <- NULL
    }

    draws[, 14] <- exp(log(draws[, 12]) -
                           qnorm(draws[, 13]) *
                               ((log(uci.obs.irr) - log(lci.obs.irr)) /
                                    (qnorm(.975) * 2)))

    corr.irr <- c(median(draws[, 12], na.rm = TRUE),
                 quantile(draws[, 12], probs = .025, na.rm = TRUE),
                 quantile(draws[, 12], probs = .975, na.rm = TRUE))
    tot.irr <- c(median(draws[, 14], na.rm = TRUE),
                quantile(draws[, 14], probs = .025, na.rm = TRUE),
                quantile(draws[, 14], probs = .975, na.rm = TRUE))

    if (is.null(rownames(tab)))
        rownames(tab) <- c("Cases", "Person-time")
    if (is.null(colnames(tab)))
        colnames(tab) <- c("Exposed", "Unexposed")
    rmat <- matrix(c(obs.irr, lci.obs.irr, uci.obs.irr), nrow = 1)
    rownames(rmat) <- " Observed Incidence Rate ratio:"
    colnames(rmat) <- c(" ",
                        paste(100 * (alpha/2), "%", sep = ""),
                        paste(100 * (1 - alpha/2), "%", sep = ""))
    rmatc <- rbind(corr.irr, tot.irr)
    rownames(rmatc) <- c("           Incidence Rate Ratio -- systematic error:",
                         "Incidence Rate Ratio -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    res <- list(obs.data = tab,
                obs.measures = rmat,
                adj.measures = rmatc,
                sim.df = as.data.frame(draws[, -13]),
                reps = reps,
                fun = "probsens.irr.conf",
                warnings = neg_warn,
                message = discard_mess)
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
