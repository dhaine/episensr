#' Probabilistic sensitivity analysis for unmeasured confounding.
#'
#' Probabilistic sensitivity analysis to correct for unknown or unmeasured confounding and random error simultaneously.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param reps Number of replications to run.
#' @param prev.exp List defining the prevalence of exposure among the exposed. The first argument provides the probability distribution function (constant, uniform, triangular, trapezoidal, logit-logistic, or logit-normal) and the second its parameters as a vector. Logit-logistic and logit-normal distributions can be shifted by providing lower and upper bounds. Avoid providing these values if a non-shifted distribution is desired.
#' \enumerate{
#' \item Constant: constant value,
#' \item Uniform: min, max,
#' \item Triangular: lower limit, upper limit, mode,
#' \item Trapezoidal: min, lower mode, upper mode, max.
#' \item Logit-logistic: location, scale, lower bound shift, upper bound shift,
#' \item Logit-normal: location, scale, lower bound shift, upper bound shift.
#' }
#' @param prev.nexp List defining the prevalence of exposure among the unexposed.
#' @param risk List defining the confounder-disease relative risk or the confounder-exposure odds ratio. The first argument provides the probability distribution function (constant, uniform, triangular, trapezoidal, log-logistic, or log-normal) and the second its parameters as a vector:
#' \enumerate{
#' \item Constant: constant value,
#' \item Uniform: min, max,
#' \item Triangular: lower limit, upper limit, mode,
#' \item Trapezoidal: min, lower mode, upper mode, max.
#' \item Log-logistic: shape, rate. Must be strictly positive,
#' \item Log-normal: meanlog, sdlog. This is the mean and standard deviation on the log scale.
#' }
#' @param corr.p Correlation between the exposure-specific confounder prevalences.
#' @param discard A logical scalar. In case of negative adjusted count, should the draws be discarded? If set to FALSE, negative counts are set to zero.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs.measures}{A table of observed relative risk and odds ratio with confidence intervals.}
#' \item{adj.measures}{A table of corrected relative risks and odds ratios.}
#' \item{sim.df}{Data frame of random parameters and computed values.}
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#'
#' @examples
#' # The data for this example come from:
#' # Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O. et al.
#' # Increased risk of infection with human immunodeficiency virus type 1 among
#' # uncircumcised men presenting with genital ulcer disease in Kenya.
#' # Clin Infect Dis 1996;23:449-53.
#' set.seed(123)
#' probsens.conf(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")), nrow = 2, byrow = TRUE),
#' reps = 20000,
#' prev.exp = list("triangular", c(.7, .9, .8)),
#' prev.nexp = list("trapezoidal", c(.03, .04, .05, .06)),
#' risk = list("triangular", c(.6, .7, .63)),
#' corr.p = .8)
#' @export
#' @importFrom stats median qnorm quantile runif rlnorm
probsens.conf <- function(case,
                          exposed,
                          reps = 1000,
                          prev.exp = list(dist = c("constant", "uniform",
                                              "triangular", "trapezoidal",
                                              "logit-logistic", "logit-normal"),
                              parms = NULL),
                          prev.nexp = list(dist = c("constant", "uniform",
                                               "triangular", "trapezoidal",
                                               "logit-logistic", "logit-normal"),
                              parms = NULL),
                          risk = list(dist = c("constant", "uniform", "triangular",
                                          "trapezoidal", "log-logistic",
                                               "log-normal"),
                              parms = NULL),
                          corr.p = NULL,
                          discard = TRUE,
                          alpha = 0.05){
    if(reps < 1)
        stop(paste("Invalid argument: reps =", reps))

    if(is.null(prev.exp) | is.null(prev.nexp))
        stop('Please provide prevalences among the exposed and unexposed.')
    if(is.null(risk))
        stop('Please provide risk of acquiring outcome.')

    if(!is.list(prev.exp))
        stop('Prevalence of exposure among the exposed should be a list.')
    else prev.exp <- prev.exp
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
        seca.parms <- list(prev.exp[[1]], c(prev.exp[[2]], c(0, 1)))
    if(prev.exp[[1]] == "logit-normal" & (length(prev.exp[[2]]) < 2 | length(prev.exp[[2]]) == 3 | length(prev.exp[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(prev.exp[[1]] == "logit-normal" & length(prev.exp[[2]]) == 4 &
       ((prev.exp[[2]][3] >= prev.exp[[2]][4]) | (!all(prev.exp[[2]][3:4] >= 0 & prev.exp[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(prev.exp[[1]] == "logit-normal" & length(prev.exp[[2]]) == 2)
        prev.exp <- list(prev.exp[[1]], c(prev.exp[[2]], c(0, 1)))
    if((prev.exp[[1]] == "constant" | prev.exp[[1]] == "uniform" | prev.exp[[1]] == "triangular" | prev.exp[[1]] == "trapezoidal") & !all(prev.exp[[2]] >= 0 & prev.exp[[2]] <= 1))
        stop('Prevalence should be between 0 and 1.')

        if(!is.list(prev.nexp))
        stop('Prevalence of exposure among the non-exposed should be a list.')
    else prev.nexp <- prev.nexp
    if(prev.nexp[[1]] == "constant" & length(prev.nexp[[2]]) != 1)
        stop('For constant value, please provide a single value.')
    if(prev.nexp[[1]] == "uniform" & length(prev.nexp[[2]]) != 2)
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
        seca.parms <- list(prev.nexp[[1]], c(prev.nexp[[2]], c(0, 1)))
    if(prev.nexp[[1]] == "logit-normal" & (length(prev.nexp[[2]]) < 2 | length(prev.nexp[[2]]) == 3 | length(prev.nexp[[2]]) > 4))
        stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
    if(prev.nexp[[1]] == "logit-normal" & length(prev.nexp[[2]]) == 4 &
       ((prev.nexp[[2]][3] >= prev.nexp[[2]][4]) | (!all(prev.nexp[[2]][3:4] >= 0 & prev.nexp[[2]][3:4] <= 1))))
        stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
    if(prev.nexp[[1]] == "logit-normal" & length(prev.nexp[[2]]) == 2)
        prev.nexp <- list(prev.nexp[[1]], c(prev.nexp[[2]], c(0, 1)))
    if((prev.nexp[[1]] == "constant" | prev.nexp[[1]] == "uniform" | prev.nexp[[1]] == "triangular" | prev.nexp[[1]] == "trapezoidal") & !all(prev.nexp[[2]] >= 0 & prev.nexp[[2]] <= 1))
        stop('Prevalence should be between 0 and 1.')

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
        stop('For log-normal distribution, please provide vector of meanlog and sdlog.')

    if(!is.null(corr.p) && (corr.p == 0 | corr.p == 1))
        stop('Correlations should be > 0 and < 1.')

    if(inherits(case, c("table", "matrix")))
        tab <- case
    else {tab.df <- table(case, exposed)
        tab <- tab.df[2:1, 2:1]
    }
    
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    draws <- matrix(NA, nrow = reps, ncol = 10)
    colnames(draws) <- c("p1", "p0", "RR.cd",
                         "A_11",
                         "A_01",
                         "RR.SMR.rr", 
                         "OR.SMR.or", 
                         "reps",
                         "tot.RRadj.smr", 
                         "tot.ORadj.smr")
    corr.draws <- matrix(NA, nrow = reps, ncol = 5)

    p1 <- c(reps, prev.exp[[2]])
    p0 <- c(reps, prev.nexp[[2]])
    rr.cd <- c(reps, risk[[2]])
    
    obs.rr <- (a/(a + c)) / (b/(b + d))
    se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
    lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
    uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)

    obs.or <- (a/b) / (c/d)
    se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
    uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

    logitlog.dstr <- function(sesp) {
        u <- runif(sesp[[1]])
        w <- sesp[[2]] + sesp[[3]] * (log(u / (1 - u)))
        p <- sesp[[4]] + (sesp[[5]] - sesp[[4]]) * exp(w) / (1 + exp(w))
        return(p)
    }
    logitnorm.dstr <- function(sesp) {
        u <- runif(sesp[[1]])
        w <- sesp[[2]] + sesp[[3]] * qnorm(u)
        p <- sesp[[4]] + (sesp[[5]] - sesp[[4]]) * exp(w) / (1 + exp(w))
        return(p)
    }

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
            draws[, 2] <- logitnorm.dstr(p1)
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
    
    draws[, 8] <- runif(reps)

    draws[, 4] <- (draws[, 3]*a*draws[, 1]*c) / ((draws[, 3]*draws[, 1]*c) + c -
                                                 (draws[, 1]*c))
    draws[, 5] <- (draws[, 3]*b*draws[, 2]*d) / ((draws[, 3]*draws[, 2]*d) + d -
                                                 (draws[, 2]*d))
    draws[, 6] <- (draws[, 4]/(draws[, 4]+draws[, 1]*c)) /
        (draws[, 5]/(draws[, 5]+(draws[, 2]*d)))
    draws[, 7] <- (draws[, 4]*draws[, 2]*d) / (draws[, 5]*draws[, 1]*c)
    
    draws[, 6] <- ifelse(draws[, 4] < 1 |
                               draws[, 5] < 1, NA, draws[, 6])
    draws[, 7] <- ifelse(draws[, 4] < 1 |
                               draws[, 5] < 1, NA, draws[, 7])
    if(all(is.na(draws[, 6])) | all(is.na(draws[, 7])))
        warning('Prior Prevalence distributions lead to all negative adjusted counts.')
    if (discard) {
        if(sum(is.na(draws[, 6])) > 0)
            message('Chosen prior Prevalence distributions lead to ',
                    sum(is.na(draws[, 6])),
                    ' negative adjusted counts which were discarded.')
    }
    else {
        if(sum(is.na(draws[, 6])) > 0) {
            message('Chosen prior Prevalence distributions lead to ',
                    sum(is.na(draws[, 6])),
                    ' negative adjusted counts which were set to zero.')
        }
        draws[, 6] <- ifelse(is.na(draws[, 6]), 0, draws[, 6])
        draws[, 7] <- ifelse(is.na(draws[, 7]), 0, draws[, 7])
        }

    draws[, 9] <- exp(log(draws[, 6]) -
                               qnorm(draws[, 8]) *
                                         ((log(uci.obs.rr) - log(lci.obs.rr)) /
                                              (qnorm(.975) * 2)))
    draws[, 10] <- exp(log(draws[, 7]) -
                               qnorm(draws[, 8]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))

    corr.rr.smr <- c(median(draws[, 6], na.rm = TRUE),
                     quantile(draws[, 6], probs = .025, na.rm = TRUE),
                     quantile(draws[, 6], probs = .975, na.rm = TRUE))
    corr.or.smr <- c(median(draws[, 7], na.rm = TRUE),
                     quantile(draws[, 7], probs = .025, na.rm = TRUE),
                     quantile(draws[, 7], probs = .975, na.rm = TRUE))
    tot.rr.smr <- c(median(draws[, 9], na.rm = TRUE),
                    quantile(draws[, 9], probs = .025, na.rm = TRUE),
                    quantile(draws[, 9], probs = .975, na.rm = TRUE))
    tot.or.smr <- c(median(draws[, 10], na.rm = TRUE),
                    quantile(draws[, 10], probs = .025, na.rm = TRUE),
                    quantile(draws[, 10], probs = .975, na.rm = TRUE))

    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr),
                  c(obs.or, lci.obs.or, uci.obs.or))
    rownames(rmat) <- c("Observed Relative Risk:",
                        "   Observed Odds Ratio:")
    colnames(rmat) <- c(" ",
                        paste(100 * (alpha/2), "%", sep = ""),
                        paste(100 * (1 - alpha/2), "%", sep = ""))
    rmatc <- rbind(corr.rr.smr, tot.rr.smr, corr.or.smr, tot.or.smr)
    rownames(rmatc) <- c("           RR (SMR) -- systematic error:",
                         "RR (SMR) -- systematic and random error:",
                         "           OR (SMR) -- systematic error:",
                         "OR (SMR) -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    res <- list(obs.data = tab,
                obs.measures = rmat, 
                adj.measures = rmatc, 
                sim.df = as.data.frame(draws[, -18]))
    class(res) <- c("episensr", "list")
    res
}

