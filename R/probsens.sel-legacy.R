#' Legacy version of `probsens.sel()`.
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
#' probsens.sel <- probsens.sel_legacy
#' ```
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param reps Number of replications to run.
#' @param or.parms List defining the selection bias odds. The first argument provides the probability distribution function (constant, uniform, triangular, trapezoidal, log-logistic or log-normal) and the second its parameters as a vector:
#' \enumerate{
#' \item constant: constant value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item log-logistic: shape, rate. Must be strictly positive,
#' \item log-normal: meanlog, sdlog. This is the mean and standard deviation on the log scale.
#' }
#' @param case.exp If or.parms not provided, defines the selection probability among case exposed. The first argument provides the probability distribution function and the second its parameters as a vector:
#' \enumerate{
#' \item constant: constant value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item logit-logistic: location, scale, lower bound shift, upper bound shift,
#' \item logit-normal: location, scale, lower bound shift, upper bound shift,
#' \item beta: alpha, beta.
#' }
#' @param case.nexp Same among cases non-exposed.
#' @param ncase.exp Same among non-cases exposed.
#' @param ncase.nexp Same among non-cases non-exposed.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs.measures}{A table of observed odds ratio with confidence intervals.}
#' \item{adj.measures}{A table of corrected odds ratios.}
#' \item{sim.df}{Data frame of random parameters and computed values.}
#' \item{reps}{Number of replications.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#'
#' @examples
#' # The data for this example come from:
#' # Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N. et al.
#' # Population-based incidence estimates of uveal melanoma in Germany.
#' # Supplementing cancer registry data by case-control data.
#' # Eur J Cancer Prev 2006;15:165-70.
#' set.seed(123)
#' probsens.sel(matrix(c(136, 107, 297, 165),
#' dimnames = list(c("Melanoma+", "Melanoma-"), c("Mobile+", "Mobile-")), nrow = 2, byrow = TRUE),
#' reps = 20000,
#' or.parms = list("triangular", c(.35, 1.1, .43)))
#' @export
#' @importFrom stats median qnorm quantile runif rlnorm rbeta
#' @rdname probsens.sel_legacy
probsens.sel_legacy <- function(case,
                                exposed,
                                reps = 1000,
                                or.parms = list(dist = c("constant", "uniform", "triangular",
                                                         "trapezoidal", "log-logistic",
                                                         "log-normal"),
                                                parms = NULL),
                                case.exp = list(dist = c("constant", "uniform", "triangular",
                                                         "trapezoidal", "logit-logistic",
                                                         "logit-normal", "beta"),
                                                parms = NULL),
                                case.nexp = list(dist = c("constant", "uniform", "triangular",
                                                          "trapezoidal", "logit-logistic",
                                                          "logit-normal", "beta"),
                                                 parms = NULL),
                                ncase.exp = list(dist = c("constant", "uniform", "triangular",
                                                          "trapezoidal", "logit-logistic",
                                                          "logit-normal", "beta"),
                                                 parms = NULL),
                                ncase.nexp = list(dist = c("constant", "uniform",
                                                           "triangular", "trapezoidal",
                                                           "logit-logistic", "logit-normal",
                                                           "beta"),
                                                  parms = NULL),
                                alpha = 0.05){
    if(reps < 1)
        stop(paste("Invalid argument: reps =", reps))

    if(is.null(or.parms) & (is.null(case.exp) | is.null(case.nexp) | is.null(ncase.exp) | is.null(ncase.nexp)))
        stop('Please provide selection probabilities.')
    if(!is.null(or.parms[[2]]) & (!is.null(case.exp[[2]]) | !is.null(case.nexp[[2]]) | !is.null(ncase.exp[[2]]) | !is.null(ncase.nexp[[2]])))
        stop('Please use either odds ratio of being selected or selection probabilities.')
    if(!is.list(or.parms))
        stop('Odds ratio for the probability of being selected should be a list.')
    else or.parms <- or.parms
    if(!is.null(or.parms[[2]])){
        if(!(or.parms[[1]] %in% c("constant", "uniform", "triangular", "trapezoidal",
                                  "log-logistic", "log-normal")))
            stop('Wrong distribution for selection odds ratio.')
        if(or.parms[[1]] == "constant" & length(or.parms[[2]]) != 1)
            stop('For constant value, please provide a single value.')
        if(or.parms[[1]] == "uniform" & length(or.parms[[2]]) != 2)
            stop('For uniform distribution, please provide vector of lower and upper limits.')
        if(or.parms[[1]] == "uniform" & or.parms[[2]][1] >= or.parms[[2]][2])
            stop('Lower limit of your uniform distribution is greater than upper limit.')
        if(or.parms[[1]] == "triangular" & length(or.parms[[2]]) != 3)
            stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
        if(or.parms[[1]] == "triangular" & ((or.parms[[2]][1] > or.parms[[2]][3]) |
                                            (or.parms[[2]][2] < or.parms[[2]][3])))
            stop('Wrong arguments for your triangular distribution.')
        if(or.parms[[1]] == "trapezoidal" & length(or.parms[[2]]) != 4)
            stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
        if(or.parms[[1]] == "trapezoidal" & ((or.parms[[2]][1] > or.parms[[2]][2]) |
                                             (or.parms[[2]][2] > or.parms[[2]][3]) |
                                             (or.parms[[2]][3] > or.parms[[2]][4])))
            stop('Wrong arguments for your trapezoidal distribution.')
        if(or.parms[[1]] == "logit-logistic" & length(or.parms[[2]]) != 2)
            stop('For log-logistic distribution, please provide vector of location and scale.')
        if(or.parms[[1]] == "log-normal" & length(or.parms[[2]]) != 2)
            stop('For log-normal distribution, please provide vector of meanlog and sdlog.')
    }

    if(!is.null(case.exp[[1]]) & !is.list(case.exp))
        stop("Please provide a list for case exposed parameters.")
    if(!is.null(case.exp[[2]])){
        if(case.exp[[1]] == "constant" & length(case.exp[[2]]) != 1)
            stop('For constant value, please provide a single value.')
        if(case.exp[[1]] == "uniform" & length(case.exp[[2]]) != 2)
            stop('For uniform distribution, please provide vector of lower and upper limits.')
        if(case.exp[[1]] == "uniform" & case.exp[[2]][1] >= case.exp[[2]][2])
            stop('Lower limit of your uniform distribution is greater than upper limit.')
        if(case.exp[[1]] == "triangular" & length(case.exp[[2]]) != 3)
            stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
        if(case.exp[[1]] == "triangular" & ((case.exp[[2]][1] > case.exp[[2]][3]) |
                                            (case.exp[[2]][2] < case.exp[[2]][3])))
            stop('Wrong arguments for your triangular distribution.')
        if(case.exp[[1]] == "trapezoidal" & length(case.exp[[2]]) != 4)
            stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
        if(case.exp[[1]] == "trapezoidal" & ((case.exp[[2]][1] > case.exp[[2]][2]) |
                                             (case.exp[[2]][2] > case.exp[[2]][3]) |
                                             (case.exp[[2]][3] > case.exp[[2]][4])))
            stop('Wrong arguments for your trapezoidal distribution.')
        if(case.exp[[1]] == "logit-logistic" & (length(case.exp[[2]]) < 2 | length(case.exp[[2]]) == 3 | length(case.exp[[2]]) > 4))
            stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(case.exp[[1]] == "logit-logistic" & length(case.exp[[2]]) == 4 &
           ((case.exp[[2]][3] >= case.exp[[2]][4]) | (!all(case.exp[[2]][3:4] >= 0 & case.exp[[2]][3:4] <= 1))))
            stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(case.exp[[1]] == "logit-logistic" & length(case.exp[[2]]) == 2)
            case.exp <- list(case.exp[[1]], c(case.exp[[2]], c(0, 1)))
        if(case.exp[[1]] == "logit-normal" & (length(case.exp[[2]]) < 2 | length(case.exp[[2]]) == 3 | length(case.exp[[2]]) > 4))
            stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(case.exp[[1]] == "logit-normal" & length(case.exp[[2]]) == 4 &
           ((case.exp[[2]][3] >= case.exp[[2]][4]) | (!all(case.exp[[2]][3:4] >= 0 & case.exp[[2]][3:4] <= 1))))
            stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(case.exp[[1]] == "logit-normal" & length(case.exp[[2]]) == 2)
            case.exp <- list(case.exp[[1]], c(case.exp[[2]], c(0, 1)))
        if((case.exp[[1]] == "constant" | case.exp[[1]] == "uniform" | case.exp[[1]] == "triangular" | case.exp[[1]] == "trapezoidal") & !all(case.exp[[2]] >= 0 & case.exp[[2]] <= 1))
            stop('Selection probability should be between 0 and 1.')
        if(case.exp[[1]] == "beta" & (case.exp[[2]][1] < 0 | case.exp[[2]][1] < 0))
            stop('Wrong arguments for your beta distribution. Alpha and Beta should be > 0')
        if(case.exp[[1]] == "beta" & length(case.exp[[2]]) != 2)
            stop('Wrong arguments for your beta distribution.')
    }

    if(!is.null(case.nexp[[2]])){
        if(case.nexp[[1]] == "constant" & length(case.nexp[[2]]) != 1)
            stop('For constant value, please provide a single value.')
        if(case.nexp[[1]] == "uniform" & length(case.nexp[[2]]) != 2)
            stop('For uniform distribution, please provide vector of lower and upper limits.')
        if(case.nexp[[1]] == "uniform" & case.nexp[[2]][1] >= case.nexp[[2]][2])
            stop('Lower limit of your uniform distribution is greater than upper limit.')
        if(case.nexp[[1]] == "triangular" & length(case.nexp[[2]]) != 3)
            stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
        if(case.nexp[[1]] == "triangular" & ((case.nexp[[2]][1] > case.nexp[[2]][3]) |
                                             (case.nexp[[2]][2] < case.nexp[[2]][3])))
            stop('Wrong arguments for your triangular distribution.')
        if(case.nexp[[1]] == "trapezoidal" & length(case.nexp[[2]]) != 4)
            stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
        if(case.nexp[[1]] == "trapezoidal" & ((case.nexp[[2]][1] > case.nexp[[2]][2]) |
                                              (case.nexp[[2]][2] > case.nexp[[2]][3]) |
                                              (case.nexp[[2]][3] > case.nexp[[2]][4])))
            stop('Wrong arguments for your trapezoidal distribution.')
        if(case.nexp[[1]] == "logit-logistic" & (length(case.nexp[[2]]) < 2 | length(case.nexp[[2]]) == 3 | length(case.nexp[[2]]) > 4))
            stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(case.nexp[[1]] == "logit-logistic" & length(case.nexp[[2]]) == 4 &
           ((case.nexp[[2]][3] >= case.nexp[[2]][4]) | (!all(case.nexp[[2]][3:4] >= 0 & case.nexp[[2]][3:4] <= 1))))
            stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(case.nexp[[1]] == "logit-logistic" & length(case.nexp[[2]]) == 2)
            case.nexp <- list(case.nexp[[1]], c(case.nexp[[2]], c(0, 1)))
        if(case.nexp[[1]] == "logit-normal" & (length(case.nexp[[2]]) < 2 | length(case.nexp[[2]]) == 3 | length(case.nexp[[2]]) > 4))
            stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(case.nexp[[1]] == "logit-normal" & length(case.nexp[[2]]) == 4 &
           ((case.nexp[[2]][3] >= case.nexp[[2]][4]) | (!all(case.nexp[[2]][3:4] >= 0 & case.nexp[[2]][3:4] <= 1))))
            stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(case.nexp[[1]] == "logit-normal" & length(case.nexp[[2]]) == 2)
            case.nexp <- list(case.nexp[[1]], c(case.nexp[[2]], c(0, 1)))
        if((case.nexp[[1]] == "constant" | case.nexp[[1]] == "uniform" | case.nexp[[1]] == "triangular" | case.nexp[[1]] == "trapezoidal") & !all(case.nexp[[2]] >= 0 & case.nexp[[2]] <= 1))
            stop('Selection probability should be between 0 and 1.')
        if(case.nexp[[1]] == "beta" & (case.nexp[[2]][1] < 0 | case.nexp[[2]][1] < 0))
            stop('Wrong arguments for your beta distribution. Alpha and Beta should be > 0')
        if(case.nexp[[1]] == "beta" & length(case.nexp[[2]]) != 2)
            stop('Wrong arguments for your beta distribution.')
    }

    if(!is.null(ncase.exp[[2]])){
        if(ncase.exp[[1]] == "constant" & length(ncase.exp[[2]]) != 1)
            stop('For constant value, please provide a single value.')
        if(ncase.exp[[1]] == "uniform" & length(ncase.exp[[2]]) != 2)
            stop('For uniform distribution, please provide vector of lower and upper limits.')
        if(ncase.exp[[1]] == "uniform" & ncase.exp[[2]][1] >= ncase.exp[[2]][2])
            stop('Lower limit of your uniform distribution is greater than upper limit.')
        if(ncase.exp[[1]] == "triangular" & length(ncase.exp[[2]]) != 3)
            stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
        if(ncase.exp[[1]] == "triangular" & ((ncase.exp[[2]][1] > ncase.exp[[2]][3]) |
                                             (ncase.exp[[2]][2] < ncase.exp[[2]][3])))
            stop('Wrong arguments for your triangular distribution.')
        if(ncase.exp[[1]] == "trapezoidal" & length(ncase.exp[[2]]) != 4)
            stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
        if(ncase.exp[[1]] == "trapezoidal" & ((ncase.exp[[2]][1] > ncase.exp[[2]][2]) |
                                              (ncase.exp[[2]][2] > ncase.exp[[2]][3]) |
                                              (ncase.exp[[2]][3] > ncase.exp[[2]][4])))
            stop('Wrong arguments for your trapezoidal distribution.')
        if(ncase.exp[[1]] == "logit-logistic" & (length(ncase.exp[[2]]) < 2 | length(ncase.exp[[2]]) == 3 | length(ncase.exp[[2]]) > 4))
            stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(ncase.exp[[1]] == "logit-logistic" & length(ncase.exp[[2]]) == 4 &
           ((ncase.exp[[2]][3] >= ncase.exp[[2]][4]) | (!all(ncase.exp[[2]][3:4] >= 0 & ncase.exp[[2]][3:4] <= 1))))
            stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(ncase.exp[[1]] == "logit-logistic" & length(ncase.exp[[2]]) == 2)
            ncase.exp <- list(ncase.exp[[1]], c(ncase.exp[[2]], c(0, 1)))
        if(ncase.exp[[1]] == "logit-normal" & (length(ncase.exp[[2]]) < 2 | length(ncase.exp[[2]]) == 3 | length(ncase.exp[[2]]) > 4))
            stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(ncase.exp[[1]] == "logit-normal" & length(ncase.exp[[2]]) == 4 &
           ((ncase.exp[[2]][3] >= ncase.exp[[2]][4]) | (!all(ncase.exp[[2]][3:4] >= 0 & ncase.exp[[2]][3:4] <= 1))))
            stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(ncase.exp[[1]] == "logit-normal" & length(ncase.exp[[2]]) == 2)
            ncase.exp <- list(ncase.exp[[1]], c(ncase.exp[[2]], c(0, 1)))
        if((ncase.exp[[1]] == "constant" | ncase.exp[[1]] == "uniform" | ncase.exp[[1]] == "triangular" | ncase.exp[[1]] == "trapezoidal") & !all(ncase.exp[[2]] >= 0 & ncase.exp[[2]] <= 1))
            stop('Selection probability should be between 0 and 1.')
        if(ncase.exp[[1]] == "beta" & (ncase.exp[[2]][1] < 0 | ncase.exp[[2]][1] < 0))
            stop('Wrong arguments for your beta distribution. Alpha and Beta should be > 0')
        if(ncase.exp[[1]] == "beta" & length(ncase.exp[[2]]) != 2)
            stop('Wrong arguments for your beta distribution.')
    }

    if(!is.null(ncase.nexp[[2]])){
        if(ncase.nexp[[1]] == "constant" & length(ncase.nexp[[2]]) != 1)
            stop('For constant value, please provide a single value.')
        if(ncase.nexp[[1]] == "uniform" & length(ncase.nexp[[2]]) != 2)
            stop('For uniform distribution, please provide vector of lower and upper limits.')
        if(ncase.nexp[[1]] == "uniform" & ncase.nexp[[2]][1] >= ncase.nexp[[2]][2])
            stop('Lower limit of your uniform distribution is greater than upper limit.')
        if(ncase.nexp[[1]] == "triangular" & length(ncase.nexp[[2]]) != 3)
            stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
        if(ncase.nexp[[1]] == "triangular" & ((ncase.nexp[[2]][1] > ncase.nexp[[2]][3]) |
                                              (ncase.nexp[[2]][2] < ncase.nexp[[2]][3])))
            stop('Wrong arguments for your triangular distribution.')
        if(ncase.nexp[[1]] == "trapezoidal" & length(ncase.nexp[[2]]) != 4)
            stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
        if(ncase.nexp[[1]] == "trapezoidal" & ((ncase.nexp[[2]][1] > ncase.nexp[[2]][2]) |
                                               (ncase.nexp[[2]][2] > ncase.nexp[[2]][3]) |
                                               (ncase.nexp[[2]][3] > ncase.nexp[[2]][4])))
            stop('Wrong arguments for your trapezoidal distribution.')
        if(ncase.nexp[[1]] == "logit-logistic" & (length(ncase.nexp[[2]]) < 2 | length(ncase.nexp[[2]]) == 3 | length(ncase.nexp[[2]]) > 4))
            stop('For logit-logistic distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(ncase.nexp[[1]] == "logit-logistic" & length(ncase.nexp[[2]]) == 4 &
           ((ncase.nexp[[2]][3] >= ncase.nexp[[2]][4]) | (!all(ncase.nexp[[2]][3:4] >= 0 & ncase.nexp[[2]][3:4] <= 1))))
            stop('For logit-logistic distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(ncase.nexp[[1]] == "logit-logistic" & length(ncase.nexp[[2]]) == 2)
            ncase.nexp <- list(ncase.nexp[[1]], c(ncase.nexp[[2]], c(0, 1)))
        if(ncase.nexp[[1]] == "logit-normal" & (length(ncase.nexp[[2]]) < 2 | length(ncase.nexp[[2]]) == 3 | length(ncase.nexp[[2]]) > 4))
            stop('For logit-normal distribution, please provide vector of location, scale, and eventually lower and upper bound limits if you want to shift and rescale the distribution.')
        if(ncase.nexp[[1]] == "logit-normal" & length(ncase.nexp[[2]]) == 4 &
           ((ncase.nexp[[2]][3] >= ncase.nexp[[2]][4]) | (!all(ncase.nexp[[2]][3:4] >= 0 & ncase.nexp[[2]][3:4] <= 1))))
            stop('For logit-normal distribution, please provide sensible values for lower and upper bound limits (between 0 and 1; lower limit < upper limit).')
        if(ncase.nexp[[1]] == "logit-normal" & length(ncase.nexp[[2]]) == 2)
            ncase.nexp <- list(ncase.nexp[[1]], c(ncase.nexp[[2]], c(0, 1)))
        if((ncase.nexp[[1]] == "constant" | ncase.nexp[[1]] == "uniform" | ncase.nexp[[1]] == "triangular" | ncase.nexp[[1]] == "trapezoidal") & !all(ncase.nexp[[2]] >= 0 & ncase.nexp[[2]] <= 1))
            stop('Selection probability should be between 0 and 1.')
        if(ncase.nexp[[1]] == "beta" & (ncase.nexp[[2]][1] < 0 | ncase.nexp[[2]][1] < 0))
            stop('Wrong arguments for your beta distribution. Alpha and Beta should be > 0')
        if(ncase.nexp[[1]] == "beta" & length(ncase.nexp[[2]]) != 2)
            stop('Wrong arguments for your beta distribution.')
    }


    if(!inherits(case, "episensr.probsens")){
        if(inherits(case, c("table", "matrix")))
            tab <- case
        else {tab.df <- table(case, exposed)
            tab <- tab.df[2:1, 2:1]
        }

        a <- as.numeric(tab[1, 1])
        b <- as.numeric(tab[1, 2])
        c <- as.numeric(tab[2, 1])
        d <- as.numeric(tab[2, 2])

        obs.or <- (a/b) / (c/d)
        se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
        uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)
    } else {
        a <- as.numeric(case[[3]][, 1])
        b <- as.numeric(case[[3]][, 2])
        c <- as.numeric(case[[3]][, 3])
        d <- as.numeric(case[[3]][, 4])

        obs.or <- (a/b) / (c/d)
        se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
        uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

        reps <- case[[4]]
    }

    draws <- matrix(NA, nrow = reps, ncol = 8)
    colnames(draws) <- c("or.sel", "corr.OR", "reps", "tot.OR",
                         "A1", "B1", "C1", "D1")

    if(!is.null(or.parms[[2]])){
        or.sel <- c(reps, or.parms[[2]])

        if (or.parms[[1]] == "constant") {
            draws[, 1] <- or.parms[[2]]
        }
        if (or.parms[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(or.sel))
        }
        if (or.parms[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(or.sel))
        }
        if (or.parms[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(or.sel))
        }
        if (or.parms[[1]] == "log-logistic") {
            draws[, 1] <- do.call(actuar::rllogis, as.list(or.sel))
        }
        if (or.parms[[1]] == "log-normal") {
            draws[, 1] <- do.call(rlnorm, as.list(or.sel))
        }
    } else {
        bias_factor <- matrix(NA, nrow = reps, ncol = 4)

        case_exp <- c(reps, case.exp[[2]])
        case_nexp <- c(reps, case.nexp[[2]])
        ncase_exp <- c(reps, ncase.exp[[2]])
        ncase_nexp <- c(reps, ncase.nexp[[2]])

        if (case.exp[[1]] == "constant") {
            bias_factor[, 1] <- case.exp[[2]]
        }
        if (case.exp[[1]] == "uniform") {
            bias_factor[, 1] <- do.call(runif, as.list(case_exp))
        }
        if (case.exp[[1]] == "triangular") {
            bias_factor[, 1] <- do.call(triangle::rtriangle, as.list(case_exp))
        }
        if (case.exp[[1]] == "trapezoidal") {
            bias_factor[, 1] <- do.call(trapezoid::rtrapezoid, as.list(case_exp))
        }
        if (case.exp[[1]] == "logit-logistic") {
            bias_factor[, 1] <- logitlog.dstr(case_exp)
        }
        if (case.exp[[1]] == "logit-normal") {
            bias_factor[, 1] <- logitnorm.dstr(case_exp)
        }
        if (case.exp[[1]] == "beta") {
            bias_factor[, 1] <- do.call(rbeta, as.list(case_exp))
        }

        if (case.nexp[[1]] == "constant") {
            bias_factor[, 2] <- case.nexp[[2]]
        }
        if (case.nexp[[1]] == "uniform") {
            bias_factor[, 2] <- do.call(runif, as.list(case_nexp))
        }
        if (case.nexp[[1]] == "triangular") {
            bias_factor[, 2] <- do.call(triangle::rtriangle, as.list(case_nexp))
        }
        if (case.nexp[[1]] == "trapezoidal") {
            bias_factor[, 2] <- do.call(trapezoid::rtrapezoid, as.list(case_nexp))
        }
        if (case.nexp[[1]] == "logit-logistic") {
            bias_factor[, 2] <- logitlog.dstr(case_nexp)
        }
        if (case.nexp[[1]] == "logit-normal") {
            bias_factor[, 2] <- logitnorm.dstr(case_nexp)
        }
        if (case.nexp[[1]] == "beta") {
            bias_factor[, 2] <- do.call(rbeta, as.list(case_nexp))
        }

        if (ncase.exp[[1]] == "constant") {
            bias_factor[, 3] <- ncase.exp[[2]]
        }
        if (ncase.exp[[1]] == "uniform") {
            bias_factor[, 3] <- do.call(runif, as.list(ncase_exp))
        }
        if (ncase.exp[[1]] == "triangular") {
            bias_factor[, 3] <- do.call(triangle::rtriangle, as.list(ncase_exp))
        }
        if (ncase.exp[[1]] == "trapezoidal") {
            bias_factor[, 3] <- do.call(trapezoid::rtrapezoid, as.list(ncase_exp))
        }
        if (ncase.exp[[1]] == "logit-logistic") {
            bias_factor[, 3] <- logitlog.dstr(ncase_exp)
        }
        if (ncase.exp[[1]] == "logit-normal") {
            bias_factor[, 3] <- logitnorm.dstr(ncase_exp)
        }
        if (ncase.exp[[1]] == "beta") {
            bias_factor[, 3] <- do.call(rbeta, as.list(ncase_exp))
        }

        if (ncase.nexp[[1]] == "constant") {
            bias_factor[, 4] <- ncase.nexp[[2]]
        }
        if (ncase.nexp[[1]] == "uniform") {
            bias_factor[, 4] <- do.call(runif, as.list(ncase_nexp))
        }
        if (ncase.nexp[[1]] == "triangular") {
            bias_factor[, 4] <- do.call(triangle::rtriangle, as.list(ncase_nexp))
        }
        if (ncase.nexp[[1]] == "trapezoidal") {
            bias_factor[, 4] <- do.call(trapezoid::rtrapezoid, as.list(ncase_nexp))
        }
        if (ncase.nexp[[1]] == "logit-logistic") {
            bias_factor[, 4] <- logitlog.dstr(ncase_nexp)
        }
        if (ncase.nexp[[1]] == "logit-normal") {
            bias_factor[, 4] <- logitnorm.dstr(ncase_nexp)
        }
        if (ncase.nexp[[1]] == "beta") {
            bias_factor[, 4] <- do.call(rbeta, as.list(ncase_nexp))
        }

        draws[, 1] <- (bias_factor[, 1]*bias_factor[, 4]) /
            (bias_factor[, 2]*bias_factor[, 3])
    }

    draws[, 3] <- runif(reps)

    draws[, 2] <- obs.or / draws[, 1]

    draws[, 4] <- exp(log(draws[, 2]) -
                               qnorm(draws[, 3]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                          (qnorm(.975) * 2)))

    draws[, 5] <- a / draws[, 1]
    draws[, 6] <- b / draws[, 1]
    draws[, 7] <- c / draws[, 1]
    draws[, 8] <- d / draws[, 1]

    corr.OR <- c(median(draws[, 2], na.rm = TRUE),
                 quantile(draws[, 2], probs = .025, na.rm = TRUE),
                 quantile(draws[, 2], probs = .975, na.rm = TRUE))
    tot.OR <- c(median(draws[, 4], na.rm = TRUE),
                quantile(draws[, 4], probs = .025, na.rm = TRUE),
                quantile(draws[, 4], probs = .975, na.rm = TRUE))

    if(!inherits(case, "episensr.probsens")){
        tab <- tab
        rmat <- matrix(c(obs.or, lci.obs.or, uci.obs.or), nrow = 1)
        rownames(rmat) <- c("Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
    } else {
        tab <- case[[1]]
        rmat <- case[[2]]
    }
    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    rmatc <- rbind(corr.OR, tot.OR)
    rownames(rmatc) <- c("           Odds Ratio -- systematic error:",
                         "Odds Ratio -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    res <- list(obs.data = tab,
                obs.measures = rmat,
                adj.measures = rmatc,
                sim.df = as.data.frame(draws[, -3]),
                reps = reps,
                fun = "probsens.sel")
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
