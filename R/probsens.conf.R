#' Probabilistic sensitivity analysis for unmeasured confounding.
#'
#' Probabilistic sensitivity analysis to correct for unknown or unmeasured
#' confounding and random error simultaneously.
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
#' computation using [probsens.conf_legacy()] as follows:
#'
#' ```
#' library(episensr)
#' probsens.conf <- probsens.conf_legacy
#' ```
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param reps Number of replications to run.
#' @param prev_exp List defining the prevalence of exposure among the exposed.
#' The first argument provides the probability distribution function (constant,
#' uniform, triangular, trapezoidal, truncated normal, or beta) and the second its
#' parameters as a vector. Lower bound of the truncated normal cannot be less than
#' zero. Upper bound is Inf by default.
#' \enumerate{
#' \item constant: constant value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item normal: lower bound, upper bound, mean, sd.
#' \item beta: alpha, beta.
#' }
#' @param prev_nexp List defining the prevalence of exposure among the unexposed.
#' @param risk List defining the confounder-disease relative risk or the confounder-exposure
#' odds ratio. The first argument provides the probability distribution function
#' (constant, uniform, triangular, trapezoidal, log-logistic, or log-normal) and
#' the second its parameters as a vector:
#' \enumerate{
#' \item constant: constant value,
#' \item uniform: min, max,
#' \item triangular: lower limit, upper limit, mode,
#' \item trapezoidal: min, lower mode, upper mode, max.
#' \item log-logistic: shape, rate. Must be strictly positive,
#' \item log-normal: meanlog, sdlog. This is the mean and standard deviation on the log scale.
#' }
#' @param corr_p Correlation between the exposure-specific confounder prevalences.
#' @param discard A logical scalar. In case of negative adjusted count, should
#' the draws be discarded? If set to FALSE, negative counts are set to zero.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs_measures}{A table of observed relative risk and odds ratio with confidence intervals.}
#' \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#' \item{sim_df}{Data frame of random parameters and computed values.}
#' \item{reps}{Number of replications.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.
#'
#' Li, S.T., Hammond, J.L., 1975. \emph{Generation of Pseudorandom Numbers
#' with Specified Univariate Distributions and Correlation Coefficients}.
#' IEEE Trans Syst Man Cybern 5:557-561.
#' @examples
#' # The data for this example come from:
#' # Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O. et al.
#' # Increased risk of infection with human immunodeficiency virus type 1 among
#' # uncircumcised men presenting with genital ulcer disease in Kenya.
#' # Clin Infect Dis 1996;23:449-53.
#' tyndall <- matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")), nrow = 2, byrow = TRUE)
#' set.seed(123)
#' probsens.conf(tyndall, reps = 20000,
#' prev_exp = list("triangular", c(.7, .9, .8)),
#' prev_nexp = list("trapezoidal", c(.03, .04, .05, .06)),
#' risk = list("triangular", c(.6, .7, .63)),
#' corr_p = .8)
#'
#' set.seed(123)
#' probsens.conf(tyndall, reps = 20000,
#' prev_exp = list("beta", c(200, 56)),
#' prev_nexp = list("beta", c(10, 16)),
#' risk = list("triangular", c(.6, .7, .63)),
#' corr_p = .8)
#'
#' set.seed(123)
#' probsens.conf(tyndall, reps = 20000,
#' prev_exp = list("normal", c(.01, .12, 0.03, 0.005)),
#' prev_nexp = list("normal", c(0, Inf, 0.01, 0.0001)),
#' risk = list("triangular", c(.6, .7, .63)), corr_p = .8)
#' @export
#' @importFrom stats median qnorm quantile runif rlnorm rbeta qbeta
probsens.conf <- function(case,
                          exposed,
                          reps = 1000,
                          prev_exp = list(dist = c("constant", "uniform",
                                                   "triangular", "trapezoidal",
                                                   "normal", "beta"),
                                          parms = NULL),
                          prev_nexp = list(dist = c("constant", "uniform",
                                                    "triangular", "trapezoidal",
                                                    "normal", "beta"),
                                           parms = NULL),
                          risk = list(dist = c("constant", "uniform",
                                               "triangular", "trapezoidal",
                                               "log-logistic", "log-normal"),
                                      parms = NULL),
                          corr_p = NULL,
                          discard = TRUE,
                          alpha = 0.05) {
    if (reps < 1)
        stop(paste("Invalid argument: reps =", reps))

    if (is.null(prev_exp) | is.null(prev_nexp))
        stop("Please provide prevalences among the exposed and unexposed.")
    if (is.null(risk))
        stop("Please provide risk of acquiring outcome.")

    if (!is.list(prev_exp))
        stop("Prevalence of exposure among the exposed should be a list.")
    else prev_exp <- prev_exp
    if ((length(prev_exp) != 2) | (length(prev_nexp) != 2) | (length(risk) != 2))
        stop("Check distribution parameters.")
    if ((length(prev_exp[[1]]) != 1) | (length(prev_nexp[[1]]) != 1) | (length(risk[[1]]) != 1))
        stop("Which distribution?")
    if (!is.null(corr_p) && (prev_exp[[1]] == "constant" | prev_nexp[[1]] == "constant"))
        stop("No correlated distributions with constant values.")
    if (prev_exp[[1]] == "constant" & length(prev_exp[[2]]) != 1)
        stop("For constant value, please provide a single value.")
    if (prev_exp[[1]] == "uniform" & length(prev_exp[[2]]) != 2)
        stop("For uniform distribution, please provide vector of lower and upper limits.")
    if (prev_exp[[1]] == "uniform" & prev_exp[[2]][1] >= prev_exp[[2]][2])
        stop("Lower limit of your uniform distribution is greater than upper limit.")
    if (prev_exp[[1]] == "triangular" & length(prev_exp[[2]]) != 3)
        stop("For triangular distribution, please provide vector of lower, upper limits, and mode.")
    if (prev_exp[[1]] == "triangular" & ((prev_exp[[2]][1] > prev_exp[[2]][3]) |
                                         (prev_exp[[2]][2] < prev_exp[[2]][3])))
        stop("Wrong arguments for your triangular distribution.")
    if (prev_exp[[1]] == "trapezoidal" & length(prev_exp[[2]]) != 4)
        stop("For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")
    if (prev_exp[[1]] == "trapezoidal" & ((prev_exp[[2]][1] > prev_exp[[2]][2]) |
                                          (prev_exp[[2]][2] > prev_exp[[2]][3]) |
                                          (prev_exp[[2]][3] > prev_exp[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (prev_exp[[1]] == "normal" & (length(prev_exp[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bound limits, mean and sd.")
    if (prev_exp[[1]] == "normal" & length(prev_exp[[2]]) == 4 &
        ((prev_exp[[2]][1] >= prev_exp[[2]][2]) | (prev_exp[[2]][1] < 0)))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (lower limit >= 0; lower limit < upper limit).")
    if ((prev_exp[[1]] == "constant" | prev_exp[[1]] == "uniform" |
         prev_exp[[1]] == "triangular" | prev_exp[[1]] == "trapezoidal") &
        !all(prev_exp[[2]] >= 0 & prev_exp[[2]] <= 1))
        stop("Prevalence should be between 0 and 1.")
    if (!is.null(prev_exp) && prev_exp[[1]] == "beta" && length(prev_exp[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (!is.null(prev_exp) && prev_exp[[1]] == "beta" &&
       (prev_exp[[2]][1] < 0 | prev_exp[[2]][2] < 0))
        stop("Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")

    if (!is.list(prev_nexp))
        stop("Prevalence of exposure among the non-exposed should be a list.")
    else prev_nexp <- prev_nexp
    if (prev_nexp[[1]] == "constant" & length(prev_nexp[[2]]) != 1)
        stop("For constant value, please provide a single value.")
    if (prev_nexp[[1]] == "uniform" & length(prev_nexp[[2]]) != 2)
        stop("For uniform distribution, please provide vector of lower and upper limits.")
    if (prev_nexp[[1]] == "uniform" & prev_nexp[[2]][1] >= prev_nexp[[2]][2])
        stop("Lower limit of your uniform distribution is greater than upper limit.")
    if (prev_nexp[[1]] == "triangular" & length(prev_nexp[[2]]) != 3)
        stop("For triangular distribution, please provide vector of lower, upper limits, and mode.")
    if (prev_nexp[[1]] == "triangular" & ((prev_nexp[[2]][1] > prev_nexp[[2]][3]) |
                                          (prev_nexp[[2]][2] < prev_nexp[[2]][3])))
        stop("Wrong arguments for your triangular distribution.")
    if (prev_nexp[[1]] == "trapezoidal" & length(prev_nexp[[2]]) != 4)
        stop("For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")
    if (prev_nexp[[1]] == "trapezoidal" & ((prev_nexp[[2]][1] > prev_nexp[[2]][2]) |
                                           (prev_nexp[[2]][2] > prev_nexp[[2]][3]) |
                                           (prev_nexp[[2]][3] > prev_nexp[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (prev_nexp[[1]] == "normal" & (length(prev_nexp[[2]]) != 4))
        stop("For truncated normal distribution, please provide vector of lower and upper bound limits, mean and sd.")
    if (prev_nexp[[1]] == "normal" & length(prev_nexp[[2]]) == 4 &
        ((prev_nexp[[2]][1] >= prev_nexp[[2]][2]) | (prev_nexp[[2]][1] < 0)))
        stop("For truncated normal distribution, please provide sensible values for lower and upper bound limits (lower limit >= 0; lower limit < upper limit).")
    if ((prev_nexp[[1]] == "constant" | prev_nexp[[1]] == "uniform" |
         prev_nexp[[1]] == "triangular" | prev_nexp[[1]] == "trapezoidal") &
        !all(prev_nexp[[2]] >= 0 & prev_nexp[[2]] <= 1))
        stop("Prevalence should be between 0 and 1.")
    if (!is.null(prev_nexp) && prev_nexp[[1]] == "beta" && length(prev_nexp[[2]]) != 2)
        stop("For beta distribution, please provide alpha and beta.")
    if (!is.null(prev_nexp) && prev_nexp[[1]] == "beta" && (prev_nexp[[2]][1] < 0 | prev_nexp[[2]][2] < 0))
        stop("Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")

    if (!is.list(risk))
        stop("Risk should be a list.")
    else risk <- risk
    if (risk[[1]] == "constant" & length(risk[[2]]) != 1)
        stop("For constant value, please provide a single value.")
    if (risk[[1]] == "uniform" & length(risk[[2]]) != 2)
        stop("For uniform distribution, please provide vector of lower and upper limits.")
    if (risk[[1]] == "uniform" & risk[[2]][1] >= risk[[2]][2])
        stop("Lower limit of your uniform distribution is greater than upper limit.")
    if (risk[[1]] == "triangular" & length(risk[[2]]) != 3)
        stop("For triangular distribution, please provide vector of lower, upper limits, and mode.")
    if (risk[[1]] == "triangular" & ((risk[[2]][1] > risk[[2]][3]) | (risk[[2]][2] < risk[[2]][3])))
        stop("Wrong arguments for your triangular distribution.")
    if (risk[[1]] == "trapezoidal" & length(risk[[2]]) != 4)
        stop("For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")
    if (risk[[1]] == "trapezoidal" & ((risk[[2]][1] > risk[[2]][2]) |
                                      (risk[[2]][2] > risk[[2]][3]) |
                                      (risk[[2]][3] > risk[[2]][4])))
        stop("Wrong arguments for your trapezoidal distribution.")
    if (risk[[1]] == "log-logistic" & length(risk[[2]]) != 2)
        stop("For log-logistic distribution, please provide vector of location and scale.")
    if (risk[[1]] == "log-normal" & length(risk[[2]]) != 2)
        stop("For log-normal distribution, please provide vector of meanlog and sdlog.")

    if (!is.null(corr_p) && (corr_p == 0 | corr_p == 1))
        stop("Correlations should be > 0 and < 1.")

    if (!inherits(case, "episensr.probsens")) {
        if (inherits(case, c("table", "matrix")))
            tab <- case
        else {
            tab_df <- table(case, exposed)
            tab <- tab_df[2:1, 2:1]
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

    obs_rr <- (a / (a + c)) / (b / (b + d))
    se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
    lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
    uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

    obs_or <- (a / b) / (c / d)
    se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
    uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

    draws <- matrix(NA, nrow = reps, ncol = 16)
    colnames(draws) <- c("p1", "p0", "RR_cd",
                         "A_11",
                         "B_11",
                         "C_11",
                         "D_11",
                         "corr_RR",
                         "corr_OR",
                         "reps",
                         "tot_RR",
                         "tot_OR",
                         "A1", "B1", "C1", "D1")
    corr_draws <- matrix(NA, nrow = reps, ncol = 2)

    p1 <- c(reps, prev_exp[[2]])
    p0 <- c(reps, prev_nexp[[2]])
    rr_cd <- c(reps, risk[[2]])

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

    draws[, 10] <- runif(reps)

    draws[, 4] <- (draws[, 3] * ((a + c) * draws[, 1]) * a) /
        (draws[, 3] * ((a + c) * draws[, 1]) + (a + c) - ((a + c) * draws[, 1]))  # A_11
    draws[, 5] <- (draws[, 3] * b * (draws[, 2] * (b + d))) /
        (draws[, 3] * (draws[, 2] * (b + d)) + (b + d) - (draws[, 2] * (b + d)))  # B_11
    draws[, 6] <- ((a + c) * draws[, 1]) - draws[, 4]  # C_11
    draws[, 7] <- ((b + d) * draws[, 2]) - draws[, 5]  # D_11
    draws[, 8] <- a / ((((a + c) * draws[, 1]) * draws[, 5] / ((b + d) * draws[, 2])) +
                       (((a + c) - ((a + c) * draws[, 1])) * (b - draws[, 5]) /
                        ((b + d) - ((b + d) * draws[, 2]))))  # RR.SMR.rr
    draws[, 9] <- a / ((draws[, 6] * draws[, 5] / draws[, 7]) +
                       ((c - draws[, 6]) * (b - draws[, 5]) / (d - draws[, 7])))  # OR.SMR.or
    draws[, 8] <- ifelse(draws[, 4] < 0 |
                         draws[, 5] < 0 |
                         draws[, 6] < 0 |
                         draws[, 7] < 0 |
                         (a - draws[, 4]) < 0 |
                         (b - draws[, 5]) < 0 |
                         (c - draws[, 6]) < 0 |
                         (d - draws[, 7]) < 0, NA, draws[, 8])
    draws[, 9] <- ifelse(draws[, 4] < 0 |
                         draws[, 5] < 0 |
                         draws[, 6] < 0 |
                         draws[, 7] < 0 |
                         (a - draws[, 4]) < 0 |
                         (b - draws[, 5]) < 0 |
                         (c - draws[, 6]) < 0 |
                         (d - draws[, 7]) < 0, NA, draws[, 9])
    if (all(is.na(draws[, 8])) | all(is.na(draws[, 9]))) {
        warning("Prior Prevalence distributions lead to all negative adjusted counts.")
        neg_warn <- "Prior Se/Sp distributions lead to all negative adjusted counts."
    } else {
        neg_warn <- NULL
    }
    if (discard) {
        if (sum(is.na(draws[, 8])) > 0) {
            message("Chosen prior Prevalence distributions lead to ",
                    sum(is.na(draws[, 8])),
                    " negative adjusted counts which were discarded.")
            discard_mess <- c(paste("Chosen prior Se/Sp distributions lead to ",
                                    sum(is.na(draws[, 8])),
                                    " negative adjusted counts which were discarded."))
        } else discard_mess <- NULL
    } else {
        if (sum(is.na(draws[, 8])) > 0) {
            message("Chosen prior Prevalence distributions lead to ",
                    sum(is.na(draws[, 8])),
                    " negative adjusted counts which were set to zero.")
            discard_mess <- c(paste("Chosen prior Se/Sp distributions lead to ",
                                    sum(is.na(draws[, 8])),
                                    " negative adjusted counts which were set to zero."))
            draws[, 8] <- ifelse(is.na(draws[, 8]), 0, draws[, 8])
            draws[, 9] <- ifelse(is.na(draws[, 9]), 0, draws[, 9])
        } else discard_mess <- NULL
    }

    draws[, 11] <- exp(log(draws[, 8]) - qnorm(draws[, 10]) * ((log(uci_obs_rr) - log(lci_obs_rr)) /
                                                               (qnorm(.975) * 2)))
    draws[, 12] <- exp(log(draws[, 9]) - qnorm(draws[, 10]) * ((log(uci_obs_or) - log(lci_obs_or)) /
                                                               (qnorm(.975) * 2)))
    draws[, 13] <- a
    draws[, 14] <- b
    draws[, 15] <- c
    draws[, 16] <- d

    corr_RR <- c(median(draws[, 8], na.rm = TRUE),
                 quantile(draws[, 8], probs = .025, na.rm = TRUE),
                 quantile(draws[, 8], probs = .975, na.rm = TRUE))
    corr_OR <- c(median(draws[, 9], na.rm = TRUE),
                 quantile(draws[, 9], probs = .025, na.rm = TRUE),
                 quantile(draws[, 9], probs = .975, na.rm = TRUE))
    tot_RR <- c(median(draws[, 11], na.rm = TRUE),
                quantile(draws[, 11], probs = .025, na.rm = TRUE),
                quantile(draws[, 11], probs = .975, na.rm = TRUE))
    tot_OR <- c(median(draws[, 12], na.rm = TRUE),
                quantile(draws[, 12], probs = .025, na.rm = TRUE),
                quantile(draws[, 12], probs = .975, na.rm = TRUE))

    if (!inherits(case, "episensr.probsens")) {
        tab <- tab
        rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr),
                      c(obs_or, lci_obs_or, uci_obs_or))
        rownames(rmat) <- c("Observed Relative Risk:",
                            "   Observed Odds Ratio:")
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
    rmatc <- rbind(corr_RR, tot_RR, corr_OR, tot_OR)
    rownames(rmatc) <- c("           RR (SMR) -- systematic error:",
                         "RR (SMR) -- systematic and random error:",
                         "           OR (SMR) -- systematic error:",
                         "OR (SMR) -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    res <- list(obs_data = tab,
                obs_measures = rmat,
                adj_measures = rmatc,
                sim_df = as.data.frame(draws[, -10]),
                reps = reps,
                fun = "probsens.conf",
                warnings = neg_warn,
                message = discard_mess)
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
