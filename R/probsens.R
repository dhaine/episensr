#' Probabilistic sensitivity analysis.
#'
#' Summary-level probabilistic sensitivity analysis to correct for exposure
#' misclassification or outcome misclassification and random error.
#' Non-differential misclassification is assumed when only the two bias
#' parameters \code{seca} and \code{spca} are provided. Adding the 2
#' parameters \code{seexp} and \code{spexp} (i.e. providing the 4 bias
#' parameters) evaluates a differential misclassification.
#'
#' Correlations between sensitivity (or specificity) of exposure classification
#' among cases and controls can be specified and use the NORmal To Anything
#' (NORTA) transformation (Li & Hammond, 1975).
#'
#' In case of negative (<=0) adjusted count in the 2-by-2 table following given
#' prior Se/Sp distribution(s), draws are discarded.
#'
#' @section Updated calculations:
#' episensr 2.0.0 introduced updated calculations of probabilistic bias analyses
#' by (1) using the NORTA transformation to define a correlation between
#' distributions, and (2) sampling true prevalences and then sampling the
#' adjusted cell counts rather than just using the expected cell counts from a
#' simple quantitative bias analysis. This updated version should be preferred
#' but if you need to run an old analysis, you can easily revert to the
#' computation using [probsens_legacy()] as follows:
#'
#' ```
#' library(episensr)
#' probsens <- probsens_legacy
#' ```
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param type Choice of correction for exposure or outcome misclassification.
#' @param reps Number of replications to run.
#' @param seca List defining sensitivity among cases:
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
#' @param seexp List defining sensitivity among controls:
#' \enumerate{
#' \item The sensitivity of exposure classification among those without the
#' outcome (when \code{type = "exposure"}), or
#' \item The sensitivity of outcome classification among those without the
#' exposure (when \code{type = "outcome"}).
#' }
#' @param spca List as above for \code{seca} but for specificity.
#' @param spexp List as above for \code{seexp} but for specificity.
#' @param corr_se Correlation between case and non-case sensitivities.
#' @param corr_sp Correlation between case and non-case specificities.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs_measures}{A table of observed relative risk and odds ratio with
#' confidence intervals.}
#' \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#' \item{sim_df}{Data frame of random parameters and computed values.}
#' \item{reps}{Number of replications.}
#'
#' @references
#' Fox, M.P, MacLehose, R.F., Lash, T.L., 2021 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.233--290, Springer.
#'
#' Li, S.T., Hammond, J.L., 1975. \emph{Generation of Pseudorandom Numbers
#' with Specified Univariate Distributions and Correlation Coefficients}.
#' IEEE Trans Syst Man Cybern 5:557-561.
#' @examples
#' # The data for this example come from:
#' # Greenland S., Salvan A., Wegman D.H., Hallock M.F., Smith T.J.
#' # A case-control study of cancer mortality at a transformer-assembly facility.
#' # Int Arch Occup Environ Health 1994; 66(1):49-54.
#' greenland <- matrix(c(45, 94, 257, 945), dimnames = list(c("BC+", "BC-"),
#' c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE)
#' set.seed(123)
#' # Exposure misclassification, non-differential
#' probsens(greenland, type = "exposure", reps = 20000,
#' seca = list("trapezoidal", c(.75, .85, .95, 1)),
#' spca = list("trapezoidal", c(.75, .85, .95, 1)))
#'
#' # Exposure misclassification, differential
#' probsens(greenland, type = "exposure", reps = 20000,
#' seca = list("trapezoidal", c(.75, .85, .95, 1)),
#' seexp = list("trapezoidal", c(.7, .8, .9, .95)),
#' spca = list("trapezoidal", c(.75, .85, .95, 1)),
#' spexp = list("trapezoidal", c(.7, .8, .9, .95)),
#' corr_se = .8,
#' corr_sp = .8)
#'
#' probsens(greenland, type = "exposure", reps = 20000,
#' seca = list("beta", c(908, 16)),
#' seexp = list("beta", c(156, 56)),
#' spca = list("beta", c(153, 6)),
#' spexp = list("beta", c(205, 18)),
#' corr_se = .8,
#' corr_sp = .8)
#'
#' probsens(matrix(c(338, 490, 17984, 32024),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "exposure",
#' reps = 1000,
#' seca = list("trapezoidal", c(.8, .9, .9, 1)),
#' spca = list("trapezoidal", c(.8, .9, .9, 1)))
#'
#' # Disease misclassification
#' probsens(matrix(c(173, 602, 134, 663),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "outcome",
#' reps = 20000,
#' seca = list("uniform", c(.8, 1)),
#' spca = list("uniform", c(.8, 1)))
#'
#' probsens(matrix(c(338, 490, 17984, 32024),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "outcome",
#' reps = 20000,
#' seca = list("uniform", c(.2, .6)),
#' seexp = list("uniform", c(.1, .5)),
#' spca = list("uniform", c(.99, 1)),
#' spexp = list("uniform", c(.99, 1)),
#' corr_se = .8,
#' corr_sp = .8)
#'
#' probsens(matrix(c(173, 602, 134, 663),
#' dimnames = list(c("BC+", "BC-"), c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE),
#' type = "outcome",
#' reps = 20000,
#' seca = list("beta", c(100, 5)),
#' seexp = list("beta", c(110, 10)),
#' spca = list("beta", c(120, 15)),
#' spexp = list("beta", c(130, 30)),
#' corr_se = .8,
#' corr_sp = .8)
#'
#' # Fox M.P., MacLehose R.F., Lash T.L.
#' # SAS and R code for probabilistic quantitative bias analysis for
#' # misclassified binary variables and binary unmeasured confounders
#' # Int J Epidemiol 2023:1624-1633.
#' fox <- matrix(c(40, 20, 60, 80),
#' dimnames = list(c("Diseased", "Non-diseased"), c("Exposed", "Unexposed")),
#' nrow = 2, byrow = TRUE)
#' set.seed(123)
#' probsens(fox, type = "exposure", reps = 20000,
#' seca = list("beta", c(25, 3)),
#' spca = list("trapezoidal", c(.9, .93, .97, 1)),
#' spexp = list("beta", c(47, 7)),
#' spexp = list("trapezoidal", c(.8, .83, .87, .9)),
#' corr_se = .8,
#' corr_sp = .8)
#' @export
#' @importFrom stats median pnorm qnorm quantile qunif runif qbeta rbeta
probsens <- function(case,
                     exposed,
                     type = c("exposure", "outcome"),
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

    if (is.null(seca[[2]]) | is.null(spca[[2]]))
        stop(cli::format_error(c("x" = "Missing argument(s) for seca or spca",
                                 "i" = "At least one Se and one Sp should be
provided through outcome parameters.")))
    if (!is.list(seca))
        stop(cli::format_error(c("i" = "Sensitivity of exposure classification among
those with the outcome should be a list.")))
    else seca <- seca
    if ((length(seca) != 2) | (length(spca) != 2))
        stop(cli::format_error(c("i" = "Check distribution parameters")))
    if ((!is.null(seexp) & length(seexp) != 2) | (!is.null(spexp) & length(spexp) != 2))
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
        stop(cli::format_error(c("i" = "Lower limit of your uniform distribution is
greater than upper limit.")))
    if (seca[[1]] == "triangular" & length(seca[[2]]) != 3)
        stop(cli::format_error(c("i" = "For triangular distribution, please provide
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
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if ((seca[[1]] == "constant" | seca[[1]] == "uniform" | seca[[1]] == "triangular" |
         seca[[1]] == "trapezoidal") & !all(seca[[2]] >= 0 & seca[[2]] <= 1))
        stop(cli::format_error(c("i" = "Sensitivity of exposure classification
among those with the outcome should be between 0 and 1.")))
    if (seca[[1]] == "beta" & length(seca[[2]]) != 2)
        stop(cli::format_error(c("x" = "For beta distribution, please provide alpha and beta.")))
    if (seca[[1]] == "beta" & (seca[[2]][1] < 0 | seca[[2]][2] < 0))
        stop(cli::format_error(c("i" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(seexp) & !is.list(seexp))
        stop(cli::format_error(c("x" = "Sensitivity of exposure classification
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
        stop(cli::format_error(c("i" = "For triangular distribution, please provide
vector of lower, upper limits, and mode.")))
    if (!is.null(seexp) && seexp[[1]] == "triangular" &&
        ((seexp[[2]][1] > seexp[[2]][3]) | (seexp[[2]][2] < seexp[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (!is.null(seexp) && seexp[[1]] == "trapezoidal" & length(seexp[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide
vector of lower limit, lower mode, upper mode, and upper limit.")))
    if (!is.null(seexp) && seexp[[1]] == "trapezoidal" &&
        ((seexp[[2]][1] > seexp[[2]][2]) | (seexp[[2]][2] > seexp[[2]][3]) |
         (seexp[[2]][3] > seexp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (!is.null(seexp) && seexp[[1]] == "normal" & (length(seexp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide vector of lower and upper bound limits, mean and SD.")))
    if (!is.null(seexp) && seexp[[1]] == "normal" &&
        ((seexp[[2]][1] >= seexp[[2]][2]) | (!all(seexp[[2]][1:2] >= 0 &
                                                  seexp[[2]][1:2] <= 1))))
        stop(cli::format_error(c("x" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if (!is.null(seexp) && (seexp[[1]] == "constant" | seexp[[1]] == "uniform" |
                            seexp[[1]] == "triangular" | seexp[[1]] == "trapezoidal") &
        !all(seexp[[2]] >= 0 & seexp[[2]] <= 1))
        stop(cli::format_error(c("x" = "Sensitivity of exposure classification among
those without the outcome should be between 0 and 1.")))
    if (!is.null(seexp) && seexp[[1]] == "beta" && length(seexp[[2]]) != 2)
        stop(cli::format_error(c("x" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(seexp) && seexp[[1]] == "beta" && (seexp[[2]][1] < 0 | seexp[[2]][2] < 0))
        stop(cli::format_error(c("i" = "Wrong arguments for your beta distribution.
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
        stop(cli::format_error(c("i" = "Lower limit of your uniform distribution
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
        stop(cli::format_error(c("i" = "For truncated normal distribution,
please provide sensible values for lower and upper bound limits (between 0
and 1; lower limit < upper limit).")))
    if ((spca[[1]] == "constant" | spca[[1]] == "uniform" |
         spca[[1]] == "triangular" |
         spca[[1]] == "trapezoidal") & !all(spca[[2]] >= 0 & spca[[2]] <= 1))
        stop(cli::format_error(c("x" = "Specificity of exposure classification
among those with the outcome should be between 0 and 1.")))
    if (spca[[1]] == "beta" & length(spca[[2]]) != 2)
        stop(cli::format_error(c("x" = "For beta distribution, please provide alpha and beta.")))
    if (spca[[1]] == "beta" & (spca[[2]][1] < 0 | spca[[2]][2] < 0))
        stop(cli::format_error(c("i" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(spexp) & !is.list(spexp))
        stop(cli::format_error(c("x" = "Specificity of exposure classification
among those without the outcome should be a list.")))
    else spexp <- spexp
    if (!is.null(spexp) && spexp[[1]] == "constant" & length(spexp[[2]]) != 1)
        stop(cli::format_error(c("x" = "For constant value, please provide a single value.")))
    if (!is.null(spexp) && spexp[[1]] == "uniform" & length(spexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide
vector of lower and upper limits.")))
    if (!is.null(spexp) && spexp[[1]] == "uniform" && spexp[[2]][1] >= spexp[[2]][2])
        stop(cli::format_error(c("i" = "Lower limit of your uniform distribution is
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
        ((spexp[[2]][1] > spexp[[2]][2]) | (spexp[[2]][2] > spexp[[2]][3]) |
         (spexp[[2]][3] > spexp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" & (length(spexp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution, please
provide vector of lower and upper bound limits, meand and SD.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" &&
        ((spexp[[2]][1] >= spexp[[2]][2]) |
         (!all(spexp[[2]][1:2] >= 0 & spexp[[2]][1:2] <= 1))))
        stop(cli::format_error(c("i" = "For truncated normal distribution, please
provide sensible values for lower and upper bound limits (between 0 and 1;
lower limit < upper limit).")))
    if (!is.null(spexp) && (spexp[[1]] == "constant" | spexp[[1]] == "uniform" |
                            spexp[[1]] == "triangular" | spexp[[1]] == "trapezoidal") &
        !all(spexp[[2]] >= 0 & spexp[[2]] <= 1))
        stop(cli::format_error(c("i" = "Specificity of exposure classification
among those without the outcome should be between 0 and 1.")))
    if (!is.null(spexp) && spexp[[1]] == "beta" && length(spexp[[2]]) != 2)
        stop(cli::format_error(c("x" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(spexp) && spexp[[1]] == "beta" && (spexp[[2]][1] < 0 | spexp[[2]][2] < 0))
        stop(cli::format_error(c("i" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(seexp) & (is.null(spca) | is.null(spexp) | is.null(corr_se) | is.null(corr_sp)))
        stop(cli::format_error(c("i" = "For differential misclassification type,
have to provide Se and Sp for among those with and without the outcome as well as
Se and Sp correlations.")))

    if (!is.null(corr_se) && (corr_se == 0 | corr_se == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))
    if (!is.null(corr_sp) && (corr_sp == 0 | corr_sp == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))

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

    cli::cli_alert_info("Calculating observed measures")
    obs_rr <- (a / (a + c)) / (b / (b + d))
    se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
    lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
    uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

    obs_or <- (a / b) / (c / d)
    se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
    uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

    draws <- matrix(NA, nrow = reps, ncol = 26)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A1", "B1", "C1", "D1",
                         "flag",
                         "prevca", "prevexp",
                         "ppvca", "ppvexp", "npvca", "npvexp",
                         "ab", "bb", "cb", "db",
                         "corr_RR", "corr_OR",
                         "rr_se_b", "or_se_b", "z",
                         "tot_RR", "tot_OR")
    corr_draws <- matrix(NA, nrow = reps, ncol = 4)

    se1 <- c(reps, seca[[2]])
    se0 <- c(reps, seexp[[2]])
    sp1 <- c(reps, spca[[2]])
    sp0 <- c(reps, spexp[[2]])

    ## Step3: Assign probability distributions to each bias parameter
    ## and Step 4a draw Se's and Sp's
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
            draws[, 1] <- do.call(triangle::qtriangle,
                                  c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "beta") {
            draws[, 1] <- do.call(qbeta, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seexp[[1]] == "uniform") {
            draws[, 2] <- do.call(qunif, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::qtriangle,
                                  c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "normal") {
            draws[, 2] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "beta") {
            draws[, 2] <- do.call(qbeta, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (spca[[1]] == "uniform") {
            draws[, 3] <- do.call(qunif, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::qtriangle,
                                  c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "normal") {
            draws[, 3] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "beta") {
            draws[, 3] <- do.call(qbeta, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spexp[[1]] == "uniform") {
            draws[, 4] <- do.call(qunif, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "triangular") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "trapezoidal") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "normal") {
            draws[, 4] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "beta") {
            draws[, 4] <- do.call(qbeta, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
    }

    type <- match.arg(type)
    if (type == "exposure") {
        ## Step 4b: Bias-adjusted cell frequencies using simple bias analysis
        ## methods and the sampled bias parameters
        cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) / (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (a + b) - draws[, 5]
        draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) / (draws[, 2] - (1 - draws[, 4]))
        draws[, 8] <- (c + d) - draws[, 7]

        ## Prevalence of exposure in cases and controls, accounting for sampling error
        suppressWarnings({
                             draws[, 10] <- rbeta(reps, draws[, 5], draws[, 6])
                             draws[, 11] <- rbeta(reps, draws[, 7], draws[, 8])
                         })
        ## PPV and NPV of exposure classification in cases and controls
        draws[, 12] <- (draws[, 1] * draws[, 10]) /
            ((draws[, 1] * draws[, 10]) + (1 - draws[, 3]) * (1 - draws[, 10]))
        draws[, 13] <- (draws[, 2] * draws[, 11]) /
            ((draws[, 2] * draws[, 11]) + (1 - draws[, 4]) * (1 - draws[, 11]))
        draws[, 14] <- (draws[, 3] * (1 - draws[, 10])) /
            ((1 - draws[, 1]) * draws[, 10] + draws[, 3] * (1 - draws[, 10]))
        draws[, 15] <- (draws[, 4] * (1 - draws[, 11])) /
            ((1 - draws[, 2]) * draws[, 11] + draws[, 4] * (1 - draws[, 11]))
        ## Expected number of exposed cases and controls
        suppressWarnings(draws[, 16] <- rbinom(reps, a, draws[, 12]) +
                             rbinom(reps, b, 1 - draws[, 14]))
        draws[, 17] <- (a + b) - draws[, 16]
        suppressWarnings(draws[, 18] <- rbinom(reps, c, draws[, 13]) +
                             rbinom(reps, d, 1 - draws[, 15]))
        draws[, 19] <- (c + d) - draws[, 18]

        ## Bias-adjusted RR and OR with second source of uncertainty
        draws[, 20] <- (draws[, 16] / (draws[, 16] + draws[, 18])) /
            (draws[, 17] / (draws[, 17] + draws[, 19]))
        draws[, 21] <- (draws[, 16] / draws[, 18]) / (draws[, 17] / draws[, 19])

        ## Step 4c: Incorporate conventional random error by sampling summary
        ## statistics
        ## Calculate bias-adjusted RR and OR, third source of uncertainty,
        ## bias-adjusted SE
        cli::cli_progress_step("Incorporating random error", spinner = TRUE)
        draws[, 22] <- sqrt(1 / draws[, 16] + 1 / draws[, 17] -
                            1 / (draws[, 16] + draws[, 18]) -
                            1 / (draws[, 17] + draws[, 19]))
        draws[, 23] <- sqrt((1 / draws[, 16]) + (1 / draws[, 17]) +
                            (1 / draws[, 18]) + (1 / draws[, 19]))
        draws[, 24] <- rnorm(reps)
        draws[, 25] <- exp(log(draws[, 20]) - (draws[, 24] * draws[, 22]))
        draws[, 26] <- exp(log(draws[, 21]) - (draws[, 24] * draws[, 23]))

        ## Clean up
        draws[, 9] <- apply(draws[, c(5:8, 16:19)], MARGIN = 1, function(x) sum(x > 0))
        draws[, 9] <- ifelse(draws[, 9] != 8 | is.na(draws[, 9]), NA, 1)
        discard <- sum(is.na(draws[, 9]))
        if (sum(is.na(draws[, 9])) > 0) {
            cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
            neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
        } else neg_warn <- NULL

        draws <- draws[draws[, 9] == 1 & !is.na(draws[, 9]), ]

        rr_syst <- c(median(draws[, 20], na.rm = TRUE),
                     quantile(draws[, 20], probs = .025, na.rm = TRUE),
                     quantile(draws[, 20], probs = .975, na.rm = TRUE))
        or_syst <- c(median(draws[, 21], na.rm = TRUE),
                     quantile(draws[, 21], probs = .025, na.rm = TRUE),
                     quantile(draws[, 21], probs = .975, na.rm = TRUE))
        rr_tot <- c(median(draws[, 25], na.rm = TRUE),
                    quantile(draws[, 25], probs = .025, na.rm = TRUE),
                    quantile(draws[, 25], probs = .975, na.rm = TRUE))
        or_tot <- c(median(draws[, 26], na.rm = TRUE),
                    quantile(draws[, 26], probs = .025, na.rm = TRUE),
                    quantile(draws[, 26], probs = .975, na.rm = TRUE))

        if (!inherits(case, "episensr.probsens")) {
            tab <- tab
            rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr),
                          c(obs_or, lci_obs_or, uci_obs_or))
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
        rmatc <- rbind(rr_syst, rr_tot, or_syst, or_tot)
        rownames(rmatc) <- c("Relative Risk -- systematic error:",
                             "                      total error:",
                             "   Odds Ratio -- systematic error:",
                             "                      total error:")
        colnames(rmatc) <- c("Median", "p2.5", "p97.5")

        cli::cli_progress_update()
    }

    if (type == "outcome") {
        ## Step 4b: Bias-adjusted cell frequencies using simple bias analysis
        ## methods and the sampled bias parameters
        cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + c)) / (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (b - (1 - draws[, 4]) * (b + d)) / (draws[, 2] - (1 - draws[, 4]))
        draws[, 7] <- (a + c) - draws[, 5]
        draws[, 8] <- (b + d) - draws[, 6]

        ## Prevalence of outcome in cases and controls, accounting for sampling error
        suppressWarnings({
                             draws[, 10] <- rbeta(reps, draws[, 5], draws[, 7])
                             draws[, 11] <- rbeta(reps, draws[, 6], draws[, 8])
                         })
        ## PPV and NPV of outcome classification in cases and controls
        draws[, 12] <- (draws[, 1] * draws[, 10]) /
            ((draws[, 1] * draws[, 10]) + (1 - draws[, 3]) * (1 - draws[, 10]))
        draws[, 13] <- (draws[, 2] * draws[, 11]) /
            ((draws[, 2] * draws[, 11]) + (1 - draws[, 4]) * (1 - draws[, 11]))
        draws[, 14] <- (draws[, 3] * (1 - draws[, 10])) /
            ((1 - draws[, 1]) * draws[, 10] + draws[, 3] * (1 - draws[, 10]))
        draws[, 15] <- (draws[, 4] * (1 - draws[, 11])) /
            ((1 - draws[, 2]) * draws[, 11] + draws[, 4] * (1 - draws[, 11]))
        ## Expected number of cases among exposed and unexposed
        suppressWarnings(draws[, 16] <- rbinom(reps, a, draws[, 12]) +
                             rbinom(reps, c, 1 - draws[, 14]))
        suppressWarnings(draws[, 17] <- rbinom(reps, b, draws[, 13]) +
                             rbinom(reps, d, 1 - draws[, 15]))
        draws[, 18] <- (a + c) - draws[, 16]
        draws[, 19] <- (b + d) - draws[, 17]

        ## Bias-adjusted RR and OR with second source of uncertainty
        draws[, 20] <- (draws[, 16] / (draws[, 16] + draws[, 18])) /
            (draws[, 17] / (draws[, 17] + draws[, 19]))
        draws[, 21] <- (draws[, 16] / draws[, 18]) / (draws[, 17] / draws[, 19])

        ## Step 4c: Incorporate conventional random error by sampling summary
        ## statistics
        ## Calculate bias-adjusted RR and OR, third source of uncertainty,
        ## bias-adjusted SE
        cli::cli_progress_step("Incorporating random error", spinner = TRUE)
        draws[, 22] <- sqrt(1 / draws[, 16] + 1 / draws[, 17] -
                            1 / (draws[, 16] + draws[, 18]) -
                            1 / (draws[, 17] + draws[, 19]))
        draws[, 23] <- sqrt((1 / draws[, 16]) + (1 / draws[, 17]) +
                            (1 / draws[, 18]) + (1 / draws[, 19]))
        draws[, 24] <- rnorm(reps)
        draws[, 25] <- exp(log(draws[, 20]) - (draws[, 24] * draws[, 22]))
        draws[, 26] <- exp(log(draws[, 21]) - (draws[, 24] * draws[, 23]))

        ## Clean up
        draws[, 9] <- apply(draws[, c(5:8, 16:19)], MARGIN = 1, function(x) sum(x > 0))
        draws[, 9] <- ifelse(draws[, 9] != 8 | is.na(draws[, 9]), NA, 1)
        discard <- sum(is.na(draws[, 9]))
        if (sum(is.na(draws[, 9])) > 0) {
            cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
            neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
        } else neg_warn <- NULL

        draws <- draws[draws[, 9] == 1 & !is.na(draws[, 9]), ]

        rr_syst <- c(median(draws[, 20], na.rm = TRUE),
                     quantile(draws[, 20], probs = .025, na.rm = TRUE),
                     quantile(draws[, 20], probs = .975, na.rm = TRUE))
        or_syst <- c(median(draws[, 21], na.rm = TRUE),
                     quantile(draws[, 21], probs = .025, na.rm = TRUE),
                     quantile(draws[, 21], probs = .975, na.rm = TRUE))
        rr_tot <- c(median(draws[, 25], na.rm = TRUE),
                    quantile(draws[, 25], probs = .025, na.rm = TRUE),
                    quantile(draws[, 25], probs = .975, na.rm = TRUE))
        or_tot <- c(median(draws[, 26], na.rm = TRUE),
                    quantile(draws[, 26], probs = .025, na.rm = TRUE),
                    quantile(draws[, 26], probs = .975, na.rm = TRUE))

        if(!inherits(case, "episensr.probsens")) {
            tab <- tab
            rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr),
                          c(obs_or, lci_obs_or, uci_obs_or))
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
        rmatc <- rbind(rr_syst, rr_tot, or_syst, or_tot)
        rownames(rmatc) <- c("Relative Risk -- systematic error:",
                             "                      total error:",
                             "   Odds Ratio -- systematic error:",
                             "                      total error:")
        colnames(rmatc) <- c("Median", "p2.5", "p97.5")

        cli::cli_progress_update()
    }
    res <- list(obs.data = tab,
                obs.measures = rmat,
                adj.measures = rmatc,
                sim.df = as.data.frame(draws[, -9]),
                reps = reps,
                fun = "probsens",
                warnings = neg_warn#,
#                message = discard_mess
                )
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
