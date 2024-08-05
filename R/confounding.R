#' Uncontrolled confounding
#'
#' `confounders()` and `probsens_conf()` allow to provide adjusted measures of
#' association corrected for unknown or unmeasured confounding without effect
#' modification.
#'
#' @section Simple bias analysis with `confounders()`:
#' `confounders()` allows you to run a simple sensitivity analysis to correct for
#' unknown or unmeasured confounding without effect modification. Implementation
#' for ratio measures (relative risk -- RR, or odds ratio -- OR) and difference
#' measures (risk difference -- RD).
#'
#' The analytic approach uses the "relative risk due to confounding" as defined by
#' Miettinen (1972), i.e. \eqn{RR_{adj} = \frac{RR_{crude}}{RR_{conf}}} where RR_adj
#' is the standardized (adjusted) risk ratio, RR_crude is the crude risk ratio, and
#' RR_conf is the relative risk component attributable to confounding by the
#' stratification factors. The output provides both RR_adj (SMR or Mantel-Haenszel)
#' and the RR_conf (i.e., RR, OR or RD due to confounding from the unmeasured
#' confounder).
#'
#' @section Probabilistic sensitivity analysis with `probsens_conf()`:
#' `probsens_conf()` performs a summary-level probabilistic sensitivity analysis
#' to correct for unknown or unmeasured confounding and random error simultaneously.
#' It returns the Mantel-Haenszel risk ratio.
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
#' @param type Choice of implementation, with no effect measure modification for
#'   ratio measures (relative risk -- RR; odds ratio -- OR) or difference measures
#'   (risk difference -- RD).
#' @param bias_parms Numeric vector defining the 3 necessary bias parameters.
#'   This vector has 3 elements, in the following order:
#'   \enumerate{
#'   \item the association between the confounder and the outcome among those who
#'   were not exposed (RR, OR, or RD according to choice of implementation),
#'   \item the prevalence of the confounder among the exposed (between 0 and 1), and
#'   \item the prevalence of the confounder among the unexposed (between 0 and 1).
#'   }
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{cfder_data}{The same table for Confounder +.}
#' \item{nocfder_data}{The same table for Confounder -.}
#' \item{obs_measures}{A table of relative risk with confidence intervals; for
#'   Total, Confounder +, and Confounder -.}
#' \item{adj_measures}{A table of Standardized Morbidity Ratio and Mantel-Haenszel
#'   estimates.}
#' \item{bias_parms}{Input bias parameters.}
#'
#' @family confounding
#'
#' @references
#' Fox, M.P, MacLehose, R.F., Lash, T.L., 2021 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.105--140, 256--262, Springer.
#'
#' Miettinen, 1971. Components of the Crude Risk Ratio. \emph{Am J Epidemiol}
#' 96(2):168-172.
#'
#' Li, S.T., Hammond, J.L., 1975. \emph{Generation of Pseudorandom Numbers
#' with Specified Univariate Distributions and Correlation Coefficients}.
#' IEEE Trans Syst Man Cybern 5:557-561.
#'
#' @examples
#' # The data for this example come from:
#' # Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O.
#' # et al.
#' # Increased risk of infection with human immunodeficiency virus type 1 among
#' # uncircumcised men presenting with genital ulcer disease in Kenya.
#' # Clin Infect Dis 1996;23:449-53.
#' confounders(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RR",
#' bias_parms = c(.63, .8, .05))
#'
#' confounders(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "OR",
#' bias_parms = c(.63, .8, .05))
#'
#' confounders(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RD",
#' bias_parms = c(-.37, .8, .05))
#' @export
#' @importFrom stats qnorm
confounders <- function(case,
                        exposed,
                        type = c("RR", "OR", "RD"),
                        bias_parms = NULL,
                        alpha = 0.05) {
    if (length(type) > 1)
        stop(cli::format_error(c("i" = "Choose between RR, OR, or RD implementation.")))

    if (is.null(bias_parms))
        bias_parms <- c(1, 0, 0)
    else bias_parms <- bias_parms
    if (length(bias_parms) != 3)
        stop(cli::format_error(c("i" = "The argument bias_parms should be made of
the following components: (1) Association between the confounder and the outcome
among those who were not exposed, (2) Prevalence of the confounder among the exposed,
and (3) Prevalence of the confounder among the unexposed.")))
    if (!all(bias_parms[-1] >= 0 & bias_parms[-1] <= 1))
        stop(cli::format_error(c("x" = "Prevalences should be between 0 and 1.")))
    if (bias_parms[1] <= 0 & type != "RD")
        stop(cli::format_error(c("x" = "Association between the confounder and the
outcome among those who were not exposed should be greater than 0.")))

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

    type <- match.arg(type)
    if (type == "RR") {
        crude_rr <- (a / (a + c)) / (b / (b + d))
        se_log_crude_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
        lci_crude_rr <- exp(log(crude_rr) - qnorm(1 - alpha / 2) * se_log_crude_rr)
        uci_crude_rr <- exp(log(crude_rr) + qnorm(1 - alpha / 2) * se_log_crude_rr)

        M1 <- (a + c) * bias_parms[2]
        N1 <- (b + d) * bias_parms[3]
        A1 <- (bias_parms[1] * M1 * a) / (bias_parms[1] * M1 + (a + c) - M1)
        B1 <- (bias_parms[1] * N1 * b) / (bias_parms[1] * N1 + (b + d) - N1)
        C1 <- M1 - A1
        D1 <- N1 - B1
        M0 <- a + c - M1
        N0 <- b + d - N1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1

        if (A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 | A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop(cli::format_error(c("x" = "Parameters chosen lead to negative cell(s)
in adjusted 2x2 table(s).")))

        tab_cfder <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab_nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRrr <- a / ((M1 * B1 / N1) + (M0 * B0 / N0))
        MHrr <- (A1 * N1 / (M1 + N1) + A0 * N0 / (M0 + N0)) / (B1 * M1 / (M1 + N1) + B0 * M0 / (M0 + N0))
        cfder_rr <- (A1 / (A1 + C1)) / (B1 / (B1 + D1))
        nocfder_rr <- (A0 / (A0 + C0)) / (B0 / (B0 + D0))
        RRadj_smr <- crude_rr / SMRrr
        RRadj_mh <- crude_rr / MHrr

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(tab_cfder) <- paste("Row", 1:2)
        } else {
            rownames(tab_cfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(tab_cfder) <- paste("Col", 1:2)
        } else {
            colnames(tab_cfder) <- colnames(tab)
        }
        if (is.null(rownames(tab))) {
            rownames(tab_nocfder) <- paste("Row", 1:2)
        } else {
            rownames(tab_nocfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(tab_nocfder) <- paste("Col", 1:2)
        } else {
            colnames(tab_nocfder) <- colnames(tab)
        }
        rmat <- rbind(c(crude_rr, lci_crude_rr, uci_crude_rr))
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
        rmatc <- rbind(c(SMRrr, RRadj_smr), c(MHrr, RRadj_mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio:",
                             "             Mantel-Haenszel:")
        colnames(rmatc) <- c(" ", "RR due to confounding")
        rmat <- rbind(rmat, c(cfder_rr, NA, NA), c(nocfder_rr, NA, NA))
        rownames(rmat) <- c("        Crude Relative Risk:",
                            "Relative Risk, Confounder +:",
                            "Relative Risk, Confounder -:")
    }

    if (type == "OR") {
        crude_or <- (a / b) / (c / d)
        se_log_crude_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        lci_crude_or <- exp(log(crude_or) - qnorm(1 - alpha / 2) * se_log_crude_or)
        uci_crude_or <- exp(log(crude_or) + qnorm(1 - alpha / 2) * se_log_crude_or)

        C1 <- c * bias_parms[2]
        D1 <- d * bias_parms[3]
        A1 <- (bias_parms[1] * C1 * a) / (bias_parms[1] * C1 + c - C1)
        B1 <- (bias_parms[1] * D1 * b) / (bias_parms[1] * D1 + d - D1)
        M1 <- A1 + C1
        N1 <- B1 + D1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1
        M0 <- A0 + C0
        N0 <- B0 + C0

        if (A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 | A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop(cli::format_error(c("x" = "Parameters chosen lead to negative cell(s)
in adjusted 2x2 table(s).")))

        tab_cfder <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab_nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRor <- a / ((C1 * B1 / D1) + (C0 * B0 / D0))
        MHor <- (A1 * D1 / (M1 + N1) + A0 * D0 / (M0 + N0)) / (B1 * C1 / (M1 + N1) + B0 * C0 / (M0 + N0))
        cfder_or <- (A1 / C1) / (B1 / D1)
        nocfder_or <- (A0 / C0) / (B0 / D0)
        ORadj_smr <- crude_or / SMRor
        ORadj_mh <- crude_or / MHor

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(tab_cfder) <- paste("Row", 1:2)
        } else {
            rownames(tab_cfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(tab_cfder) <- paste("Col", 1:2)
        } else {
            colnames(tab_cfder) <- colnames(tab)
        }
        if (is.null(rownames(tab))) {
            rownames(tab_nocfder) <- paste("Row", 1:2)
        } else {
            rownames(tab_nocfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(tab_nocfder) <- paste("Col", 1:2)
        } else {
            colnames(tab_nocfder) <- colnames(tab)
        }
        rmat <- rbind(c(crude_or, lci_crude_or, uci_crude_or))
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
        rmatc <- rbind(c(SMRor, ORadj_smr), c(MHor, ORadj_mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio:",
                             "             Mantel-Haenszel:")
        colnames(rmatc) <- c(" ", "OR due to confounding")
        rmat <- rbind(rmat, c(cfder_or, NA, NA), c(nocfder_or, NA, NA))
        rownames(rmat) <- c("        Crude Odds Ratio:",
                            "Odds Ratio, Confounder +:",
                            "Odds Ratio, Confounder -:")
    }

    if (type == "RD") {
        crude_rd <- (a / (a + c)) - (b / (b + d))
        se_log_crude_rd <- sqrt((a * c) / (a + c)^3 + (b * d) / (b + d)^3)
        lci_crude_rd <- crude_rd - qnorm(1 - alpha / 2) * se_log_crude_rd
        uci_crude_rd <- crude_rd + qnorm(1 - alpha / 2) * se_log_crude_rd

        M1 <- (a + c) * bias_parms[2]
        N1 <- (b + d) * bias_parms[3]
        M0 <- (a + c) - M1
        N0 <- (b + d) - N1
        A1 <- (bias_parms[1] * M1 * M0 + M1 * a) / (a + c)
        B1 <- (bias_parms[1] * N1 * N0 + N1 * b) / (b + d)
        C1 <- M1 - A1
        D1 <- N1 - B1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1

        if (A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 | A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop(cli::format_error(c("x" = "Parameters chosen lead to negative
cell(s) in adjusted 2x2 table(s).")))

        tab_cfder <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab_nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        MHrd <- (((A1 * N1 - B1 * M1) / (M1 + N1)) + ((A0 * N0 - B0 * M0) / (M0 + N0))) /
            ((M1 * N1 / (M1 + N1)) + (M0 * N0 / (M0 + N0)))
        cfder_rd <- (A1 / M1) - (B1 / N1)
        nocfder_rd <- (A0 / M0) - (B0 / N0)
        RDadj_mh <- crude_rd - MHrd

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(tab_cfder) <- paste("Row", 1:2)
        } else {
            rownames(tab_cfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(tab_cfder) <- paste("Col", 1:2)
        } else {
            colnames(tab_cfder) <- colnames(tab)
        }
        if (is.null(rownames(tab))) {
            rownames(tab_nocfder) <- paste("Row", 1:2)
        } else {
            rownames(tab_nocfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(tab_nocfder) <- paste("Col", 1:2)
        } else {
            colnames(tab_nocfder) <- colnames(tab)
        }
        rmat <- rbind(c(crude_rd, lci_crude_rd, uci_crude_rd))
        colnames(rmat) <- c("     ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
        rmatc <- rbind(c(MHrd, RDadj_mh))
        rownames(rmatc) <- "Mantel-Haenszel:"
        colnames(rmatc) <- c(" ", "RD due to unmeasured confounder")
        rmat <- rbind(rmat, c(cfder_rd, NA, NA), c(nocfder_rd, NA, NA))
        rownames(rmat) <- c("        Crude Risk Difference:",
                            "Risk Difference, Confounder +:",
                            "Risk Difference, Confounder -:")
    }
    res <- list(obs_data = tab,
                cfder_data = tab_cfder,
                nocfder_data = tab_nocfder,
                obs_measures = rmat,
                adj_measures = rmatc,
                bias_parms = bias_parms)
    class(res) <- c("episensr", "list")
    res
}

#' @rdname confounders
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param reps Number of replications to run.
#' @param prev_exp List defining the prevalence of exposure among the exposed.
#'   The first argument provides the probability distribution function (constant,
#'   uniform, triangular, trapezoidal, truncated normal, or beta) and the second
#'   its parameters as a vector. Lower bound of the truncated normal cannot be
#'   less than zero. Upper bound is Inf by default.
#'   \enumerate{
#'   \item constant: constant value,
#'   \item uniform: min, max,
#'   \item triangular: lower limit, upper limit, mode,
#'   \item trapezoidal: min, lower mode, upper mode, max.
#'   \item normal: lower bound, upper bound, mean, sd.
#'   \item beta: alpha, beta.
#'   }
#' @param prev_nexp List defining the prevalence of exposure among the unexposed.
#' @param risk List defining the confounder-disease relative risk or the
#'   confounder-exposure odds ratio. The first argument provides the probability
#'   distribution function (constant, uniform, triangular, trapezoidal,
#'   log-logistic, or log-normal) and the second its parameters as a vector:
#'   \enumerate{
#'   \item constant: constant value,
#'   \item uniform: min, max,
#'   \item triangular: lower limit, upper limit, mode,
#'   \item trapezoidal: min, lower mode, upper mode, max.
#'   \item log-logistic: shape, rate. Must be strictly positive,
#'   \item log-normal: meanlog, sdlog. This is the mean and standard deviation on the log scale.
#'   }
#' @param corr_p Correlation between the exposure-specific confounder prevalences.
#' @param alpha Significance level.
#'
#' @return A list with elements (for `probsens_conf()`):
#'   \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#'   \item{obs_measures}{A table of observed relative risk and odds ratio with confidence intervals.}
#'   \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#'   \item{sim_df}{Data frame of random parameters and computed values.}
#'   \item{reps}{Number of replications.}
#'
#'@examples
#' # The data for this example come from:
#' # Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O. et al.
#' # Increased risk of infection with human immunodeficiency virus type 1 among
#' # uncircumcised men presenting with genital ulcer disease in Kenya.
#' # Clin Infect Dis 1996;23:449-53.
#' tyndall <- matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")), nrow = 2, byrow = TRUE)
#' set.seed(1234)
#' probsens_conf(tyndall, reps = 100000,
#' prev_exp = list("trapezoidal", c(.7, .75, .85, .9)),
#' prev_nexp = list("trapezoidal", c(.03, .04, .07, .1)),
#' risk = list("trapezoidal", c(.5, .6, .7, .8)))
#'
#' set.seed(123)
#' probsens_conf(tyndall, reps = 20000,
#' prev_exp = list("beta", c(200, 56)),
#' prev_nexp = list("beta", c(10, 16)),
#' risk = list("triangular", c(.6, .7, .63)),
#' corr_p = .8)
#'
#' set.seed(123)
#' probsens_conf(tyndall, reps = 20000,
#' prev_exp = list("normal", c(.01, .12, 0.03, 0.005)),
#' prev_nexp = list("normal", c(0, Inf, 0.01, 0.0001)),
#' risk = list("triangular", c(.6, .7, .63)), corr_p = .8)
#'
#' # Fox M.P., MacLehose R.F., Lash T.L.
#' # SAS and R code for probabilistic quantitative bias analysis for
#' # misclassified binary variables and binary unmeasured confounders
#' # Int J Epidemiol 2023:1624-1633.
#' fox <- matrix(c(40, 20, 60, 80),
#' dimnames = list(c("Diseased", "Non-diseased"), c("Exposed", "Unexposed")),
#' nrow = 2, byrow = TRUE)
#' set.seed(1234)
#' probsens_conf(fox, reps = 10^5,
#' prev_exp = list("beta", c(10, 20)),
#' prev_nexp = list("beta", c(5, 20)),
#' risk = list("trapezoidal", c(1.5, 1.7, 2.3, 2.5)))
#'
#' set.seed(1234)
#' probsens_conf(fox, reps = 20000,
#' prev_exp = list("beta", c(10, 20)),
#' prev_nexp = list("beta", c(5, 20)),
#' risk = list("log-normal", c(log(2), .23)))
#' @export
#' @importFrom stats median qnorm quantile runif rlnorm rbeta qbeta
probsens_conf <- function(case,
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
                          alpha = 0.05) {
    if (reps < 1)
        stop(cli::format_error(c("x" = "Wrong number of replications: reps = {reps}",
                                 "i" = "reps must be >= 1")))

    if (is.null(prev_exp) | is.null(prev_nexp))
        stop(cli::format_error(c("x" = "Missing argument(s) for prev_exp or prev_nexp",
                                 "i" = "Please provide prevalences among the exposed
and unexposed.")))
    if (is.null(risk))
        stop(cli::format_error(c("x" = "Please provide risk of acquiring outcome.")))

    if (!is.list(prev_exp))
        stop(cli::format_error(c("i" = "Prevalence of exposure among the exposed
should be a list.")))
    else prev_exp <- prev_exp
    if ((length(prev_exp) != 2) | (length(prev_nexp) != 2) | (length(risk) != 2))
        stop(cli::format_error(c("i" = "Check distribution parameters.")))
    if ((length(prev_exp[[1]]) != 1) | (length(prev_nexp[[1]]) != 1) | (length(risk[[1]]) != 1))
        stop(cli::format_error(c("x" = "Which distribution?")))
    if (!is.null(corr_p) && (prev_exp[[1]] == "constant" | prev_nexp[[1]] == "constant"))
        stop(cli::format_error(c("x" = "No correlated distributions with constant values.")))
    if (prev_exp[[1]] == "constant" & length(prev_exp[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (prev_exp[[1]] == "uniform" & length(prev_exp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide vector of
lower and upper limits.")))
    if (prev_exp[[1]] == "uniform" & prev_exp[[2]][1] >= prev_exp[[2]][2])
        stop(cli::format_error(c("i" = "Lower limit of your uniform distribution is greater than
upper limit.")))
    if (prev_exp[[1]] == "triangular" & length(prev_exp[[2]]) != 3)
        stop(cli::format_error(c("i" = "For triangular distribution, please provide vector of
lower, upper limits, and mode.")))
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
        stop(cli::format_error(c("i" = "For truncated normal distribution, please
provide vector of lower and upper bound limits, mean and sd.")))
    if (prev_exp[[1]] == "normal" & length(prev_exp[[2]]) == 4 &
        ((prev_exp[[2]][1] >= prev_exp[[2]][2]) | (prev_exp[[2]][1] < 0)))
        stop(cli::format_error(c("i" = "For truncated normal distribution, please provide
sensible values for lower and upper bound limits (lower limit >= 0; lower limit < upper limit).")))
    if ((prev_exp[[1]] == "constant" | prev_exp[[1]] == "uniform" |
         prev_exp[[1]] == "triangular" | prev_exp[[1]] == "trapezoidal") &
        !all(prev_exp[[2]] >= 0 & prev_exp[[2]] <= 1))
        stop(cli::format_error(c("x" = "Prevalence should be between 0 and 1.")))
    if (!is.null(prev_exp) && prev_exp[[1]] == "beta" && length(prev_exp[[2]]) != 2)
        stop(cli::format_error(c("x" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(prev_exp) && prev_exp[[1]] == "beta" &&
       (prev_exp[[2]][1] < 0 | prev_exp[[2]][2] < 0))
        stop(cli::format_error(c("i" = "Wrong arguments for your beta distribution. Alpha and
Beta should be > 0.")))

    if (!is.list(prev_nexp))
        stop(cli::format_error(c("x" = "Prevalence of exposure among the non-exposed
should be a list.")))
    else prev_nexp <- prev_nexp
    if (prev_nexp[[1]] == "constant" & length(prev_nexp[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (prev_nexp[[1]] == "uniform" & length(prev_nexp[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide vector of
lower and upper limits.")))
    if (prev_nexp[[1]] == "uniform" & prev_nexp[[2]][1] >= prev_nexp[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is greater
than upper limit.")))
    if (prev_nexp[[1]] == "triangular" & length(prev_nexp[[2]]) != 3)
        stop(cli::format_error(c("i" = "For triangular distribution, please provide vector of
lower, upper limits, and mode.")))
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
        stop(cli::format_error(c("i" = "For truncated normal distribution, please provide
vector of lower and upper bound limits, mean and sd.")))
    if (prev_nexp[[1]] == "normal" & length(prev_nexp[[2]]) == 4 &
        ((prev_nexp[[2]][1] >= prev_nexp[[2]][2]) | (prev_nexp[[2]][1] < 0)))
        stop(cli::format_error(c("i" = "For truncated normal distribution, please provide
sensible values for lower and upper bound limits (lower limit >= 0; lower limit < upper limit).")))
    if ((prev_nexp[[1]] == "constant" | prev_nexp[[1]] == "uniform" |
         prev_nexp[[1]] == "triangular" | prev_nexp[[1]] == "trapezoidal") &
        !all(prev_nexp[[2]] >= 0 & prev_nexp[[2]] <= 1))
        stop(cli::format_error(c("x" = "Prevalence should be between 0 and 1.")))
    if (!is.null(prev_nexp) && prev_nexp[[1]] == "beta" && length(prev_nexp[[2]]) != 2)
        stop(cli::format_error(c("x" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(prev_nexp) && prev_nexp[[1]] == "beta" && (prev_nexp[[2]][1] < 0 | prev_nexp[[2]][2] < 0))
        stop(cli::format_error(c("i" = "Wrong arguments for your beta distribution. Alpha
and Beta should be > 0.")))

    if (!is.list(risk))
        stop(cli::format_error(c("x" = "Risk should be a list.")))
    else risk <- risk
    if (risk[[1]] == "constant" & length(risk[[2]]) != 1)
        stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
    if (risk[[1]] == "uniform" & length(risk[[2]]) != 2)
        stop(cli::format_error(c("i" = "For uniform distribution, please provide vector of
lower and upper limits.")))
    if (risk[[1]] == "uniform" & risk[[2]][1] >= risk[[2]][2])
        stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is greater
than upper limit.")))
    if (risk[[1]] == "triangular" & length(risk[[2]]) != 3)
        stop(cli::format_error(c("x" = "For triangular distribution, please provide vector of
lower, upper limits, and mode.")))
    if (risk[[1]] == "triangular" & ((risk[[2]][1] > risk[[2]][3]) | (risk[[2]][2] < risk[[2]][3])))
        stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
    if (risk[[1]] == "trapezoidal" & length(risk[[2]]) != 4)
        stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide vector of
lower limit, lower mode, upper mode, and upper limit.")))
    if (risk[[1]] == "trapezoidal" & ((risk[[2]][1] > risk[[2]][2]) |
                                      (risk[[2]][2] > risk[[2]][3]) |
                                      (risk[[2]][3] > risk[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (risk[[1]] == "log-logistic" & length(risk[[2]]) != 2)
        stop(cli::format_error(c("i" = "For log-logistic distribution, please provide vector
of location and scale.")))
    if (risk[[1]] == "log-normal" & length(risk[[2]]) != 2)
        stop(cli::format_error(c("i" = "For log-normal distribution, please provide vector of
meanlog and sdlog.")))

    if (!is.null(corr_p) && (corr_p == 0 | corr_p == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))

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

    cli::cli_alert_info("Calculating observed measures")
    obs_rr <- (a / (a + c)) / (b / (b + d))
    se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
    lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
    uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

    obs_or <- (a / b) / (c / d)
    se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
    uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

    draws <- matrix(NA, nrow = reps, ncol = 21)
    colnames(draws) <- c("p1", "p0", "RR_cd",
                         "A1", "B1", "C1", "D1",
                         "A1_2", "B1_2", "C1_2", "D1_2",
                         "A0_2", "B0_2", "C0_2", "D0_2",
                         "flag",
                         "corr_RR", "corr_OR",
                         "var",
                         "tot_RR", "tot_OR")
    corr_draws <- matrix(NA, nrow = reps, ncol = 2)

    p1 <- c(reps, prev_exp[[2]])
    p0 <- c(reps, prev_nexp[[2]])
    rr_cd <- c(reps, risk[[2]])

    ## Step3: Assign probability distributions to each bias parameter
    ## and Step 4a draw p_1, p_0, and RR_CD
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

    ## Step 4b: Bias-adjusted cell date using simple bias analysis methods and
    ## the sampled bias parameters
    cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
    draws[, 4] <- (draws[, 3] * ((a + c) * draws[, 1]) * a) /
        (draws[, 3] * ((a + c) * draws[, 1]) + (a + c) - ((a + c) * draws[, 1]))
    draws[, 5] <- (draws[, 3] * b * (draws[, 2] * (b + d))) /
        (draws[, 3] * (draws[, 2] * (b + d)) + (b + d) - (draws[, 2] * (b + d)))
    draws[, 6] <- ((a + c) * draws[, 1]) - draws[, 4]
    draws[, 7] <- ((b + d) * draws[, 2]) - draws[, 5]

    ## Randomly impute confounder status
    suppressWarnings(draws[, 8] <- rbinom(reps, a, draws[, 4] / a))
    suppressWarnings(draws[, 9] <- rbinom(reps, b, draws[, 5] / b))
    suppressWarnings(draws[, 10] <- rbinom(reps, c, draws[, 6] / c))
    suppressWarnings(draws[, 11] <- rbinom(reps, d, draws[, 7] / d))
    draws[, 12] <- a - draws[, 8]
    draws[, 13] <- b - draws[, 9]
    draws[, 14] <- c - draws[, 10]
    draws[, 15] <- d - draws[, 11]

    ## Clean up
    draws[, 16] <- apply(draws[, c(4:15)], MARGIN = 1, function(x) sum(x > 0))
    draws[, 16] <- ifelse(draws[, 16] != 12 | is.na(draws[, 16]), NA, 1)
    discard <- sum(is.na(draws[, 16]))
    if (sum(is.na(draws[, 16])) > 0) {
        cli::cli_alert_warning("Samplings lead to {discard} instance{?s} in which
sampled cell counts w{?as/ere} zero and discarded.")
        neg_warn <- paste("Samplings lead to",  discard, "non-sensible value(s).")
    } else neg_warn <- NULL

    draws <- draws[draws[, 16] == 1 & !is.na(draws[, 16]), ]
    reps <- dim(draws)[1]

    draws[, 17] <- (draws[, 4] * (draws[, 5] + draws[, 7]) /
                    ((draws[, 4] + draws[, 6]) + (draws[, 5] + draws[, 7])) +
                    draws[, 4] * (draws[, 5] + draws[, 7]) /
                    ((draws[, 4] + draws[, 6]) + (draws[, 5] + draws[, 7]))) /
        (draws[, 5] * (draws[, 4] + draws[, 6]) /
         ((draws[, 4] + draws[, 6]) + (draws[, 5] + draws[, 7])) +
         draws[, 5] * (draws[, 4] + draws[, 6]) /
         ((draws[, 4] + draws[, 6]) + (draws[, 5] + draws[, 7])))  ## RR_MH
    draws[, 18] <- (draws[, 4] * draws[, 7] / (draws[, 6] + draws[, 7]) +
                    draws[, 8] * draws[, 7] / (draws[, 6] + draws[, 7])) /
        (draws[, 5] * draws[, 6] / (draws[, 6] + draws[, 7]) + draws[, 5] * draws[, 6] /
         (draws[, 6] + draws[, 7]))  ## OR_MH

    ## Step 4c: Sample the bias-adjusted effect estimate
    cli::cli_progress_step("Incorporating random error", spinner = TRUE)
    draws[, 19] <- (((draws[, 12] + draws[, 13]) * (draws[, 12] + draws[, 14]) * (draws[, 13] + draws[, 15]) / ((draws[, 12] + draws[, 14]) + (draws[, 13] + draws[, 15]))^2 - draws[, 12] * draws[, 13] / ((draws[, 12] + draws[, 14]) + (draws[, 13] + draws[, 15]))) +
                    ((draws[, 8] + draws[, 9]) * (draws[, 8] + draws[, 10]) * (draws[, 9] + draws[, 11]) / ((draws[, 8] + draws[, 10]) + (draws[, 9] + draws[, 11]))^2 - draws[, 8] * draws[, 9] / ((draws[, 8] + draws[, 10]) + (draws[, 9] + draws[, 11])))) /
        (((draws[, 12] * (draws[, 13] + draws[, 15]) / ((draws[, 12] + draws[, 14]) + (draws[, 13] + draws[, 15])) + draws[, 8] * (draws[, 9] + draws[, 11]) / ((draws[, 8] + draws[, 10]) + (draws[, 9] + draws[, 11])))) * (draws[, 13] * (draws[, 12] + draws[, 14]) / ((draws[, 12] + draws[, 14]) + (draws[, 13] + draws[, 15])) + draws[, 9] * (draws[, 8] + draws[, 10]) / ((draws[, 8] + draws[, 10]) + (draws[, 9] + draws[, 11]))))

    draws[, 20] <- exp(rnorm(reps, log(draws[, 17]), sqrt(draws[, 19])))
    draws[, 21] <- exp(rnorm(reps, log(draws[, 18]), sqrt(draws[, 19])))

    corr_RR <- c(median(draws[, 17], na.rm = TRUE),
                 quantile(draws[, 17], probs = .025, na.rm = TRUE),
                 quantile(draws[, 17], probs = .975, na.rm = TRUE))
    corr_OR <- c(median(draws[, 18], na.rm = TRUE),
                 quantile(draws[, 18], probs = .025, na.rm = TRUE),
                 quantile(draws[, 18], probs = .975, na.rm = TRUE))
    tot_RR <- c(median(draws[, 20], na.rm = TRUE),
                quantile(draws[, 20], probs = .025, na.rm = TRUE),
                quantile(draws[, 20], probs = .975, na.rm = TRUE))
    tot_OR <- c(median(draws[, 21], na.rm = TRUE),
                quantile(draws[, 21], probs = .025, na.rm = TRUE),
                quantile(draws[, 21], probs = .975, na.rm = TRUE))

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
    rownames(rmatc) <- c("           RR (MH) -- systematic error:",
                         "RR (MH) -- systematic and random error:",
                         "           OR (SMR) -- systematic error:",
                         "OR (SMR) -- systematic and random error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    res <- list(obs_data = tab,
                obs_measures = rmat,
                adj_measures = rmatc,
                sim_df = as.data.frame(draws[, -16]),
                reps = reps,
                fun = "probsens_conf",
                warnings = neg_warn)
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
