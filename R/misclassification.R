#' Misclassification of exposure or outcome
#'
#' `misclass()`, `probsens()` and `probcase()` allow to provide adjusted measures
#' of association corrected for misclassification of the exposure or the outcome.
#'
#' @section Simple bias analysis with `misclass()`:
#' `misclass()` allows you to run a simple sensitivity analysis for disease or
#' exposure misclassification. Confidence interval for odds ratio adjusted using
#' sensitivity and specificity is computed as in Chu et al. (2006), for exposure
#' misclassification.
#'
#' For exposure misclassification, bias-adjusted measures are available using
#' sensitivity and specificity, or using predictive values.
#'
#' @section Probabilistic sensitivity analysis with `probsens()`:
#' `probsens()` performs a summary-level probabilistic sensitivity analysis to
#' correct for exposure misclassification or outcome misclassification and random
#' error. Non-differential misclassification is assumed when only the two bias
#' parameters \code{seca} and \code{spca} are provided. Adding the 2 parameters
#' \code{seexp} and \code{spexp} (i.e. providing the 4 bias parameters) evaluates
#' a differential misclassification.
#'
#' For exposure misclassification, bias-adjusted measures are available using
#' sensitivity and specificity, or using predictive values. However, only a beta
#' distribution is available for predictive values.
#'
#' Correlations between sensitivity (or specificity) of exposure classification
#' among cases and controls can be specified and use the NORmal To Anything
#' (NORTA) transformation (Li & Hammond, 1975).
#'
#' In case of negative (<=0) adjusted count in the 2-by-2 table following given
#' prior Se/Sp distribution(s), draws are discarded.
#'
#' @section Probabilistic sensitivity analysis with `probcase()`:
#' `probcase()` performs a record-level probabilistic sensitivity analysis to
#' correct for exposure misclassification or outcome misclassification and random
#' error. Non-differential misclassification is assumed when only the two bias
#' parameters \code{seca} and \code{spca} are provided. Adding the 2 parameters
#' \code{seexp} and \code{spexp} (i.e. providing the 4 bias parameters) evaluates
#' a differential misclassification.
#'
#' ORs are estimated with a logistic regression. RRs are estimated with a Poisson
#' regression (Barros & Hirakata, 2003; McNutt et al., 2003) with "robust"
#' standard errors using the sandwich estimator (Greenland, 2004; Zhou, 2004).
#'
#' Note that the exposure and the outcome have to be numeric variables coded 0
#' (absence of the event) or 1 (presence of the event).
#'
#' Correlations between sensitivity (or specificity) of exposure classification
#' among cases and controls can be specified and use the NORmal To Anything
#' (NORTA) transformation (Li & Hammond, 1975).
#'
#' In case of negative (<=0) adjusted count in the 2-by-2 table following given
#' prior Se/Sp distribution(s), draws are discarded.
#'
#' @section Updated calculations, probabilistic bias analysis:
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
#' @param type Choice of misclassification:
#'   \enumerate{
#'   \item exposure: bias analysis for exposure misclassification; corrections
#'   using sensitivity and specificity: nondifferential and independent errors,
#'   \item exposure_pv: bias analysis for exposure misclassification; corrections
#'   using PPV/NPV: nondifferential and independent errors,
#'   \item outcome: bias analysis for outcome misclassification.
#'   }
#' @param bias_parms Vector defining the bias parameters. This vector has 4
#'   elements between 0 and 1, in the following order:
#'   \enumerate{
#'   \item Sensitivity of exposure (when \code{type = "exposure"}) or outcome
#'   (when \code{type = "outcome"}) classification among those with the outcome
#'   (when \code{type = "exposure"}) or exposure (when \code{type = "outcome"}),
#'   \item Sensitivity of exposure (or outcome) classification among those without
#'   the outcome (or exposure),
#'   \item Specificity of exposure (or outcome) classification among those with the
#'   outcome (or exposure), and
#'   \item Specificity of exposure (or outcome) classification among those without
#'   the outcome (or exposure).
#'   }
#'   If PPV/NPV is chosen in case of exposure misclassification, this vector is the
#'   following:
#'   \enumerate{
#'   \item Positive predictive value among those with the outcome,
#'   \item Positive predictive value among those without the outcome,
#'   \item Negative predictive value among those with the outcome,
#'   \item Negative predictive value among those without the outcome.
#'   }
#' @param alpha Significance level.
#'
#' @return A list with elements (for `misclass()`):
#'   \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#'   \item{corr_data}{The expected observed data given the true data assuming misclassification.}
#'   \item{obs_measures}{A table of observed relative risk and odds ratio with
#'   confidence intervals.}
#'   \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#'   \item{bias_parms}{Input bias parameters.}
#'
#' @family misclassification
#'
#' @references
#' Fox, M.P, MacLehose, R.F., Lash, T.L., 2021 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.141--176, 233--256, 293--308, Springer.
#'
#' Li, S.T., Hammond, J.L., 1975. \emph{Generation of Pseudorandom Numbers
#' with Specified Univariate Distributions and Correlation Coefficients}.
#' IEEE Trans Syst Man Cybern 5:557-561.
#'
#' Chu, H., Zhaojie, W., Cole, S.R., Greenland, S., \emph{Sensitivity analysis of
#' misclassification: A graphical and a Bayesian approach},
#' Annals of Epidemiology 2006;16:834-841.
#'
#' Barros, A. & Hirakata, V.N., 2003. Alternatives for Logistic Regression in
#' Cross-sectional Studies: An Empirical Comparison of Models that Directly
#' Estimate the Prevalence Ratio. BMC Medical Research Methodology 3:21.
#'
#' McNutt, L-A, Wu, C., Xue, X., Hafner J.P., 2003. \emph{Estimating the Relative
#' Risk in Cohort Studies and Clinical Trials of Common Outcomes}. American
#' Journal of Epidemiology 157(10):940-943.
#'
#' Greenland, S. (2004). Model-based Estimation of Relative Risks and Other
#' Epidemiologic Measures in Studies of Common Outcomes and in Case-Control
#' Studies. American Journal of Epidemiology 160(4):301-305.
#'
#' Zhou, G. (2004). A Modified Poisson Regression Approach to Prospective Studies
#' with Binary Data. American Journal of Epidemiology 159(7):702-706.
#'
#' @examples
#' # The data for this example come from:
#' # Fink, A.K., Lash,  T.L. A null association between smoking during pregnancy
#' # and breast cancer using Massachusetts registry data (United States).
#' # Cancer Causes Control 2003;14:497-503.
#' misclass(matrix(c(215, 1449, 668, 4296),
#' dimnames = list(c("Breast cancer+", "Breast cancer-"),
#' c("Smoker+", "Smoker-")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure",
#' bias_parms = c(.78, .78, .99, .99))
#'
#' misclass(matrix(c(4558, 3428, 46305, 46085),
#' dimnames = list(c("AMI death+", "AMI death-"),
#' c("Male+", "Male-")),
#' nrow = 2, byrow = TRUE),
#' type = "outcome",
#' bias_parms = c(.53, .53, .99, .99))
#'
#' # The following example comes from Chu et al. Sensitivity analysis of
#' # misclassification: A graphical and a Bayesian approach.
#' # Annals of Epidemiology 2006;16:834-841.
#' misclass(matrix(c(126, 92, 71, 224),
#' dimnames = list(c("Case", "Control"), c("Smoker +", "Smoker -")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure",
#' bias_parms = c(.94, .94, .97, .97))
#'
#' # The next example, using PPV/NPV, comes from Bodnar et al. Validity of birth
#' # certificate-derived maternal weight data.
#' # Paediatric and Perinatal Epidemiology 2014;28:203-212.
#' misclass(matrix(c(599, 4978, 31175, 391851),
#' dimnames = list(c("Preterm", "Term"), c("Underweight", "Normal weight")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure_pv",
#' bias_parms = c(0.65, 0.74, 1, 0.98))
#' @export
#' @importFrom stats qnorm
misclass <- function(case,
                     exposed,
                     type = c("exposure", "exposure_pv", "outcome"),
                     bias_parms = NULL,
                     alpha = 0.05) {
    if (is.null(bias_parms))
        bias_parms <- c(1, 1, 1, 1)
    else bias_parms <- bias_parms
    if (length(bias_parms) != 4)
        stop(cli::format_error(c("i" = "The argument bias_parms should be made of
the following components: (1) Sensitivity of classification among those with the
outcome, (2) Sensitivity of classification among those without the outcome,
(3) Specificity of classification among those with the outcome, and (4) Specificity
of classification among those without the outcome.")))
    if (!all(bias_parms >= 0 & bias_parms <=1))
        stop(cli::format_error(c("x" = "Bias parameters should be between 0 and 1.")))

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
    if (type == "exposure") {
        obs_rr <- (a / (a + c)) / (b / (b + d))
        se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
        lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
        uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

        obs_or <- (a / b) / (c / d)
        se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
        uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

        A <- (a - (1 - bias_parms[3]) * (a + b)) / (bias_parms[1] - (1 - bias_parms[3]))
        C <- (c - (1 - bias_parms[4]) * (c + d)) / (bias_parms[2] - (1 - bias_parms[4]))
        B <- (a + b) - A
        D <- (c + d) - C

        if (A < 1 | B < 1 | C < 1 | D < 1)
            stop(cli::format_error(c("x" = "Parameters chosen lead to negative
cell(s) in adjusted 2x2 table.")))

        corr_tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr_rr <- (A / (A + C)) / (B / (B + D))
        corr_or <- (A / B) / (C / D)

        mle_corr_or <- ((a + ((a + b) * (bias_parms[3] - 1))) * (((c + d) * bias_parms[2]) - c)) /
            ((c + ((c + d) * (bias_parms[4] - 1))) * (((a + b) * bias_parms[1]) - a))

        se_corr_or <- sqrt((((a + b) * a * b * ((bias_parms[1] + bias_parms[3] - 1)^2)) /
                            (((((a + b) * bias_parms[1]) - a)^2) * ((((a + b) * bias_parms[3]) - b)^2))) +
                           (((c + d) * c * d * ((bias_parms[2] + bias_parms[4] - 1)^2)) /
                            (((((c + d) * bias_parms[2]) - c)^2) * ((((c + d) * bias_parms[4]) - d)^2))))
        lci_corr_or <- exp(log(mle_corr_or) - qnorm(1 - alpha / 2) * se_corr_or)
        uci_corr_or <- exp(log(mle_corr_or) + qnorm(1 - alpha / 2) * se_corr_or)

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(corr_tab) <- paste("Row", 1:2)
        } else {
            rownames(corr_tab) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(corr_tab) <- paste("Col", 1:2)
        } else {
           colnames(corr_tab) <- colnames(tab)
        }
        rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr), c(obs_or, lci_obs_or, uci_obs_or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
        rmatc <- rbind(c(corr_rr, NA, NA), c(corr_or, lci_corr_or, uci_corr_or))
        rownames(rmatc) <- c("Misclassification Bias Corrected Relative Risk:",
                             "   Misclassification Bias Corrected Odds Ratio:")
        colnames(rmatc) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
    }

    if (type == "exposure_pv") {
        obs_rr <- (a / (a + c)) / (b / (b + d))
        se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
        lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
        uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

        obs_or <- (a / b) / (c / d)
        se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
        uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

        A <- ((a * bias_parms[1]) + (b * (1 - bias_parms[3])))
        B <- a + b - A
        C <- ((c * bias_parms[2]) + (d * (1 - bias_parms[4])))
        D <- c + d - C

        if (A < 1 | B < 1 | C < 1 | D < 1)
            stop(cli::format_error(c("x" = "Parameters chosen lead to negative
cell(s) in adjusted 2x2 table.")))

        corr_tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr_rr <- (A / (A + C)) / (B / (B + D))
        corr_or <- (A / B) / (C / D)

        lci_corr_or <- NA
        uci_corr_or <- NA

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(corr_tab) <- paste("Row", 1:2)
        } else {
            rownames(corr_tab) <- rownames(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(corr_tab) <- paste("Col", 1:2)
        } else {
            colnames(corr_tab) <- colnames(tab)
        }
        rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr), c(obs_or, lci_obs_or, uci_obs_or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
        rmatc <- rbind(c(corr_rr, NA, NA), c(corr_or, lci_corr_or, uci_corr_or))
        rownames(rmatc) <- c("Misclassification Bias Corrected Relative Risk:",
                             "   Misclassification Bias Corrected Odds Ratio:")
        colnames(rmatc) <- c(" ",
                             paste(100 * (alpha / 2), "%", sep = ""),
                             paste(100 * (1 - alpha / 2), "%", sep = ""))
    }

    if (type == "outcome") {
        obs_rr <- (a / (a + c)) / (b / (b + d))
        se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
        lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
        uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

        obs_or <- (a / b) / (c / d)
        se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
        uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

        A <- (a - (1 - bias_parms[3]) * (a + c)) / (bias_parms[1] - (1 - bias_parms[3]))
        B <- (b - (1 - bias_parms[4]) * (b + d)) / (bias_parms[2] - (1 - bias_parms[4]))
        C <- (a + c) - A
        D <- (b + d) - B

        if (A < 1 | B < 1 | C < 1 | D < 1)
            stop(cli::format_error(c("x" = "Parameters chosen lead to negative cell(s)
in adjusted 2x2 table.")))

        corr_tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr_rr <- (A / (A + C)) / (B / (B + D))
        corr_or <- (A / B) / (C / D)

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(corr_tab) <- paste("Row", 1:2)
        } else {
            rownames(corr_tab) <- rownames(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(corr_tab) <- paste("Col", 1:2)
        } else {
            colnames(corr_tab) <- colnames(tab)
        }
        rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr), c(obs_or, lci_obs_or, uci_obs_or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
        rmatc <- rbind(corr_rr, corr_or)
        rownames(rmatc) <- c("Misclassification Bias Corrected Relative Risk:",
                             "   Misclassification Bias Corrected Odds Ratio:")
        colnames(rmatc) <- " "
    }

    res <- list(model = "misclassification",
                type = type,
                obs_data = tab,
                corr_data = corr_tab,
                obs_measures = rmat,
                adj_measures = rmatc,
                bias_parms = bias_parms)
    class(res) <- c("episensr", "episensr.boot", "list")
    res
}

#' @rdname misclass
#' @param reps Number of replications to run.
#' @param seca List defining sensitivity among cases:
#'   \enumerate{
#'   \item The sensitivity of exposure classification among those with the outcome
#'   (when \code{type = "exposure"}), or
#'   \item The sensitivity of outcome classification among those with the exposure
#'   (when \code{type = "outcome"}).
#'   }
#'   The first argument provides the probability distribution function (constant,
#'   uniform, triangular, trapezoidal, truncated normal, or beta) and the second
#'   its parameters as a vector. Lower and upper bounds of the truncated normal
#'   have to be between 0 and 1.
#'   \enumerate{
#'   \item constant: constant value,
#'   \item uniform: min, max,
#'   \item triangular: lower limit, upper limit, mode,
#'   \item trapezoidal: min, lower mode, upper mode, max,
#'   \item normal: lower bound, upper bound, mean, sd.
#'   \item beta: alpha, beta.
#'   }
#'   If PPV/NPV is chosen in case of exposure misclassification, the same four (4)
#'   parameters `seca`, `seexp`, `spca`, `spexp` as for Se/Sp have to be used but
#'   with the following meaning, and only for a beta distributions and no
#'   correlation between distributions:
#'   \enumerate{
#'   \item Positive predictive value among those with the outcome,
#'   \item Positive predictive value among those without the outcome,
#'   \item Negative predictive value among those with the outcome,
#'   \item Negative predictive value among those without the outcome.
#'   }
#' @param seexp List defining sensitivity among controls:
#'   \enumerate{
#'   \item The sensitivity of exposure classification among those without the
#'   outcome (when \code{type = "exposure"}), or
#'   \item The sensitivity of outcome classification among those without the
#'   exposure (when \code{type = "outcome"}).
#'   }
#' @param spca List as above for \code{seca} but for specificity.
#' @param spexp List as above for \code{seexp} but for specificity.
#' @param corr_se Correlation between case and non-case sensitivities. If PPV/NPV is
#' chosen in case of exposure misclassification, correlations are set to NULL.
#' @param corr_sp Correlation between case and non-case specificities.
#'
#' @return A list with elements (for `probsens()`):
#'   \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#'   \item{obs_measures}{A table of observed relative risk and odds ratio with
#'   confidence intervals.}
#'   \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#'   \item{sim_df}{Data frame of random parameters and computed values.}
#'   \item{reps}{Number of replications.}
#'
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
#' set.seed(1234)
#' probsens(fox, type = "exposure", reps = 10^6,
#' seca = list("beta", c(25, 3)),
#' spca = list("trapezoidal", c(.9, .93, .97, 1)),
#' seexp = list("beta", c(47, 7)),
#' spexp = list("trapezoidal", c(.8, .83, .87, .9)),
#' corr_se = .8,
#' corr_sp = .8)
#'
#' # Using PPV/NPV, from Bodnar et al. Validity of birth certificate-derived maternal
#' # weight data. Paediatric and Perinatal Epidemiology 2014;28:203-212.
#' set.seed(1234)
#' probsens(matrix(c(599, 4978, 31175, 391851),
#' dimnames = list(c("Preterm", "Term"), c("Underweight", "Normal weight")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure_pv", reps = 10^6,
#' seca = list("beta", c(50, 27)),  ## PPV_case
#' spca = list("beta", c(120, .5)),  ## NPV_case
#' seexp = list("beta", c(132, 47)),  ## PPV_ctrl
#' spexp = list("beta", c(115, 2)))  ## NPV_ctrl
#' @export
#' @importFrom stats median pnorm qnorm quantile qunif runif rnorm rbinom qbeta rbeta
probsens <- function(case,
                     exposed,
                     type = c("exposure", "exposure_pv", "outcome"),
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
    if ((seca[[1]] == "constant" | seca[[1]] == "uniform" | seca[[1]] == "triangular" |
         seca[[1]] == "trapezoidal") & !all(seca[[2]] >= 0 & seca[[2]] <= 1))
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
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(seexp) && seexp[[1]] == "beta" && (seexp[[2]][1] < 0 | seexp[[2]][2] < 0))
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
         spca[[1]] == "triangular" |
         spca[[1]] == "trapezoidal") & !all(spca[[2]] >= 0 & spca[[2]] <= 1))
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
        ((spexp[[2]][1] > spexp[[2]][2]) | (spexp[[2]][2] > spexp[[2]][3]) |
         (spexp[[2]][3] > spexp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" & (length(spexp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution, please
provide vector of lower and upper bound limits, meand and SD.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" &&
        ((spexp[[2]][1] >= spexp[[2]][2]) |
         (!all(spexp[[2]][1:2] >= 0 & spexp[[2]][1:2] <= 1))))
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
    if (!is.null(spexp) && spexp[[1]] == "beta" && (spexp[[2]][1] < 0 | spexp[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(seexp) & (type != "exposure_pv") &
        (is.null(spca) | is.null(spexp) | is.null(corr_se) | is.null(corr_sp)))
        stop(cli::format_error(c("i" = "For differential misclassification type,
have to provide Se and Sp for among those with and without the outcome as well as
Se and Sp correlations.")))

    if (type == "exposure_pv" &&
        (!(seca[[1]] %in% c("beta")) | !(seexp[[1]] %in% c("beta")) | !(spca[[1]] %in% c("beta")) |
         !(spexp[[1]] %in% c("beta"))))
        stop(cli::format_error(c("x" = "Wrong distributions provided for exposure misclassification using PPV/NPV.")))
    if (type == "exposure_pv" &&
        ((length(seca[[2]]) != 2) | (length(seexp[[2]]) != 2) | (length(spca[[2]]) != 2) |
         (length(spexp[[2]]) != 2)))
        stop(cli::format_error(c("x" = "Wrong distributions provided for exposure misclassification using PPV/NPV.")))
    if (type == "exposure_pv" && ((seca[[2]][1] < 0 | seca[[2]][2] < 0)))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution(s).")))
    if (type == "exposure_pv" && ((seexp[[2]][1] < 0 | seexp[[2]][2] < 0)))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution(s).")))
    if (type == "exposure_pv" && ((spca[[2]][1] < 0 | spca[[2]][2] < 0)))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution(s).")))
    if (type == "exposure_pv" && ((spexp[[2]][1] < 0 | spexp[[2]][2] < 0)))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution(s).")))

    if (!is.null(corr_se) && (corr_se == 0 | corr_se == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))
    if (!is.null(corr_sp) && (corr_sp == 0 | corr_sp == 1))
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

    draws <- matrix(NA, nrow = reps, ncol = 28)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A1", "B1", "C1", "D1",
                         "flag",
                         "prevca", "prevexp",
                         "ppvca", "ppvexp", "npvca", "npvexp",
                         "ab", "bb", "cb", "db",
                         "corr_RR", "corr_OR",
                         "rr_se_b", "or_se_b", "z",
                         "tot_RR", "tot_OR", "syst_RR", "syst_OR")
    corr_draws <- matrix(NA, nrow = reps, ncol = 4)

    se1 <- c(reps, seca[[2]])
    se0 <- c(reps, seexp[[2]])
    sp1 <- c(reps, spca[[2]])
    sp0 <- c(reps, spexp[[2]])

    if (type == "exposure_pv") {
        corr_se <- NULL
        corr_sp <- NULL
    }

    ## Step3: Assign probability distributions to each bias parameter
    ## and Step 4a draw Se's and Sp's
    cli::cli_progress_step("Assign probability distributions", spinner = TRUE)
    if (type == "exposure_pv") {
        draws[, 1] <- do.call(rbeta, as.list(se1))
        draws[, 2] <- do.call(rbeta, as.list(se0))
        draws[, 3] <- do.call(rbeta, as.list(sp1))
        draws[, 4] <- do.call(rbeta, as.list(sp0))
    } else if (is.null(seexp) & !is.null(spca) &
               is.null(spexp) & is.null(corr_se) & is.null(corr_sp)) {
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

        ## Systematic error
        draws[, 27] <- (draws[, 5] / (draws[, 5] + draws[, 7])) / (draws[, 6] / (draws[, 6] + draws[, 8]))
        draws[, 28] <- (draws[, 5] / draws[, 7]) / (draws[, 6] / draws[, 8])

        ## Clean up
        draws[, 9] <- apply(draws[, c(5:8, 16:19)], MARGIN = 1, function(x) sum(x > 0))
        draws[, 9] <- ifelse(draws[, 9] != 8 | is.na(draws[, 9]), NA, 1)
        discard <- sum(is.na(draws[, 9]))
        if (sum(is.na(draws[, 9])) > 0) {
            cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
            neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
        } else neg_warn <- NULL

        draws <- draws[draws[, 9] == 1 & !is.na(draws[, 9]), ]

        rr_syst <- c(median(draws[, 27], na.rm = TRUE),
                     quantile(draws[, 27], probs = .025, na.rm = TRUE),
                     quantile(draws[, 27], probs = .975, na.rm = TRUE))
        or_syst <- c(median(draws[, 28], na.rm = TRUE),
                     quantile(draws[, 28], probs = .025, na.rm = TRUE),
                     quantile(draws[, 28], probs = .975, na.rm = TRUE))
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

    if (type == "exposure_pv") {
        cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
        draws <- draws[, 1:12]
        colnames(draws) <- c("PPV_case", "PPV_ctrl", "NPV_case", "NPV_ctrl",
                             "r0", "r1",
                             "q1", "q0", "flag",
                             "psi",
                             "corr_RR", "corr_OR")
        draws[, 5] <- rbeta(reps, c, d)
        draws[, 6] <- rbeta(reps, a, b)

        draws[, 7] <- draws[, 1] * draws[, 6] + (1 - draws[, 3]) * (1 - draws[, 6])
        draws[, 8] <- draws[, 2] * draws[, 5] + (1 - draws[, 4]) * (1 - draws[, 5])
        draws[, 10] <- logit(draws[, 7]) - logit(draws[, 8])
        draws[, 12] <- exp(draws[, 10])

        or_syst <- c(median(draws[, 12], na.rm = TRUE),
                     quantile(draws[, 12], probs = .025, na.rm = TRUE),
                     quantile(draws[, 12], probs = .975, na.rm = TRUE))

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
        rmatc <- matrix(or_syst, ncol = 3, byrow = TRUE)
        rownames(rmatc) <- "Odds Ratio -- systematic error:"
        colnames(rmatc) <- c("Median", "p2.5", "p97.5")
        neg_warn <- NULL

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

        ## Systematic error
        draws[, 27] <- (draws[, 5] / (draws[, 5] + draws[, 7])) / (draws[, 6] / (draws[, 6] + draws[, 8]))
        draws[, 28] <- (draws[, 5] / draws[, 7]) / (draws[, 6] / draws[, 8])

        ## Clean up
        draws[, 9] <- apply(draws[, c(5:8, 16:19)], MARGIN = 1, function(x) sum(x > 0))
        draws[, 9] <- ifelse(draws[, 9] != 8 | is.na(draws[, 9]), NA, 1)
        discard <- sum(is.na(draws[, 9]))
        if (sum(is.na(draws[, 9])) > 0) {
            cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
            neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
        } else neg_warn <- NULL

        draws <- draws[draws[, 9] == 1 & !is.na(draws[, 9]), ]

        rr_syst <- c(median(draws[, 27], na.rm = TRUE),
                     quantile(draws[, 27], probs = .025, na.rm = TRUE),
                     quantile(draws[, 27], probs = .975, na.rm = TRUE))
        or_syst <- c(median(draws[, 28], na.rm = TRUE),
                     quantile(draws[, 28], probs = .025, na.rm = TRUE),
                     quantile(draws[, 28], probs = .975, na.rm = TRUE))
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
    res <- list(obs_data = tab,
                obs_measures = rmat,
                adj_measures = rmatc,
                sim_df = as.data.frame(draws[, -9]),
                reps = reps,
                fun = "probsens",
                warnings = neg_warn
                )
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}


#' @rdname misclass
#' @param df Default dataset, a data frame, to use for analysis.
#' @param x,y,... <[`data-masking`][rlang::topic-data-mask]> List of name-value
#' pairs describing which variables in the data should be used. The expression
#' `variable` is evaluated within `df`, so there is no need to refer to the
#' original dataset (i.e., use `probcase(df, variable)` instead of
#' `probcase(df, df$variable)`). Note that `x` and `y` should be numeric
#' dichotomous 0,1 variables.
#' @param type Choice of misclassification:
#'   \enumerate{
#'   \item exposure: bias analysis for exposure misclassification; corrections
#'   using sensitivity and specificity: nondifferential and independent errors,
#'   \item outcome: bias analysis for outcome misclassification.
#'   }
#' @param reps Number of replications to run.
#' @param seca List defining sensitivity among cases:
#'   \enumerate{
#'   \item The sensitivity of exposure classification among those with the outcome
#'   (when \code{type = "exposure"}), or
#'   \item The sensitivity of outcome classification among those with the exposure
#'   (when \code{type = "outcome"}).
#'   }
#'   The first argument provides the probability distribution function (constant,
#'   uniform, triangular, trapezoidal, truncated normal, or beta) and the second
#'   its parameters as a vector. Lower and upper bounds of the truncated normal
#'   have to be between 0 and 1.
#'   \enumerate{
#'   \item constant: constant value,
#'   \item uniform: min, max,
#'   \item triangular: lower limit, upper limit, mode,
#'   \item trapezoidal: min, lower mode, upper mode, max,
#'   \item normal: lower bound, upper bound, mean, sd.
#'   \item beta: alpha, beta.
#'   }
#' @param seexp List defining sensitivity among controls:
#'   \enumerate{
#'   \item The sensitivity of exposure classification among those without the
#'   outcome (when \code{type = "exposure"}), or
#'   \item The sensitivity of outcome classification among those without the
#'   exposure (when \code{type = "outcome"}).
#'   }
#' @param spca List as above for \code{seca} but for specificity.
#' @param spexp List as above for \code{seexp} but for specificity.
#' @param corr_se Correlation between case and non-case sensitivities. If PPV/NPV is
#' chosen in case of exposure misclassification, correlations are set to NULL.
#' @param corr_sp Correlation between case and non-case specificities.
#' @param alpha Significance level.
#' @return A list with elements (for `probcase()`):
#'   \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#'   \item{obs_measures}{A table of observed relative risk and odds ratio with
#'   confidence intervals.}
#'   \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#'   \item{sim_df}{Data frame of random parameters and computed values at each run.}
#'   \item{reps}{Number of replications.}
#'
#' @examples
#' # Fox M.P., MacLehose R.F., Lash T.L.
#' # SAS and R code for probabilistic quantitative bias analysis for
#' # misclassified binary variables and binary unmeasured confounders
#' # Int J Epidemiol 2023:1624-1633.
#' \dontrun{
#' set.seed(1234)
#' probcase(D, x = e_obs, y = d, reps = 10^3,
#' type = "exposure",
#' seca = list("beta", c(25, 3)),
#' spca = list("trapezoidal", c(.9, .93, .97, 1)),
#' seexp = list("beta", c(45, 7)),
#' spexp = list("trapezoidal", c(.8, .83, .87, .9)),
#' corr_se = .8,
#' corr_sp = .8)
#'
#' set.seed(1234)
#' probcase(D, x = e, y = d_obs, reps = 10^3,
#' type = "outcome",
#' seca = list("beta", c(254, 24)),
#' spca = list("trapezoidal", c(.94, .96, .98, 1)),
#' seexp = list("beta", c(450, 67)),
#' spexp = list("trapezoidal", c(.9, .92, .93, .95)),
#' corr_se = .8, corr_sp = .8)
#'
#' library(aplore3)  # to get ICU data
#' data(icu)
#' icu$sta2 <- as.numeric(icu$sta) - 1
#' icu$inf2 <- as.numeric(icu$inf) - 1
#' probcase(icu, x = inf2, y = sta2, reps = 1000,
#' seca = list("beta", c(20, 3)), spca = list("trapezoidal", c(.8, .83, .87, .9)),
#' seexp = list("beta", c(60, 9)), spexp = list("trapezoidal", c(.88, .93, .97, 1)),
#' corr_se = .8, corr_sp = .8)
#' }
#' @export
#' @importFrom stats median pnorm qnorm quantile qunif runif rnorm rbinom qbeta rbeta reformulate glm binomial poisson coef confint
#' @importFrom sandwich sandwich
#' @importFrom lmtest coeftest
probcase <- function(df,
                     x,
                     y,
                     ...,
                     type = c("exposure", "outcome"),
                     reps = 100,
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
    x <- deparse(substitute(x))
    y <- deparse(substitute(y))
    argList <- as.character(match.call(expand.dots = FALSE)$...)

    if (!is.data.frame(df)) {
        stop(cli::format_error(c("x" = "Argument {deparse(substitute(df))} must be a data.frame.")))
    }
    if (reps < 1)
        stop(cli::format_error(c("x" = "Wrong number of replications: reps = {reps}",
                                 "i" = "reps must be >= 1")))
    if (reps > 10^5) {
        stop(cli::format_error(c("x" = "Are you sure you want {reps} replications?")))
    }
    if (dim(df)[2] < 2) {
        stop(cli::format_error(c("x" = "{df} does not have enough variables!")))
    }
    if (!(x %in% colnames(df)) | !(y %in% colnames(df)) | !(all(c(argList) %in% colnames(df)))) {
        stop(cli::format_error(c("x" = "{deparse(substitute(df))} does not contain x, y and/or confounders")))
    }
    if (!(is.numeric(df[[x]])) | !(is.numeric(df[[y]]))) {
        stop(cli::format_error(c("x" = "Exposure and outcome should be numeric variables.")))
    }
    if (!(all(df[[x]] %in% c(0, 1))) | !(all(df[[y]] %in% c(0, 1)))) {
        stop(cli::format_error(c("x" = "Exposure and outcome should be dichotmous 0,1 variables.")))
    }

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
    if ((seca[[1]] == "constant" | seca[[1]] == "uniform" | seca[[1]] == "triangular" |
         seca[[1]] == "trapezoidal") & !all(seca[[2]] >= 0 & seca[[2]] <= 1))
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
        stop(cli::format_error(c("i" = "For beta distribution, please provide alpha and beta.")))
    if (!is.null(seexp) && seexp[[1]] == "beta" && (seexp[[2]][1] < 0 | seexp[[2]][2] < 0))
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
         spca[[1]] == "triangular" |
         spca[[1]] == "trapezoidal") & !all(spca[[2]] >= 0 & spca[[2]] <= 1))
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
        ((spexp[[2]][1] > spexp[[2]][2]) | (spexp[[2]][2] > spexp[[2]][3]) |
         (spexp[[2]][3] > spexp[[2]][4])))
        stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" & (length(spexp[[2]]) != 4))
        stop(cli::format_error(c("i" = "For truncated normal distribution, please
provide vector of lower and upper bound limits, meand and SD.")))
    if (!is.null(spexp) && spexp[[1]] == "normal" &&
        ((spexp[[2]][1] >= spexp[[2]][2]) |
         (!all(spexp[[2]][1:2] >= 0 & spexp[[2]][1:2] <= 1))))
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
    if (!is.null(spexp) && spexp[[1]] == "beta" && (spexp[[2]][1] < 0 | spexp[[2]][2] < 0))
        stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!is.null(seexp) &
        (is.null(spca) | is.null(spexp) | is.null(corr_se) | is.null(corr_sp)))
        stop(cli::format_error(c("i" = "For differential misclassification type,
have to provide Se and Sp for among those with and without the outcome as well as
Se and Sp correlations.")))

    if (!is.null(corr_se) && (corr_se == 0 | corr_se == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))
    if (!is.null(corr_sp) && (corr_sp == 0 | corr_sp == 1))
        stop(cli::format_error(c("x" = "Correlations should be > 0 and < 1.")))

    obs_data <- df[, c(x, y, c(argList))]
    outcome <- colnames(obs_data[1])
    disease <- colnames(obs_data[2])
    confounder_names <- colnames(obs_data[-c(1:2)])
    tab_df <- table(obs_data[[x]], obs_data[[y]],
                    dnn = list(colnames(obs_data[1]), colnames(obs_data[2])))
    tab <- t(tab_df[2:1, 2:1])
    a <- as.numeric(tab[1, 1])
    b <- as.numeric(tab[1, 2])
    c <- as.numeric(tab[2, 1])
    d <- as.numeric(tab[2, 2])

    cli::cli_alert_info("Conventional analysis assuming no bias")
    formula_conv <- reformulate(c(outcome, confounder_names), response = disease)
    convrr_mod <- glm(formula_conv, data = obs_data, family = poisson(link = "log"))
    obs_rr <- exp(coef(convrr_mod))[2]
    obsci_rr <- exp(confint(convrr_mod, level = 1 - alpha))[2, ]

    convor_mod <- glm(formula_conv, data = obs_data, family = binomial)
    obs_or <- exp(coef(convor_mod))[2]
    obsci_or <- exp(confint(convor_mod, level = 1 - alpha))[2, ]

    draws <- matrix(NA, nrow = reps, ncol = 21)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A0", "B0", "C0", "D0",
                         "flag",
                         "prevca", "prevexp",
                         "ppvca", "ppvexp", "npvca", "npvexp",
                         "ped1", "ped0",
                         "PPV_d1", "PPV_d0", "NPV_d1", "NPV_d0")
    corr_draws <- matrix(NA, nrow = reps, ncol = 4)

    se1 <- c(reps, seca[[2]])
    se0 <- c(reps, seexp[[2]])
    sp1 <- c(reps, spca[[2]])
    sp0 <- c(reps, spexp[[2]])

    ## Step3: Assign probability distributions to each bias parameter
    ## and Step 4a draw Se's and Sp's
    cli::cli_progress_step("Assign probability distributions", spinner = TRUE)
    if (is.null(seexp) & !is.null(spca) &
        is.null(spexp) & is.null(corr_se) & is.null(corr_sp)) {
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
        ## Step 4b: Using simple bias analysis (A0, B0, C0, D0)
        cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) / (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (a + b) - draws[, 5]
        draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) / (draws[, 2] - (1 - draws[, 4]))
        draws[, 8] <- (c + d) - draws[, 7]

        ## Clean up
        draws[, 9] <- apply(draws[, 5:8], MARGIN = 1, function(x) sum(x > 0))
        draws[, 9] <- ifelse(draws[, 9] != 4 | is.na(draws[, 9]), NA, 1)
        discard <- sum(is.na(draws[, 9]))
        if (sum(is.na(draws[, 9])) > 0) {
            cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
            neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
        } else neg_warn <- NULL
        reps2 <- nrow(draws)

        ## Prevalence of exposure in cases and controls, accounting for sampling error
        ## For systematic error only
        draws[, 16] <- draws[, 5] / (draws[, 5] + draws[, 6])
        draws[, 17] <- draws[, 7] / (draws[, 7] + draws[, 8])
        ## Computing predictive values for imputation (PPV_d1, PPV_d0, NPV_d1, NPV_d0)
        ## For systematic error only
        draws[, 18] <- draws[, 1] * draws[, 16] /
            (draws[, 1] * draws[, 16] + (1 - draws[, 16]) * (1 - draws[, 3]))
        draws[, 19] <- draws[, 2] * draws[, 17] /
            (draws[, 2] * draws[, 17] + (1 - draws[, 17]) * (1 - draws[, 4]))
        draws[, 20] <- (draws[, 3] * (1 - draws[, 16])) /
            ((draws[, 3] * (1 - draws[, 16])) + (1 - draws[, 1]) * (draws[, 16]))
        draws[, 21] <- (draws[, 4] * (1 - draws[, 17])) /
            ((draws[, 4] * (1 - draws[, 17])) + (1 - draws[, 2]) * (draws[, 17]))
        ## For systematic and random error
        suppressWarnings({
                             draws[, 10] <- rbeta(reps2, draws[, 5], draws[, 6])
                             draws[, 11] <- rbeta(reps2, draws[, 7], draws[, 8])
                         })
        ## Computing predictive values for imputation (PPV_d1, PPV_d0, NPV_d1, NPV_d0)
        ## For systematic and random error
        draws[, 12] <- (draws[, 1] * draws[, 10]) /
            ((draws[, 1] * draws[, 10]) + (1 - draws[, 3]) * (1 - draws[, 10]))
        draws[, 13] <- (draws[, 2] * draws[, 11]) /
            ((draws[, 2] * draws[, 11]) + (1 - draws[, 4]) * (1 - draws[, 11]))
        draws[, 14] <- (draws[, 3] * (1 - draws[, 10])) /
            ((1 - draws[, 1]) * draws[, 10] + draws[, 3] * (1 - draws[, 10]))
        draws[, 15] <- (draws[, 4] * (1 - draws[, 11])) /
            ((1 - draws[, 2]) * draws[, 11] + draws[, 4] * (1 - draws[, 11]))

        ## Loop through draws and impute new exposures at each step
        cli::cli_progress_step("Assign probability of exposure for the record", spinner = TRUE)
        names(obs_data)[1] <- "e_obs"
        names(obs_data)[2] <- "d"
        obs_data[, "ppvca"] <- NA
        obs_data[, "ppvexp"] <- NA
        obs_data[, "npvca"] <- NA
        obs_data[, "npvexp"] <- NA
        obs_data[, "PPV_d1"] <- NA
        obs_data[, "PPV_d0"] <- NA
        obs_data[, "NPV_d1"] <- NA
        obs_data[, "NPV_d0"] <- NA
        obs_data[, "p"] <- NA
        obs_data[, "e_syst"] <- NA
        obs_data[, "e"] <- NA
        res_mat <- matrix(NA, nrow = reps2, ncol = 11)
        colnames(res_mat) <- c("rr_adj", "rr_tot", "rr_se",
                               "or_adj", "or_tot", "or_se",
                               "or_adj_syst",
                               "se_D1", "se_D0", "sp_D1", "sp_D0")
        formula <- reformulate(c("e", confounder_names), response = "d")
        cli::cli_progress_bar("Processing...", total = reps2)
        for (i in 1:reps2) {
            obs_data[, "ppvca"] <- draws[i, 12]
            obs_data[, "ppvexp"] <- draws[i, 13]
            obs_data[, "npvca"] <- draws[i, 14]
            obs_data[, "npvexp"] <- draws[i, 15]
            obs_data[, "PPV_d1"] <- draws[i, 18]
            obs_data[, "PPV_d0"] <- draws[i, 19]
            obs_data[, "NPV_d1"] <- draws[i, 20]
            obs_data[, "NPV_d0"] <- draws[i, 21]
            obs_data[, "p"] <- obs_data[, "e_obs"] * obs_data[, "d"] * obs_data[, "ppvca"] +
                obs_data[, "e_obs"] * (1 - obs_data[, "d"]) * obs_data[, "ppvexp"] +
                (1 - obs_data[, "e_obs"]) * obs_data[, "d"] * (1 - obs_data[, "npvca"]) +
                (1 - obs_data[, "e_obs"]) * (1 - obs_data[, "d"]) * (1 - obs_data[, "npvexp"])
            obs_data[, "e"] <- suppressWarnings({
                                                    rbinom(nrow(obs_data), 1, obs_data[, "p"])
                                                })
            obs_data[, "e_syst"] <- obs_data[, "e_obs"] * obs_data[, "d"] * obs_data[, "PPV_d1"] +
                obs_data[, "e_obs"] * (1 - obs_data[, "d"]) * obs_data[, "PPV_d0"] +
                (1 - obs_data[, "e_obs"]) * obs_data[, "d"] * (1 - obs_data[, "NPV_d1"]) +
                (1 - obs_data[, "e_obs"]) * (1 - obs_data[, "d"]) * (1 - obs_data[, "NPV_d0"])
            ## Logistic regression, systematic and random error
            if (all(is.na(obs_data[, "e"]))) {
                modrr_coef <- NA
                modrr_se <- NA
                modor_coef <- NA
                modor_se <- NA
            } else {
                modrr_pois <- glm(formula, data = obs_data,
                                  family = poisson(link = "log"))
                modrr_coef <- coef(summary(modrr_pois))[2, 1]
                modrr_se <- lmtest::coeftest(modrr_pois,
                                             vcov = sandwich::sandwich)[2, 2]
                modor_log <- glm(formula, data = obs_data,
                                 family = binomial(link = "log"))
                modor_coef <- coef(summary(modor_log))[2, 1]
                modor_se <- coef(summary(modor_log))[2, 2]
            }
            z <- rnorm(1)
            ## For systematic error only:
            if (all(is.na(obs_data[, "e_syst"]))) {
                rr_coef <- NA
                or_coef <- NA
            } else {
                at <- sum(obs_data$e_syst[obs_data$d == 1])
                ct <- sum(obs_data$e_syst[obs_data$d == 0])
                dt <- ct - sum(obs_data$d == 0)
                bt <- at - sum(obs_data$d == 1)
                or_coef <- log(at * dt / bt / ct)
            }

            res_mat[i, 1] <- exp(modrr_coef)
            res_mat[i, 2] <- exp(modrr_coef + z * modrr_se)
            res_mat[i, 3] <- modrr_se
            res_mat[i, 4] <- exp(modor_coef)
            res_mat[i, 5] <- exp(modor_coef + z * modor_se)
            res_mat[i, 6] <- modor_se
            res_mat[i, 7] <- exp(or_coef)
            res_mat[i, 8] <- draws[i, 1]
            res_mat[i, 9] <- draws[i, 2]
            res_mat[i, 10] <- draws[i, 3]
            res_mat[i, 11] <- draws[i, 4]
            cli::cli_progress_update()
        }

        rr_syst <- c(median(res_mat[, 1], na.rm = TRUE),
                     quantile(res_mat[, 1], probs = .025, na.rm = TRUE),
                     quantile(res_mat[, 1], probs = .975, na.rm = TRUE))
        rr_tot <- c(median(res_mat[, 2], na.rm = TRUE),
                    quantile(res_mat[, 2], probs = .025, na.rm = TRUE),
                    quantile(res_mat[, 2], probs = .975, na.rm = TRUE))
        or_syst <- c(median(res_mat[, 4], na.rm = TRUE),
                     quantile(res_mat[, 4], probs = .025, na.rm = TRUE),
                     quantile(res_mat[, 4], probs = .975, na.rm = TRUE))
        or_tot <- c(median(res_mat[, 5], na.rm = TRUE),
                    quantile(res_mat[, 5], probs = .025, na.rm = TRUE),
                    quantile(res_mat[, 5], probs = .975, na.rm = TRUE))
        or_syst2 <- c(median(res_mat[, 7], na.rm = TRUE),
                    quantile(res_mat[, 7], probs = .025, na.rm = TRUE),
                    quantile(res_mat[, 7], probs = .975, na.rm = TRUE))

        rmat <- rbind(c(obs_rr, obsci_rr[1], obsci_rr[2]),
                      c(obs_or, obsci_or[1], obsci_or[2]))
        rownames(rmat) <- c(" Observed Relative Risk:",
                            "    Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        rmatc <- rbind(rr_syst, rr_tot, or_syst, or_tot, or_syst2)
        rownames(rmatc) <- c("Relative Risk -- systematic error:",
                             "                      total error:",
                             "   Odds Ratio -- systematic error:",
                             "                      total error:",
                             "   Odds Ratio -- systematic err.2:")
        colnames(rmatc) <- c("Median", "p2.5", "p97.5")
    }

    if (type == "outcome") {
        ## Step 4b: Using simple bias analysis (A0, B0, C0, D0)
        cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + c)) / (draws[, 1] - (1 - draws[, 3]))
        draws[, 6] <- (b - (1 - draws[, 4]) * (b + d)) / (draws[, 2] - (1 - draws[, 4]))
        draws[, 7] <- (a + c) - draws[, 5]
        draws[, 8] <- (b + d) - draws[, 6]

        ## Clean up
        draws[, 9] <- apply(draws[, 5:8], MARGIN = 1, function(x) sum(x > 0))
        draws[, 9] <- ifelse(draws[, 9] != 4 | is.na(draws[, 9]), NA, 1)
        discard <- sum(is.na(draws[, 9]))
        if (sum(is.na(draws[, 9])) > 0) {
            cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
            neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
        } else neg_warn <- NULL
        reps2 <- nrow(draws)

        ## Prevalence of exposure in cases and controls, accounting for sampling error
        suppressWarnings({
                             draws[, 10] <- rbeta(reps2, draws[, 5], draws[, 7])
                             draws[, 11] <- rbeta(reps2, draws[, 6], draws[, 8])
                         })
        ## Computing predictive values for imputation (PPV_d1, PPV_d0, NPV_d1, NPV_d0)
        draws[, 12] <- (draws[, 1] * draws[, 10]) /
            ((draws[, 1] * draws[, 10]) + (1 - draws[, 3]) * (1 - draws[, 10]))
        draws[, 13] <- (draws[, 2] * draws[, 11]) /
            ((draws[, 2] * draws[, 11]) + (1 - draws[, 4]) * (1 - draws[, 11]))
        draws[, 14] <- (draws[, 3] * (1 - draws[, 10])) /
            ((1 - draws[, 1]) * draws[, 10] + draws[, 3] * (1 - draws[, 10]))
        draws[, 15] <- (draws[, 4] * (1 - draws[, 11])) /
            ((1 - draws[, 2]) * draws[, 11] + draws[, 4] * (1 - draws[, 11]))

        ## Loop through draws and impute new exposures at each step
        cli::cli_progress_step("Assign probability of outcome for the record", spinner = TRUE)
        names(obs_data)[1] <- "e"
        names(obs_data)[2] <- "d_obs"
        obs_data[, "ppvca"] <- NA
        obs_data[, "ppvexp"] <- NA
        obs_data[, "npvca"] <- NA
        obs_data[, "npvexp"] <- NA
        obs_data[, "p"] <- NA
        obs_data[, "d"] <- NA
        res_mat <- matrix(NA, nrow = reps2, ncol = 10)
        colnames(res_mat) <- c("rr_adj", "rr_tot", "rr_se",
                               "or_adj", "or_tot", "or_se",
                               "se_D1", "se_D0", "sp_D1", "sp_D0")
        formula <- reformulate(c("e", confounder_names), response = "d")
        cli::cli_progress_bar("Processing...", total = reps)
        for (i in 1:reps2) {
            obs_data[, "ppvca"] <- draws[i, 12]
            obs_data[, "ppvexp"] <- draws[i, 13]
            obs_data[, "npvca"] <- draws[i, 14]
            obs_data[, "npvexp"] <- draws[i, 15]
            obs_data[, "p"] <- obs_data[, "e"] * obs_data[, "d_obs"] * obs_data[, "ppvca"] +
                obs_data[, "e"] * (1 - obs_data[, "d_obs"]) * (1 - obs_data[, "npvca"]) +
                (1 - obs_data[, "e"]) * obs_data[, "d_obs"] * obs_data[, "ppvexp"] +
                (1 - obs_data[, "e"]) * (1 - obs_data[, "d_obs"]) * (1 - obs_data[, "npvexp"])
            obs_data[, "d"] <- suppressWarnings({
                                                    rbinom(nrow(obs_data), 1, obs_data[, "p"])
                                                })
            ## Logistic regression
            if (all(is.na(obs_data[, "d"]))) {
                modrr_coef <- NA
                modrr_se <- NA
                modor_coef <- NA
                modor_se <- NA
            } else {
                modrr_pois <- glm(formula, data = obs_data,
                                  family = poisson(link = "log"))
                modrr_coef <- coef(summary(modrr_pois))[2, 1]
                modrr_se <- lmtest::coeftest(modrr_pois,
                                             vcov = sandwich::sandwich)[2, 2]
                modor_log <- glm(formula, data = obs_data,
                                 family = binomial(link = "log"))
                modor_coef <- coef(summary(modor_log))[2, 1]
                modor_se <- coef(summary(modor_log))[2, 2]
            }
            z <- rnorm(1)

            res_mat[i, 1] <- exp(modrr_coef)
            res_mat[i, 2] <- exp(modrr_coef + z * modrr_se)
            res_mat[i, 3] <- modrr_se
            res_mat[i, 4] <- exp(modor_coef)
            res_mat[i, 5] <- exp(modor_coef + z * modor_se)
            res_mat[i, 6] <- modor_se
            res_mat[i, 7] <- draws[i, 1]
            res_mat[i, 8] <- draws[i, 2]
            res_mat[i, 9] <- draws[i, 3]
            res_mat[i, 10] <- draws[i, 4]
            cli::cli_progress_update()
        }

        rr_syst <- c(median(res_mat[, 1], na.rm = TRUE),
                     quantile(res_mat[, 1], probs = .025, na.rm = TRUE),
                     quantile(res_mat[, 1], probs = .975, na.rm = TRUE))
        rr_tot <- c(median(res_mat[, 2], na.rm = TRUE),
                    quantile(res_mat[, 2], probs = .025, na.rm = TRUE),
                    quantile(res_mat[, 2], probs = .975, na.rm = TRUE))
        or_syst <- c(median(res_mat[, 4], na.rm = TRUE),
                     quantile(res_mat[, 4], probs = .025, na.rm = TRUE),
                     quantile(res_mat[, 4], probs = .975, na.rm = TRUE))
        or_tot <- c(median(res_mat[, 5], na.rm = TRUE),
                    quantile(res_mat[, 5], probs = .025, na.rm = TRUE),
                    quantile(res_mat[, 5], probs = .975, na.rm = TRUE))

        rmat <- rbind(c(obs_rr, obsci_rr[1], obsci_rr[2]),
                      c(obs_or, obsci_or[1], obsci_or[2]))
        rownames(rmat) <- c(" Observed Relative Risk:",
                            "    Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
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
    }

    res <- list(obs_data = tab,
                obs_measures = rmat,
                adj_measures = rmatc,
                sim_df = as.data.frame(res_mat),
                reps = reps,
                fun = "probsens",
                warnings = neg_warn
                )
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
