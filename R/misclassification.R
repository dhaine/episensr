#' Sensitivity analysis for disease or exposure misclassification.
#'
#' Simple sensitivity analysis for disease or exposure misclassification.
#' Confidence interval for odds ratio adjusted using sensitivity and specificity is
#' computed as in Chu et al. (2006), for exposure misclassification.
#'
#' For exposure misclassification, bias-adjusted measures are available using sensitivity
#' and specificity, or using predictive values.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param type Choice of misclassification:
#' \enumerate{
#' \item exposure: bias analysis for exposure misclassification; corrections using
#' sensitivity and specificity: nondifferential and independent errors,
#' \item exposure_pv: bias analysis for exposure misclassification; corrections using
#' PPV/NPV: nondifferential and independent errors,
#' \item outcome: bias analysis for outcome misclassification.
#' }
#' @param bias_parms Vector defining the bias parameters. This vector has 4 elements
#' between 0 and 1, in the following order:
#' \enumerate{
#' \item Sensitivity of exposure (when \code{type = "exposure"}) or outcome (when \code{type = "outcome"}) classification among those with the outcome (when \code{type = "exposure"}) or exposure (when \code{type = "outcome"}),
#' \item Sensitivity of exposure (or outcome) classification among those without the outcome (or exposure),
#' \item Specificity of exposure (or outcome) classification among those with the outcome (or exposure), and
#' \item Specificity of exposure (or outcome) classification among those without the outcome (or exposure).
#' }
#' If PPV/NPV is chosen in case of exposure misclassification, this vector is the following:
#' \enumerate{
#' \item Positive predictive value among those with the outcome,
#' \item Positive predictive value among those without the outcome,
#' \item Negative predictive value among those with the outcome,
#' \item Negative predictive value among those without the outcome.
#' }
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{corr.data}{The expected observed data given the true data assuming
#' misclassification.}
#' \item{obs.measures}{A table of observed relative risk and odds ratio with
#' confidence intervals.}
#' \item{adj.measures}{A table of adjusted relative risk and odds ratio with confidence
#' interval for odds ratio.}
#' \item{bias.parms}{Input bias parameters.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.79--108, Springer.
#'
#' Chu, H., Zhaojie, W., Cole, S.R., Greenland, S., \emph{Sensitivity analysis of
#' misclassification: A graphical and a Bayesian approach},
#' Annals of Epidemiology 2006;16:834-841.
#'
#' @examples
#' # The data for this example come from:
#' # Fink, A.K., Lash,  T.L. A null association between smoking during pregnancy
#' # and breast cancer using Massachusetts registry data (United States).
#' # Cancer Causes Control 2003;14:497-503.
#' misclassification(matrix(c(215, 1449, 668, 4296),
#' dimnames = list(c("Breast cancer+", "Breast cancer-"),
#' c("Smoker+", "Smoker-")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure",
#' bias_parms = c(.78, .78, .99, .99))
#'
#' misclassification(matrix(c(4558, 3428, 46305, 46085),
#' dimnames = list(c("AMI death+", "AMI death-"),
#' c("Male+", "Male-")),
#' nrow = 2, byrow = TRUE),
#' type = "outcome",
#' bias_parms = c(.53, .53, .99, .99))
#'
#' # The following example comes from Chu et al. Sensitivity analysis of
#' # misclassification: A graphical and a Bayesian approach.
#' # Annals of Epidemiology 2006;16:834-841.
#' misclassification(matrix(c(126, 92, 71, 224),
#' dimnames = list(c("Case", "Control"), c("Smoker +", "Smoker -")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure",
#' bias_parms = c(.94, .94, .97, .97))
#'
#' # The next example, using PPV/NPV, comes from Bodnar et al. Validity of birth
#' # certificate-derived maternal weight data.
#' # Paediatric and Perinatal Epidemiology 2014;28:203-212.
#' misclassification(matrix(c(599, 4978, 31175, 391851),
#' dimnames = list(c("Preterm", "Term"), c("Underweight", "Normal weight")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure_pv",
#' bias_parms = c(0.65, 0.74, 1, 0.98))
#' @export
#' @importFrom stats qnorm
misclassification <- function(case,
                              exposed,
                              type = c("exposure", "exposure_pv", "outcome"),
                              bias_parms = NULL,
                              alpha = 0.05) {
    if(is.null(bias_parms))
        bias_parms <- c(1, 1, 1, 1)
   else bias_parms <- bias_parms
    if(length(bias_parms) != 4)
        stop('The argument bias_parms should be made of the following components: (1) Sensitivity of classification among those with the outcome, (2) Sensitivity of classification among those without the outcome, (3) Specificity of classification among those with the outcome, and (4) Specificity of classification among those without the outcome.')
    if(!all(bias_parms >= 0 & bias_parms <=1))
        stop('Bias parameters should be between 0 and 1.')

    if(inherits(case, c("table", "matrix")))
        tab <- case
    else {
        tab.df <- table(case, exposed)
        tab <- tab.df[2:1, 2:1]
        }

    a <- as.numeric(tab[1, 1])
    b <- as.numeric(tab[1, 2])
    c <- as.numeric(tab[2, 1])
    d <- as.numeric(tab[2, 2])

    type <- match.arg(type)
    if (type == "exposure") {
        obs.rr <- (a/(a + c)) / (b/(b + d))
        se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
        uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)

        obs.or <- (a/b) / (c/d)
        se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
        uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

        A <- (a - (1 - bias_parms[3]) * (a + b)) / (bias_parms[1] - (1 - bias_parms[3]))
        C <- (c - (1 - bias_parms[4]) * (c + d)) / (bias_parms[2] - (1 - bias_parms[4]))
        B <- (a + b) - A
        D <- (c + d) - C

        if(A < 1 | B < 1 | C < 1 | D < 1)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table.')

        corr.tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr.rr <- (A/(A + C)) / (B/(B + D))
        corr.or <- (A/B) / (C/D)

        mle.corr.or <- ((a+((a+b)*(bias_parms[3]-1))) * (((c+d)*bias_parms[2])-c)) /
            ((c+((c+d)*(bias_parms[4]-1))) * (((a+b)*bias_parms[1])-a))

        se.corr.or <- sqrt(
                           (
                               (
                                   (a+b)*a*b*((bias_parms[1]+bias_parms[3]-1)^2)
                               ) /
                               (
                                   ((((a+b)*bias_parms[1])-a)^2) * ((((a+b)*bias_parms[3])-b)^2)
                               )
                           ) +
                           (
                               (
                                   (c+d)*c*d*((bias_parms[2]+bias_parms[4]-1)^2)
                               ) /
                               (
                                   ((((c+d)*bias_parms[2])-c)^2) * ((((c+d)*bias_parms[4])-d)^2)
                               )
                           )
        )
        lci.corr.or <- exp(log(mle.corr.or) - qnorm(1 - alpha/2) * se.corr.or)
        uci.corr.or <- exp(log(mle.corr.or) + qnorm(1 - alpha/2) * se.corr.or)

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(corr.tab) <- paste("Row", 1:2)
        } else {
            rownames(corr.tab) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(corr.tab) <- paste("Col", 1:2)
        } else {
           colnames(corr.tab) <- colnames(tab)
        }
        rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr), c(obs.or, lci.obs.or, uci.obs.or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(corr.rr, NA, NA), c(corr.or, lci.corr.or, uci.corr.or))
        rownames(rmatc) <- c("Misclassification Bias Corrected Relative Risk:",
                             "   Misclassification Bias Corrected Odds Ratio:")
        colnames(rmatc) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
    }

        if (type == "exposure_pv") {
        obs.rr <- (a/(a + c)) / (b/(b + d))
        se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
        uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)

        obs.or <- (a/b) / (c/d)
        se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
        uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

        A <- ((a * bias_parms[1]) + (b * (1 - bias_parms[3])))
        B <- a + b - A
        C <- ((c * bias_parms[2]) + (d * (1 - bias_parms[4])))
        D <- c + d - C

        if(A < 1 | B < 1 | C < 1 | D < 1)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table.')

        corr.tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr.rr <- (A/(A + C)) / (B/(B + D))
        corr.or <- (A/B) / (C/D)

        lci.corr.or <- NA
        uci.corr.or <- NA

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(corr.tab) <- paste("Row", 1:2)
        } else {
            rownames(corr.tab) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(corr.tab) <- paste("Col", 1:2)
        } else {
            colnames(corr.tab) <- colnames(tab)
        }
        rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr), c(obs.or, lci.obs.or, uci.obs.or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(corr.rr, NA, NA), c(corr.or, lci.corr.or, uci.corr.or))
        rownames(rmatc) <- c("Misclassification Bias Corrected Relative Risk:",
                             "   Misclassification Bias Corrected Odds Ratio:")
        colnames(rmatc) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
    }

    if (type == "outcome") {
        obs.rr <- (a/(a + c)) / (b/(b + d))
        se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
        uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)

        obs.or <- (a/b) / (c/d)
        se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
        uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

        A <- (a - (1 - bias_parms[3]) * (a + c)) / (bias_parms[1] - (1 - bias_parms[3]))
        B <- (b - (1 - bias_parms[4]) * (b + d)) / (bias_parms[2] - (1 - bias_parms[4]))
        C <- (a + c) - A
        D <- (b + d) - B

        if(A < 1 | B < 1 | C < 1 | D < 1)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table.')

        corr.tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr.rr <- (A/(A + C)) / (B/(B + D))
        corr.or <- (A/B) / (C/D)

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))) {
            rownames(corr.tab) <- paste("Row", 1:2)
        } else {
            rownames(corr.tab) <- row.names(tab)
        }
        if (is.null(colnames(tab))) {
            colnames(corr.tab) <- paste("Col", 1:2)
        } else {
            colnames(corr.tab) <- colnames(tab)
        }
        rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr), c(obs.or, lci.obs.or, uci.obs.or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(corr.rr, corr.or)
        rownames(rmatc) <- c("Misclassification Bias Corrected Relative Risk:",
                             "   Misclassification Bias Corrected Odds Ratio:")
        colnames(rmatc) <- " "
    }

    res <- list(model = "misclassification",
                type = type,
                obs.data = tab,
                corr.data = corr.tab,
                obs.measures = rmat,
                adj.measures = rmatc,
                bias.parms = bias_parms)
    class(res) <- c("episensr", "episensr.boot", "list")
    res
}
