#' Covariate misclassification
#'
#' `misclass_cov()` allows to provide adjusted measures of association corrected
#' for misclassification of a covariate (a potential confounder or effect measure
#' modifier).
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param covariate Covariate to stratify on.
#' @param bias_parms Vector defining the bias parameters. This vector has 4 elements
#'   between 0 and 1, in the following order:
#'   \enumerate{
#'   \item Sensitivity of confounder classification among those with the outcome,
#'   \item Sensitivity of confounder classification among those without the outcome,
#'   \item Specificity of confounder classification among those with the outcome,and
#'   \item Specificity of confounder classification among those without the outcome.
#'   }
#' @param alpha Significance level.
#'
#' @return A list with elements (for `misclass_cov()`):
#'   \item{obs_data}{The analyzed stratified 2 x 2 tables from the observed data.}
#'   \item{corr_data}{The expected stratified observed data given the true data
#'   assuming misclassification.}
#'   \item{obs_measures}{A table of observed relative risk and odds ratio with
#'   confidence intervals.}
#'   \item{adj_measures}{A table of adjusted relative risk and odds ratio.}
#'   \item{bias_parms}{Input bias parameters.}
#'
#' @family misclassification
#'
#' @references
#' Fox, M.P, MacLehose, R.F., Lash, T.L., 2021 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.176--179, Springer.
#'
#' @examples
#' # The data for this example come from:
#' # Berry, R.J., Kihlberg, R., and Devine, O. Impact of misclassification of in vitro
#' # fertilisation in studies of folic acid and twinning: modelling using population
#' # based Swedish vital records.
#' # BMJ, doi:10.1136/bmj.38369.437789.82 (published 17 March 2004)
#' misclass_cov(array(c(1319, 38054, 5641, 405546, 565, 3583, 781, 21958,
#' 754, 34471, 4860, 383588),
#' dimnames = list(c("Twins+", "Twins-"),
#' c("Folic acid+", "Folic acid-"), c("Total", "IVF+", "IVF-")),
#' dim = c(2, 2, 3)),
#' bias_parms = c(.6, .6, .95, .95))
#' @export
#' @importFrom stats qnorm
misclass_cov <- function(case,
                         exposed,
                         covariate,
                         bias_parms = NULL,
                         alpha = 0.05) {
    if (is.null(bias_parms))
        bias_parms <- c(1, 1, 1, 1)
    else bias_parms <- bias_parms
    if (length(bias_parms) != 4)
        stop(cli::format_error(c("i" = "The argument bias_parms should be made of
the following components: (1) Sensitivity of classification among those with the
outcome, (2) Sensitivity of classification among those without the outcome, (3)
Specificity of classification among those with the outcome, and (4) Specificity
of classification among those without the outcome.")))
    if (!all(bias_parms >= 0 & bias_parms <=1))
        stop(cli::format_error(c("x" = "Bias parameters should be between 0 and 1.")))

    if (inherits(case, c("table", "array")))
        tab <- case
    else {
        tab_df <- table(case, exposed, covariate)
        tab <- tab_df[2:1, 2:1, ]
    }

    a <- as.numeric(tab[1, 1, 1])
    b <- as.numeric(tab[1, 2, 1])
    c <- as.numeric(tab[2, 1, 1])
    d <- as.numeric(tab[2, 2, 1])

    Ac1 <- as.numeric(tab[1, 1, 2])
    Bc1 <- as.numeric(tab[1, 2, 2])
    Cc1 <- as.numeric(tab[2, 1, 2])
    Dc1 <- as.numeric(tab[2, 2, 2])

    Ac0 <- as.numeric(tab[1, 1, 3])
    Bc0 <- as.numeric(tab[1, 2, 3])
    Cc0 <- as.numeric(tab[2, 1, 3])
    Dc0 <- as.numeric(tab[2, 2, 3])

    obs_rr <- (a / (a + c)) / (b / (b + d))
    se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
    lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
    uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

    obs_or <- (a / c) / (b / d)
    se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
    uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

    SMR_RR <- a / (((Ac1 + Cc1) * Bc1 / (Bc1 + Dc1)) + ((Ac0 + Cc0) * Bc0 / (Bc0 + Dc0)))
    SMR_OR <- a / ((Cc1 * Bc1 / Dc1) + (Cc0 * Bc0 / Dc0))
    SMR_RR_C <- obs_rr / SMR_RR
    SMR_OR_C <- obs_or / SMR_OR
    MH_RR <- (Ac1 * (Bc1 + Dc1) / sum(tab[, , 2]) + Ac0 * (Bc0 + Dc0) / sum(tab[, , 3])) /
        (Bc1 * (Ac1 + Cc1) / sum(tab[, , 2]) + Bc0 * (Ac0 + Cc0) / sum(tab[, , 3]))
    MH_OR <- (Ac1 * Dc1 / sum(tab[, , 2]) + Ac0 * Dc0 / sum(tab[, , 3])) /
        (Bc1 * Cc1 / sum(tab[, , 2]) + Bc0 * Cc0 / sum(tab[, , 3]))
    MH_RR_C <- obs_rr / MH_RR
    MH_OR_C <- obs_or / MH_OR

    Ac1_Cr <- (Ac1 - (1 - bias_parms[3]) * a) / (bias_parms[1] - (1 - bias_parms[3]))
    Bc1_Cr <- (Bc1 - (1 - bias_parms[3]) * b) / (bias_parms[1] - (1 - bias_parms[3]))
    Cc1_Cr <- (Cc1 - (1 - bias_parms[4]) * c) / (bias_parms[2] - (1 - bias_parms[4]))
    Dc1_Cr <- (Dc1 - (1 - bias_parms[4]) * d) / (bias_parms[2] - (1 - bias_parms[4]))

    Ac0_Cr <- a - Ac1_Cr
    Bc0_Cr <- b - Bc1_Cr
    Cc0_Cr <- c - Cc1_Cr
    Dc0_Cr <- d - Dc1_Cr

    if (Ac1_Cr < 1 | Ac0_Cr < 1 | Bc1_Cr < 1 | Bc0_Cr < 1 | Cc1_Cr < 1 | Cc0_Cr < 1 |
        Dc1_Cr < 1 | Dc0_Cr < 1)
        stop(cli::format_error(c("x" = "Parameters chosen lead to negative cell(s)
in adjusted stratified 2x2 table.")))

    corr_tab1 <- matrix(c(Ac1_Cr, Bc1_Cr, Cc1_Cr, Dc1_Cr), nrow = 2, byrow = TRUE)
    corr_tab0 <- matrix(c(Ac0_Cr, Bc0_Cr, Cc0_Cr, Dc0_Cr), nrow = 2, byrow = TRUE)

    corr_SMR_RR <- a / (((Ac1_Cr + Cc1_Cr) * Bc1_Cr / (Bc1_Cr + Dc1_Cr)) +
                        ((Ac0_Cr + Cc0_Cr) * Bc0_Cr / (Bc0_Cr + Dc0_Cr)))
    corr_SMR_OR <- a / ((Cc1_Cr * Bc1_Cr / Dc1_Cr) + (Cc0_Cr * Bc0_Cr / Dc0_Cr))
    corr_RR_C <- obs_rr / corr_SMR_RR
    corr_OR_C <- obs_or / corr_SMR_OR
    corr_MH_RR <- (Ac1_Cr * (Bc1_Cr + Dc1_Cr) / sum(corr_tab1) + Ac0_Cr *
                   (Bc0_Cr + Dc0_Cr) / sum(corr_tab0)) /
        (Bc1_Cr * (Ac1_Cr + Cc1_Cr) / sum(corr_tab1) + Bc0_Cr * (Ac0_Cr + Cc0_Cr) /
         sum(corr_tab0))
    corr_MH_OR <- (Ac1_Cr * Dc1_Cr / sum(corr_tab1) + Ac0_Cr * Dc0_Cr / sum(corr_tab0)) /
        (Bc1_Cr * Cc1_Cr / sum(corr_tab1) + Bc0_Cr * Cc0_Cr / sum(corr_tab0))
    corr_MH_RR_C <- obs_rr / corr_MH_RR
    corr_MH_OR_C <- obs_or / corr_MH_OR

    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    if (is.null(rownames(tab))) {
        rownames(corr_tab1) <- paste("Row", 1:2)
        rownames(corr_tab0) <- paste("Row", 1:2)
        } else {
            rownames(corr_tab1) <- row.names(tab)
            rownames(corr_tab0) <- row.names(tab)
        }
    if (is.null(colnames(tab))) {
        colnames(corr_tab1) <- paste("Col", 1:2)
        colnames(corr_tab0) <- paste("Col", 1:2)
        } else {
            colnames(corr_tab1) <- colnames(tab)
            colnames(corr_tab0) <- colnames(tab)
        }
    rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr),
                  c(obs_or, lci_obs_or, uci_obs_or))
    rownames(rmat) <- c("Observed Relative Risk:",
                        "   Observed Odds Ratio:")
    colnames(rmat) <- c(" ",
                        paste(100 * (alpha / 2), "%", sep = ""),
                        paste(100 * (1 - alpha / 2), "%", sep = ""))
    rmata <- rbind(SMR_RR, SMR_RR_C, MH_RR, MH_RR_C,
                   SMR_OR, SMR_OR_C, MH_OR, MH_OR_C)
    rmatar <- rbind(corr_SMR_RR, corr_RR_C, corr_MH_RR, corr_MH_RR_C,
                    corr_SMR_OR, corr_OR_C, corr_MH_OR, corr_MH_OR_C)
    rmatc <- cbind(rmata, rmatar)
    rownames(rmatc) <- c("                      SMR RR adjusted for confounder:",
                         "   RR due to confounding by misclassified confounder:",
                         "          Mantel-Haenszel RR adjusted for confounder:",
                         "MH RR due to confounding by misclassified confounder:",
                         "                      SMR OR adjusted for confounder:",
                         "   OR due to confounding by misclassified confounder:",
                         "          Mantel-Haenszel OR adjusted for confounder:",
                         "MH OR due to confounding by misclassified confounder:")
    colnames(rmatc) <- c("Observed", "Corrected")

    res <- list(model = "misclass_cov",
                obs_data = tab,
                corr_data = c(corr_tab1, corr_tab0),
                obs_measures = rmat,
                adj_measures = rmatc,
                bias_parms = bias_parms)
    class(res) <- c("episensr", "list")
    res
}
