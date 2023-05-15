#' Sensitivity analysis for unmeasured confounders based on external adjustment
#'
#' Sensitivity analysis to explore effect of residual confounding using simple
#' algebraic transformation. It provides the relative risk adjusted for unmeasured
#' confounders based on available external information (i.e. from the literature) on
#' the relation between confounders and outcome.
#'
#' @param RR "True" or fully adjusted exposure relative risk.
#' @param bias_parms Numeric vector defining the necessary bias parameters. This
#' vector has 4 elements, in the following order:
#' \enumerate{
#' \item the association between the confounder and the outcome (RR, relative risk),
#' \item the association between exposure category and the confounder (OR, odds ratio),
#' \item the prevalence of the confounder (between 0 and 1), and
#' \item the prevalence of the exposure (between 0 and 1).
#' }
#'
#' @return A list with elements:
#' \item{model}{Bias analysis performed.}
#' \item{bias.parms}{Input bias parameters.}
#' \item{adj.measures}{Output results, with bias as a percentage: (crude.RR - RR)/RR * 100.}
#'
#' @references Schneeweiss, S., 2006. Sensitivity analysis and external adjustment for
#' unmeasured confounders in epidemiologic database studies of therapeutics.
#' \emph{Pharmacoepidemiol Drug Safety} 15: 291-303.
#'
#' @examples
#' # Schneeweiss, S, Glynn, R.J., Tsai, E.H., Avorn, J., Solomon, D.H. Adjusting for
#' # unmeasured confounders in pharmacoepidemiologic claims data using external
#' # information. Epidemiology 2005; 16: 17-24.
#' confounders.ext(RR = 1, bias_parms = c(0.1, 0.9, 0.1, 0.4))
#' @export
confounders.ext <- function(RR,
                            bias_parms = NULL) {
    if (is.null(bias_parms))
        bias_parms <- c(1, 1, 1, 1)
    else bias_parms <- bias_parms
    if (length(bias_parms) != 4)
        stop('The argument bias_parms should be made of the following components: (1) Association between the confounder and the outcome, (2) Association between exposure category and confounder, (3) Prevalence of the confounder, and (4) Prevalence of the exposure.')
    if (!all(bias_parms[-c(1, 2)] >= 0 & bias_parms[-c(1, 2)] <= 1))
        stop('Prevalences should be between 0 and 1.')
    if (bias_parms[1] < 0)
        stop('Association between the confounder and the outcome should be >= 0.')
    if (bias_parms[2] < 0)
        stop('Association between exposure category and confounder should be >= 0.')
    if (RR < 0)
        stop('True relative risk should be greater than 0.')

    a <- bias_parms[2] - 1
    b <- (-bias_parms[3] * bias_parms[2]) - (bias_parms[4] * bias_parms[2]) +
        bias_parms[4] + bias_parms[3] - 1
    c <- bias_parms[3] * bias_parms[2] * bias_parms[4]
    P_C1 <- (-b - (sqrt(b^2 - (4 * a * c)))) / (2 * a)

    crude_RR <- ((P_C1 * (bias_parms[1] - 1) + bias_parms[4]) /
                 (((bias_parms[3] - P_C1) * (bias_parms[1] - 1)) -
                  bias_parms[4] + 1)) *
        ((1 - bias_parms[4]) / bias_parms[4])
    bias_perc <- (crude_RR - RR) / RR * 100

    rmatc <- rbind(crude_RR, bias_perc)
    rownames(rmatc) <- c("Crude RR", "Percent bias")
    colnames(rmatc) <- " "

    bias.parms <- matrix(bias_parms)
    colnames(bias.parms) <- " "
    rownames(bias.parms) <- c("RR(Confounder-Disease):",
                              "OR(Exposure category-Confounder):",
                              "p(Confounder):",
                              "p(Exposure):")

    res <- list(model = "confounder",
                bias.parms = bias.parms,
                adj.measures = rmatc)
    class(res) <- c("episensr.confounder", "episensr", "list")
    res
}
