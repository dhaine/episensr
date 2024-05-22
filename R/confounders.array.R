#' Sensitivity analysis for unmeasured confounders based on confounding imbalance
#' among exposed and unexposed
#'
#' Sensitivity analysis to explore effect of residual confounding using simple
#' algebraic transformation (array approach). It indicates the strength of an
#' unmeasured confounder and the necessary imbalance among exposure categories to
#' affect the observed (crude) relative risk.
#'
#' @param crude.risk Crude (apparent or observed) relative risk between the exposure
#' and the outcome. If type `RD`, this is the crude (observed) risk difference.
#' @param type Choice of implementation, for binary covariates, continuous
#' covariates, or on risk difference scale.
#' @param bias_parms Numeric vector defining the necessary bias parameters. This
#' vector has 3 elements, in the following order:
#' \enumerate{
#' \item the association between the confounder and the outcome (RR, relative risk),
#' \item the prevalence of the confounder among the exposed (between 0 and 1, if
#' type `binary`), or mean value of the confounder among the exposed (if type
#' `continuous` or `RD`), and
#' \item the prevalence of the confounder among the unexposed (between 0 and 1, if
#' type `binary`), or mean value of the confounder among the unexposed (if type
#' `continuous` or `RD`).
#' }
#'
#' @return A list with elements:
#' \item{model}{Bias analysis performed.}
#' \item{bias.parms}{Input bias parameters.}
#' \item{adj.measures}{Output results, with bias as a percentage: (crude.RR - risk_adj)/risk_adj * 100.}
#'
#' @references Schneeweiss, S., 2006. Sensitivity analysis and external adjustment
#' for unmeasured confounders in epidemiologic database studies of therapeutics.
#' \emph{Pharmacoepidemiol Drug Safety} 15: 291-303.
#'
#' @examples
#' # Example from Schneeweiss, S. Sensitivity analysis and external adjustment for
#' # unmeasured confounders in epidemiologic database studies of therapeutics.
#' # Pharmacoepidemiol Drug Safety 2006; 15: 291-303.
#' confounders.array(crude.risk = 1.5, type = "binary",
#' bias_parms = c(5.5, 0.5, 0.1))
#' # Examples from Patorno E., Gopalakrishnan, C., Franklin, J.M., Brodovicz, K.G.,
#' # Masso-Gonzalez, E., Bartels, D.B., Liu, J., and Schneeweiss, S. Claims-based
#' # studies of oral glucose-lowering medications can achieve balance in critical
#' # clinical variables only observed in electronic health records 2017; 20(4): 974-
#' # 984.
#' confounders.array(crude.risk = 1.5, type = "binary",
#' bias_parms = c(3.25, 0.333, 0.384))
#' confounders.array(crude.risk = 1.5, type = "continuous",
#' bias_parms = c(1.009, 7.8, 7.9))
#' confounders.array(crude.risk = 0.05, type = "RD", bias_parms = c(0.009, 8.5, 8))
#' @export
confounders.array <- function(crude.risk,
                              type = c("binary", "continuous", "RD"),
                              bias_parms = NULL) {
    if (length(type) > 1)
        stop("Choose between binary, continuous, or RD implementation.")

    if (is.null(bias_parms))
        bias_parms <- c(1, 1, 1)
    else bias_parms <- bias_parms
    if (length(bias_parms) != 3)
        stop("The argument bias_parms should be made of the following components: (1) Association between the confounder and the outcome, (2) Prevalence of the confounder among the exposed, and (3) Prevalence of the confounder among the unexposed.")
    if (!all(bias_parms[-1] >= 0 & bias_parms[-1] <= 1) & type == "binary")
        stop("Prevalences should be between 0 and 1.")
    if (bias_parms[1] < 0)
        stop("Association between the confounder and the outcome should be >= 0.")
    if (crude.risk < 0 & type != "RD")
        stop("Crude risk should be greater than 0.")
    if ((crude.risk > 1 | crude.risk < -1) & type == "RD")
        stop("Crude risk should be between -1 and 1.")

    type <- match.arg(type)
    if (type == "binary") {
        risk_adj <- crude.risk /
            (
                (bias_parms[2] * (bias_parms[1] - 1) + 1) /
                (bias_parms[3] * (bias_parms[1] - 1) + 1)
            )
        bias_perc <- (crude.risk - risk_adj) / risk_adj * 100
    }

    if (type == "continuous") {
        risk_adj <- exp(log(crude.risk) - (log(bias_parms[1]) * (bias_parms[2] - bias_parms[3])))
        bias_perc <- (crude.risk - risk_adj) / risk_adj * 100
    }

    if (type == "RD") {
        risk_adj <- crude.risk - (bias_parms[1] * (bias_parms[2] - bias_parms[3]))
        bias_perc <- (crude.risk - risk_adj) / risk_adj * 100
    }

    if (type != "RD") {
        rmatc <- rbind(risk_adj, bias_perc)
        rownames(rmatc) <- c("Adjusted RR", "Percent bias")
        colnames(rmatc) <- " "
    } else {
        rmatc <- rbind(risk_adj, bias_perc)
        rownames(rmatc) <- c("Adjusted RD", "Percent bias")
        colnames(rmatc) <- " "
    }

    bias.parms <- matrix(bias_parms)
    colnames(bias.parms) <- " "
    if (type != "RD") {
        rownames(bias.parms) <- c("RR(Confounder-Disease):",
                                  "p(Confounder+|Exposure+):",
                                  "p(Confounder+|Exposure-):")
    } else {
        rownames(bias.parms) <- c("RR(Confounder-Disease):",
                                  "mean(Confounder+|Exposure+):",
                                  "mean(Confounder+|Exposure-):")
    }

    res <- list(model = "confounder",
                bias.parms = bias.parms,
                adj.measures = rmatc)
    class(res) <- c("episensr.confounder", "episensr", "list")
    res
}
