#' Bounding the bias limits of unmeasured confounding.
#'
#' Function to elicit the limits on measures of effect corrected for an unmeasured
#' confounder when only some of the bias parameters are known. Crude relative
#' risk between exposure and outcome has minimally to be provided. Up to 3 other
#' parameters can be entered.
#'
#' @param p Proportion with the confounder among the unexposed group.
#' @param RR Relative risk between the confounder and the outcome.
#' @param OR Odds ratio between the confounder and the outcome.
#' @param crude_RR Crude relative risk between the exposure and the outcome.
#'
#' @return A list with elements:
#' \item{model}{Bias analysis performed.}
#' \item{bias_parms}{Input bias parameters.}
#' \item{adj_measures}{Output results.}
#'
#' @family confounding
#'
#' @references
#' Fox, M.P, MacLehose, R.F., Lash, T.L., 2021 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.129--131, Springer.
#'
#' Flanders, W. Dana, Khoury, Muin J., 1990. Indirect Assessment of
#' Confounding: Graphic Description and Limits on Effect of Adjusting for
#' Covariates. \emph{Epidemiology} 1(3): 239--246.
#'
#' @examples
#' confounders.limit(OR = 1.65, crude_RR = 1.5)
#' @export
confounders.limit <- function(p = NA,
                              RR = NA,
                              OR = NA,
                              crude_RR = NULL) {
    if (is.null(crude_RR))
        stop(cli::format_error(c("x" = "Please provide crude relative risk.")))
    if(is.na(p) & is.na(RR) & is.na(OR))
        stop(cli::format_error(c("x" = "Not enough information.")))

    q <- ifelse(is.null(p), NA, 1 - p)
    lower_bound <- crude_RR / min(RR, OR, 1 / p, RR / (q + RR * p), OR / (q + OR * p), na.rm = TRUE)
    upper_bound <- crude_RR

    rmatc <- cbind(lower_bound, upper_bound)
    rownames(rmatc) <- "Limits on adjusted RR:"
    colnames(rmatc) <- c("Lower bound", "Upper bound")

    bias_parms <- matrix(c(p, RR, OR, crude_RR))
    colnames(bias_parms) <- " "
    rownames(bias_parms) <- c("p(Confounder+|Exposure-):",
                              "RR(Confounder-Disease):",
                              "OR(Confounder-Exposure):",
                              "Crude RR(Exposure-Disease):")

    res <- list(model = "confounder",
                bias_parms = bias_parms,
                adj_measures = rmatc)
    class(res) <- c("episensr.confounder", "episensr", "list")
    res
}
