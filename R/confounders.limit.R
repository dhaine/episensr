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
#' @param crude.RR Crude relative risk between the exposure and the outcome.
#'
#' @return A list with elements:
#' \item{model}{Bias analysis performed.}
#' \item{bias.parms}{Input bias parameters.}
#' \item{adj.measures}{Output results.}
#'
#' @references
#' Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.59--78, Springer.
#'
#' Flanders, W. Dana, Khoury, Muin J., 1990. Indirect Assessment of
#' Confounding: Graphic Description and Limits on Effect of Adjusting for
#' Covariates. \emph{Epidemiology} 1(3): 239--246.
#'
#' @examples
#' confounders.limit(OR = 1.65, crude.RR = 1.5)
#' @export
confounders.limit <- function(p = NA,
                              RR = NA,
                              OR = NA,
                              crude.RR = NULL) {
    if(is.null(crude.RR))
        stop('Please provide crude relative risk.')
    if(is.na(p) & is.na(RR) & is.na(OR))
        stop('Not enough information.')

    q <- ifelse(is.null(p), NA, 1 - p)
    lower.bound <- crude.RR /
        min(RR, OR, 1/p, RR/(q+RR*p), OR/(q+OR*p), na.rm = TRUE)
    upper.bound <- crude.RR

    rmatc <- cbind(lower.bound, upper.bound)
    rownames(rmatc) <- "Limits on adjusted RR:"
    colnames(rmatc) <- c("Lower bound", "Upper bound")

    bias.parms <- matrix(c(p, RR, OR, crude.RR))
    colnames(bias.parms) <- " "
    rownames(bias.parms) <- c("p(Confounder+|Exposure-):",
                              "RR(Confounder-Disease):",
                              "OR(Confounder-Exposure):",
                              "Crude RR(Exposure-Disease):")

    res <- list(model = "confounder",
                bias.parms = bias.parms,
                adj.measures = rmatc)
    class(res) <- c("episensr.confounder", "episensr", "list")
    res
}
