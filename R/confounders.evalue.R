#' Compute E-value to assess bias due to unmeasured confounder.
#'
#' Help to quantify the evidence strength for causality in presence of unmeasured
#' confounding. The E-value is the minimum strength of association that an unmeasured
#' counfounder would need to have with both the exposure and the outcome, conditional
#' on the measured covariates, to fully explain away a specific exposure-outcome
#' association.
#'
#' @param est Point estimate for the effect measure.
#' @param lower_ci Lower limit of the confidence interval for the association (relative
#' risk, odds ratio, hazard ratio, incidence rate ratio, risk differece).
#' @param upper_ci Upper limit of the confidence interval for the association (relative
#' risk, odds ratio, hazard ratio, incidence rate ratio, risk differece).
#' @param sd 
#' @param type Choice of effect measure (relative risk, and odds ratio or hazard ratio
#' for rare outcomes i.e. < 15% at end of follow-up -- RR; incidence rate ratio for count
#' and continuous outcomes --- IRR; odds ratio for common outcome -- ORc; hazard ratio
#' for common outcome i.e. > 15% at end of follow-up -- HRc; difference in continuous
#' outcomes -- diff; and risk difference -- RD).
#' @param true_est True estimate to assess E-value for. Default to 1 to assess
#' against null value. Set to a different value to assess for non-null hypotheses.
#' 
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{cfder.data}{The same table for Confounder +.}
#' \item{nocfder.data}{The same table for Confounder -.}
#' \item{obs.measures}{A table of relative risk with confidence intervals; for
#' Total, Confounder +, and Confounder -.}
#' \item{adj.measures}{A table of Standardized Morbidity Ratio and Mantel-Haenszel
#' estimates.}
#' \item{bias.parms}{Input bias parameters.}
#'
#' @references VanderWeele T.J and Ding P. Sensitivity analysis in observational research:
#' Introducing the E-value. Annals of Internal Medicine 2017;167:268-274.
#' 
#' @examples
#' # The data for this example come from:
#' # Victoria C.G., Smith P.G., Vaughan J.P., Nobre L.C., Lombardi C., Teixeira A.M.
#' # et al.
#' # Evidence for protection by breast-feeding against infant deaths from infectious
#' # diseases in Brazil.
#' # Lancet 1987;2:319-22.
#' confounders.evalue(est = 3.9, type = "RR")
#' @export
#' @importFrom stats qnorm
confounders.evalue <- function(est,
                               lower_ci = NULL,
                               upper_ci = NULL,
                               sd = NULL,
                               type = c("RR", "IRR", "ORc", "HRc", "diff", "RD"),
                               true_est = 1){
    if (length(type) > 1)
        stop('Choose between RR, OR, IRR, ORc, HRc, diff, or RD implementation.')

    if ((type != "RD" | type != "diff") & est < 0)
        stop('Association cannot be negative.')

    if ((type != "RD" | type != "diff") & est == 1)
        stop('Association = 1! No E-value.')

    if(!is.null(lower_ci) & !is.null(upper_ci)){
        if (!is.na(lower_ci) & !is.na(upper_ci) & lower_ci > upper_ci)
            stop("Lower confidence interval should be less than the upper limit.")
        if ((type != "diff") & (!is.na(lower_ci) & est < lower_ci) |
            (!is.na(upper_ci) & est > upper_ci))
            stop("Provided association is not within provided confidence interval.")        
        if ((type != "diff") & true_est == 1 &
            (!is.na(lower_ci) & lower_ci < 1) &
            (!is.na(upper_ci) & 1 < upper_ci))
            stop("True association within provided confidence interval: E-value = 1.")
    }

    if ((type != "RD" | type != "diff") & true_est < 0)
        stop("True association should not be negative.")
    
    if (!inherits(est, "numeric"))
        stop("Please provide a valid value for association measure.")



    .closest <- function(x, y) {
        x[which(abs(x - y) == min(abs(x - y)))]
    }

    if (!is.null(lower_ci) | !is.null(upper_ci)) {
        closest_ci <- .closest(c(lower_ci, upper_ci), true_est)
    } else {
        closest_ci <- NA
    }

    tab <- c(est, closest_ci)

    .compute_evalue <- function(x, true_x) {
        if (is.na(x)) return(NA)
        
        if (true_x != 1) {
            x = true_x / x
        }
        
        if (x < 1) {
            x <- 1 / x
        }
        return(x + sqrt(x * (x - 1)))
    }

    type <- match.arg(type)
    if (type == "RR") {
        e.value <- sapply(tab, function(x) .compute_evalue(x, true_x = true_est))
    }

    if (type == "ORc") {
        tab <- sqrt(tab)
        e.value <- sapply(tab, function(x) .compute_evalue(x, true_x = sqrt(true_est)))
    }

    if (type == "diff") {
        
    }

#    if (type == "HRc") {
#        assoc <- (1 - 0.5^sqrt(assoc)) / (1 - 0.5^sqrt(1 / assoc))
#        lower__ci <- (1 - 0.5^sqrt(lower_ci)) / (1 - 0.5^sqrt(1 / lower_ci))
#        upper_ci <- (1 - 0.5^sqrt(upper_ci)) / (1 - 0.5^sqrt(1 / upper_ci))
#        RR_true <- (1 - 0.5^sqrt(RR_true)) / (1 - 0.5^sqrt(1/RR_true))
#        assoc <- sqrt(assoc)
#        lower_ci <- sqrt(lower_ci)
#        upper_ci <- sqrt(upper_ci)
#        RR_true <- sqrt(RR_true)
#        
#        if (!is.na(lower_ci) & !is.na(upper_ci) & (lower_ci < 1) & (upper_ci > 1)) {
#            lci.e.value <- 1
#            uci.e.value <- 1
#        }
#        if ((!is.na(lower_ci) | !is.na(upper_ci))) {
#            if (assoc > RR_true) uci.e.value <- NA
#            if (assoc < RR_true) lci.e.value <- NA
#            if (assoc == RR_true) {
#                lci.e.value <- 1
#                uci.e.value <- NA
#            }
#        }
#    }

e.value
#    print(paste("With an observed risk ratio of", assoc, "an unmeasured confounder that was associated with both the outcome and the exposure by a risk ratioof", round(e.value, 2), "-fold each, above and beyond the measured confounders, could explain away the estimate, but weaker confounding could not."))
}
