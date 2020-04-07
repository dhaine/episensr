#' Compute E-value to assess bias due to unmeasured confounder.
#'
#' Help to quantify the evidence strength for causality in presence of unmeasured
#' confounding. The E-value is the minimum strength of association that an unmeasured
#' confounder would need to have with both the exposure and the outcome, conditional
#' on the measured covariates, to fully explain away a specific exposure-outcome
#' association.
#'
#' @param est Point estimate for the effect measure. For difference in continuous outcomes,
#' it is the standardized effect size (i.e. mean of the outcome divided by its standard
#' deviation).
#' @param lower_ci Lower limit of the confidence interval for the association (relative
#' risk, odds ratio, hazard ratio, incidence rate ratio, risk difference).
#' @param upper_ci Upper limit of the confidence interval for the association (relative
#' risk, odds ratio, hazard ratio, incidence rate ratio, risk difference).
#' @param sd For difference in continuous outcomes, the standard error of the outcome
#' divided by its standard deviation.
#' @param type Choice of effect measure (relative risk, and odds ratio or hazard ratio
#' for rare outcomes i.e. < 15% at end of follow-up -- RR; odds ratio for common
#' outcome -- ORc; hazard ratio for common outcome i.e. > 15% at end of follow-up -- HRc;
#' difference in continuous outcomes, RR approximation -- diff_RR; difference in
#' continuous outcomes, OR approximation -- diff_OR).
#' @param true_est True estimate to assess E-value for. Default to 1 on risk scale
#' to assess against null value. Set to a different value to assess for non-null
#' hypotheses.
#' 
#' @return A matrix with the observed point estimate and closest confidence interval to
#' the null hypothesis, expressed as a relative risk, and their corresponding E-value.
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
#'
#' # The data for this example come from:
#' # Oddy W.H, Smith G.J., Jacony P. 
#' # A possible strategy for developing a model to account for attrition bias in a
#' # longitudinal cohort to investigate associations between exclusive breastfeeding and
#' # overweight and obesity at 20 years.
#' # Annals of Nutrition and Metabolism 2014;65:234-235.
#' confounders.evalue(est = 1.47, lower_ci = 1.12, upper_ci = 1.93, type = "ORc")
#'
#' # The data for this example come from:
#' # Reinisch J., Sanders S., Mortensen E., Rubin D.B.
#' # In-utero exposure to phenobarbital and intelligence deficits in adult men.
#' # Journal of the American Medical Association 1995;274:1518-1525
#' confounders.evalue(est = -0.42, sd = 0.14, type = "diff_RR")
#' @export
#' @importFrom stats qnorm
confounders.evalue <- function(est,
                               lower_ci = NULL,
                               upper_ci = NULL,
                               sd = NA,
                               type = c("RR", "ORc", "HRc", "diff_RR", "diff_OR"),
                               true_est = 1){
    if (length(type) > 1)
        stop('Choose between RR, ORc, HRc, diff_RR, or diff_OR implementation.')

    if (!(type %in% c("diff_RR", "diff_OR")) & est < 0)
        stop('Association cannot be negative.')

    if ((type != "diff_RR" | type != "diff_OR") & est == 1)
        stop('Association = 1! No E-value.')

    if(!is.null(lower_ci) & !is.null(upper_ci)){
        if (!is.na(lower_ci) & !is.na(upper_ci) & lower_ci > upper_ci)
            stop("Lower confidence interval should be less than the upper limit.")
        if ((type != "diff") & (!is.na(lower_ci) & est < lower_ci) |
            (!is.na(upper_ci) & est > upper_ci))
            stop("Provided association is not within provided confidence interval.")        
        if ((type != "diff_RR" | type != "diff_OR") & true_est == 1 &
            (!is.na(lower_ci) & lower_ci < 1) &
            (!is.na(upper_ci) & 1 < upper_ci))
            stop("True association within provided confidence interval: E-value = 1.")
    }

    if ((type != "diff_RR" | type != "diff_OR") & true_est < 0)
        stop("True association should not be negative.")
    
    if (!inherits(est, "numeric"))
        stop("Please provide a valid value for association measure.")

    if ((type == "diff_RR" | type == "diff_OR") & is.na(sd))
        stop("Please provide sd.")

    if (((type == "diff_RR" | type == "diff_OR") & !is.na(sd)) & sd < 0)
        stop("Standardized SE cannot be negative.")

    .closest <- function(x, y) {
        x[which(abs(x - y) == min(abs(x - y)))]
    }

    if ((type != "diff_RR" | type != "diff_OR") &
        (!is.null(lower_ci) | !is.null(upper_ci))) {
        closest_ci <- .closest(c(lower_ci, upper_ci), true_est)
    } else {
        closest_ci <- NA
    }

    if (type != "diff_RR" & type != "diff_OR") {
        tab <- c(est, closest_ci)
    } else if (type == "diff_RR") {
        tab <- c(exp(0.91*est),
                 .closest(c(exp(0.91*est - 1.78*sd), exp(0.91*est + 1.78*sd)),
                          true_est))
    } else if (type == "diff_OR") {
        tab <- c(exp(1.81*est),
                 .closest(c(exp(1.81*est - 3.55*sd), exp(1.81*est + 3.55*sd)),
                          true_est))
    }

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

    if (type == "HRc") {
        tab <- sapply(tab, function(x) (1 - 0.5^sqrt(x)) / (1 - 0.5^sqrt(1/x)))
        e.value <- sapply(tab,
                          function(x) .compute_evalue(x,
                                                      true_x = (1 - 0.5^sqrt(true_est)) /
                                                          (1 - 0.5^sqrt(1/true_est))))
    }

    if (type == "diff_RR") {
        if (true_est == 1) true_est <- 0
        e.value <- sapply(tab,
                          function(x) .compute_evalue(x, true_x = sqrt(exp(0.91*true_est))))
    }

    if (type == "diff_OR") {
        e.value <- sapply(tab,
                          function(x) .compute_evalue(x, true_x = sqrt(exp(1.81*true_est))))
    }

    res <- rbind(tab, e.value)
    colnames(res) <- c("Point estimate", "  CI closest to H_0")
    rownames(res) <- c("     RR:",
                       "E-value:")
    class(res) <- c("episensr.evalue", "episensr", "matrix")
    res
}
