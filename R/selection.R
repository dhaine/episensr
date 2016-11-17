#' Sensitivity analysis to correct for selection bias.
#'
#' Simple sensitivity analysis to correct for selection bias using estimates of
#' the selection proportions.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated
#' against.
#' @param exposed Exposure variable.
#' @param bias_parms Numeric vector defining the selection probabilities. This vector has 4 elements between 0 and 1, in the following order:
#' \enumerate{
#' \item Selection probability among cases exposed,
#' \item Selection probability among cases unexposed,
#' \item Selection probabillity among noncases exposed, and
#' \item Selection probability among noncases unexposed.
#' }
#' @param alpha Significance level.
#' 
#' @return A list with elements:
#' \item{obs.data}{The analysed 2 x 2 table from the observed data.}
#' \item{corr.data}{The same table corrected for  selection proportions.}
#' \item{obs.measures}{A table of odds ratios and relative risk with confidence intervals.}
#' \item{adj.measures}{Selection bias corrected measures of outcome-exposure relationship.}
#' \item{bias.parms}{Input bias parameters: selection probabilities.}
#'
#' @examples
#' # The data for this example come from:
#' # Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N.
#' # et al.
#' # Population-based incidence estimates of uveal melanoma in Germany. Supplementing
#' # cancer registry data by case-control data.
#' # Eur J Cancer Prev 2006;15:165-70.
#' selection(matrix(c(136, 107, 297, 165),
#' dimnames = list(c("UM+", "UM-"), c("Mobile+", "Mobile-")),
#' nrow = 2, byrow = TRUE),
#' bias_parms = c(.94, .85, .64, .25))
#' @export
#' @importFrom stats qnorm
selection <- function(case,
                      exposed,
                      selprob = NULL,
                      bias_parms = NULL,
                      alpha = 0.05) {
    if (!missing(selprob)) {
        warning("argument selprob is deprecated; please use bias_parms instead.",
                call. = FALSE)
        bias_parms <- selprob
    }
    if(is.null(bias_parms))
        bias_parms <- c(1, 1, 1, 1)
    else bias_pamrs <- bias_parms
    if(!is.vector(bias_parms))
        stop('The argument bias_parms should be a vector of length 4')
    if(length(bias_parms) != 4)
        stop('The argument bias_parms should be made of 4 components in the following order: (1) Selection probability among cases exposed, (2) Selection probability among cases unexposed, (3) Selection probability among noncases exposed, and (4) Selection probability among noncases unexposed')
    if(!all(bias_parms >= 0 & bias_parms <=1))
        stop('Selection probabilities should be between 0 and 1')

    if(inherits(case, c("table", "matrix")))
        tab <- case
    else {tab.df <- table(case, exposed)
          tab <- tab.df[2:1, 2:1]}
    tab <- tab[1:2, 1:2]

    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    rr <- (a/(a + c)) / (b/(b + d))
    se.log.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
    lci.rr <- exp(log(rr) - qnorm(1 - alpha/2) * se.log.rr)
    uci.rr <- exp(log(rr) + qnorm(1 - alpha/2) * se.log.rr)

    or <- (a/b) / (c/d)
    se.log.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.or <- exp(log(or) - qnorm(1 - alpha/2) * se.log.or)
    uci.or <- exp(log(or) + qnorm(1 - alpha/2) * se.log.or)

    A0 <- a / bias_parms[1]
    B0 <- b / bias_parms[2]
    C0 <- c / bias_parms[3]
    D0 <- d / bias_parms[4]

    tab.corr <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)
    rr.corr <- (A0/(A0 + C0)) / (B0/(B0 + D0))
    or.corr <- (A0/B0) / (C0/D0)
   
    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    if (is.null(rownames(tab))){
        rownames(tab.corr) <- paste("Row", 1:2)
        } else {
        rownames(tab.corr) <- row.names(tab)
    }
    if (is.null(colnames(tab))){ 
        colnames(tab.corr) <- paste("Col", 1:2)
        } else {
        colnames(tab.corr) <- colnames(tab)
    }
    rmat <- rbind(c(rr, lci.rr, uci.rr), c(or, lci.or, uci.or))
    rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
    colnames(rmat) <- c("",
                        paste(100 * (alpha/2), "%", sep = ""),
                        paste(100 * (1 - alpha/2), "%", sep = ""))
    rmatc <- rbind(rr.corr, or.corr)
    rownames(rmatc) <- c("Selection Bias Corrected Relative Risk:",
                         "   Selection Bias Corrected Odds Ratio:")
    colnames(rmatc) <- ""
    bias <- rbind(bias_parms[1], bias_parms[2], bias_parms[3], bias_parms[4])
    rownames(bias) <- c("Selection probability among cases exposed:",
                        "Selection probability among cases unexposed:",
                        "Selection probability among noncases exposed:",
                        "Selection probability among noncases unexposed:")
    colnames(bias) <- ""
    res <- list(obs.data = tab,
                corr.data = tab.corr,
                obs.measures = rmat,
                adj.measures = rmatc,
                bias.parms = bias)
    class(res) <- c("episensr", "list")
    res
}
