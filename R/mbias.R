#' Sensitivity analysis to correct for selection bias caused by M bias.
#'
#' Simple sensitivity analysis to correct for selection bias caused by M bias using estimates of the odds ratios relating the variables.
#'
#' @param or Vector defining the input bias parameters, in the following order:
#' \enumerate{
#' \item Odds ratio between A and the exposure,
#' \item Odds ratio between A and P,
#' \item Odds ratio between B and P,
#' \item Odds ratio between B and the outcome,
#' \item Odds ratio observed between the exposure and the outcome.
#' }
#' @param var Vector defining variable names, in the following order:
#' \enumerate{
#' \item Outcome,
#' \item Exposure,
#' \item A,
#' \item B,
#' \item P.
#' }
#' 
#' @return A list with elements:
#' \item{mbias.parms}{Maximum bias parameters.}
#' \item{adj.measures}{Selection bias corrected measures.}
#' \item{bias.parms}{Input bias parameters.}
#'
#' @examples
#' # The data for this example come from:
#' # Greenland S.
#' # Quantifying biases in causal models: classical confounding vs. collider-stratification bias.
#' # Epidemiology 2003;14:300-6.
#' mbias(or = c(2, 5.4, 2.5, 1.5, 1),
#' var = c("HIV", "Circumcision", "Muslim", "Low CD4", "Participation"))
#' @export
mbias <- function(or,
                  var) {
    if(is.null(or))
        stop('Missing input bias parameters.')
    else or <- or
    if(!is.vector(or))
        stop('The argument or should be a vector of length 5')
    if(length(or) != 5)
        stop('The argument or should be made of 5 components in the following order: (1) Odds ratio between A and the exposure, (2) Odds ratio between A and P, (3) Odds ratio between B and P, (4) Odds ratio between B and the outcome, and (5) Odds ratio observed between the exposure and the outcome.')
    if(!all(or >= 0))
        stop('Odds ratios should be greater than 0.')

    or.ae <- or[1]
    or.ap <- or[2]
    or.bp <- or[3]
    or.bd <- or[4]
    or.ed <- or[5]

    mbias.pe <- (or.ae * or.ap * (1 / (1 + sqrt(or.ae * or.ap))) +
                     1 - (1 / (1 + sqrt(or.ae * or.ap)))) /
        ((or.ae * (1 / (1 + sqrt(or.ae * or.ap))) +
              1 - (1 / (1 + sqrt(or.ae * or.ap)))) *
             (or.ap * (1 / (1 + sqrt(or.ae * or.ap))) +
                  1 - (1 / (1 + sqrt(or.ae * or.ap)))))
    mbias.pd <- (or.bp * or.bd * (1 / (1 + sqrt(or.bp * or.bd))) +
                     1 - (1 / (1 + sqrt(or.bp * or.bd)))) /
        ((or.bp * (1 / (1 + sqrt(or.bp * or.bd))) +
              1 - (1 / (1 + sqrt(or.bp * or.bd)))) *
             (or.bd * (1 / (1 + sqrt(or.bp * or.bd))) +
                  1 - (1 / (1 + sqrt(or.bp * or.bd)))))
    mbias.ed <- (mbias.pe * mbias.pd * (1 / (1 + sqrt(mbias.pe * mbias.pd))) +
                     1 - (1 / (1 + sqrt(mbias.pe * mbias.pd)))) /
        ((mbias.pe * (1 / (1 + sqrt(mbias.pe * mbias.pd))) +
              1 - (1 / (1 + sqrt(mbias.pe * mbias.pd)))) *
             (mbias.pd * (1 / (1 + sqrt(mbias.pe * mbias.pd))) +
                  1 - (1 / (1 + sqrt(mbias.pe * mbias.pd)))))

    or.corr <- or.ed / mbias.ed
   
    res <- list(mbias.parms = c(mbias.pe, mbias.pd, mbias.ed),
                adj.measures = or.corr,
                bias.parms = or,
                labels = var)
    class(res) <- c("mbias", "list")
    res
}
