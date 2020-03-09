#' Sensitivity analysis to correct for selection bias caused by M bias.
#'
#' Simple sensitivity analysis to correct for selection bias caused by M bias using
#' estimates of the odds ratios relating the variables.
#'
#' @param or Vector defining the input bias parameters, in the following order:
#' \enumerate{
#' \item Odds ratio between A and the exposure E,
#' \item Odds ratio between A and the collider M,
#' \item Odds ratio between B and the collider M,
#' \item Odds ratio between B and the outcome D,
#' \item Odds ratio observed between the exposure E and the outcome D.
#' }
#' @param var Vector defining variable names, in the following order:
#' \enumerate{
#' \item Outcome,
#' \item Exposure,
#' \item A,
#' \item B,
#' \item Collider.
#' }
#' 
#' @return A list with elements:
#' \item{model}{Bias analysis performed.}
#' \item{mbias.parms}{Three maximum bias parameters: in collider-exposure relationship
#' created by conditioning on the collider, in collider-outcome relationship created by
#' conditioning on the collider, and in exposure-outcome relationship created by
#' conditioning on the collider.}
#' \item{adj.measures}{Selection bias corrected odds ratio.}
#' \item{bias.parms}{Input bias parameters.}
#' \item{labels}{Variables' labels.}
#'
#' @references Greenland S. Quantifying biases in causal models: classical
#' confounding vs. collider-stratification bias. Epidemiology 2003;14:300-6. 
#' @examples
#' mbias(or = c(2, 5.4, 2.5, 1.5, 1),
#' var = c("HIV", "Circumcision", "Muslim", "Low CD4", "Participation"))
#' @export
mbias <- function(or,
                  var = c("y", "x", "a", "b", "m")) {
    if(is.null(or))
        stop('Missing input bias parameters.')
    else or <- or
    if(!is.vector(or))
        stop('The argument or should be a vector of length 5')
    if(length(or) != 5)
        stop('The argument or should be made of 5 components in the following order: (1) Odds ratio between A and the exposure E, (2) Odds ratio between A and the collider M, (3) Odds ratio between B and the collider M, (4) Odds ratio between B and the outcome Y, and (5) Odds ratio observed between the exposure E and the outcome Y.')
    if(!all(or >= 0))
        stop('Odds ratios should be greater than 0.')

    or.ae <- or[1]
    or.am <- or[2]
    or.bm <- or[3]
    or.by <- or[4]
    or.ey <- or[5]

    mbias.me <- (or.ae * or.am * (1 / (1 + sqrt(or.ae * or.am))) +
                     1 - (1 / (1 + sqrt(or.ae * or.am)))) /
        ((or.ae * (1 / (1 + sqrt(or.ae * or.am))) +
              1 - (1 / (1 + sqrt(or.ae * or.am)))) *
             (or.am * (1 / (1 + sqrt(or.ae * or.am))) +
                  1 - (1 / (1 + sqrt(or.ae * or.am)))))
    mbias.my <- (or.bm * or.by * (1 / (1 + sqrt(or.bm * or.by))) +
                     1 - (1/(1 + sqrt(or.bm * or.by)))) /
        ((or.bm * (1 / (1 + sqrt(or.bm * or.by))) +
              1 -(1 / (1 + sqrt(or.bm * or.by)))) *
             (or.by * (1 / (1 + sqrt(or.bm * or.by))) +
                  1 - (1 / (1 + sqrt(or.bm * or.by)))))
    mbias.ey <- (mbias.me * mbias.my * (1 / (1 + sqrt(mbias.me * mbias.my))) +
                     1 - (1 / (1 + sqrt(mbias.me * mbias.my)))) /
        ((mbias.me * (1 / (1 + sqrt(mbias.me * mbias.my))) +
              1 - (1 / (1 + sqrt(mbias.me * mbias.my)))) *
             (mbias.my * (1 / (1 + sqrt(mbias.me * mbias.my))) +
                  1 - (1 / (1 + sqrt(mbias.me * mbias.my)))))

    or.corr <- or.ey / mbias.ey
   
    res <- list(model = "mbias",
                mbias.parms = c(mbias.me, mbias.my, mbias.ey),
                adj.measures = or.corr,
                bias.parms = or,
                labels = var)
    class(res) <- c("mbias", "list")
    res
}
