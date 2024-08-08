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
#' @family selection
#'
#' @references Greenland S. Quantifying biases in causal models: classical
#' confounding vs. collider-stratification bias. Epidemiology 2003;14:300-6.
#' @examples
#' mbias(or = c(2, 5.4, 2.5, 1.5, 1),
#' var = c("HIV", "Circumcision", "Muslim", "Low CD4", "Participation"))
#' @export
mbias <- function(or,
                  var = c("y", "x", "a", "b", "m")) {
    if (is.null(or))
        stop(cli::format_error(c("x" = "Missing input bias parameters.")))
    else or <- or
    if (!is.vector(or))
        stop(cli::format_error(c("i" = "The argument or should be a vector of length 5.")))
    if (length(or) != 5)
        stop(cli::format_error(c("i" = "The argument or should be made of 5 components in the following order: (1) Odds ratio between A and the exposure E, (2) Odds ratio between A and the collider M, (3) Odds ratio between B and the collider M, (4) Odds ratio between B and the outcome Y, and (5) Odds ratio observed between the exposure E and the outcome Y.")))
    if (!all(or >= 0))
        stop(cli::format_error(c("x" = "Odds ratios should be greater than 0.")))

    or_ae <- or[1]
    or_am <- or[2]
    or_bm <- or[3]
    or_by <- or[4]
    or_ey <- or[5]

    mbias_me <- (or_ae * or_am * (1 / (1 + sqrt(or_ae * or_am))) +
                     1 - (1 / (1 + sqrt(or_ae * or_am)))) /
        ((or_ae * (1 / (1 + sqrt(or_ae * or_am))) +
              1 - (1 / (1 + sqrt(or_ae * or_am)))) *
             (or_am * (1 / (1 + sqrt(or_ae * or_am))) +
                  1 - (1 / (1 + sqrt(or_ae * or_am)))))
    mbias_my <- (or_bm * or_by * (1 / (1 + sqrt(or_bm * or_by))) +
                 1 - (1 / (1 + sqrt(or_bm * or_by)))) /
        ((or_bm * (1 / (1 + sqrt(or_bm * or_by))) +
          1 - (1 / (1 + sqrt(or_bm * or_by)))) *
         (or_by * (1 / (1 + sqrt(or_bm * or_by))) +
          1 - (1 / (1 + sqrt(or_bm * or_by)))))
    mbias_ey <- (mbias_me * mbias_my * (1 / (1 + sqrt(mbias_me * mbias_my))) +
                 1 - (1 / (1 + sqrt(mbias_me * mbias_my)))) /
        ((mbias_me * (1 / (1 + sqrt(mbias_me * mbias_my))) +
          1 - (1 / (1 + sqrt(mbias_me * mbias_my)))) *
         (mbias_my * (1 / (1 + sqrt(mbias_me * mbias_my))) +
          1 - (1 / (1 + sqrt(mbias_me * mbias_my)))))

    or_corr <- or_ey / mbias_ey

    res <- list(model = "mbias",
                mbias_parms = c(mbias_me, mbias_my, mbias_ey),
                adj_measures = or_corr,
                bias_parms = or,
                labels = var)
    class(res) <- c("mbias", "list")
    res
}
