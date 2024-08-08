#' Multidimensional sensitivity analysis for different sources of bias
#'
#' Multidimensional sensitivity analysis for different sources of bias,
#' where the bias analysis is repeated within a range of values for the
#' bias parameter(s).
#'
#' @param case Outcome variable. If a variable, this variable is tabulated
#' against.
#' @param exposed Exposure variable.
#' @param type Implement analysis for exposure misclassification, outcome
#' misclassification, unmeasured confounder, or selection bias.
#' @param se Numeric vector of sensitivities. Parameter used with exposure or outcome
#' misclassification.
#' @param sp Numeric vector of specificities. Parameter used with exposure or outcome
#' misclassification. Should be the same length as `se`.
#' @param bias_parms List of bias parameters used with unmeasured confounder. The list
#' is made of 3 vectors of the same length:
#' \enumerate{
#' \item Prevalence of Confounder in Exposure+ population,
#' \item Prevalence of Confounder in Exposure- population, and
#' \item Relative risk between Confounder and Outcome.
#' }
#' @param OR_sel Selection odds ratios, for selection bias implementation.
#' @param alpha Significance level.
#' @param dec Number of decimals in the printout.
#' @param print A logical scalar. Should the results be printed?
#'
#' @return A list with elements:
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs_measures}{A table of odds ratios and relative risk with confidence
#' intervals.}
#' \item{adj_measures}{Multidimensional corrected relative risk and/or odds ratio
#' data.}
#' \item{bias_parms}{Bias parameters.}
#'
#' @examples
#' multidimBias(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure",
#' se = c(1, 1, 1, .9, .9, .9, .8, .8, .8),
#' sp = c(1, .9, .8, 1, .9, .8, 1, .9, .8))
#'
#' multidimBias(matrix(c(45, 94, 257, 945),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "outcome",
#' se = c(1, 1, 1, .9, .9, .9, .8, .8, .8),
#' sp = c(1, .9, .8, 1, .9, .8, 1, .9, .8))
#'
#' multidimBias(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "confounder",
#' bias_parms = list(seq(.72, .92, by = .02),
#' seq(.01, .11, by = .01), seq(.13, 1.13, by = .1)))
#'
#' multidimBias(matrix(c(136, 107, 297, 165),
#' dimnames = list(c("Uveal Melanoma+", "Uveal Melanoma-"),
#' c("Mobile Use+", "Mobile Use -")),
#' nrow = 2, byrow = TRUE),
#' type = "selection",
#' OR_sel = seq(1.5, 6.5, by = .5))
#' @export
#' @importFrom stats qnorm setNames
multidimBias <- function(case,
                         exposed,
                         type = c("exposure", "outcome", "confounder", "selection"),
                         se = NULL,
                         sp = NULL,
                         bias_parms = NULL,
                         OR_sel = NULL,
                         alpha = 0.05,
                         dec = 4,
                         print = TRUE) {
    if (type %in% c("exposure", "outcome") && (!is.null(bias_parms) | !is.null(OR_sel)))
        stop(cli::format_error(c("x" = "Please provide adequate bias parameters (se and sp).")))
    if (type == "confounder" && (!is.null(se) | !is.null(sp) | !is.null(OR_sel)))
        stop(cli::format_error(c("x" = "Please provide adequate bias parameters (bias_parms).")))
    if (type == "selection" && (!is.null(se) | !is.null(sp) | !is.null(bias_parms)))
        stop(cli::format_error(c("x" = "Please provide adequate bias parameters (OR_sel).")))

    if (is.null(se))
        se <- c(1, 1)
    else se <- se
    if (is.null(sp))
        sp <- c(1, 1)
    else sp <- sp
    if (!is.vector(se))
        stop(cli::format_error(c("x" = "Sensitivity should be a vector.")))
    if (!is.vector(sp))
        stop(cli::format_error(c("x" = "Specificity should be a vector.")))
    if (!all(se >= 0 & se <=1))
        stop(cli::format_error(c("i" = "Sensitivity should be between 0 and 1.")))
    if (!all(sp >= 0 & sp <=1))
        stop(cli::format_error(c("i" = "Specificity should be between 0 and 1.")))
    if (length(se) != length(sp))
        stop(cli::format_error(c("x" = "Sensitivity and specificity should be vectors of the same length.")))

    if (is.null(bias_parms))
        bias_parms <- list(1, 1, 1)
    else bias_parms <- bias_parms
    if (!is.list(bias_parms))
        stop(cli::format_error(c("x" = "Bias parameters for the impact of unmeasured confounder should be provided as a list made of 3 elements.")))
    if (length(bias_parms) != 3)
        stop(cli::format_error(c("i" = "The argument bias should be made of the following vectors: (1) Prevalence of Confounder in Exposure+ population, (2) Prevalence of Confounder in Exposure- population, and (3) Relative risk between Confounder and Outcome.")))
    if (!all((bias_parms[[1]] >= 0 & bias_parms[[1]] <= 1) | (bias_parms[[2]] >= 0 & bias_parms[[2]] <= 1)))
        stop(cli::format_error(c("i" = "Prevalences should be between 0 and 1.")))
    if (length(bias_parms[[1]]) != length(bias_parms[[2]]) | length(bias_parms[[2]]) != length(bias_parms[[3]]))
        stop(cli::format_error(c("i" = "Prevalences of Confounder in Exposure+ and Exposure- populations and Relative risk between Confounder and Outcome should be of the same length.")))

    if (is.null(OR_sel))
        OR_sel <- 1
    else OR_sel <- OR_sel
    if (!is.vector(OR_sel))
        stop(cli::format_error(c("x" = "Selection odds ratios should be a vector.")))
    if (!all(OR_sel > 0))
        stop(cli::format_error(c("i" = "Selection odds ratios should be positive.")))

    if (inherits(case, c("table", "matrix")))
        tab <- case
    else {tab_df <- table(case, exposed)
          tab <- tab_df[2:1, 2:1]
         }

    a <- as.numeric(tab[1, 1])
    b <- as.numeric(tab[1, 2])
    c <- as.numeric(tab[2, 1])
    d <- as.numeric(tab[2, 2])

    rr <- (a / (a + c)) / (b / (b + d))
    se_log_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
    lci_rr <- exp(log(rr) - qnorm(1 - alpha / 2) * se_log_rr)
    uci_rr <- exp(log(rr) + qnorm(1 - alpha / 2) * se_log_rr)

    or <- (a / b) / (c / d)
    se_log_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lci_or <- exp(log(or) - qnorm(1 - alpha / 2) * se_log_or)
    uci_or <- exp(log(or) + qnorm(1 - alpha / 2) * se_log_or)

    rr_mat <- matrix(NA, nrow = length(se), ncol = length(se))
    or_mat <- matrix(NA, nrow = length(se), ncol = length(se))
    rrc_mat <- matrix(NA, nrow = length(bias_parms[[1]]), ncol = length(bias_parms[[1]]))
    ors_mat <- matrix(NA, nrow = length(OR_sel), ncol = 2)

    type <- match.arg(type)
    if (type == "exposure") {
        for (i in 1:nrow(rr_mat)) {
            for (j in 1:nrow(rr_mat)) {
                rr_mat[i, j] <- (((a - (1 - sp[j]) * (a + b)) / (se[j] - (1 - sp[j]))) /
                                 (((a - (1 - sp[j]) * (a + b)) /
                                   (se[j] - (1 - sp[j]))) +
                                  ((c - (1 - sp[i]) * (c + d))) /
                                  (se[i] - (1 - sp[i])))) /
                    (((a + b) - ((a - (1 - sp[j]) * (a + b)) /
                                 (se[j] - (1 - sp[j])))) /
                     (((a + b) - ((a - (1 - sp[j]) * (a + b)) /
                                  (se[j] - (1 - sp[j])))) +
                      ((c + d) - ((c - (1 - sp[i]) * (c + d)) /
                                  (se[i] - (1 - sp[i]))))))
            }
        }

        for (i in 1:nrow(or_mat)) {
            for (j in 1:nrow(or_mat)) {
                or_mat[i, j] <- (((a - (1 - sp[j]) * (a + b)) / (se[j] - (1 - sp[j]))) /
                                 (((c - (1 - sp[i]) * (c + d))) /
                                  (se[i] - (1 - sp[i])))) /
                    (((a + b) - ((a - (1 - sp[j]) * (a + b)) /
                                 (se[j] - (1 - sp[j])))) /
                     ((c + d) - ((c - (1 - sp[i]) * (c + d)) /
                                 (se[i] - (1 - sp[i])))))
            }
        }

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        rownames(rr_mat) <- paste("Se:", se, "Sp:", sp)
        colnames(rr_mat) <- paste("Se:", se, "Sp:", sp)
        rownames(or_mat) <- paste("Se:", se, "Sp:", sp)
        colnames(or_mat) <- paste("Se:", se, "Sp:", sp)
        if (print)
            cli::cli_h1("Multidimensional Exposure Misclassification")
        cli::cli_h2("Observed Data:")
        cli::cli_par()
        cli::cli_ul(c("Outcome: {rownames(tab)[1]}",
                      "Comparing: {colnames(tab)[1]} vs. {colnames(tab)[2]}"))
        cli::cli_end()
        print(round(tab, dec))
        rmat <- rbind(c(rr, lci_rr, uci_rr), c(or, lci_or, uci_or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                           sep = ""), "interval")
        rmatc <- list("Multidimensional Corrected Relative Risk Data" = rr_mat,
                      "Multidimensional Corrected Odds Ratio Data" = or_mat)
        cli::cli_par()
        print(round(rmat, dec))
        cli::cli_end()
        cli::cli_h2("Multidimensional Corrected Relative Risk Data:")
        cat("           Outcome + ->", "\nOutcome - |", "\n          V\n")
        print(rr_mat)
        cli::cli_h2("Multidimensional Corrected Odds Ratio Data:")
        cat("          Cases ->", "\nControls |", "\n         V\n")
        print(or_mat)
        bias_parms <- rbind(se, sp)
        rownames(bias_parms) <- c("Sensitivities:",
                                  "Specificities:")
        print(bias_parms)
        }

    if (type == "outcome") {
        for (i in 1:nrow(rr_mat)) {
            for (j in 1:nrow(rr_mat)) {
                rr_mat[i, j] <- (((a - (1 - sp[j]) * (a + c)) / (se[j] - (1 - sp[j]))) /
                                 ((a + c))) / (((b - (1 - sp[i]) * (b + d)) /
                                                (se[i] - (1 - sp[i]))) /
                                               ((b + d)))
            }
        }

        for (i in 1:nrow(or_mat)) {
            for (j in 1:nrow(or_mat)) {
                or_mat[i, j] <- (((a - (1 - sp[j]) * (a + c)) / (se[j] - (1 - sp[j]))) /
                                 ((a + c) - ((a - (1 - sp[j]) * (a + c)) /
                                             (se[j] - (1 - sp[j]))))) /
                    (((b - (1 - sp[i]) * (b + d)) / (se[i] - (1 - sp[i]))) /
                     ((b + d) - ((b - (1 - sp[i]) * (b + d)) /
                                 (se[i] - (1 - sp[i])))))
            }
        }

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        rownames(rr_mat) <- paste("Se:", se, "Sp:", sp)
        colnames(rr_mat) <- paste("Se:", se, "Sp:", sp)
        rownames(or_mat) <- paste("Se:", se, "Sp:", sp)
        colnames(or_mat) <- paste("Se:", se, "Sp:", sp)
        if (print)
            cli::cli_h1("Multidimensional Outcome Misclassification")
        cli::cli_h2("Observed Data:")
        cli::cli_par()
        cli::cli_ul(c("Outcome: {rownames(tab)[1]}",
                      "Comparing: {colnames(tab)[1]} vs. {colnames(tab)[2]}"))
        cli::cli_end()
        print(round(tab, dec))
        rmat <- rbind(c(rr, lci_rr, uci_rr), c(or, lci_or, uci_or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                           sep = ""), "interval")
        rmatc <- list("Multidimensional Corrected Relative Risk Data" = rr_mat,
                      "Multidimensional Corrected Odds Ratio Data" = or_mat)
        cli::cli_par()
        print(round(rmat, dec))
        cli::cli_end()
        cli::cli_h2("Multidimensional Corrected Relative Risk Data:")
        cat("           Outcome + ->", "\nOutcome - |", "\n          V\n")
        print(rr_mat)
        cli::cli_h2("Multidimensional Corrected Odds Ratio Data:")
        cat("          Cases ->", "\nControls |", "\n         V\n")
        print(or_mat)
        bias_parms <- rbind(se, sp)
        rownames(bias_parms) <- c("Sensitivities:",
                                  "Specificities:")
        print(bias_parms)
    }

    if (type == "confounder") {
        for (i in 1:nrow(rrc_mat)) {
            for (j in 1:nrow(rrc_mat)) {
                rrc_mat[i, j] <- rr / ((bias_parms[[1]][i] * (bias_parms[[3]][j] - 1) + 1) /
                                       (bias_parms[[2]][i] * (bias_parms[[3]][j] - 1) + 1))
            }
        }

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        rownames(rrc_mat) <- paste("p(Conf+|Exp+):", bias_parms[[1]],
                                   "p(Conf+|Exp-):", bias_parms[[2]])
        colnames(rrc_mat) <- paste("RR(Conf-Outc):", bias_parms[[3]])
        if (print)
            cli::cli_h1("Multidimensional Unmeasured Confounding")
        cli::cli_h2("Observed Data:")
        cli::cli_par()
        cli::cli_ul(c("Outcome: {rownames(tab)[1]}",
                      "Comparing: {colnames(tab)[1]} vs. {colnames(tab)[2]}"))
        cli::cli_end()
        print(round(tab, dec))
        rmat <- rbind(c(rr, lci_rr, uci_rr), c(or, lci_or, uci_or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                           sep = ""), "interval")
        rmatc <- rrc_mat
        cli::cli_par()
        print(round(rmat, dec))
        cli::cli_end()
        cli::cli_h2("Multidimensional Relative Risk Exposure-Data Relationship Adjusted for Confounder:")
        print(rrc_mat)
        bias_parms <- matrix(unlist(bias_parms),
                             dimnames = list(c("p(Confounder+|Exposure+):",
                                  "p(Confounder+|Exposure-):",
                                  "RR(Confounder-Outcome)")),
                             nrow = 3, byrow = TRUE)
        print(bias_parms)
        bias_parms <- setNames(bias_parms, c("p(Confounder+|Exposure+):",
                                             "p(Confounder+|Exposure-):",
                                             "RR(Confounder-Outcome)"))
    }

    if (type == "selection") {
        ors_mat[, 1] <- OR_sel
        for (i in 1:nrow(ors_mat)) {
            ors_mat[i, 2] <- or / OR_sel[i]
        }

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        colnames(ors_mat) <- paste(c("OR selection:", "OR corrected:"))
        if (print)
            cli::cli_h1("Multidimensional Selection Bias")
        cli::cli_h2("Observed Data:")
        cli::cli_par()
        cli::cli_ul(c("Outcome: {rownames(tab)[1]}",
                      "Comparing: {colnames(tab)[1]} vs. {colnames(tab)[2]}"))
        cli::cli_end()
        print(round(tab, dec))
        rmat <- rbind(c(rr, lci_rr, uci_rr), c(or, lci_or, uci_or))
        rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                           sep = ""), "interval")
        rmatc <- ors_mat
        cli::cli_par()
        print(round(rmat, dec))
        cli::cli_end()
        cli::cli_h2("Observed and Selection Bias Corrected Measures:")
        print(ors_mat)
        bias_parms <- ors_mat[, 1]
        }
    invisible(list(obs_data = tab,
                   obs_measures = rmat,
                   adj_measures = rmatc,
                   bias_parms = bias_parms))
}
