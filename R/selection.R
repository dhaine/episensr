#' Selection bias.
#'
#' `selection()` and `probsens.sel()` allow to provide adjusted measures of
#' association corrected for selection bias.
#'
#' @section Simple bias analysis with `selection()`:
#' `selection()` allows you to run a simple sensitivity analysis to correct for
#' selection bias using estimates of the selection proportions.
#'
#' @section Probabilistic sensitivity analysis with `probsens.sel()`:
#' `probsens.sel()` performs a summary-level probabilistic sensitivity analysis to
#' correct for selection bias.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated
#' against.
#' @param exposed Exposure variable.
#' @param bias_parms Selection probabilities. Either a vector of 4 elements
#' between 0 and 1 defining the following probabilities in this order can be
#' provided:
#' \enumerate{
#' \item Selection probability among cases exposed (1),
#' \item Selection probability among cases unexposed (2),
#' \item Selection probability among noncases exposed (3), and
#' \item Selection probability among noncases unexposed (4).
#' }
#' or a single positive selection-bias factor which is the ratio of the exposed
#' versus unexposed selection probabilities comparing cases and noncases [(1*4)/(2*3)
#' from above].
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{model}{Bias analysis performed.}
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{corr_data}{The same table corrected for  selection proportions.}
#' \item{obs_measures}{A table of odds ratios and relative risk with confidence intervals.}
#' \item{adj_measures}{Selection bias corrected measures of outcome-exposure relationship.}
#' \item{bias_parms}{Input bias parameters: selection probabilities.}
#' \item{selbias_or}{Selection bias odds ratio based on the bias parameters chosen.}
#'
#' @family selection
#'
#' @references
#' Fox, M.P, MacLehose, R.F., Lash, T.L., 2021 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.90--91, 274--279, Springer.
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
#'
#'
#' selection(matrix(c(136, 107, 297, 165),
#' dimnames = list(c("UM+", "UM-"), c("Mobile+", "Mobile-")),
#' nrow = 2, byrow = TRUE),
#' bias_parms = 0.43)
#' @export
#' @importFrom stats qnorm
selection <- function(case,
                      exposed,
                      bias_parms = NULL,
                      alpha = 0.05) {
    if (is.null(bias_parms))
        bias_parms <- c(1, 1, 1, 1)
    else bias_parms <- bias_parms
    if (!is.vector(bias_parms))
        stop(cli::format_error(c("x" = "The argument bias_parms should be a vector of length 4.")))
    if (length(bias_parms) != 4 & length(bias_parms) != 1)
        stop(cli::format_error(c("x" = "The argument bias_parms should be made of either a) 4 components in the following order: (1) Selection probability among cases exposed, (2) Selection probability among cases unexposed, (3) Selection probability among noncases exposed, and (4) Selection probability among noncases unexposed; or b) the selection probability.")))
    if (length(bias_parms) == 4 & !all(bias_parms >= 0 & bias_parms <=1))
        stop(cli::format_error(c("x" = "Selection probabilities should be between 0 and 1.")))
    if (length(bias_parms) == 1 & !all(bias_parms >= 0))
        stop(cli::format_error(c("x" = "Selection probability should be positive.")))

    if (inherits(case, c("table", "matrix")))
        tab <- case
    else {
        tab_df <- table(case, exposed)
        tab <- tab_df[2:1, 2:1]
    }
    tab <- tab[1:2, 1:2]

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

    if (length(bias_parms) == 4) {
        A0 <- a / bias_parms[1]
        B0 <- b / bias_parms[2]
        C0 <- c / bias_parms[3]
        D0 <- d / bias_parms[4]

        tab_corr <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)
        rr_corr <- (A0 / (A0 + C0)) / (B0 / (B0 + D0))
        or_corr <- (A0 / B0) / (C0 / D0)
    } else {
        tab_corr <- matrix(c(NA, NA, NA, NA), nrow = 2, byrow = TRUE)
        rr_corr <- rr / bias_parms
        or_corr <- or / bias_parms
    }


    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    if (is.null(rownames(tab))) {
        rownames(tab_corr) <- paste("Row", 1:2)
        } else {
        rownames(tab_corr) <- rownames(tab)
    }
    if (is.null(colnames(tab))) {
        colnames(tab_corr) <- paste("Col", 1:2)
        } else {
        colnames(tab_corr) <- colnames(tab)
    }
    rmat <- rbind(c(rr, lci_rr, uci_rr), c(or, lci_or, uci_or))
    rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
    colnames(rmat) <- c(" ",
                        paste(100 * (alpha / 2), "%", sep = ""),
                        paste(100 * (1 - alpha / 2), "%", sep = ""))
    rmatc <- rbind(rr_corr, or_corr)
    rownames(rmatc) <- c("Selection Bias Corrected Relative Risk:",
                         "   Selection Bias Corrected Odds Ratio:")
    colnames(rmatc) <- " "
    if (length(bias_parms) == 4) {
        selbias <- (bias_parms[1] * bias_parms[4]) / (bias_parms[2] * bias_parms[3])
    } else {
        selbias <- bias_parms
    }
    res <- list(model = "selection",
                obs_data = tab,
                corr_data = tab_corr,
                obs_measures = rmat,
                adj_measures = rmatc,
                bias_parms = bias_parms,
                selbias_or = selbias)
    class(res) <- c("episensr", "episensr.boot", "list")
    res
}


#' @rdname selection
#' @param reps Number of replications to run.
#' @param case_exp If or_parms not provided, defines the selection probability
#' among case exposed. The first argument provides the probability distribution
#' function and the second its parameters as a vector:
#'   \enumerate{
#'   \item constant: constant value,
#'   \item uniform: min, max,
#'   \item triangular: lower limit, upper limit, mode,
#'   \item trapezoidal: min, lower mode, upper mode, max.
#'   \item normal: truncated normal with lower bound, upper bound, mean, sd,
#'   \item beta: alpha, beta.
#'   }
#' @param case_nexp Same among cases non-exposed.
#' @param ncase_exp Same among non-cases exposed.
#' @param ncase_nexp Same among non-cases non-exposed.
#' @param alpha Significance level.
#'
#' @return A list with elements (for `probsens.sel()`):
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs_measures}{A table of observed odds ratio with confidence intervals.}
#' \item{adj_measures}{A table of corrected odds ratios.}
#' \item{sim_df}{Data frame of random parameters and computed values.}
#' \item{reps}{Number of replications.}
#'
#' @examples
#' # The data for this example come from:
#' # Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N. et al.
#' # Population-based incidence estimates of uveal melanoma in Germany.
#' # Supplementing cancer registry data by case-control data.
#' # Eur J Cancer Prev 2006;15:165-70.
#' set.seed(1234)
#' probsens.sel(matrix(c(139, 114, 369, 377),
#' dimnames = list(c("Melanoma+", "Melanoma-"), c("Mobile+", "Mobile-")), nrow = 2, byrow = TRUE),
#' reps = 5000,
#' case_exp = list("beta", c(139, 5.1)),
#' case_nexp = list("beta", c(114, 11.9)),
#' ncase_exp = list("beta", c(369, 96.1)),
#' ncase_nexp = list("beta", c(377, 282.9)))
#' @export
#' @importFrom stats median qnorm quantile runif rlnorm rbeta
probsens.sel <- function(case,
                         exposed,
                         reps = 1000,
                         case_exp = list(dist = c("constant", "uniform",
                                                  "triangular", "trapezoidal",
                                                  "normal", "beta"),
                                         parms = NULL),
                         case_nexp = list(dist = c("constant", "uniform",
                                                   "triangular", "trapezoidal",
                                                   "normal", "beta"),
                                          parms = NULL),
                         ncase_exp = list(dist = c("constant", "uniform",
                                                   "triangular", "trapezoidal",
                                                   "normal", "beta"),
                                          parms = NULL),
                         ncase_nexp = list(dist = c("constant", "uniform",
                                                    "triangular", "trapezoidal",
                                                    "normal", "beta"),
                                           parms = NULL),
                         alpha = 0.05) {
    if (reps < 1)
        stop(cli::format_error(c("x" = "Wrong number of replications: reps = {reps}",
                                 "i" = "reps must be >= 1")))

    if (is.null(case_exp) | is.null(case_nexp) | is.null(ncase_exp) | is.null(ncase_nexp))
        stop(cli::format_error(c("x" = "Please provide selection probabilities.")))
    if (!is.null(case_exp[[1]]) & !is.list(case_exp))
        stop(cli::format_error(c("x" = "Please provide a list for case exposed parameters.")))
    if (!is.null(case_exp[[2]])) {
        if (case_exp[[1]] == "constant" & length(case_exp[[2]]) != 1)
            stop(cli::format_error(c("i" = "For constant value, please provide a single value.")))
        if (case_exp[[1]] == "uniform" & length(case_exp[[2]]) != 2)
            stop(cli::format_error(c("i" = "For uniform distribution, please provide vector of lower and upper limits.")))
        if (case_exp[[1]] == "uniform" & case_exp[[2]][1] >= case_exp[[2]][2])
            stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is greater than upper limit.")))
        if (case_exp[[1]] == "triangular" & length(case_exp[[2]]) != 3)
            stop(cli::format_error(c("x" = "For triangular distribution, please provide vector of lower, upper limits, and mode.")))
        if (case_exp[[1]] == "triangular" & ((case_exp[[2]][1] > case_exp[[2]][3]) |
                                             (case_exp[[2]][2] < case_exp[[2]][3])))
            stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
        if (case_exp[[1]] == "trapezoidal" & length(case_exp[[2]]) != 4)
            stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide
vector of lower limit, lower mode, upper mode, and upper limit.")))
        if (case_exp[[1]] == "trapezoidal" & ((case_exp[[2]][1] > case_exp[[2]][2]) |
                                              (case_exp[[2]][2] > case_exp[[2]][3]) |
                                              (case_exp[[2]][3] > case_exp[[2]][4])))
            stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
        if (case_exp[[1]] == "normal" & (length(case_exp[[2]]) != 4))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")))
        if (case_exp[[1]] == "normal" & length(case_exp[[2]]) == 4 &
            ((case_exp[[2]][1] >= case_exp[[2]][2]) | (case_exp[[2]][1] < 0)))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide sensible values for lower and upper bounds (lower bound >= 0; lower limit < upper limit).")))
        if ((case_exp[[1]] == "constant" | case_exp[[1]] == "uniform" |
             case_exp[[1]] == "triangular" | case_exp[[1]] == "trapezoidal") &
            !all(case_exp[[2]] >= 0 & case_exp[[2]] <= 1))
            stop(cli::format_error(c("x" = "Selection probability should be between 0 and 1.")))
        if (case_exp[[1]] == "beta" & (case_exp[[2]][1] < 0 | case_exp[[2]][1] < 0))
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")))
        if (case_exp[[1]] == "beta" & length(case_exp[[2]]) != 2)
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.")))
    }

    if (!is.null(case_nexp[[2]])) {
        if (case_nexp[[1]] == "constant" & length(case_nexp[[2]]) != 1)
            stop(cli::format_error(c("x" = "For constant value, please provide a single value.")))
        if (case_nexp[[1]] == "uniform" & length(case_nexp[[2]]) != 2)
            stop(cli::format_error(c("i" = "For uniform distribution, please provide vector of lower and upper limits.")))
        if (case_nexp[[1]] == "uniform" & case_nexp[[2]][1] >= case_nexp[[2]][2])
            stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is greater than upper limit.")))
        if (case_nexp[[1]] == "triangular" & length(case_nexp[[2]]) != 3)
            stop(cli::format_error(c("i" = "For triangular distribution, please provide vector of lower, upper limits, and mode.")))
        if (case_nexp[[1]] == "triangular" & ((case_nexp[[2]][1] > case_nexp[[2]][3]) |
                                              (case_nexp[[2]][2] < case_nexp[[2]][3])))
            stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
        if (case_nexp[[1]] == "trapezoidal" & length(case_nexp[[2]]) != 4)
            stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")))
        if (case_nexp[[1]] == "trapezoidal" & ((case_nexp[[2]][1] > case_nexp[[2]][2]) |
                                               (case_nexp[[2]][2] > case_nexp[[2]][3]) |
                                               (case_nexp[[2]][3] > case_nexp[[2]][4])))
            stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
        if (case_nexp[[1]] == "normal" & (length(case_nexp[[2]]) != 4))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")))
        if (case_nexp[[1]] == "normal" & length(case_nexp[[2]]) == 4 &
            ((case_nexp[[2]][1] >= case_nexp[[2]][2]) | (case_nexp[[2]][1] < 0)))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide sensible values for lower and upper bounds (lower bound >=0; lower limit < upper limit).")))
        if ((case_nexp[[1]] == "constant" | case_nexp[[1]] == "uniform" |
             case_nexp[[1]] == "triangular" | case_nexp[[1]] == "trapezoidal") &
            !all(case_nexp[[2]] >= 0 & case_nexp[[2]] <= 1))
            stop(cli::format_error(c("x" = "Selection probability should be between 0 and 1.")))
        if (case_nexp[[1]] == "beta" & (case_nexp[[2]][1] < 0 | case_nexp[[2]][1] < 0))
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")))
        if (case_nexp[[1]] == "beta" & length(case_nexp[[2]]) != 2)
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.")))
    }

    if (!is.null(ncase_exp[[2]])) {
        if (ncase_exp[[1]] == "constant" & length(ncase_exp[[2]]) != 1)
            stop(cli::format_error(c("x" = "For constant value, please provide a single value.")))
        if (ncase_exp[[1]] == "uniform" & length(ncase_exp[[2]]) != 2)
            stop(cli::format_error(c("i" = "For uniform distribution, please provide vector of lower and upper limits.")))
        if (ncase_exp[[1]] == "uniform" & ncase_exp[[2]][1] >= ncase_exp[[2]][2])
            stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is greater than upper limit.")))
        if (ncase_exp[[1]] == "triangular" & length(ncase_exp[[2]]) != 3)
            stop(cli::format_error(c("i" = "For triangular distribution, please provide vector of lower, upper limits, and mode.")))
        if (ncase_exp[[1]] == "triangular" & ((ncase_exp[[2]][1] > ncase_exp[[2]][3]) |
                                              (ncase_exp[[2]][2] < ncase_exp[[2]][3])))
            stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
        if (ncase_exp[[1]] == "trapezoidal" & length(ncase_exp[[2]]) != 4)
            stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")))
        if (ncase_exp[[1]] == "trapezoidal" & ((ncase_exp[[2]][1] > ncase_exp[[2]][2]) |
                                               (ncase_exp[[2]][2] > ncase_exp[[2]][3]) |
                                               (ncase_exp[[2]][3] > ncase_exp[[2]][4])))
            stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
        if (ncase_exp[[1]] == "normal" & (length(ncase_exp[[2]]) != 4))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")))
        if (ncase_exp[[1]] == "normal" & length(ncase_exp[[2]]) == 4 &
            ((ncase_exp[[2]][1] >= ncase_exp[[2]][2]) | (ncase_exp[[2]][1] < 0)))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide sensible values for lower and upper bounds (lower bound >=0; lower limit < upper limit).")))
        if ((ncase_exp[[1]] == "constant" | ncase_exp[[1]] == "uniform" |
             ncase_exp[[1]] == "triangular" | ncase_exp[[1]] == "trapezoidal") &
            !all(ncase_exp[[2]] >= 0 & ncase_exp[[2]] <= 1))
            stop(cli::format_error(c("x" = "Selection probability should be between 0 and 1.")))
        if (ncase_exp[[1]] == "beta" & (ncase_exp[[2]][1] < 0 | ncase_exp[[2]][1] < 0))
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")))
        if (ncase_exp[[1]] == "beta" & length(ncase_exp[[2]]) != 2)
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.")))
    }

    if (!is.null(ncase_nexp[[2]])) {
        if (ncase_nexp[[1]] == "constant" & length(ncase_nexp[[2]]) != 1)
            stop(cli::format_error(c("x" = "For constant value, please provide a single value.")))
        if (ncase_nexp[[1]] == "uniform" & length(ncase_nexp[[2]]) != 2)
            stop(cli::format_error(c("i" = "For uniform distribution, please provide vector of lower and upper limits.")))
        if (ncase_nexp[[1]] == "uniform" & ncase_nexp[[2]][1] >= ncase_nexp[[2]][2])
            stop(cli::format_error(c("x" = "Lower limit of your uniform distribution is greater than upper limit.")))
        if (ncase_nexp[[1]] == "triangular" & length(ncase_nexp[[2]]) != 3)
            stop(cli::format_error(c("i" = "For triangular distribution, please provide vector of lower, upper limits, and mode.")))
        if (ncase_nexp[[1]] == "triangular" & ((ncase_nexp[[2]][1] > ncase_nexp[[2]][3]) |
                                               (ncase_nexp[[2]][2] < ncase_nexp[[2]][3])))
            stop(cli::format_error(c("x" = "Wrong arguments for your triangular distribution.")))
        if (ncase_nexp[[1]] == "trapezoidal" & length(ncase_nexp[[2]]) != 4)
            stop(cli::format_error(c("i" = "For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.")))
        if (ncase_nexp[[1]] == "trapezoidal" & ((ncase_nexp[[2]][1] > ncase_nexp[[2]][2]) |
                                                (ncase_nexp[[2]][2] > ncase_nexp[[2]][3]) |
                                                (ncase_nexp[[2]][3] > ncase_nexp[[2]][4])))
            stop(cli::format_error(c("x" = "Wrong arguments for your trapezoidal distribution.")))
        if (ncase_nexp[[1]] == "normal" & (length(ncase_nexp[[2]]) != 4))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide vector of lower and upper bounds, mean and sd.")))
        if (ncase_nexp[[1]] == "normal" & length(ncase_nexp[[2]]) == 4 &
            ((ncase_nexp[[2]][1] >= ncase_nexp[[2]][2]) | (ncase_nexp[[2]][1] < 0)))
            stop(cli::format_error(c("i" = "For truncated normal distribution, please provide sensible values for lower and upper bounds (lower bound >=0; lower limit < upper limit).")))
        if ((ncase_nexp[[1]] == "constant" | ncase_nexp[[1]] == "uniform" |
             ncase_nexp[[1]] == "triangular" | ncase_nexp[[1]] == "trapezoidal") &
            !all(ncase_nexp[[2]] >= 0 & ncase_nexp[[2]] <= 1))
            stop(cli::format_error(c("x" = "Selection probability should be between 0 and 1.")))
        if (ncase_nexp[[1]] == "beta" & (ncase_nexp[[2]][1] < 0 | ncase_nexp[[2]][1] < 0))
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution. Alpha and Beta should be > 0.")))
        if (ncase_nexp[[1]] == "beta" & length(ncase_nexp[[2]]) != 2)
            stop(cli::format_error(c("x" = "Wrong arguments for your beta distribution.")))
    }


    if (!inherits(case, "episensr.probsens")) {
        if (inherits(case, c("table", "matrix")))
            tab <- case
        else {
            tab_df <- table(case, exposed)
            tab <- tab_df[2:1, 2:1]
        }

        a <- as.numeric(tab[1, 1])
        b <- as.numeric(tab[1, 2])
        c <- as.numeric(tab[2, 1])
        d <- as.numeric(tab[2, 2])
    } else {
        a <- as.numeric(case[[3]][, 1])
        b <- as.numeric(case[[3]][, 2])
        c <- as.numeric(case[[3]][, 3])
        d <- as.numeric(case[[3]][, 4])

        reps <- case[[4]]
    }

    cli::cli_alert_info("Calculating observed measures")
    obs_rr <- (a / (a + c)) / (b / (b + d))
    se_log_obs_rr <- sqrt(1/a + 1/b - 1/(a+c) - 1/(b+d))
    lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
    uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

    obs_or <- (a / b) / (c / d)
    se_log_obs_or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
    uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

    draws <- matrix(NA, nrow = reps, ncol = 18)
    colnames(draws) <- c("S1_1", "S0_1", "S1_0", "S0_0",
                         "A1", "B1", "C1", "D1",
                         "A0", "B0", "C0", "D0",
                         "flag", "corr_RR", "cor_OR", "tot_RR", "tot_OR", "reps")

    case1 <- c(reps, case_exp[[2]])
    case0 <- c(reps, case_nexp[[2]])
    ctrl1 <- c(reps, ncase_exp[[2]])
    ctrl0 <- c(reps, ncase_nexp[[2]])

    cli::cli_progress_step("Assign probability distributions", spinner = TRUE)
    if (case_exp[[1]] == "constant") {
        draws[, 1] <- case_exp[[2]]
    }
    if (case_exp[[1]] == "uniform") {
        draws[, 1] <- do.call(runif, as.list(case1))
    }
    if (case_exp[[1]] == "triangular") {
        draws[, 1] <- do.call(triangle::rtriangle, as.list(case1))
    }
    if (case_exp[[1]] == "trapezoidal") {
        draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(case1))
    }
    if (case_exp[[1]] == "normal") {
        draws[, 1] <- do.call(truncnorm::rtruncnorm, as.list(case1))
    }
    if (case_exp[[1]] == "beta") {
        draws[, 1] <- do.call(rbeta, as.list(case1))
    }

    if (case_nexp[[1]] == "constant") {
        draws[, 2] <- case_nexp[[2]]
    }
    if (case_nexp[[1]] == "uniform") {
        draws[, 2] <- do.call(runif, as.list(case0))
    }
    if (case_nexp[[1]] == "triangular") {
        draws[, 2] <- do.call(triangle::rtriangle, as.list(case0))
    }
    if (case_nexp[[1]] == "trapezoidal") {
        draws[, 2] <- do.call(trapezoid::rtrapezoid, as.list(case0))
    }
    if (case_nexp[[1]] == "normal") {
        draws[, 2] <- do.call(truncnorm::rtruncnorm, as.list(case0))
    }
    if (case_nexp[[1]] == "beta") {
        draws[, 2] <- do.call(rbeta, as.list(case0))
    }

    if (ncase_exp[[1]] == "constant") {
        draws[, 3] <- ncase_exp[[2]]
    }
    if (ncase_exp[[1]] == "uniform") {
        draws[, 3] <- do.call(runif, as.list(ctrl1))
    }
    if (ncase_exp[[1]] == "triangular") {
        draws[, 3] <- do.call(triangle::rtriangle, as.list(ctrl1))
    }
    if (ncase_exp[[1]] == "trapezoidal") {
        draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(ctrl1))
    }
    if (ncase_exp[[1]] == "normal") {
        draws[, 3] <- do.call(truncnorm::rtruncnorm, as.list(ctrl1))
    }
    if (ncase_exp[[1]] == "beta") {
        draws[, 3] <- do.call(rbeta, as.list(ctrl1))
    }

    if (ncase_nexp[[1]] == "constant") {
        draws[, 4] <- ncase_nexp[[2]]
    }
    if (ncase_nexp[[1]] == "uniform") {
        draws[, 4] <- do.call(runif, as.list(ctrl0))
    }
    if (ncase_nexp[[1]] == "triangular") {
        draws[, 4] <- do.call(triangle::rtriangle, as.list(ctrl0))
    }
    if (ncase_nexp[[1]] == "trapezoidal") {
        draws[, 4] <- do.call(trapezoid::rtrapezoid, as.list(ctrl0))
    }
    if (ncase_nexp[[1]] == "normal") {
        draws[, 4] <- do.call(truncnorm::rtruncnorm, as.list(ctrl0))
    }
    if (ncase_nexp[[1]] == "beta") {
        draws[, 4] <- do.call(rbeta, as.list(ctrl0))
    }

    cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
    draws[, 9] <- a / draws[, 1]
    draws[, 10] <- b / draws[, 2]
    draws[, 11] <- c / draws[, 3]
    draws[, 12] <- d / draws[, 4]

    draws[, 5] <- draws[, 9] - a
    draws[, 6] <- draws[, 10] - b
    draws[, 7] <- draws[, 11] - c
    draws[, 8] <- draws[, 12] - d

    draws[, 14] <- (draws[, 9] / (draws[, 9] + draws[, 11])) /
        (draws[, 10] / (draws[, 10] + draws[, 12]))
    draws[, 15] <- (draws[, 9] / draws[, 11]) / (draws[, 10] / draws[, 12])

    draws[, 18] <- runif(reps)

    cli::cli_progress_step("Incorporating random error", spinner = TRUE)
    draws[, 16] <- exp(log(draws[, 14]) - draws[, 18] * se_log_obs_rr)
    draws[, 17] <- exp(log(draws[, 15]) - draws[, 18] * se_log_obs_or)

    ## Clean up
    draws[, 13] <- apply(draws[, 5:8], MARGIN = 1, function(x) sum(x > 0))
    draws[, 13] <- ifelse(draws[, 13] != 4 | is.na(draws[, 13]), NA, 1)
    discard <- sum(is.na(draws[, 13]))
    if (sum(is.na(draws[, 13])) > 0) {
        cli::cli_alert_warning("Chosen distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
        neg_warn <- paste("Prior distributions lead to",  discard, "impossible value(s).")
    } else neg_warn <- NULL

    draws <- draws[draws[, 13] == 1 & !is.na(draws[, 13]), ]

    corr_RR <- c(median(draws[, 14], na.rm = TRUE),
                 quantile(draws[, 14], probs = .025, na.rm = TRUE),
                 quantile(draws[, 14], probs = .975, na.rm = TRUE))
    tot_RR <- c(median(draws[, 16], na.rm = TRUE),
                quantile(draws[, 16], probs = .025, na.rm = TRUE),
                quantile(draws[, 16], probs = .975, na.rm = TRUE))
    corr_OR <- c(median(draws[, 15], na.rm = TRUE),
                 quantile(draws[, 15], probs = .025, na.rm = TRUE),
                 quantile(draws[, 15], probs = .975, na.rm = TRUE))
    tot_OR <- c(median(draws[, 17], na.rm = TRUE),
                quantile(draws[, 17], probs = .025, na.rm = TRUE),
                quantile(draws[, 17], probs = .975, na.rm = TRUE))

    if (!inherits(case, "episensr.probsens")) {
        tab <- tab
        rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr),
                      c(obs_or, lci_obs_or, uci_obs_or))
        rownames(rmat) <- c("Observed Relative Risk:",
                            "   Observed Odds Ratio:")
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha / 2), "%", sep = ""),
                            paste(100 * (1 - alpha / 2), "%", sep = ""))
    } else {
        tab <- case[[1]]
        rmat <- case[[2]]
    }
    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    rmatc <- rbind(corr_RR, tot_RR, corr_OR, tot_OR)
    rownames(rmatc) <- c("Relative Risk -- systematic error:",
                         "                      total error:",
                         "   Odds Ratio -- systematic error:",
                         "                      total error:")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")

    cli::cli_progress_update()

    res <- list(obs_data = tab,
                obs_measures = rmat,
                adj_measures = rmatc,
                sim_df = as.data.frame(draws[, -c(13, 18)]),
                reps = reps,
                fun = "probsens.sel")
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
