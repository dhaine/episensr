#' Probabilistic sensitivity analysis using predictive values.
#'
#' Probabilistic sensitivity analysis to correct for exposure misclassification
#' or outcome misclassification and random error, using predictive values.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param reps Number of replications to run.
#' @param ppvca Vector providing beta distribution parameters (alpha and beta)
#' for the positive predictive values among cases.
#' @param ppvexp Vector providing beta distribution parameters (alpha and beta)
#' for the positive predictive values among controls.
#' @param npvca Vector as above for \code{ppvca} but for negative predictive value.
#' @param npvexp Vector as above for \code{ppvexp} but for negative predictive value.
#' @param alpha Significance level.
#'
#' @return A list with elements:
#' \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#' \item{obs_measures}{A table of observed relative risk and odds ratio with
#' confidence intervals.}
#' \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#' \item{sim_df}{Data frame of random parameters and computed values.}
#' \item{reps}{Number of replications.}
#'
#' @references
#' Fox, M.P, MacLehose, R.F, Lash, T.L., 2021 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.253--255, Springer.
#' @examples
#' # The data for this example come from:
#' # Greenland S., Salvan A., Wegman D.H., Hallock M.F., Smith T.J.
#' # A case-control study of cancer mortality at a transformer-assembly facility.
#' # Int Arch Occup Environ Health 1994; 66(1):49-54.
#' greenland <- matrix(c(45, 94, 257, 945), dimnames = list(c("BC+", "BC-"),
#' c("Smoke+", "Smoke-")), nrow = 2, byrow = TRUE)
#' set.seed(123)
#' # Exposure misclassification, non-differential
#' probsens(greenland, type = "exposure", reps = 20000,
#' seca = list("trapezoidal", c(.75, .85, .95, 1)),
#' spca = list("trapezoidal", c(.75, .85, .95, 1)))
#' @export
#' @importFrom stats median pnorm qnorm quantile qunif runif qbeta rbeta
probsens_pv <- function(case,
                        exposed,
                        reps = 1000,
                        ppvca = c(shape1, shape2),
                        ppvexp = c(shape1, shape2),
                        npvca = c(shape1, shape2),
                        npvexp = c(shape1, shape2),
                        alpha = 0.05) {
    if (reps < 1)
        stop(cli::format_error(c("x" = "Wrong number of replications: reps = {reps}",
                                 "i" = "reps must be >= 1")))

    if (is.null(ppvca) | is.null(ppvexp) | is.null(npvca) | is.null(npvexp))
        stop(cli::format_error(c("x" = "Missing argument(s) for PPV or NPV",
                                 "i" = "PPV and NPV for cases and controls have
to be provided.")))
    if ((length(ppvca) != 2) | (length(ppvexp) != 2) | (length(npvca) != 2) |
        length(npvexp) != 2)
        stop(cli::format_error(c("i" = "Check distribution parameters")))
    if (ppvca[1] < 0 | ppvca[2] < 0 | ppvexp[1] < 0 | ppvexp[2] < 0 | npvca[1] < 0 |
        npvca[2] < 0 | npvexp[1] < 0 | npvexp[2] < 0)
        stop(cli::format_error(c("i" = "Wrong arguments for your beta distribution.
Alpha and Beta should be > 0.")))

    if (!inherits(case, "episensr.probsens")) {
        if (inherits(case, c("table", "matrix")))
            tab <- case
        else {
            tab.df <- table(case, exposed)
            tab <- tab.df[2:1, 2:1]
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
    se_log_obs_rr <- sqrt((c / a) / (a + c) + (d / b) / (b + d))
    lci_obs_rr <- exp(log(obs_rr) - qnorm(1 - alpha / 2) * se_log_obs_rr)
    uci_obs_rr <- exp(log(obs_rr) + qnorm(1 - alpha / 2) * se_log_obs_rr)

    obs_or <- (a / b) / (c / d)
    se_log_obs_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lci_obs_or <- exp(log(obs_or) - qnorm(1 - alpha / 2) * se_log_obs_or)
    uci_obs_or <- exp(log(obs_or) + qnorm(1 - alpha / 2) * se_log_obs_or)

    draws <- matrix(NA, nrow = reps, ncol = 11)
    colnames(draws) <- c("PPV_case", "NPV_case", "PPV_ctrl", "NPV_ctrl",
                         "r0", "r1",
                         "q1", "q0", "psi",
                         "corr_RR", "corr_OR")

    cli::cli_progress_step("Assign probability distributions", spinner = TRUE)
    draws[, 1] <- rbeta(reps, ppvca[1], ppvca[2])
    draws[, 2] <- rbeta(reps, npvca[1], npvca[2])
    draws[, 3] <- rbeta(reps, ppvexp[1], ppvexp[2])
    draws[, 4] <- rbeta(reps, npvexp[1], npvexp[2])
    draws[, 5] <- rbeta(reps, c, d)
    draws[, 6] <- rbeta(reps, a, b)

    draws[, 7] <- draws[, 1] * draws[, 6] + (1 - draws[, 2]) * (1 - draws[, 6])
    draws[, 8] <- draws[, 3] * draws[, 5] + (1 - draws[, 4]) * (1 - draws[, 5])
    draws[, 9] <- logit(draws[, 7]) - logit(draws[, 8])
    draws[, 11] <- exp(draws[, 9])

    or_syst <- c(median(draws[, 11], na.rm = TRUE),
                 quantile(draws[, 11], probs = .025, na.rm = TRUE),
                 quantile(draws[, 11], probs = .975, na.rm = TRUE))

    if (!inherits(case, "episensr.probsens")) {
        tab <- tab
        rmat <- rbind(c(obs_rr, lci_obs_rr, uci_obs_rr),
                      c(obs_or, lci_obs_or, uci_obs_or))
        rownames(rmat) <- c(" Observed Relative Risk:",
                            "    Observed Odds Ratio:")
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
    rmatc <- or_syst
    neg_warn <- NULL
#    rownames(rmatc) <- c(#"Relative Risk -- systematic error:",
#                         #    "                      total error:",
#                             "   Odds Ratio -- systematic error:",
    #    "                      total error:"
#    )
#    colnames(rmatc) <- c("Median", "p2.5", "p97.5")

    cli::cli_progress_update()

    res <- list(obs.data = tab,
                obs.measures = rmat,
                adj.measures = rmatc,
                sim.df = as.data.frame(draws),
                reps = reps,
                fun = "probsens",
                warnings = neg_warn#,
#                message = discard_mess
                )
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
