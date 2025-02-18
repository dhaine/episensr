#' @rdname misclass
#' @param df Default dataset, a data frame, to use for analysis.
#' @param x,y,... <[`data-masking`][rlang::topic-data-mask]> List of name-value
#' pairs describing which variables in the data should be used. The expression
#' `variable` is evaluated within `df`, so there is no need to refer to the
#' original dataset (i.e., use `probcase(df, variable)` instead of
#' `probcase(df, df$variable)`). Note that `x` and `y` should be numeric
#' dichotomous 0,1 variables.
#' @param measure Risk ratio or odds ratio.
#' @param type Choice of misclassification:
#'   \enumerate{
#'   \item exposure: bias analysis for exposure misclassification; corrections
#'   using sensitivity and specificity: nondifferential and independent errors,
#'   \item outcome: bias analysis for outcome misclassification.
#'   }
#' @param reps Number of replications to run.
#' @param seca List defining sensitivity among cases:
#'   \enumerate{
#'   \item The sensitivity of exposure classification among those with the outcome
#'   (when \code{type = "exposure"}), or
#'   \item The sensitivity of outcome classification among those with the exposure
#'   (when \code{type = "outcome"}).
#'   }
#'   The first argument provides the probability distribution function (constant,
#'   uniform, triangular, trapezoidal, truncated normal, or beta) and the second
#'   its parameters as a vector. Lower and upper bounds of the truncated normal
#'   have to be between 0 and 1.
#'   \enumerate{
#'   \item constant: constant value,
#'   \item uniform: min, max,
#'   \item triangular: lower limit, upper limit, mode,
#'   \item trapezoidal: min, lower mode, upper mode, max,
#'   \item normal: lower bound, upper bound, mean, sd.
#'   \item beta: alpha, beta.
#'   }
#' @param seexp List defining sensitivity among controls:
#'   \enumerate{
#'   \item The sensitivity of exposure classification among those without the
#'   outcome (when \code{type = "exposure"}), or
#'   \item The sensitivity of outcome classification among those without the
#'   exposure (when \code{type = "outcome"}).
#'   }
#' @param spca List as above for \code{seca} but for specificity.
#' @param spexp List as above for \code{seexp} but for specificity.
#' @param corr_se Correlation between case and non-case sensitivities. If PPV/NPV is
#' chosen in case of exposure misclassification, correlations are set to NULL.
#' @param corr_sp Correlation between case and non-case specificities.
#' @param alpha Significance level.
#' @return A list with elements (for `probcase()`):
#'   \item{obs_data}{The analyzed 2 x 2 table from the observed data.}
#'   \item{obs_measures}{A table of observed relative risk and odds ratio with
#'   confidence intervals.}
#'   \item{adj_measures}{A table of corrected relative risks and odds ratios.}
#'   \item{sim_df}{Data frame of random parameters and computed values at each run.}
#'   \item{reps}{Number of replications.}
#' @examples
#' # Fox M.P., MacLehose R.F., Lash T.L.
#' # SAS and R code for probabilistic quantitative bias analysis for
#' # misclassified binary variables and binary unmeasured confounders
#' # Int J Epidemiol 2023:1624-1633.
#' \dontrun{
#' a <- 40; b <- 20; c <- 60; d <- 80
#' D <- data.frame(e_obs = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
#' d = c(rep(1, a), rep(1, b), rep(0, c), rep(0, d)))
#' set.seed(1234)
#' probcase(D, x = e_obs, y = d, reps = 10^3,
#' measure = "RR", type = "exposure",
#' seca = list("beta", c(25, 3)),
#' spca = list("trapezoidal", c(.9, .93, .97, 1)),
#' seexp = list("beta", c(45, 7)),
#' spexp = list("trapezoidal", c(.8, .83, .87, .9)),
#' corr_se = .8,
#' corr_sp = .8)
#'
#' set.seed(1234)
#' probcase(D, x = e, y = d_obs, reps = 10^3,
#' meausre = "RR", type = "outcome",
#' seca = list("beta", c(254, 24)),
#' spca = list("trapezoidal", c(.94, .96, .98, 1)),
#' seexp = list("beta", c(450, 67)),
#' spexp = list("trapezoidal", c(.9, .92, .93, .95)),
#' corr_se = .8, corr_sp = .8)
#'
#' library(aplore3)  # to get ICU data
#' data(icu)
#' icu$sta2 <- as.numeric(icu$sta) - 1
#' icu$inf2 <- as.numeric(icu$inf) - 1
#' probcase(icu, x = inf2, y = sta2, reps = 1000,
#' measure = "OR", type = "exposure",
#' seca = list("beta", c(20, 3)), spca = list("trapezoidal", c(.8, .83, .87, .9)),
#' seexp = list("beta", c(60, 9)), spexp = list("trapezoidal", c(.88, .93, .97, 1)),
#' corr_se = .8, corr_sp = .8)
#' }
#' @export
#' @importFrom stats median pnorm qnorm quantile qunif runif rnorm rbinom qbeta rbeta reformulate glm binomial poisson coef confint
#' @importFrom sandwich sandwich
probrcpp <- function(df,
                     x,
                     y,
                     ...,
                     measure = c("RR", "OR"),
                     type = c("exposure", "outcome"),
                     reps = 100,
                     seca = list(dist = c("constant", "uniform",
                                          "triangular", "trapezoidal",
                                          "normal", "beta"),
                                 parms = NULL),
                     seexp = NULL,
                     spca = list(dist = c("constant", "uniform",
                                          "triangular", "trapezoidal",
                                          "normal", "beta"),
                                 parms = NULL),
                     spexp = NULL,
                     corr_se = NULL,
                     corr_sp = NULL,
                     alpha = 0.05) {
    x <- deparse(substitute(x))
    y <- deparse(substitute(y))
    argList <- as.character(match.call(expand.dots = FALSE)$...)

        obs_data <- df[, c(x, y, c(argList))]
    disease <- colnames(obs_data[1])
    outcome <- colnames(obs_data[2])
    confounder_names <- colnames(obs_data[-c(1:2)])
    tab_df <- table(obs_data[[x]], obs_data[[y]],
                    dnn = list(colnames(obs_data[1]), colnames(obs_data[2])))
    tab <- t(tab_df[2:1, 2:1])
    a <- as.numeric(tab[1, 1])
    b <- as.numeric(tab[1, 2])
    c <- as.numeric(tab[2, 1])
    d <- as.numeric(tab[2, 2])

    cli::cli_alert_info("Ran conventional analysis assuming no bias")
    formula_conv <- reformulate(c(disease, confounder_names), response = outcome)
    convrr_mod <- glm(formula_conv, data = obs_data, family = poisson(link = "log"))
    obs_rr <- exp(coef(convrr_mod))[2]
    convrr_cov <- sandwich::vcovHC(convrr_mod, type = "HC0")
    convrr_se <- sqrt(diag(convrr_cov))[2]
    obsci_rr <- c(exp(coef(convrr_mod)[2] - qnorm(1 - alpha / 2) * convrr_se),
                  exp(coef(convrr_mod)[2] + qnorm(1 - alpha / 2) * convrr_se))

    convor_mod <- glm(formula_conv, data = obs_data, family = binomial(link = "logit"))
    obs_or <- exp(coef(convor_mod))[2]
    convor_cov <- sandwich::vcovHC(convor_mod, type = "HC0")
    convor_se <- sqrt(diag(convor_cov))[2]
    obsci_or <- c(exp(coef(convor_mod)[2] - qnorm(1 - alpha / 2) * convor_se),
                  exp(coef(convor_mod)[2] + qnorm(1 - alpha / 2) * convor_se))

    draws <- matrix(NA, nrow = reps, ncol = 21)
    colnames(draws) <- c("seca", "seexp", "spca", "spexp",
                         "A0", "B0", "C0", "D0",
                         "flag",
                         "prevca", "prevexp",
                         "ppvca", "ppvexp", "npvca", "npvexp",
                         "ped1", "ped0",
                         "PPV_d1", "PPV_d0", "NPV_d1", "NPV_d0")
    corr_draws <- matrix(NA, nrow = reps, ncol = 4)

    se1 <- c(reps, seca[[2]])
    se0 <- c(reps, seexp[[2]])
    sp1 <- c(reps, spca[[2]])
    sp0 <- c(reps, spexp[[2]])

    ## Step3: Assign probability distributions to each bias parameter
    ## and Step 4a draw Se's and Sp's
    cli::cli_progress_step("Assigned probability distributions", spinner = TRUE)
    if (is.null(seexp) & !is.null(spca) &
        is.null(spexp) & is.null(corr_se) & is.null(corr_sp)) {
        if (seca[[1]] == "constant") {
            draws[, 1] <- seca[[2]]
        }
        if (seca[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(se1))
        }
        if (seca[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(se1))
        }
        if (seca[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(se1))
        }
        if (seca[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::rtruncnorm, as.list(se1))
        }
        if (seca[[1]] == "beta") {
            draws[, 1] <- do.call(rbeta, as.list(se1))
        }
        draws[, 2] <- draws[, 1]
        if (spca[[1]] == "constant") {
            draws[, 3] <- spca[[2]]
        }
        if (spca[[1]] == "uniform") {
            draws[, 3] <- do.call(runif, as.list(sp1))
        }
        if (spca[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::rtriangle, as.list(sp1))
        }
        if (spca[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(sp1))
        }
        if (spca[[1]] == "normal") {
            draws[, 3] <- do.call(truncnorm::rtruncnorm, as.list(sp1))
        }
        if (spca[[1]] == "beta") {
            draws[, 3] <- do.call(rbeta, as.list(sp1))
        }
        draws[, 4] <- draws[, 3]
    } else {
        norta_se <- matrix(c(1, corr_se, corr_se, 1), ncol = 2)
        norta_sp <- matrix(c(1, corr_sp, corr_sp, 1), ncol = 2)
        corr_draws[, 1:2] <- MASS::mvrnorm(reps, c(0, 0), norta_se)
        corr_draws[, 3:4] <- MASS::mvrnorm(reps, c(0, 0), norta_sp)
        corr_draws <- pnorm(corr_draws)

        if (seca[[1]] == "uniform") {
            draws[, 1] <- do.call(qunif, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::qtriangle,
                                  c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "normal") {
            draws[, 1] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seca[[1]] == "beta") {
            draws[, 1] <- do.call(qbeta, c(list(corr_draws[, 1]), as.list(se1[-1])))
        }
        if (seexp[[1]] == "uniform") {
            draws[, 2] <- do.call(qunif, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::qtriangle,
                                  c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "normal") {
            draws[, 2] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (seexp[[1]] == "beta") {
            draws[, 2] <- do.call(qbeta, c(list(corr_draws[, 2]), as.list(se0[-1])))
        }
        if (spca[[1]] == "uniform") {
            draws[, 3] <- do.call(qunif, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::qtriangle,
                                  c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "normal") {
            draws[, 3] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spca[[1]] == "beta") {
            draws[, 3] <- do.call(qbeta, c(list(corr_draws[, 3]), as.list(sp1[-1])))
        }
        if (spexp[[1]] == "uniform") {
            draws[, 4] <- do.call(qunif, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "triangular") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "trapezoidal") {
            draws[, 4] <- do.call(trapezoid::qtrapezoid,
                                  c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "normal") {
            draws[, 4] <- do.call(truncnorm::qtruncnorm,
                                  c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
        if (spexp[[1]] == "beta") {
            draws[, 4] <- do.call(qbeta, c(list(corr_draws[, 4]), as.list(sp0[-1])))
        }
    }

    measure <- match.arg(measure)
    type <- match.arg(type)
    if (type == "exposure") {
        ## Step 4b: Using simple bias analysis (A0, B0, C0, D0)
        cli::cli_progress_step("Simple bias analysis", spinner = TRUE)
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) / (draws[, 1] - (1 - draws[, 3]))  # A0
        draws[, 6] <- (a + b) - draws[, 5]  # B0
        draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) / (draws[, 2] - (1 - draws[, 4]))  #C0
        draws[, 8] <- (c + d) - draws[, 7]  # D0

        ## Clean up
        draws[, 9] <- apply(draws[, 5:8], MARGIN = 1, function(x) sum(x > 0))
        draws[, 9] <- ifelse(draws[, 9] != 4 | is.na(draws[, 9]), NA, 1)
        discard <- sum(is.na(draws[, 9]))
        if (sum(is.na(draws[, 9])) > 0) {
            cli::cli_alert_warning("Chosen Se/Sp distributions lead to {discard} impossible value{?s} which w{?as/ere} discarded.")
            neg_warn <- paste("Prior Se/Sp distributions lead to",  discard, "impossible value(s).")
        } else neg_warn <- NULL
        reps2 <- nrow(draws)

        ## Prevalence of exposure in cases and controls, accounting for sampling error
        ## For systematic error only
        draws[, 16] <- draws[, 5] / (draws[, 5] + draws[, 6])  # ped1 syst. err
        draws[, 17] <- draws[, 7] / (draws[, 7] + draws[, 8])  # ped0 syst. err.
        ## Computing predictive values for imputation (PPV_d1, PPV_d0, NPV_d1, NPV_d0)
        ## For systematic error only
        draws[, 18] <- draws[, 1] * draws[, 16] /
            (draws[, 1] * draws[, 16] + (1 - draws[, 16]) * (1 - draws[, 3]))  # PPV_d1
        draws[, 19] <- draws[, 2] * draws[, 17] /
            (draws[, 2] * draws[, 17] + (1 - draws[, 17]) * (1 - draws[, 4]))  # PPV_d0
        draws[, 20] <- (draws[, 3] * (1 - draws[, 16])) /
            ((draws[, 3] * (1 - draws[, 16])) + (1 - draws[, 1]) * (draws[, 16]))  # NPV_d1
        draws[, 21] <- (draws[, 4] * (1 - draws[, 17])) /
            ((draws[, 4] * (1 - draws[, 17])) + (1 - draws[, 2]) * (draws[, 17]))  # NPV_d0
        ## For systematic and random error
        suppressWarnings({
                             draws[, 10] <- rbeta(reps2, draws[, 5], draws[, 6])  # ped1 syst. + rdm err.
                             draws[, 11] <- rbeta(reps2, draws[, 7], draws[, 8])  # ped0 syst. + rdm err.
                         })
        ## Computing predictive values based on sampled prevalences (PPV_d1, PPV_d0, NPV_d1, NPV_d0)
        ## For systematic and random error
        draws[, 12] <- (draws[, 1] * draws[, 10]) /
            ((draws[, 1] * draws[, 10]) + (1 - draws[, 3]) * (1 - draws[, 10]))  # ppvca
        draws[, 13] <- (draws[, 2] * draws[, 11]) /
            ((draws[, 2] * draws[, 11]) + (1 - draws[, 4]) * (1 - draws[, 11]))  # ppvexp
        draws[, 14] <- (draws[, 3] * (1 - draws[, 10])) /
            ((1 - draws[, 1]) * draws[, 10] + draws[, 3] * (1 - draws[, 10]))  # npvca
        draws[, 15] <- (draws[, 4] * (1 - draws[, 11])) /
            ((1 - draws[, 2]) * draws[, 11] + draws[, 4] * (1 - draws[, 11]))  #npvexp

        ## Loop through draws and impute new exposures at each step
        cli::cli_progress_step("Rcpp!", spinner = TRUE)
        nrow_obs <- nrow(obs_data)
        names(obs_data)[1] <- "e_obs"
        names(obs_data)[2] <- "d"
        obs_mat <- as.matrix(obs_data[, c("e_obs", "d")])
        obs_mat <- cbind(obs_mat, obs_mat[, "e_obs"] * obs_mat[, "d"],  ## e_obs * d
                         obs_mat[, "e_obs"] * (1 - obs_mat[, "d"]),  ## e * (1 - d)
                         (1 - obs_mat[, "e_obs"]) * obs_mat[, "d"],  ## (1 - e_obs) * d
                         (1 - obs_mat[, "e_obs"]) * (1 - obs_mat[, "d"])  ## (1 - e_obs) * (1 - d)
                         )
        colnames(obs_mat) <- c("e_obs", "d", "e_d", "e_1d", "e1_d", "e1_1d")

        cli::cli_alert_info("Compute systematic and total error")
        if (measure == "RR") adj_risk <- calc_RRexpo(reps2, obs_mat, draws)
        if (measure == "OR") adj_risk <- calc_ORexpo(reps2, obs_mat, draws)

        ptest <- adj_risk[[1]]
        etest <- adj_risk[[2]]

#        res_mat2 <- cbind(exp(adj_risk[[2]][, 2]),
#                          exp(adj_risk[[1]][, 1]),
#                          exp(adj_risk[[1]][, 1] + adj_risk[[2]][, 1] * adj_risk[[1]][, 2]))

#        meas_syst <- c(median(res_mat[, 1], na.rm = TRUE),
#                       quantile(res_mat[, 1], probs = .025, na.rm = TRUE),
#                       quantile(res_mat[, 1], probs = .975, na.rm = TRUE))
#        meas_tot <- c(median(res_mat[, 3], na.rm = TRUE),
#                    quantile(res_mat[, 3], probs = .025, na.rm = TRUE),
#                    quantile(res_mat[, 3], probs = .975, na.rm = TRUE))

#        rmat <- rbind(c(obs_rr, obsci_rr[1], obsci_rr[2]),
#                      c(obs_or, obsci_or[1], obsci_or[2]))
#        rownames(rmat) <- c(" Observed Relative Risk:",
#                            "    Observed Odds Ratio:")
#        colnames(rmat) <- c(" ",
#                            paste(100 * (alpha / 2), "%", sep = ""),
#                            paste(100 * (1 - alpha / 2), "%", sep = ""))
#        if (is.null(rownames(tab)))
#            rownames(tab) <- paste("Row", 1:2)
#        if (is.null(colnames(tab)))
#            colnames(tab) <- paste("Col", 1:2)
#        rmatc <- rbind(meas_syst, meas_tot)
#        if (measure == "RR") {
#            rownames(rmatc) <- c("Relative Risk -- systematic error:",
#                                 "                      total error:")
#        }
#        if (measure == "OR") {
#            rownames(rmatc) <- c("Odds Ratio -- systematic error:",
#                                 "                   total error:")
#        }
#        colnames(rmatc) <- c("Median", "p2.5", "p97.5")
    }

    res <- list(#obs_data = tab,
                ptest = ptest,
                etest = etest
#                obs_measures = rmat,
#                adj_measures = rmatc,
#                sim_df = as.data.frame(res_mat2),
#                reps = reps,
#                fun = "probsens",
#                warnings = neg_warn
                )
    class(res) <- c("episensr", "episensr.probsens", "list")
    res
}
