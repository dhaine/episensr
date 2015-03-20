probsens <- function(exposed, case,
                     reps = 1000,
                     prev.exp = list(dist = c("uniform", "triangular",
                                         "trapezoidal"),
                                     parms = NULL),
                     prev.nexp = list(dist = c("uniform", "triangular",
                                               "trapezoidal"),
                                      parms = NULL),
                     risk = list(dist = c("uniform", "triangular",
                                        "trapezoidal"),
                               parms = NULL),
                     corr.p = NULL,
                     alpha = 0.05,
                     dec = 4,
                     print = TRUE){
    if(is.null(prev.exp) | is.null(prev.nexp))
        stop('Please provide prevalences among the exposed and unexposed.')
    if(is.null(risk))
        stop('Please provide risk of acquiring outcome.')
    if(!is.null(corr.p) && (corr.p == 0 | corr.p == 1))
        stop('Correlations should be > 0 and < 1.')
    if(!all(prev.exp >= 0 & prev.exp <=1))
        stop('Prevalences should be between 0 and 1.')
    if(!all(prev.nexp >= 0 & prev.nexp <=1))
        stop('Prevalences should be between 0 and 1.')
    
    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, cased)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    draws <- matrix(NA, nrow = reps, ncol = 43)
    corr.draws <- matrix(NA, nrow = reps, ncol = 4)

    p1 <- c(reps, prev.exp[[2]])
    p0 <- c(reps, prev.nexp[[2]])
    rr.cd <- c(reps, risk[[2]])
    
    obs.rr <- (a/(a + c)) / (b/(b + d))
    se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
    lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
    uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)

    obs.or <- (a/b) / (c/d)
    se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
    uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

    if (is.null(corr.sp)) {
        if (prev.exp[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(p1))
            }
        if (prev.exp[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(p1))
            }
        if (prev.exp[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(p1))
            }
        if (prev.nexp[[1]] == "uniform") {
            draws[, 2] <- do.call(runif, as.list(p0))
            }
        if (prev.nexp[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::rtriangle, as.list(p0))
            }
        if (prev.nexp[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::rtrapezoid, as.list(p0))
            }
    } else {
        corr.draws[, 1:2] <- apply(corr.draws[, 1:2],
                                   2,
                                   function(x) x = runif(reps))
        corr.draws[, 1:2] <- apply(corr.draws[, 1:2],
                                   2,
                                   function(x) log(x / (1 - x)))
        corr.draws[, 3] <- exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 2]) /
            (1 + (exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 2])))
        corr.draws[, 4] <- exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 2]) /
            (1 + (exp(sqrt(corr.p) * corr.draws[, 1] + sqrt(1 - corr.p) * corr.draws[, 2])))

    if (prev.exp[[1]] == "uniform") {
        draws[, 1] <- prev.exp[[2]][2] -
            (prev.exp[[2]][2] - prev.exp[[2]][1]) * corr.draws[, 3]
    }
    if (prev.exp[[1]] == "triangular" | prev.exp[[1]] == "trapezoidal") {
        draws[, 1] <- (corr.draws[, 3] *
            (prev.exp[[2]][4] + prev.exp[[2]][3] - prev.exp[[2]][1] - prev.exp[[2]][2]) + (prev.exp[[2]][1] + prev.exp[[2]][2])) / 2
        draws[, 1] <- ifelse(draws[, 1] < prev.exp[[2]][2],
                             prev.exp[[2]][1] + sqrt((prev.exp[[2]][2] - prev.exp[[2]][1]) * (2 * draws[, 1] - prev.exp[[2]][1] - prev.exp[[2]][2])),
                             draws[, 1])
        draws[, 1] <- ifelse(draws[, 1] > prev.exp[[2]][3],
                             prev.exp[[2]][4] - sqrt(2 * (prev.exp[[2]][4] - prev.exp[[2]][3]) * (draws[, 1] - prev.exp[[2]][3])),
                             draws[, 1])
    }
    if (prev.nexp[[1]] == "uniform") {
        draws[, 2] <- prev.nexp[[2]][2] -
            (prev.nexp[[2]][2] - prev.nexp[[2]][1]) * corr.draws[, 4]
    }
    if (prev.nexp[[1]] == "triangular" | prev.nexp[[1]] == "trapezoidal") {
        draws[, 2] <- (corr.draws[, 4] *
                           (prev.nexp[[2]][4] + prev.nexp[[2]][3] - prev.nexp[[2]][1] - prev.nexp[[2]][2]) + (prev.nexp[[2]][1] + prev.nexp[[2]][2])) / 2
        draws[, 2] <- ifelse(draws[, 2] < prev.nexp[[2]][2],
                             prev.nexp[[2]][1] + sqrt((prev.nexp[[2]][2] - prev.nexp[[2]][1]) * (2 * draws[, 2] - prev.nexp[[2]][1] - prev.nexp[[2]][2])),
                             draws[, 2])
        draws[, 2] <- ifelse(draws[, 2] > prev.nexp[[2]][3],
                             prev.nexp[[2]][4] - sqrt(2 * (prev.nexp[[2]][4] - prev.nexp[[2]][3]) * (draws[, 2] - prev.nexp[[2]][3])),
                             draws[, 2])
    }
    }
    if (rr.cd[[1]] == "uniform") {
        draws[, 3] <- do.call(runif, as.list(rr.cd))
    }
    if (rr.cd[[1]] == "triangular") {
        draws[, 3] <- do.call(triangle::rtriangle, as.list(rr.cd))
    }
    if (rr.cd[[1]] == "trapezoidal") {
        draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(rr.cd))
    }
    
    draws[, 40] <- runif(reps)

    draws[, 4] <- (a + c) * draws[, 1]
    draws[, 5] <- (b + d) * draws[, 2]
    draws[, 6] <- (draws[, 3] * draws[, 4] * a) /
        (draws[, 3] * draws[, 4] + (a + c) - draws[, 4])
    draws[, 7] <- (draws[, 3] * draws[, 5] * b) /
        (draws[, 3] * draws[, 5] + (b + d) - draws[, 5])
    draws[, 8] <- draws[, 4] - draws[, 6]
    draws[, 9] <- draws[, 5] - draws[, 7]
    draws[, 10] <- a + c - draws[, 4]
    draws[, 11] <- b + d - draws[, 5]
    draws[, 12] <- a - draws[, 6]
    draws[, 13] <- b - draws[, 7]
    draws[, 14] <- c - draws[, 8]
    draws[, 15] <- d - draws[, 9]

    draws[, 16] <- a /
        ((draws[, 4] * draws[, 7]/draws[, 5]) +
             (draws[, 10] * draws[, 13]/draws[, 11]))
    draws[, 17] <- (draws[, 6] * draws[, 5]/(draws[, 4] + draws[, 5]) +
                        draws[, 12] * draws[, 11]/(draws[, 10] + draws[, 11])) /
                   (draws[, 7] * draws[, 4]/(draws[, 4] + draws[, 5]) +
                        draws[, 13] * draws[, 10]/(draws[, 10] + draws[, 11])) 
    draws[, 18] <- (draws[, 6]/(draws[, 6] + draws[, 8])) /
        (draws[, 7]/(draws[, 7] + draws[, 9]))
    draws[, 19] <- (draws[, 12]/(draws[, 12] + draws[, 14])) /
        (draws[, 13]/(draws[, 13] + draws[, 15]))
    draws[, 20] <- obs.rr / draws[, 16]
    draws[, 21] <- obs.rr / draws[, 17]

    draws[, 22] <- c * draws[, 1] 
    draws[, 23] <- d * draws[, 2]
    draws[, 24] <- (draws[, 3] * draws[, 22] * a) /
        (draws[, 3] * draws[, 22] + c - draws[, 22])
    draws[, 25] <- (draws[, 3] * draws[, 23] * b) /
        (draws[, 3] * draws[, 23] + d - draws[, 23])
    draws[, 26] <- draws[, 24] + draws[, 22]
    draws[, 27] <- draws[, 25] + draws[, 23]
    draws[, 28] <- a - draws[, 24]
    draws[, 29] <- b - draws[, 25]
    draws[, 30] <- c - draws[, 22]
    draws[, 31] <- d - draws[, 23]
    draws[, 32] <- draws[, 28] + draws[, 30]
    draws[, 33] <- draws[, 29] + draws[, 30]

    draws[, 34] <- a /
        ((draws[, 22] * draws[, 25]/draws[, 23]) +
             (draws[, 30] * draws[, 29]/draws[, 31]))
    draws[, 35] <- (draws[, 24] * draws[, 23]/(draws[, 26] + draws[, 27]) +
                        draws[, 28] * draws[, 31]/(draws[, 32] + draws[, 33])) /
                   (draws[, 25] * draws[, 22]/(draws[, 26] + draws[, 27]) +
                        draws[, 29] * draws[, 30]/(draws[, 32] + draws[, 33])) 
    draws[, 36] <- (draws[, 24] / draws[, 22]) / (draws[, 25] / draws[, 23])
    draws[, 37] <- (draws[, 28] / draws[, 30]) / (draws[, 29] / draws[, 31])
    draws[, 38] <- obs.or / draws[, 34]
    draws[, 39] <- obs.or / draws[, 35]
    
    draws[, 20] <- ifelse(draws[, 6] < 1 |
                               draws[, 7] < 1 |
                                 draws[, 8] < 1 |
                                   draws[, 9] < 1 |
                          draws[, 12] < 1 |
                            draws[, 13] < 1 |
                              draws[, 14] < 1 |
                                draws[, 15] < 1, NA, draws[, 20])
    draws[, 21] <- ifelse(draws[, 6] < 1 |
                            draws[, 7] < 1 |
                              draws[, 8] < 1 |
                                draws[, 9] < 1 |
                          draws[, 12] < 1 |
                            draws[, 13] < 1 |
                              draws[, 14] < 1 |
                                draws[, 15] < 1, NA, draws[, 21])
    draws[, 38] <- ifelse(draws[, 24] < 1 |
                               draws[, 25] < 1 |
                                 draws[, 22] < 1 |
                                   draws[, 23] < 1 |
                          draws[, 28] < 1 |
                            draws[, 29] < 1 |
                              draws[, 30] < 1 |
                                draws[, 31] < 1, NA, draws[, 38])
    draws[, 39] <- ifelse(draws[, 24] < 1 |
                               draws[, 25] < 1 |
                                 draws[, 22] < 1 |
                                   draws[, 23] < 1 |
                          draws[, 28] < 1 |
                            draws[, 29] < 1 |
                              draws[, 30] < 1 |
                                draws[, 31] < 1, NA, draws[, 39])

    draws[, 41] <- exp(log(draws[, 20]) -
                               qnorm(draws[, 40]) *
                                         ((log(uci.obs.rr) - log(lci.obs.rr)) /
                                              (qnorm(.975) * 2)))
    draws[, 42] <- exp(log(draws[, 21]) -
                               qnorm(draws[, 40]) *
                                         ((log(uci.obs.rr) - log(lci.obs.rr)) /
                                              (qnorm(.975) * 2)))
    draws[, 43] <- exp(log(draws[, 38]) -
                               qnorm(draws[, 40]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))
    draws[, 44] <- exp(log(draws[, 39]) -
                               qnorm(draws[, 40]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))

    corr.rr.smr <- c(median(draws[, 20], na.rm = TRUE),
                     quantile(draws[, 20], probs = .025, na.rm = TRUE),
                     quantile(draws[, 20], probs = .975, na.rm = TRUE))
    corr.rr.mh <- c(median(draws[, 21], na.rm = TRUE),
                    quantile(draws[, 21], probs = .025, na.rm = TRUE),
                    quantile(draws[, 21], probs = .975, na.rm = TRUE))
    corr.or.smr <- c(median(draws[, 38], na.rm = TRUE),
                     quantile(draws[, 38], probs = .025, na.rm = TRUE),
                     quantile(draws[, 38], probs = .975, na.rm = TRUE))
    corr.or.mh <- c(median(draws[, 39], na.rm = TRUE),
                    quantile(draws[, 39], probs = .025, na.rm = TRUE),
                    quantile(draws[, 39], probs = .975, na.rm = TRUE))
    tot.rr.smr <- c(median(draws[, 40], na.rm = TRUE),
                    quantile(draws[, 40], probs = .025, na.rm = TRUE),
                    quantile(draws[, 40], probs = .975, na.rm = TRUE))
    tot.rr.mh <- c(median(draws[, 41], na.rm = TRUE),
                   quantile(draws[, 41], probs = .025, na.rm = TRUE),
                   quantile(draws[, 41], probs = .975, na.rm = TRUE))
    tot.or.smr <- c(median(draws[, 42], na.rm = TRUE),
                    quantile(draws[, 42], probs = .025, na.rm = TRUE),
                    quantile(draws[, 42], probs = .975, na.rm = TRUE))
    tot.or.mh <- c(median(draws[, 43], na.rm = TRUE),
                   quantile(draws[, 43], probs = .025, na.rm = TRUE),
                   quantile(draws[, 43], probs = .975, na.rm = TRUE))

    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    if (print)
        cat("Observed Data:",
            "\n--------------", 
            "nOutcome   :", rownames(tab)[1],
            "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
    if (print) 
        print(round(tab, dec))
    if (print) 
        cat("\n")
    rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr),
                  c(obs.or, lci.obs.or, uci.obs.or))
    rownames(rmat) <- c(" Observed Relative Risk:", "    Observed Odds Ratio:")
    colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                               sep = ""), "interval")
    if (print)
        cat("Observed Measures of Exposure-Outcome Relationship:",
            "\n-----------------------------------------------------------------------------------\n\n")
    if (print) 
        print(round(rmat, dec))
    if (print)
        cat("\n")
    rmatc <- rbind(corr.rr.smr, corr.rr.mh, tot.rr.smr, tot.rr.mh, corr.or.smr,
                   corr.or.mh, tot.or.smr, tot.or.mh)
    rownames(rmatc) <- c("RR adjusted using SMR estimate -- systematic error:", "RR adjusted using MH estimate -- systematic error", "RR adjusted using SMR estimate -- systematic error:", "RR adjusted using MH estimate -- systematic error", "RR adjusted using SMR estimate -- systematic and random error:", "RR adjusted using MH estimate -- systematic and random error", "OR adjusted using SMR estimate -- systematic error:", "OR adjusted using MH estimate -- systematic error", "OR adjusted using SMR estimate -- systematic error:", "OR adjusted using MH estimate -- systematic error", "OR adjusted using SMR estimate -- systematic and random error:", "OR adjusted using MH estimate -- systematic and random error")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    if (print)
        print(round(rmatc, dec))
    if (print)
        cat("\nBias Parameters:",
            "\n----------------\n\n")
    if (print)
        cat("p(Confounder+|Exposure+):", prev.exp[[1]], "(", prev.exp[[2]], ")",
            "\np(Confounder+|Exposure-):", prev.nexp[[1]], "(", prev.nexp[[2]], ")",
            "\nRisk(Confounder-Outcome):", risk[[1]], "(", risk[[2]], ")",
            "\n")
    invisible(list(obs.data = tab,
                   obs.measures = rmat, 
                   corr.rr.smr = corr.rr.smr, corr.rr.mh = corr.rr.mh,
                   corr.or.smr = corr.or.smr, corr.or.mh = corr.or.mh,
                   sim.df = as.data.frame(draws)))
}

