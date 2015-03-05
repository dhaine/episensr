probsens <- function(exposed, case,
                     implement = c("exposure", "outcome"),
                     reps = 1000,
                     seca.parms = list(dist = c("uniform", "triangular",
                                           "trapezoidal"),
                                       parms = NULL),
                     seexp.parms = NULL,
                     spca.parms = list(dist = c("uniform", "triangular",
                                           "trapezoidal"),
                                       parms = NULL),
                     spexp.parms = NULL,
                     corr.se = NULL,
                     corr.sp = NULL,
                     alpha = 0.05,
                     dec = 4,
                     print = TRUE){
    if(is.null(seca.parms) | is.null(spca.parms))
        stop('At least one Se and one Sp should be provided through outcome parameters.')
    if(!is.list(seca.parms))
        stop('Sensitivity of exposure classification among those with the outcome should be a list.')
    else seca.parms <- seca.parms
    if(!is.null(seexp.parms) & !is.list(seexp.parms))
        stop('Sensitivity of exposure classification among those without the outcome should be a list.')
    else seexp.parms <- seexp.parms
    if(!is.list(spca.parms))
        stop('Specificity of exposure classification among those with the outcome should be a list.')
    else spca.parms <- spca.parms
    if(!is.null(spexp.parms) & !is.list(spexp.parms))
        stop('Specificity of exposure classification among those without the outcome should be a list.')
    else spexp.parms <- spexp.parms
    if(!is.null(seexp.parms) & (is.null(spca.parms) | is.null(spexp.parms) |
                                is.null(corr.se) | is.null(corr.sp)))
        stop('For non-differential misclassification type, have to provide Se and Sp for among those with and without the outcome as well as Se and Sp correlations.')
    if(!is.null(corr.se) && (corr.se == 0 | corr.se == 1))
        stop('Correlations should be > 0 and < 1.')
    if(!is.null(corr.sp) && (corr.sp == 0 | corr.sp == 1))
        stop('Correlations should be > 0 and < 1.')
    
    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, cased)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    draws <- matrix(NA, nrow = reps, ncol = 10)

    seca <- c(reps, seca.parms[[2]])
    seexp <- c(reps, seexp.parms[[2]])
    spca <- c(reps, spca.parms[[2]])
    spexp <- c(reps, spexp.parms[[2]])

    implement <- match.arg(implement)
    if (implement == "exposure") {
        obs.rr <- (a/(a + c)) / (b/(b + d))
        se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
        uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)

        obs.or <- (a/b) / (c/d)
        se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
        uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

        if (is.null(seexp.parms) & !is.null(spca.parms) & is.null(spexp.parms) &
            is.null(corr.se) & is.null(corr.sp)) {
            if (seca.parms[[1]] == "uniform") {
                draws[, 1] <- do.call(runif, as.list(seca))
            }
            if (seca.parms[[1]] == "triangular") {
                draws[, 1] <- do.call(triangle::rtriangle, as.list(seca))
            }
            if (seca.parms[[1]] == "trapezoidal") {
                draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(seca))
            }
            draws[, 2] <- draws[, 1]
            if (spca.parms[[1]] == "uniform") {
                draws[, 3] <- do.call(runif, as.list(spca))
            }
            if (spca.parms[[1]] == "triangular") {
                draws[, 3] <- do.call(triangle::rtriangle, as.list(spca))
            }
            if (spca.parms[[1]] == "trapezoidal") {
                draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(spca))
            }
            draws[, 4] <- draws[, 3]
        } else {
            corr.draws <- matrix(NA, nrow = reps, ncol = 10)
            corr.draws[, 1:6] <- apply(corr.draws[, 1:6],
                                       2,
                                       function(x) x = runif(reps))
            corr.draws[, 1:6] <- apply(corr.draws[, 1:6],
                                       2,
                                       function(x) log(x / (1 - x)))
            corr.draws[, 7] <- exp((sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 2]) / (1 + (sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 2])))
            corr.draws[, 8] <- exp((sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 3]) / (1 + (sqrt(corr.se) * corr.draws[, 1] + sqrt(1 - corr.se) * corr.draws[, 3])))
            corr.draws[, 9] <- exp((sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 5]) / (1 + (sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 5])))
            corr.draws[, 10] <- exp((sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 6]) / (1 + (sqrt(corr.sp) * corr.draws[, 4] + sqrt(1 - corr.sp) * corr.draws[, 6])))

            if (seca.parms[[1]] == "uniform") {
                draws[, 1] <- seca.parms[[2]][2] -
                    (seca.parms[[2]][2] - seca.parms[[2]][1]) * corr.draws[, 7]
            }
            if (seca.parms[[1]] == "triangular" | seca.parms[[1]] == "trapezoidal") {
                draws[, 1] <- (corr.draws[, 1] *
                                   (seca.parms[[2]][4] + seca.parms[[2]][3] - seca.parms[[2]][1] - seca.parms[[2]][2]) + (seca.parms[[2]][1] + seca.parms[[2]][2])) / 2
                draws[, 1] <- ifelse(draws[, 1] < seca.parms[[2]][2],
                                     seca.parms[[2]][1] + sqrt((seca.parms[[2]][2] - seca.parms[[2]][1]) * (2 * draws[, 1] - seca.parms[[2]][1] - seca.parms[[2]][2])),
                                     draws[, 1])
                draws[, 1] <- ifelse(draws[, 1] > seca.parms[[2]][3],
                                     seca.parms[[2]][4] - sqrt(2 * (seca.parms[[2]][4] - seca.parms[[2]][3]) * (draws[, 1] - seca.parms[[2]][3])),
                                     draws[, 1])
            }
            if (seexp.parms[[1]] == "uniform") {
                draws[, 2] <- seexp.parms[[2]][2] -
                    (seexp.parms[[2]][2] - seexp.parms[[2]][1]) * corr.draws[, 7]
            }
            if (seexp.parms[[1]] == "triangular" | seexp.parms[[1]] == "trapezoidal") {
                draws[, 2] <- (corr.draws[, 1] *
                                   (seexp.parms[[2]][4] + seexp.parms[[2]][3] - seexp.parms[[2]][1] - seexp.parms[[2]][2]) + (seexp.parms[[2]][1] + seexp.parms[[2]][2])) / 2
                draws[, 2] <- ifelse(draws[, 2] < seexp.parms[[2]][2],
                                     seexp.parms[[2]][1] + sqrt((seexp.parms[[2]][2] - seexp.parms[[2]][1]) * (2 * draws[, 2] - seexp.parms[[2]][1] - seexp.parms[[2]][2])),
                                     draws[, 2])
                draws[, 2] <- ifelse(draws[, 2] > seexp.parms[[2]][3],
                                     seexp.parms[[2]][4] - sqrt(2 * (seexp.parms[[2]][4] - seexp.parms[[2]][3]) * (draws[, 2] - seexp.parms[[2]][3])),
                                     draws[, 2])
            }
            if (spca.parms[[1]] == "uniform") {
                draws[, 3] <- spca.parms[[2]][2] -
                    (spca.parms[[2]][2] - spca.parms[[2]][1]) * corr.draws[, 7]
            }
            if (spca.parms[[1]] == "triangular" | spca.parms[[1]] == "trapezoidal") {
                draws[, 3] <- (corr.draws[, 1] *
                                   (spca.parms[[2]][4] + spca.parms[[2]][3] - spca.parms[[2]][1] - spca.parms[[2]][2]) + (spca.parms[[2]][1] + spca.parms[[2]][2])) / 2
                draws[, 3] <- ifelse(draws[, 3] < spca.parms[[2]][2],
                                     spca.parms[[2]][1] + sqrt((spca.parms[[2]][2] - spca.parms[[2]][1]) * (2 * draws[, 3] - spca.parms[[2]][1] - spca.parms[[2]][2])),
                                     draws[, 3])
                draws[, 3] <- ifelse(draws[, 3] > spca.parms[[2]][3],
                                     spca.parms[[2]][4] - sqrt(2 * (spca.parms[[2]][4] - spca.parms[[2]][3]) * (draws[, 3] - spca.parms[[2]][3])),
                                     draws[, 3])
            }
            if (spexp.parms[[1]] == "uniform") {
                draws[, 4] <- spexp.parms[[2]][2] -
                    (spexp.parms[[2]][2] - spexp.parms[[2]][1]) * corr.draws[, 7]
            }
            if (spexp.parms[[1]] == "triangular" | spexp.parms[[1]] == "trapezoidal") {
                draws[, 4] <- (corr.draws[, 1] *
                                   (spexp.parms[[2]][4] + spexp.parms[[2]][3] - spexp.parms[[2]][1] - spexp.parms[[2]][2]) + (spexp.parms[[2]][1] + spexp.parms[[2]][2])) / 2
                draws[, 4] <- ifelse(draws[, 4] < spexp.parms[[2]][2],
                                     spexp.parms[[2]][1] + sqrt((spexp.parms[[2]][2] - spexp.parms[[2]][1]) * (2 * draws[, 4] - spexp.parms[[2]][1] - spexp.parms[[2]][2])),
                                     draws[, 4])
                draws[, 4] <- ifelse(draws[, 4] > spexp.parms[[2]][3],
                                     spexp.parms[[2]][4] - sqrt(2 * (spexp.parms[[2]][4] - spexp.parms[[2]][3]) * (draws[, 4] - spexp.parms[[2]][3])),
                                     draws[, 4])
            }
        }
            
        draws[, 5] <- (a - (1 - draws[, 3]) * (a + b)) /
            (draws[, 1] - (1 - draws[, 3]))
        draws[, 7] <- (c - (1 - draws[, 4]) * (c + d)) /
            (draws[, 2] - (1 - draws[, 4]))
        draws[, 6] <- (a + b) - draws[, 5]
        draws[, 8] <- (c + d) - draws[, 7]

        draws[, 9] <- (draws[, 5]/(draws[, 5] + draws[, 7])) /
            (draws[, 6]/(draws[, 6] + draws[, 8]))
        draws[, 10] <- (draws[, 5]/draws[, 6]) / (draws[, 7]/draws[, 8])

        draws[, 9] <- ifelse(draws[, 5] < 0 |
                               draws[, 6] < 0 |
                                 draws[, 7] < 0 |
                                   draws[, 8] < 0 |
                                     draws[, 9] < 0, NA, draws[, 9])
        draws[, 10] <- ifelse(draws[, 5] < 0 |
                               draws[, 6] < 0 |
                                 draws[, 7] < 0 |
                                   draws[, 8] < 0 |
                                       draws[, 10] < 0, NA, draws[, 10])

        corr.rr <- c(median(draws[, 9], na.rm = TRUE),
                     quantile(draws[, 9], probs = .025, na.rm = TRUE),
                     quantile(draws[, 9], probs = .975, na.rm = TRUE))
        corr.or <- c(median(draws[, 10], na.rm = TRUE),
                     quantile(draws[, 10], probs = .025, na.rm = TRUE),
                     quantile(draws[, 10], probs = .975, na.rm = TRUE))
        
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
        rmatc <- rbind(corr.rr, corr.or)
        rownames(rmatc) <- c(" Corrected Relative Risk:", "    Corrected Odds Ratio:")
        colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
        if (print)
            print(round(rmatc, dec))
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("Se|Cases:", seca.parms[[1]], "(", seca.parms[[2]], ")",
                "\nSp|Cases:", spca.parms[[1]], "(", spca.parms[[2]], ")",
                "\nSe|No-cases:", seexp.parms[[1]], "(", seexp.parms[[2]], ")",
                "\nSp|No-cases:", spexp.parms[[1]], "(", spexp.parms[[2]], ")",
                "\n")
        invisible(list(obs.data = tab,
                       obs.measures = rmat, 
                       corr.rr = corr.rr, corr.or = corr.or,
                       sim.df = as.data.frame(draws)))
        }
}
