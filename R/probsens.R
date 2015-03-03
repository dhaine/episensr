probsens <- function(exposed, case,
                     implement = c("exposure", "outcome"),
                     reps = 1000,
                     seca.parms = list(dist = c("uniform", "triangular",
                                           "trapezoidal", "normal",
                                           "beta", "lognormal", "logistic",
                                           "logitnormal", "loglogistic"),
                                       parms = NULL),
                     seexp.parms = NULL,
                     spca.parms = list(dist = c("uniform", "triangular",
                                           "trapezoidal", "normal",
                                           "beta", "lognormal", "logistic",
                                           "logitnormal", "loglogistic"),
                                       parms = NULL),
                     spexp.parms = NULL,
                     type = c("nondiff", "diff"),
                     alpha = 0.05,
                     dec = 4,
                     print = TRUE){
#    if(is.null(bias))
#        bias <- c(1, 1, 1, 1)
#    else bias <- bias
#    if(length(bias) != 4)
#        stop('The argument bias should be made of the following components: (1) Sensitivity of exposure classification among those with the outcome, (2) Sensitivity of exposure classification among those without the outcome, (3) Specificity of exposure classification among those with the outcome, and (4) Specificity of exposure classification among those without the outcome.')
#    if(!all(bias >= 0 & bias <=1))
#        stop('Bias parameters should be between 0 and 1.')

#    if(!is.list(seca.parms))
#        stop('Sensitivity of exposure classification among those with the outcome should be a list')
#    else seca.parms <- seca.parms
#    if(!is.list(seexp.parms))
#        stop('Sensitivity of exposure classification among those without the outcome should be a list')
#    else seexp.parms <- seexp.parms
#    if(!is.list(spca.parms))
#        stop('Specificity of exposure classification among those with the outcome should be a list')
#    else spca.parms <- spca.parms
#    if(!is.list(spexp.parms))
#        stop('Specificity of exposure classification among those without the outcome should be a list')
#    else spexp.parms <- spexp.parms
    
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

        if (seca.parms[[1]] == "uniform") {
            draws[, 1] <- do.call(runif, as.list(seca))
            }
        if (seca.parms[[1]] == "triangular") {
            draws[, 1] <- do.call(triangle::rtriangle, as.list(seca))
            }
        if (seca.parms[[1]] == "trapezoidal") {
            draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(seca))
            }
        if (seca.parms[[1]] == "normal") {
            draws[, 1] <- do.call(rnorm, as.list(seca))
            }
        if (seca.parms[[1]] == "beta") {
            draws[, 1] <- do.call(rbeta, as.list(seca))
            }
        if (seca.parms[[1]] == "lognormal") {
            draws[, 1] <- do.call(rlnorm, as.list(seca))
            }
        if (seca.parms[[1]] == "logistic") {
            draws[, 1] <- do.call(rlogis, as.list(seca))
            }
        if (seca.parms[[1]] == "logitnormal") {
            draws[, 1] <- do.call(logitnorm::rlogitnorm, as.list(seca))
            }
        if (seca.parms[[1]] == "loglogistic") {
            draws[, 1] <- do.call(actuar::rllogis, as.list(seca))
            }

        if(is.null(seexp.parms)) {
            draws[, 2] <- draws[, 1]
        } else {
        if (seexp.parms[[1]] == "uniform") {
            draws[, 2] <- do.call(runif, as.list(seexp))
            }
        if (seexp.parms[[1]] == "triangular") {
            draws[, 2] <- do.call(triangle::rtriangle, as.list(seexp))
            }
        if (seexp.parms[[1]] == "trapezoidal") {
            draws[, 2] <- do.call(trapezoid::rtrapezoid, as.list(seexp))
            }
        if (seexp.parms[[1]] == "normal") {
            draws[, 2] <- do.call(rnorm, as.list(seexp))
            }
        if (seexp.parms[[1]] == "beta") {
            draws[, 2] <- do.call(rbeta, as.list(seexp))
            }
        if (seexp.parms[[1]] == "lognormal") {
            draws[, 2] <- do.call(rlnorm, as.list(seexp))
            }
        if (seexp.parms[[1]] == "logistic") {
            draws[, 2] <- do.call(rlogis, as.list(seexp))
            }
        if (seexp.parms[[1]] == "logitnormal") {
            draws[, 2] <- do.call(logitnorm::rlogitnorm, as.list(seexp))
            }
        if (seexp.parms[[1]] == "loglogistic") {
            draws[, 2] <- do.call(actuar::rllogis, as.list(seexp))
            }
        }

        if (spca.parms[[1]] == "uniform") {
            draws[, 3] <- do.call(runif, as.list(spca))
            }
        if (spca.parms[[1]] == "triangular") {
            draws[, 3] <- do.call(triangle::rtriangle, as.list(spca))
            }
        if (spca.parms[[1]] == "trapezoidal") {
            draws[, 3] <- do.call(trapezoid::rtrapezoid, as.list(spca))
            }
        if (spca.parms[[1]] == "normal") {
            draws[, 3] <- do.call(rnorm, as.list(spca))
            }
        if (spca.parms[[1]] == "beta") {
            draws[, 3] <- do.call(rbeta, as.list(spca))
            }
        if (spca.parms[[1]] == "lognormal") {
            draws[, 3] <- do.call(rlnorm, as.list(spca))
            }
        if (spca.parms[[1]] == "logistic") {
            draws[, 3] <- do.call(rlogis, as.list(spca))
            }
        if (spca.parms[[1]] == "logitnormal") {
            draws[, 3] <- do.call(logitnorm::rlogitnorm, as.list(spca))
            }
        if (spca.parms[[1]] == "loglogistic") {
            draws[, 3] <- do.call(actuar::rllogis, as.list(spca))
            }

        if (is.null(spexp.parms)) {
            draws[, 4] <- draws[, 3]
        } else {
        if (spexp.parms[[1]] == "uniform") {
            draws[, 4] <- do.call(runif, as.list(spexp))
            }
        if (spexp.parms[[1]] == "triangular") {
            draws[, 4] <- do.call(triangle::rtriangle, as.list(spexp))
            }
        if (spexp.parms[[1]] == "trapezoidal") {
            draws[, 4] <- do.call(trapezoid::rtrapezoid, as.list(spexp))
            }
        if (spexp.parms[[1]] == "normal") {
            draws[, 4] <- do.call(rnorm, as.list(spexp))
            }
        if (spexp.parms[[1]] == "beta") {
            draws[, 4] <- do.call(rbeta, as.list(spexp))
            }
        if (spexp.parms[[1]] == "lognormal") {
            draws[, 4] <- do.call(rlnorm, as.list(spexp))
            }
        if (spexp.parms[[1]] == "logistic") {
            draws[, 4] <- do.call(rlogis, as.list(spexp))
            }
        if (spexp.parms[[1]] == "logitnormal") {
            draws[, 4] <- do.call(logitnorm::rlogitnorm, as.list(spexp))
            }
        if (spexp.parms[[1]] == "loglogistic") {
            draws[, 4] <- do.call(actuar::rllogis, as.list(spexp))
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
                "\nOutcome   :", rownames(tab)[1],
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
