probsens <- function(exposed,
                     case,
                     reps = 1000,
                     or.parms = list(dist = c("uniform", "triangular",
                                           "trapezoidal"),
                                       parms = NULL),
                     alpha = 0.05,
                     dec = 4,
                     print = TRUE){
    if(reps < 1)
        stop(paste("Invalid argument: reps =", reps))
    
    if(is.null(or.parms))
        stop('Please provide odds ratio for the probability of being selected.')
    if(!is.list(or.parms))
        stop('Odds ratio for the probability of being selected should be a list.')
    else or.parms <- or.parms
    if(or.parms[[1]] == "uniform" & length(or.parms[[2]]) != 2)
        stop('For uniform distribution, please provide vector of lower and upper limits.')
    if(or.parms[[1]] == "uniform" & or.parms[[2]][1] >= or.parms[[2]][2])
        stop('Lower limit of your uniform distribution is greater than upper limit.')
    if(or.parms[[1]] == "triangular" & length(or.parms[[2]]) != 3)
        stop('For triangular distribution, please provide vector of lower, upper limits, and mode.')
    if(or.parms[[1]] == "triangular" & ((or.parms[[2]][1] > or.parms[[2]][3]) |
                                        (or.parms[[2]][2] < or.parms[[2]][3])))
        stop('Wrong arguments for your triangular distribution.')
    if(or.parms[[1]] == "trapezoidal" & length(or.parms[[2]]) != 4)
        stop('For trapezoidal distribution, please provide vector of lower limit, lower mode, upper mode, and upper limit.')
    if(or.parms[[1]] == "trapezoidal" & ((or.parms[[2]][1] > or.parms[[2]][2]) |
                                         (or.parms[[2]][2] > or.parms[[2]][3]) |
                                         (or.parms[[2]][3] > or.parms[[2]][4])))
        stop('Wrong arguments for your trapezoidal distribution.')    
    
    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, cased)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    draws <- matrix(NA, nrow = reps, ncol = 4)
    colnames(draws) <- c("or.sel", "corr.or", "reps", "tot.or")

    or.sel <- c(reps, or.parms[[2]])

    obs.or <- (a/b) / (c/d)
    se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
    uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

    if (or.parms[[1]] == "uniform") {
        draws[, 1] <- do.call(runif, as.list(or.sel))
    }
    if (or.parms[[1]] == "triangular") {
        draws[, 1] <- do.call(triangle::rtriangle, as.list(or.sel))
    }
    if (or.parms[[1]] == "trapezoidal") {
        draws[, 1] <- do.call(trapezoid::rtrapezoid, as.list(or.sel))
    }

    draws[, 3] <- runif(reps)

    draws[, 2] <- obs.or / draws[, 1]

    draws[, 4] <- exp(log(draws[, 2]) -
                               qnorm(draws[, 3]) *
                                         ((log(uci.obs.or) - log(lci.obs.or)) /
                                              (qnorm(.975) * 2)))

    corr.or <- c(median(draws[, 2], na.rm = TRUE),
                 quantile(draws[, 2], probs = .025, na.rm = TRUE),
                 quantile(draws[, 2], probs = .975, na.rm = TRUE))
    tot.or <- c(median(draws[, 4], na.rm = TRUE),
                quantile(draws[, 4], probs = .025, na.rm = TRUE),
                quantile(draws[, 4], probs = .975, na.rm = TRUE))        

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
    rmat <- data.frame(obs.or, lci.obs.or, uci.obs.or)
    rownames(rmat) <- "Observed Odds Ratio:"
    colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                       sep = ""), "interval")
    if (print)
        cat("Observed Measures of Exposure-Outcome Relationship:",
            "\n-----------------------------------------------------------------------------------\n\n")
    if (print) 
        print(round(rmat, dec))
    if (print)
        cat("\n")
    rmatc <- rbind(corr.or, tot.or)
    rownames(rmatc) <- c("Odds Ratio -- systematic error:", "Odds Ratio -- systematic and random error")
    colnames(rmatc) <- c("Median", "2.5th percentile", "97.5th percentile")
    if (print)
        print(round(rmatc, dec))
    if (print)
        cat("\nBias Parameters:",
            "\n----------------\n\n")
    if (print)
        cat("OR selection:", or.parms[[1]], "(", or.parms[[2]], ")",
            "\n")
    invisible(list(obs.data = tab,
                   obs.measures = rmat, 
                   corr.or = corr.or, tot.or = tot.or,
                   sim.df = as.data.frame(draws[, -3])))
}
