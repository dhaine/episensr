misclassification <- function(exposed,
                              case,
                              implement = c("exposure", "outcome"),
                              bias = NULL,
                              alpha = 0.05,
                              dec = 4,
                              print = TRUE){
    if(is.null(bias))
        bias <- c(1, 1, 1, 1)
    else bias <- bias
    if(length(bias) != 4)
        stop('The argument bias should be made of the following components: (1) Sensitivity of exposure classification among those with the outcome, (2) Sensitivity of exposure classification among those without the outcome, (3) Specificity of exposure classification among those with the outcome, and (4) Specificity of exposure classification among those without the outcome.')
    if(!all(bias >= 0 & bias <=1))
        stop('Bias parameters should be between 0 and 1.')

    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, cased)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

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

        A <- (a - (1 - bias[3]) * (a + b)) / (bias[1] - (1 - bias[3]))
        C <- (c - (1 - bias[4]) * (c + d)) / (bias[2] - (1 - bias[4]))
        B <- (a + b) - A
        D <- (c + d) - C

        if(A < 1 | B < 1 | C < 1 | D < 1)
            stop('Negative cell.')
        
        corr.tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr.rr <- (A/(A + C)) / (B/(B + D))
        corr.or <- (A/B) / (C/D)

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(corr.tab) <- paste("Row", 1:2)
        } else {
            rownames(corr.tab) <- row.names(tab)
        }
        if (is.null(colnames(tab))){
            colnames(corr.tab) <- paste("Col", 1:2)
        } else {
            colnames(corr.tab) <- colnames(tab)
        }
        if (print)
            cat("Observed Data:",
                "\n--------------", 
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
        if (print) 
            print(round(tab, dec))
        if (print)
            cat("\nCorrected Data:",
                "\n--------------------\n\n")
        if (print)
            print(round(corr.tab, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr), c(obs.or, lci.obs.or, uci.obs.or))
        rownames(rmat) <- c(" Observed Relative Risk:", "    Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                               sep = ""), "interval")
        rmatc <- rbind(corr.rr, corr.or)
        rownames(rmatc) <- c("Corrected Relative Risk:",
                             "   Corrected Odds Ratio:")
        if (print)
            cat("Observed Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Corrected Relative Risk:", round(corr.rr, dec), "\n   Corrected Odds Ratio:", round(corr.or, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("Se(Outcome+):", bias[1],
                "\nSe(Outcome-):", bias[2],
                "\nSp(Outcome+):", bias[3],
                "\nSp(Outcome-):", bias[4],
                "\n")
    }
    
    if (implement == "outcome"){
        obs.rr <- (a/(a + c)) / (b/(b + d))
        se.log.obs.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.obs.rr <- exp(log(obs.rr) - qnorm(1 - alpha/2) * se.log.obs.rr)
        uci.obs.rr <- exp(log(obs.rr) + qnorm(1 - alpha/2) * se.log.obs.rr)
        
        obs.or <- (a/b) / (c/d)
        se.log.obs.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.obs.or <- exp(log(obs.or) - qnorm(1 - alpha/2) * se.log.obs.or)
        uci.obs.or <- exp(log(obs.or) + qnorm(1 - alpha/2) * se.log.obs.or)

        A <- (a - (1 - bias[3]) * (a + c)) / (bias[1] - (1 - bias[3]))
        B <- (b - (1 - bias[4]) * (b + d)) / (bias[2] - (1 - bias[4]))
        C <- (a + c) - A
        D <- (b + d) - B

        if(A < 1 | B < 1 | C < 1 | D < 1)
            stop('Negative cell.')
        
        corr.tab <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)

        corr.rr <- (A/(A + C)) / (B/(B + D))
        corr.or <- (A/B) / (C/D)

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(corr.tab) <- paste("Row", 1:2)
        } else {
            rownames(corr.tab) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(corr.tab) <- paste("Col", 1:2)
        } else {
            colnames(corr.tab) <- colnames(tab)
        }
        if (print) 
            cat("Observed Data:",
                "\n--------------", 
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
        if (print) 
            print(round(tab, dec))
        if (print)
            cat("\nCorrected Data:",
                "\n--------------------\n\n")
        if (print)
            print(round(corr.tab, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(obs.rr, lci.obs.rr, uci.obs.rr), c(obs.or, lci.obs.or, uci.obs.or))
        rownames(rmat) <- c(" Observed Relative Risk:", "    Observed Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
        rmatc <- rbind(corr.rr, corr.or)
        rownames(rmatc) <- c("Corrected Relative Risk:",
                             "   Corrected Odds Ratio:")
        if (print)
            cat("Observed Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Corrected Relative Risk:", round(corr.rr, dec), "\n   Corrected Odds Ratio:", round(corr.or, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("Se(Exposure+):", bias[1],
                "\nSe(Exposure-):", bias[2],
                "\nSp(Exposure+):", bias[3],
                "\nSp(Exposure-):", bias[4],
                "\n")
    }
    invisible(list(obs.data = tab,
                   corr.data = corr.tab,
                   obs.measures = rmat, 
                   adj.measures = rmatc,
                   bias.params = bias))
}
