multidimBias <- function(exposed, case,
                         se = NULL, sp = NULL,
                         alpha = 0.05, dec = 4, print = TRUE) {
    if(is.null(se))
        stop('Please provide argument for sensitivity.')
    if(is.null(sp))
        stop('Please provide argument for specificity.')
    if(!is.vector(se))
        stop('Sensitivity should be a vector.')
    if(!is.vector(sp))
        stop('Specificity should be a vector.')
    if(!all(se >= 0 & se <=1))
        stop('Sensitivity should be between 0 and 1.')
    if(!all(sp >= 0 & sp <=1))
        stop('Specificity should be between 0 and 1.')
    if(length(se) != length(sp))
        stop('Sensitivity and specificity should be of the same length.')
    
    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, cased)

    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    rr <- (a/(a + c)) / (b/(b + d))
    se.log.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
    lci.rr <- exp(log(rr) - qnorm(1 - alpha/2) * se.log.rr)
    uci.rr <- exp(log(rr) + qnorm(1 - alpha/2) * se.log.rr)

    or <- (a/b) / (c/d)
    se.log.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.or <- exp(log(or) - qnorm(1 - alpha/2) * se.log.or)
    uci.or <- exp(log(or) + qnorm(1 - alpha/2) * se.log.or)

    rr.mat <- matrix(NA, nrow = length(se), ncol = length(se))
    or.mat <- matrix(NA, nrow = length(se), ncol = length(se))
    
    for (i in 1:nrow(rr.mat)) {
        for (j in 1:nrow(rr.mat)) {
        rr.mat[i, j] <- (((a - (1 - sp[j]) * (a + b)) / (se[j] - (1 - sp[j]))) /
                          (((a - (1 - sp[j]) * (a + b)) / (se[j] - (1 - sp[j]))) +
                     ((c - (1 - sp[i]) * (c + d))) / (se[i] - (1 - sp[i])))) /
                         (((a + b) - ((a - (1 - sp[j]) * (a + b)) /
                         (se[j] - (1 - sp[j])))) /
                             (((a + b) - ((a - (1 - sp[j]) * (a + b)) /
                                              (se[j] - (1 - sp[j])))) +
                             ((c + d) - ((c - (1 - sp[i]) * (c + d)) /
                                             (se[i] - (1 - sp[i]))))))
        }
    }

    for (i in 1:nrow(or.mat)) {
        for (j in 1:nrow(or.mat)) {
        or.mat[i, j] <- (((a - (1 - sp[j]) * (a + b)) / (se[j] - (1 - sp[j]))) /
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
    rownames(rr.mat) <- paste("Se:", se, "Sp:", sp)
    colnames(rr.mat) <- paste("Se:", se, "Sp:", sp)
    rownames(or.mat) <- paste("Se:", se, "Sp:", sp)
    colnames(or.mat) <- paste("Se:", se, "Sp:", sp)
    if (print) 
        cat("Observed Data:", "\n---------------------------------------------------", 
            "\nOutcome   :", rownames(tab)[1],
            "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
    if (print) 
        print(round(tab, dec))
    if (print) 
        cat("\n")
    rmat <- rbind(c(rr, lci.rr, uci.rr), c(or, lci.or, uci.or))
    rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
    colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
        sep = ""), "interval")
    if (print) 
        print(round(rmat, dec))
    if (print)
        cat("\nMultidimensional Corrected Relative Risk Data:",
            "\n----------------------------------------------",
            "\n           Outcome + -->",
            "\nOutcome - |",
            "\n          V\n")
    if (print)
        print(rr.mat)
    if (print)
        cat("\nMultidimensional Corrected Odds Ratio Data:",
            "\n-------------------------------------------",
            "\n          Cases -->",
            "\nControls |",
            "\n         V\n")
    if (print)
        print(or.mat)
    if (print)
        cat("\n")
    sesp <- rbind(se, sp)
    rownames(sesp) <- c("Sensitivities:",
                        "Specificities:")
    if (print)
        print(sesp)
    invisible(list(obs.data = tab, 
                   obs.measures = rmat,
                   rr.mat = rr.mat,
                   or.mat = or.mat,
                   sesp = sesp))
}
