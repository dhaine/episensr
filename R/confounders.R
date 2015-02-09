confounders <- function(exposed, case, prev.cfder = NULL, cfder.dis.RR = NULL, alpha = 0.05, dec = 4, print = TRUE){
    if(is.null(prev.cfder))
        prev.cfder <- c(0, 0)
    else prev.cfder <- prev.cfder
    if(!is.vector(prev.cfder))
        stop('The argument prev.cfder should be a vector of length 2.')
    if(length(prev.cfder) != 2)
        stop('The argument prev.cfder should be made of 2 components: (1) Prevalence of the confounder among the exposed, and (2) Prevalence of the confounder among the unexposed.')
    if(is.null(cfder.dis.RR))
        cfder.dis.RR <- 1
    else cfder.dis.RR <- cfder.dis.RR
    if(cfder.dis.RR <= 0)
        stop('Confounder-disease relative risk should be greater than 0.')
    if(!all(prev.cfder >= 0 & prev.cfder <=1))
        stop('Prevalences should be between 0 and 1')
    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, cased)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]
    crude.rr <- (a/(a + c)) / (b/(b + d))
    se.log.crude.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
    lci.crude.rr <- exp(log(crude.rr) - qnorm(1 - alpha/2) * se.log.crude.rr)
    uci.crude.rr <- exp(log(crude.rr) + qnorm(1 - alpha/2) * se.log.crude.rr)

    M1 <- (a + c) * prev.cfder[1]
    N1 <- (b + d) * prev.cfder[2]
    A1 <- (cfder.dis.RR * M1 * a) / (cfder.dis.RR * M1 + (a + c) - M1)
    B1 <- (cfder.dis.RR * N1 * b) / (cfder.dis.RR * N1 + (b + d) - N1)
    C1 <- M1 - A1
    D1 <- N1 - B1
    M0 <- a + c - M1
    N0 <- b + d - N1
    A0 <- a - A1
    B0 <- b - B1
    C0 <- c - C1
    D0 <- d - D1
    tab.cfder <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
    tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

    SMR <- a / ((M1 * B1/N1) + (M0 * B0/N0))
    MH <- (A1 * N1/(M1 + N1) + A0 * N0/(M0 + N0)) / (B1 * M1/(M1 + N1) + B0 * M0/(M0 + N0))
    cfder.rr <- (A1/(A1 + C1)) / (B1/(B1 + D1))
    nocfder.rr <- (A0/(A0 + C0)) / (B0/(B0 + D0))
    RR0 <- crude.rr / SMR
    RRc <- crude.rr / MH

    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    if (is.null(rownames(tab))){
        rownames(tab.cfder) <- paste("Row", 1:2)
        } else {
        rownames(tab.cfder) <- row.names(tab)
    }
    if (is.null(colnames(tab))){ 
        colnames(tab.cfder) <- paste("Col", 1:2)
        } else {
        colnames(tab.cfder) <- colnames(tab)
    }
    if (is.null(rownames(tab))){
        rownames(tab.nocfder) <- paste("Row", 1:2)
        } else {
        rownames(tab.nocfder) <- row.names(tab)
    }
    if (is.null(colnames(tab))){ 
        colnames(tab.nocfder) <- paste("Col", 1:2)
        } else {
        colnames(tab.nocfder) <- colnames(tab)
    }
    if (print) 
        cat("Observed Data:",
            "\n--------------", 
            "\nOutcome   :", colnames(tab)[1],
            "\nComparing :", rownames(tab)[1], "vs.", rownames(tab)[2], "\n\n")
    if (print) 
        print(round(tab, dec))
    if (print)
        cat("\nData, Counfounder +:",
            "\n--------------------\n\n")
    if (print)
        print(round(tab.cfder, dec))
    if (print)
        cat("\nData, Counfounder -:",
            "\n--------------------\n\n")
    if (print)
        print(round(tab.nocfder, dec))
    if (print) 
        cat("\n")
    rmat <- rbind(c(crude.rr, lci.crude.rr, uci.crude.rr))
    rownames(rmat) <- c("        Crude Relative Risk:")
    colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
        sep = ""), "interval")
    if (print)
        cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
            "\n-----------------------------------------------------------------------------------\n\n")
    if (print) 
        print(round(rmat, dec))
    if (print)
        cat("Relative Risk, Confounder +:", round(cfder.rr, dec), "\nRelative Risk, Confounder -:", round(nocfder.rr, dec), "\n")
    if (print)
        cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
            "\n------------------------------------------------------\n\n")
    if (print)
        cat("Standardized Morbidity Ratio", "    SMRrr:", round(SMR, dec), "   RR0:", round(RR0, dec),
            "\nMantel-Haenszel", "                 MHrr:", round(MH, dec), "   RRc:", round(RRc, dec), "\n")
    if (print)
        cat("\nBias Parameters:",
            "\n----------------\n\n")
    if (print)
        cat("p(Confounder+|Exposure+):", prev.cfder[1],
            "\np(Confounder+|Exposure-):", prev.cfder[2],
            "\n  RR(Confounder-Outcome):", cfder.dis.RR, "\n")
    invisible(list(obs.data = tab, cfder.data = tab.cfder, nocfder.data = tab.nocfder,
                   obs.measures = rmat, SMR = SMR, MH = MH, cfder.rr = cfder.rr,
                   nocfder.rr = nocfder.rr, RR0 = RR0, RRc = RRc,
                   bias.params = c(prev.cfder, cfder.dis.RR)))
}
