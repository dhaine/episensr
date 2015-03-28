confounders.emm <- function(exposed,
                            case,
                            implement = c("RR", "OR", "RD"),
                            p = NULL,
                            RR.cd = NULL,
                            OR.cd = NULL,
                            RD.cd = NULL,
                            alpha = 0.05,
                            dec = 4,
                            print = TRUE){
    if(length(implement) > 1)
        stop('Choose between RR, OR, or RD implementation.')
    if(implement == "RR" & (!is.null(OR.cd) | !is.null(RD.cd)))
        stop('Mismatch between implementation type and confounder risk.')
    if(implement == "OR" & (!is.null(RR.cd) | !is.null(RD.cd)))
        stop('Mismatch between implementation type and confounder risk.')
    if(implement == "RD" & (!is.null(RR.cd) | !is.null(OR.cd)))
        stop('Mismatch between implementation type and confounder risk.')
    
    if(is.null(p))
        p <- c(0, 0)
    else p <- p
    if(length(p) != 2)
        stop('The argument p should be made of the following components: (1) Prevalence of the confounder among the exposed, and (2) Prevalence of the confounder among the unexposed.')
    if(!all(p >= 0 & p <=1))
        stop('Prevalences should be between 0 and 1.')

    if(is.null(RR.cd))
        RR.cd <- c(1, 1)
    else RR.cd <- RR.cd
    if(length(RR.cd) != 2)
        stop('Confounder-disease relative risk needs two arguments.')
    if(!all(RR.cd > 0))
        stop('Confounder-disease relative risk should be greater than 0.')

    if(is.null(OR.cd))
        OR.cd <- c(1, 1)
    else OR.cd <- OR.cd
    if(length(OR.cd) != 2)
        stop('Confounder-disease odds ratio needs two arguments.')
    if(!all(OR.cd > 0))
        stop('Confounder-disease odds ratio should be greater than 0.')

    if(is.null(RD.cd))
        RD.cd <- c(1, 1)
    else RD.cd <- RD.cd
    if(length(RD.cd) != 2)
        stop('Confounder-disease risk difference needs two arguments.')    
    
    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, case)
    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    implement <- match.arg(implement)
    if (implement == "RR") {
        crude.rr <- (a/(a + c)) / (b/(b + d))
        se.log.crude.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.crude.rr <- exp(log(crude.rr) - qnorm(1 - alpha/2) * se.log.crude.rr)
        uci.crude.rr <- exp(log(crude.rr) + qnorm(1 - alpha/2) * se.log.crude.rr)

        M1 <- (a + c) * p[1]
        N1 <- (b + d) * p[2]
        A1 <- (RR.cd[1] * M1 * a) / (RR.cd[1] * M1 + (a + c) - M1)
        B1 <- (RR.cd[2] * N1 * b) / (RR.cd[2] * N1 + (b + d) - N1)
        C1 <- M1 - A1
        D1 <- N1 - B1
        M0 <- a + c - M1
        N0 <- b + d - N1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1

        if(A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Negative cell.')
        
        tab.cfder <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRrr <- a / ((M1 * B1/N1) + (M0 * B0/N0))
        MHrr <- (A1 * N1/(M1 + N1) + A0 * N0/(M0 + N0)) /
            (B1 * M1/(M1 + N1) + B0 * M0/(M0 + N0)) 
        cfder.rr <- (A1/(A1 + C1)) / (B1/(B1 + D1))
        nocfder.rr <- (A0/(A0 + C0)) / (B0/(B0 + D0))
        RRadj.smr <- crude.rr / SMRrr
        RRadj.mh <- crude.rr / MHrr

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
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
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
        rmatc <- rbind(c(SMRrr, RRadj.smr), c(MHrr, RRadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio", "Mantel-Haenszel")
        colnames(rmatc) <- c("SMR_RR/MH_RR", "RRc")
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
            cat("Standardized Morbidity Ratio", "    SMRrr:", round(SMRrr, dec), "   RR adjusted using SMR estimate:", round(RRadj.smr, dec),
                "\nMantel-Haenszel", "                  MHrr:", round(MHrr, dec), "    RR adjusted using MH estimate:", round(RRadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("           p(Confounder+|Exposure+):", p[1],
                "\n           p(Confounder+|Exposure-):", p[2],
                "\n  RR(Confounder-Outcome in Exposed):", RR.cd[1],
                "\nRR(Confounder-Outcome in Unexposed):", RR.cd[2],
                "\n")
        rmat <- rbind(rmat, c(cfder.rr, NA, NA), c(nocfder.rr, NA, NA))
        rownames(rmat) <- c("        Crude Relative Risk:",
                            "Relative Risk, Confounder +:",
                            "Relative Risk, Confounder -:")
    }
    
    if (implement == "OR"){
        crude.or <- (a/b) / (c/d)
        se.log.crude.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.crude.or <- exp(log(crude.or) - qnorm(1 - alpha/2) * se.log.crude.or)
        uci.crude.or <- exp(log(crude.or) + qnorm(1 - alpha/2) * se.log.crude.or)

        C1 <- c * p[1] 
        D1 <- d * p[2]
        A1 <- (OR.cd[1] * C1 * a) / (OR.cd[1] * C1 + c - C1)
        B1 <- (OR.cd[2] * D1 * b) / (OR.cd[2] * D1 + d - D1)
        M1 <- A1 + C1
        N1 <- B1 + D1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1
        M0 <- A0 + C0
        N0 <- B0 + C0

        if(A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Negative cell.')
        
        tab.cfder <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRor <- a / ((C1 * B1/D1) + (C0 * B0/D0))
        MHor <- (A1 * D1/(M1 + N1) + A0 * D0/(M0 + N0)) /
            (B1 * C1/(M1 + N1) + B0 * C0/(M0 + N0)) 
        cfder.or <- (A1 / C1) / (B1 / D1)
        nocfder.or <- (A0 / C0) / (B0 / D0)
        ORadj.smr <- crude.or / SMRor
        ORadj.mh <- crude.or / MHor

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
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
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
        rmat <- rbind(c(crude.or, lci.crude.or, uci.crude.or))
        rownames(rmat) <- c("        Crude Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
        rmatc <- rbind(c(SMRor, ORadj.smr), c(MHor, ORadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio", "Mantel-Haenszel")
        colnames(rmatc) <- c("SMR_OR/MH_OR", "ORc")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Odds Ratio, Confounder +:", round(cfder.or, dec), "\nOdds Ratio, Confounder -:", round(nocfder.or, dec), "\n")
        if (print)
            cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
                "\n------------------------------------------------------\n\n")
        if (print)
            cat("Standardized Morbidity Ratio", "    SMRor:", round(SMRor, dec), "   OR adjusted using SMR estimate:", round(ORadj.smr, dec),
                "\nMantel-Haenszel", "                  MHor:", round(MHor, dec), "    OR adjusted using MH estimate:", round(ORadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("           p(Confounder+|Exposure+):", p[1],
                "\n           p(Confounder+|Exposure-):", p[2],
                "\n  OR(Confounder-Outcome in Exposed):", OR.cd[1],
                "\nOR(Confounder-Outcome in Unexposed):", OR.cd[2],
                "\n")
        rmat <- rbind(rmat, c(cfder.or, NA, NA), c(nocfder.or, NA, NA))
        rownames(rmat) <- c("        Crude Odds Ratio:",
                            "Odds Ratio, Confounder +:",
                            "Odds Ratio, Confounder -:")
    }
    
    if (implement == "RD"){
        crude.rd <- (a / (a + c)) - (b / (b + d))
        se.log.crude.rd <- sqrt((a * c) / (a + c)^3 + (b * d) / (b + d)^3)
        lci.crude.rd <- crude.rd - qnorm(1 - alpha/2) * se.log.crude.rd
        uci.crude.rd <- crude.rd + qnorm(1 - alpha/2) * se.log.crude.rd

        M1 <- (a + c) * p[1]
        N1 <- (b + d) * p[2]
        M0 <- (a + c) - M1
        N0 <- (b + d) - N1
        A1 <- (RD.cd[1] * M1 * M0 + M1 * a) / (a + c)
        B1 <- (RD.cd[2] * N1 * N0 + N1 * b) / (b + d)
        C1 <- M1 - A1
        D1 <- N1 - B1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1

        if(A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Negative cell.')
        
        tab.cfder <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        MHrd <- (((A1 * N1 - B1 * M1) / (M1 + N1)) +
                     ((A0 * N0 - B0 * M0) / (M0 + N0))) /
                         ((M1 * N1 / (M1 + N1)) +
                              (M0 * N0 / (M0 + N0)))
        cfder.rd <- (A1 / M1) - (B1 / N1)
        nocfder.rd <- (A0 / M0) - (B0 / N0)
        RDadj.mh <- crude.rd - MHrd

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
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
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
        rmat <- rbind(c(crude.rd, lci.crude.rd, uci.crude.rd))
        rownames(rmat) <- c("        Crude Risk Difference:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
        rmatc <- rbind(c(MHrd, RDadj.mh))
        rownames(rmatc) <- "Mantel-Haenszel"
        colnames(rmatc) <- c("MH_RD", "RDc")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Risk Difference, Confounder +:", round(cfder.rd, dec), "\nRisk Difference, Confounder -:", round(nocfder.rd, dec), "\n")
        if (print)
            cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
                "\n------------------------------------------------------\n\n")
        if (print)
            cat("\nMantel-Haenszel", "                  MHrd:", round(MHrd, dec), "   RD adjusted using MH estimate:", round(RDadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("           p(Confounder+|Exposure+):", p[1],
                "\n           p(Confounder+|Exposure-):", p[2],
                "\n  RD(Confounder-Outcome in Exposed):", RD.cd[1],
                "\nRD(Confounder-Outcome in Unexposed):", RD.cd[2],
                "\n")
    rmat <- rbind(rmat, c(cfder.rd, NA, NA), c(nocfder.rd, NA, NA))
    rownames(rmat) <- c("        Crude Risk Difference:",
                        "Risk Difference, Confounder +:",
                        "Risk Difference, Confounder -:")
    }
    invisible(list(obs.data = tab,
                   cfder.data = tab.cfder, nocfder.data = tab.nocfder,
                   obs.measures = rmat,
                   adj.measures = rmatc,
                   bias.parms = c(p, RR.cd)))
}
