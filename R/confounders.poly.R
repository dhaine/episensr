confounders.poly <- function(exposed, case, implement = c("RR", "OR", "RD"), prev.cfder = NULL, cfder.dis.RR = NULL, cfder.dis.OR = NULL, cfder.dis.RD = NULL, alpha = 0.05, dec = 4, print = TRUE){
    if(is.null(prev.cfder))
        prev.cfder <- c(0, 0)
    else prev.cfder <- prev.cfder
    if(length(prev.cfder) != 4)
        stop('The argument prev.cfder should be made of the following components: (1) Prevalence of the confounder (highest level) among the exposed, (2) Prevalence of the confounder (highest level) among the unexposed, (3) Prevalence of the confounder (mid-level) among the exposed, and (4) Prevalence of the confounder (mid-level) among the unexposed.')
    if(!all(prev.cfder >= 0 & prev.cfder <=1))
        stop('Prevalences should be between 0 and 1.')
    if(!all(prev.cfder >= 0 & prev.cfder <=1))
        stop('Prevalences should be between 0 and 1')
    
    if(is.null(cfder.dis.RR))
        cfder.dis.RR <- c(1, 1)
    else cfder.dis.RR <- cfder.dis.RR
    if(length(cfder.dis.RR) > 2)
        stop('Confounder-disease relative risk: more than 2 arguments.')
    if(!all(cfder.dis.RR > 0))
        stop('Confounder-disease relative risks should be greater than 0.')
    if(is.null(cfder.dis.OR))
        cfder.dis.OR <- c(1, 1)
    else cfder.dis.OR <- cfder.dis.OR
    if(length(cfder.dis.OR) > 2)
        stop('Confounder-disease odds ratio: more than 2 arguments.')
    if(!all(cfder.dis.OR > 0))
        stop('Confounder-disease odds ratios should be greater than 0.')
    if(is.null(cfder.dis.RD))
        cfder.dis.RD <- c(1, 1)
    else cfder.dis.RD <- cfder.dis.RD
    if(length(cfder.dis.RD) > 2)
        stop('Confounder-disease risk difference: more than 2 arguments.')    

    if(inherits(exposed, c("table", "matrix")))
        tab <- exposed
    else tab <- table(exposed, cased)
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

        M2 <- (a + c) * prev.cfder[1]
        M1 <- (a + c) * prev.cfder[3]
        N2 <- (b + d) * prev.cfder[2]
        N1 <- (b + d) * prev.cfder[4]
        M0 <- a + c - M2 - M1
        N0 <- b + d - N2 - N1
        A0 <- (M0 * a) / (cfder.dis.RR[2] * M1 + M0 + cfder.dis.RR[1] * M2)
        B0 <- (N0 * b) / (cfder.dis.RR[2] * N1 + N0 + cfder.dis.RR[1] * N2)
        A1 <- cfder.dis.RR[2] * A0 * M1 / M0
        B1 <- cfder.dis.RR[2] * B0 * N1 / N0
        A2 <- cfder.dis.RR[1] * A0 * M2 / M0
        B2 <- cfder.dis.RR[1] * B0 * N2 / N0
        C2 <- M2 - A2
        D2 <- N2 - B2
        C1 <- M1 - A1
        D1 <- N1 - B1
        C0 <- M0 - A0
        D0 <- N0 - B0
        tab.cfder2 <- matrix(c(A2, B2, C2, D2), nrow = 2, byrow = TRUE)
        tab.cfder1 <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRrr <- a / ((M2 * B2/N2) + (M1 * B1/N1) + (M0 * B0/N0))
        MHrr <- (A2 * N2 / (M2 + N2) + A1 * N1 / (M1 + N1) + A0 * N0 / (M0 + N0)) /
            (B2 * M2 / (M2 + N2) + B1 * M1 / (M1 + N1) + B0 * M0 / (M0 + N0))
        cfder2.rr <- (A2/(A2 + C2)) / (B2/(B2 + D2))
        cfder1.rr <- (A1/(A1 + C1)) / (B1/(B1 + D1))
        nocfder.rr <- (A0/(A0 + C0)) / (B0/(B0 + D0))
        RRadj.smr <- crude.rr / SMRrr
        RRadj.mh <- crude.rr / MHrr

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(tab.cfder2) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder2) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder2) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder2) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.cfder1) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder1) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder1) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder1) <- colnames(tab)
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
            cat("\nData, Counfounder +, Highest Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder2, dec))
        if (print)
            cat("\nData, Counfounder +, Mid-Level Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder1, dec))
        if (print)
            cat("\nData, Counfounder -:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.nocfder, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(crude.rr, lci.crude.rr, uci.crude.rr))
        rownames(rmat) <- c("                       Crude Relative Risk:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                           sep = ""), "interval")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Relative Risk, Confounder +, Highest Level:", round(cfder2.rr, dec), "\n    Relative Risk, Confounder +, Mid-Level:", round(cfder1.rr, dec), "\n               Relative Risk, Confounder -:", round(nocfder.rr, dec), "\n")
        if (print)
            cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
                "\n------------------------------------------------------\n\n")
        if (print)
            cat("Standardized Morbidity Ratio", "    SMRrr:", round(SMRrr, dec), "   RR adjusted using SMR estimate:", round(RRadj.smr, dec),
                "\nMantel-Haenszel", "                 MHrr:", round(MHrr, dec), "   RR adjusted using MH estimate:", round(RRadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("p(Confounder+HighestLevel|Exposure+):", prev.cfder[1],
                "\np(Confounder+HighestLevel|Exposure-):", prev.cfder[2],
                "\n    p(Confounder+MidLevel|Exposure+):", prev.cfder[3],
                "\n    p(Confounder+MidLevel|Exposure-):", prev.cfder[4],
                "\n  RR(ConfounderHighestLevel-Outcome):", cfder.dis.RR[1],
                "\n      RR(ConfounderMidLevel-Outcome):", cfder.dis.RR[2],
                "\n")
        invisible(list(obs.data = tab, cfder1.data = tab.cfder1,
                       cfder2.data = tab.cfder2, nocfder.data = tab.nocfder,
                       obs.measures = rmat, SMR = SMRrr, MH = MHrr,
                       cfder1.rr = cfder1.rr, cfder2.rr = cfder2.rr,
                       nocfder.rr = nocfder.rr,
                       RRadj.smr = RRadj.smr, RRadj.mh = RRadj.mh,
                       bias.params = c(prev.cfder, cfder.dis.RR)))
    }
    if (implement == "OR"){
        crude.or <- (a/b) / (c/d)
        se.log.crude.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.crude.or <- exp(log(crude.or) - qnorm(1 - alpha/2) * se.log.crude.or)
        uci.crude.or <- exp(log(crude.or) + qnorm(1 - alpha/2) * se.log.crude.or)

        C2 <- c * prev.cfder[1]
        C1 <- c * prev.cfder[3]
        D2 <- d * prev.cfder[2]
        D1 <- d * prev.cfder[4]
        C0 <- c - C2 - C1
        D0 <- d - D2 - D1
        A0 <- (C0 * a) / (cfder.dis.OR[2] * C1 + cfder.dis.OR[1] * C2 + C0)
        B0 <- (D0 * b) / (cfder.dis.OR[2] * D1 + cfder.dis.OR[1] * D2 + D0)
        A1 <- cfder.dis.OR[2] * A0 * C1 / C0
        B1 <- cfder.dis.OR[2] * B0 * D1 / D0
        A2 <- cfder.dis.OR[1] * A0 * C2 / C0
        B2 <- cfder.dis.OR[1] * B0 * D2 / D0
        M2 <- A2 + C2
        N2 <- B2 + C2
        M1 <- A1 + C1
        N1 <- B1 + D1
        M0 <- A0 + C0
        N0 <- B0 + C0
        tab.cfder2 <- matrix(c(A2, B2, C2, D2), nrow = 2, byrow = TRUE)
        tab.cfder1 <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRor <- a / ((C2 * B2/D2) + (C1 * B1/D1) + (C0 * B0/D0))
        MHor <- (A2 * D2 / N2 + A1 * D1 / N1 + A0 * D0 / N0) /
            (B2 * C2 / N2 + B1 * C1 / N1 + B0 * C0 / N0)
        cfder2.or <- (A2 / C2) / (B2 / D2)
        cfder1.or <- (A1 / C1) / (B1 / D1)
        nocfder.or <- (A0 / C0) / (B0 / D0)
        ORadj.smr <- crude.or / SMRor
        ORadj.mh <- crude.or / MHor

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(tab.cfder2) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder2) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder2) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder2) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.cfder1) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder1) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder1) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder1) <- colnames(tab)
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
            cat("\nData, Counfounder +, Highest Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder2, dec))
        if (print)
            cat("\nData, Counfounder +, Mid-Level Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder1, dec))
        if (print)
            cat("\nData, Counfounder -:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.nocfder, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(crude.or, lci.crude.or, uci.crude.or))
        rownames(rmat) <- c("                       Crude Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Odds Ratio, Confounder +, Highest Level:", round(cfder2.or, dec), "\n    Odds Ratio, Confounder +, Mid-Level:", round(cfder1.or, dec), "\n               Odds Ratio, Confounder -:", round(nocfder.or, dec), "\n")
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
            cat("p(Confounder+HighestLevel|Exposure+):", prev.cfder[1],
                "\np(Confounder+HighestLevel|Exposure-):", prev.cfder[2],
                "\n    p(Confounder+MidLevel|Exposure+):", prev.cfder[3],
                "\n    p(Confounder+MidLevel|Exposure-):", prev.cfder[4],
                "\n  OR(ConfounderHighestLevel-Outcome):", cfder.dis.OR[1],
                "\n      OR(ConfounderMidLevel-Outcome):", cfder.dis.OR[2],
                "\n")
        invisible(list(obs.data = tab, cfder1.data = tab.cfder1,
                       cfder2.data = tab.cfder2, nocfder.data = tab.nocfder,
                       obs.measures = rmat, SMR = SMRor, MH = MHor,
                       cfder1.or = cfder1.or, cfder2.or = cfder2.or, 
                       nocfder.or = nocfder.or,
                       ORadj.smr = ORadj.smr, ORadj.mh = ORadj.mh,
                       bias.params = c(prev.cfder, cfder.dis.OR)))
    }
    if (implement == "RD"){
        crude.rd <- (a / (a + c)) - (b / (b + d))
        se.log.crude.rd <- sqrt((a * c) / (a + c)^3 + (b * d) / (b + d)^3)
        lci.crude.rd <- crude.rd - qnorm(1 - alpha/2) * se.log.crude.rd
        uci.crude.rd <- crude.rd + qnorm(1 - alpha/2) * se.log.crude.rd

        M2 <- (a + c) * prev.cfder[1]
        M1 <- (a + c) * prev.cfder[3]
        N2 <- (b + d) * prev.cfder[2]
        N1 <- (b + d) * prev.cfder[4]
        M0 <- a + c - M2 - M1
        N0 <- b + d - N2 - N1
        A0 <- M0 * (a - M2 * cfder.dis.RD[1] - M1 * cfder.dis.RD[2]) / (a + c)
        B0 <- N0 * (b - N2 * cfder.dis.RD[1] - N1 * cfder.dis.RD[2]) / (b + d)
        A1 <- M1 * cfder.dis.RD[2] + A0 * M1 / M0
        B1 <- N1 * cfder.dis.RD[2] + B0 * N1 / N0
        A2 <- M2 * cfder.dis.RD[1] + A0 * M2 / M0
        B2 <- N2 * cfder.dis.RD[1] + B0 * N2 / N0
        C2 <- M2 - A2
        D2 <- N2 - B2
        C1 <- M1 - A1
        D1 <- N1 - B1
        C0 <- M0 - A0
        D0 <- N0 - B0
        tab.cfder2 <- matrix(c(A2, B2, C2, D2), nrow = 2, byrow = TRUE)
        tab.cfder1 <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        MHrd <- (((A2 * N2 - B2 * M2) / (M2 + N2)) +
                     ((A1 * N1 - B1 * M1) / (M1 + N1)) +
                         (((A0 * N0 - B0 * M0) / (M0 + N0)))) /
                             ((M2 * N2 / (M2 + N2)) + (M1 * N1 / (M1 + N1)) +
                                  (M0 * N0 / (M0 + N0)))
        cfder2.rd <- (A2 / M2) - (B2 / N2)
        cfder1.rd <- (A1 / M1) - (B1 / N1)
        nocfder.rd <- (A0 / M0) - (B0 / N0)
        RDadj.mh <- crude.rd - MHrd

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(tab.cfder2) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder2) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder2) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder2) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.cfder1) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder1) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder1) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder1) <- colnames(tab)
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
            cat("\nData, Counfounder +, Highest Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder2, dec))
        if (print)
            cat("\nData, Counfounder +, Mid-Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder1, dec))
        if (print)
            cat("\nData, Counfounder -:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.nocfder, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(crude.rd, lci.crude.rd, uci.crude.rd))
        rownames(rmat) <- c("                       Crude Risk Difference:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Risk Difference, Confounder +, Highest Level:", round(cfder2.rd, dec), "\n    Risk Difference, Confounder +, Mid-Level:", round(cfder1.rd, dec),  "\n               Risk Difference, Confounder -:", round(nocfder.rd, dec), "\n")
        if (print)
            cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
                "\n------------------------------------------------------\n\n")
        if (print)
            cat("\nMantel-Haenszel", "               MHrd:", round(MHrd, dec), "   RD adjusted using MH estimate:", round(RDadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("p(Confounder+HighestLevel|Exposure+):", prev.cfder[1],
                "\np(Confounder+HighestLevel|Exposure-):", prev.cfder[2],
                "\n    p(Confounder+MidLevel|Exposure+):", prev.cfder[3],
                "\n    p(Confounder+MidLevel|Exposure-):", prev.cfder[4],
                "\n  RD(ConfounderHighestLevel-Outcome):", cfder.dis.RD[1],
                "\n      RD(ConfounderMidLevel-Outcome):", cfder.dis.RD[2],
                "\n")
        invisible(list(obs.data = tab, cfder1.data = tab.cfder1,
                       cfder2.data = tab.cfder2, nocfder.data = tab.nocfder,
                       obs.measures = rmat, MH = MHrd,
                       cfder1.rd = cfder1.rd, cfder2.rd = cfder2.rd,
                       nocfder.rd = nocfder.rd, RDadj.mh = RDadj.mh,
                       bias.params = c(prev.cfder, cfder.dis.RD)))
    }    
}
