misclassification <- function(exposed, case, implement = c("exposure", "outcome", "confounder"), bias = NULL, alpha = 0.05, dec = 4, print = TRUE){
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

        corr.a <- (a - (1 - bias[3]) * (a + b)) / (bias[1] - (1 - bias[3]))
        corr.c <- (c - (1 - bias[4]) * (c + d)) / (bias[2] - (1 - bias[4]))
        corr.b <- (a + b) - corr.a
        corr.d <- (c + d) - corr.c
        corr.tab <- matrix(c(corr.a, corr.b, corr.c, corr.d), nrow = 2, byrow = TRUE)

        corr.rr <- (corr.a/(corr.a + corr.c)) / (corr.b/(corr.b + corr.d))
        corr.or <- (corr.a/corr.b) / (corr.c/corr.d)

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
        invisible(list(obs.data = tab, corr.data = corr.tab,
                       obs.measures = rmat, 
                       corr.rr = corr.rr, corr.or = corr.or,
                       bias.params = bias))
    }
    if (implement == "OR"){
        crude.or <- (a/b) / (c/d)
        se.log.crude.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.crude.or <- exp(log(crude.or) - qnorm(1 - alpha/2) * se.log.crude.or)
        uci.crude.or <- exp(log(crude.or) + qnorm(1 - alpha/2) * se.log.crude.or)

        C1 <- c * prev.cfder[1] 
        D1 <- d * prev.cfder[2]
        A1 <- (cfder.dis.OR * C1 * a) / (cfder.dis.OR * C1 + c - C1)
        B1 <- (cfder.dis.OR * D1 * b) / (cfder.dis.OR * D1 + d - D1)
        M1 <- A1 + C1
        N1 <- B1 + D1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1
        M0 <- A0 + C0
        N0 <- B0 + C0
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
        rmat <- rbind(c(crude.or, lci.crude.or, uci.crude.or))
        rownames(rmat) <- c("        Crude Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
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
            cat("p(Confounder+|Exposure+):", prev.cfder[1],
                "\np(Confounder+|Exposure-):", prev.cfder[2],
                "\n  OR(Confounder-Outcome):", cfder.dis.OR, "\n")
        invisible(list(obs.data = tab, cfder.data = tab.cfder,
                       nocfder.data = tab.nocfder,
                       obs.measures = rmat, SMR = SMRor, MH = MHor,
                       cfder.or = cfder.or, nocfder.or = nocfder.or,
                       ORadj.smr = ORadj.smr, RRadj.mh = ORadj.mh,
                       bias.params = c(prev.cfder, cfder.dis.OR)))
    }
    if (implement == "RD"){
        crude.rd <- (a / (a + c)) - (b / (b + d))
        se.log.crude.rd <- sqrt((a * c) / (a + c)^3 + (b * d) / (b + d)^3)
        lci.crude.rd <- crude.rd - qnorm(1 - alpha/2) * se.log.crude.rd
        uci.crude.rd <- crude.rd + qnorm(1 - alpha/2) * se.log.crude.rd

        M1 <- (a + c) * prev.cfder[1]
        N1 <- (b + d) * prev.cfder[2]
        M0 <- (a + c) - M1
        N0 <- (b + d) - N1
        A1 <- (cfder.dis.RD * M1 * M0 + M1 * a) / (a + c)
        B1 <- (cfder.dis.RD * N1 * N0 + N1 * b) / (b + d)
        C1 <- M1 - A1
        D1 <- N1 - B1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1
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
        rmat <- rbind(c(crude.rd, lci.crude.rd, uci.crude.rd))
        rownames(rmat) <- c("        Crude Risk Difference:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
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
            cat("p(Confounder+|Exposure+):", prev.cfder[1],
                "\np(Confounder+|Exposure-):", prev.cfder[2],
                "\n  RD(Confounder-Outcome):", cfder.dis.RD, "\n")
        invisible(list(obs.data = tab, cfder.data = tab.cfder,
                       nocfder.data = tab.nocfder,
                       obs.measures = rmat, MH = MHrd,
                       cfder.rd = cfder.rd, nocfder.rd = nocfder.rd,
                       RDadj.mh = RDadj.mh,
                       bias.params = c(prev.cfder, cfder.dis.RD)))
    }    
}
