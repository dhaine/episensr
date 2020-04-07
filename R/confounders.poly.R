#' Sensitivity analysis to correct for unknown or unmeasured polychotomous confounding without effect modification
#'
#' Simple sensitivity analysis to correct for unknown or unmeasured polychotomous
#' (3-level) confounding without effect modification. Implementation for ratio
#' measures (relative risk -- RR, or odds ratio -- OR) and difference measures
#' (risk difference -- RD).
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param type Choice of implementation, with no effect measure modification for
#' ratio measures (relative risk -- RR; odds ratio -- OR) or difference measures
#' (risk difference -- RD).
#' @param bias_parms Numeric vector defining the bias parameters. This vector
#' has 6 elements, in the following order:
#' \enumerate{
#' \item the association between the highest level confounder and the outcome,
#' \item the association between the mid-level confounder and the outcome,
#' \item the prevalence of the highest level confounder among the exposed (between 0 and 1),
#' \item the prevalence of the highest level confounder among the unexposed (between 0 and 1),
#' \item the prevalence of the mid-level confounder among the exposed (between 0 and 1), and
#' \item the prevalence of the mid-level confounder among the unexposed (between 0 and 1).
#' }
#' @param alpha Significance level.
#' 
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{cfder1.data}{The same table for Mid-level Confounder +.}
#' \item{cfder2.data}{The same table for Highest-level Confounder +.}
#' \item{nocfder.data}{The same table for Confounder -.}
#' \item{obs.measures}{A table of relative risk with confidence intervals; Total
#' and by confounders.}
#' \item{adj.measures}{A table of Standardized Morbidity Ratio and Mantel-Haenszel
#' estimates.}
#' \item{bias.parms}{Input bias parameters.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.59--78, Springer.
#' 
#' @examples
#' # The data for this example come from:
#' # Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O.
#' # et al.
#' # Increased risk of infection with human immunodeficiency virus type 1 among
#' # uncircumcised men presenting with genital ulcer disease in Kenya.
#' # Clin Infect Dis 1996;23:449-53.
#' confounders.poly(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RR",
#' bias_parms = c(.4, .8, .6, .05, .2, .2))
#' 
#' confounders.poly(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "OR",
#' bias_parms = c(.4, .8, .6, .05, .2, .2))
#' 
#' confounders.poly(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RD",
#' bias_parms = c(-.4, -.2, .6, .05, .2, .2))
#' @export
#' @importFrom stats qnorm
confounders.poly <- function(case,
                             exposed,
                             type = c("RR", "OR", "RD"),
                             bias_parms = NULL,
                             alpha = 0.05){
    if(length(type) != 1)
        stop('Choose between RR, OR, or RD implementation.')
    
    if(is.null(bias_parms))
        bias_parms <- c(1, 1, 0, 0, 0, 0)
    else bias_parms <- bias_parms
    if(length(bias_parms) != 6)
        stop('The argument bias_parms should be made of the following components: (1) Association between highest level confounder and outcome, (2) Association between mid-level confounder and outcome, (3) Prevalence of the confounder (highest level) among the exposed, (4) Prevalence of the confounder (highest level) among the unexposed, (5) Prevalence of the confounder (mid-level) among the exposed, and (6) Prevalence of the confounder (mid-level) among the unexposed.')
    if(!all(bias_parms[3:6] >= 0 & bias_parms[3:6] <=1))
        stop('Prevalences should be between 0 and 1.')
    if(bias_parms[3] + bias_parms[5] >= 1)
        stop('Sum of prevalences among the exposed >= 1.')
    if(bias_parms[4] + bias_parms[6] >= 1)
        stop('Sum of prevalences among the unexposed >= 1.')
    if(bias_parms[1] <= 0 & type != "RD")
        stop('Association between confounder and outcome should be greater than 0.')
    if(bias_parms[2] <= 0 & type != "RD")
        stop('Association between confounder and outcome should be greater than 0.')

    if(inherits(case, c("table", "matrix")))
        tab <- case
    else {
        tab.df <- table(case, exposed)
        tab <- tab.df[2:1, 2:1]
    }

    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    type <- match.arg(type)
    if (type == "RR") {
        crude.rr <- (a/(a + c)) / (b/(b + d))
        se.log.crude.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.crude.rr <- exp(log(crude.rr) - qnorm(1 - alpha/2) * se.log.crude.rr)
        uci.crude.rr <- exp(log(crude.rr) + qnorm(1 - alpha/2) * se.log.crude.rr)

        M2 <- (a + c) * bias_parms[3]
        M1 <- (a + c) * bias_parms[5]
        N2 <- (b + d) * bias_parms[4]
        N1 <- (b + d) * bias_parms[6]
        M0 <- a + c - M2 - M1
        N0 <- b + d - N2 - N1
        A0 <- (M0 * a) / (bias_parms[2] * M1 + M0 + bias_parms[1] * M2)
        B0 <- (N0 * b) / (bias_parms[2] * N1 + N0 + bias_parms[1] * N2)
        A1 <- bias_parms[2] * A0 * M1 / M0
        B1 <- bias_parms[2] * B0 * N1 / N0
        A2 <- bias_parms[1] * A0 * M2 / M0
        B2 <- bias_parms[1] * B0 * N2 / N0
        C2 <- M2 - A2
        D2 <- N2 - B2
        C1 <- M1 - A1
        D1 <- N1 - B1
        C0 <- M0 - A0
        D0 <- N0 - B0

        if(A2 < 0 | B2 < 0 | C2 < 0 | D2 < 0 |
           A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')
        
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
        rmat <- rbind(c(crude.rr, lci.crude.rr, uci.crude.rr))
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(SMRrr, RRadj.smr), c(MHrr, RRadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio:",
                             "             Mantel-Haenszel:")
        colnames(rmatc) <- c(" ", "Adjusted RR")
        rmat <- rbind(rmat, c(cfder1.rr, NA, NA), c(cfder2.rr, NA, NA),
                      c(nocfder.rr, NA, NA))
        rownames(rmat) <- c("                       Crude Relative Risk:",
                            "Relative Risk, Confounder +, Highest Level:",
                            "    Relative Risk, Confounder +, Mid-Level:",
                            "               Relative Risk, Confounder -:")
    }

    if (type == "OR"){
        crude.or <- (a/b) / (c/d)
        se.log.crude.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.crude.or <- exp(log(crude.or) - qnorm(1 - alpha/2) * se.log.crude.or)
        uci.crude.or <- exp(log(crude.or) + qnorm(1 - alpha/2) * se.log.crude.or)

        C2 <- c * bias_parms[3]
        C1 <- c * bias_parms[5]
        D2 <- d * bias_parms[4]
        D1 <- d * bias_parms[6]
        C0 <- c - C2 - C1
        D0 <- d - D2 - D1
        A0 <- (C0 * a) / (bias_parms[2] * C1 + bias_parms[1] * C2 + C0)
        B0 <- (D0 * b) / (bias_parms[2] * D1 + bias_parms[1] * D2 + D0)
        A1 <- bias_parms[2] * A0 * C1 / C0
        B1 <- bias_parms[2] * B0 * D1 / D0
        A2 <- bias_parms[1] * A0 * C2 / C0
        B2 <- bias_parms[1] * B0 * D2 / D0
        M2 <- A2 + C2
        N2 <- B2 + C2
        M1 <- A1 + C1
        N1 <- B1 + D1
        M0 <- A0 + C0
        N0 <- B0 + C0

        if(A2 < 0 | B2 < 0 | C2 < 0 | D2 < 0 |
           A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')
        
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
        rmat <- rbind(c(crude.or, lci.crude.or, uci.crude.or))
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(SMRor, ORadj.smr), c(MHor, ORadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio:",
                             "             Mantel-Haenszel:")
        colnames(rmatc) <- c(" ", "Adjusted OR")
        rmat <- rbind(rmat, c(cfder1.or, NA, NA), c(cfder2.or, NA, NA),
                      c(nocfder.or, NA, NA))
        rownames(rmat) <- c("                       Crude Odds Ratio:",
                            "Odds Ratio, Confounder +, Highest Level:",
                            "    Odds Ratio, Confounder +, Mid-Level:",
                            "               Odds Ratio, Confounder -:")
    }

    if (type == "RD"){
        crude.rd <- (a / (a + c)) - (b / (b + d))
        se.log.crude.rd <- sqrt((a * c) / (a + c)^3 + (b * d) / (b + d)^3)
        lci.crude.rd <- crude.rd - qnorm(1 - alpha/2) * se.log.crude.rd
        uci.crude.rd <- crude.rd + qnorm(1 - alpha/2) * se.log.crude.rd

        M2 <- (a + c) * bias_parms[3]
        M1 <- (a + c) * bias_parms[5]
        N2 <- (b + d) * bias_parms[4]
        N1 <- (b + d) * bias_parms[6]
        M0 <- a + c - M2 - M1
        N0 <- b + d - N2 - N1
        A0 <- M0 * (a - M2 * bias_parms[1] - M1 * bias_parms[2]) / (a + c)
        B0 <- N0 * (b - N2 * bias_parms[1] - N1 * bias_parms[2]) / (b + d)
        A1 <- M1 * bias_parms[2] + A0 * M1 / M0
        B1 <- N1 * bias_parms[2] + B0 * N1 / N0
        A2 <- M2 * bias_parms[1] + A0 * M2 / M0
        B2 <- N2 * bias_parms[1] + B0 * N2 / N0
        C2 <- M2 - A2
        D2 <- N2 - B2
        C1 <- M1 - A1
        D1 <- N1 - B1
        C0 <- M0 - A0
        D0 <- N0 - B0

        if(A2 < 0 | B2 < 0 | C2 < 0 | D2 < 0 |
           A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')
        
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
        rmat <- rbind(c(crude.rd, lci.crude.rd, uci.crude.rd))
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(MHrd, RDadj.mh))
        rownames(rmatc) <- "Mantel-Haenszel:"
        colnames(rmatc) <- c(" ", "Adjusted RD")
        rmat <- rbind(rmat, c(cfder1.rd, NA, NA), c(cfder2.rd, NA, NA),
                      c(nocfder.rd, NA, NA))
        rownames(rmat) <- c("                       Crude Risk Difference:",
                            "Risk Difference, Confounder +, Highest Level:",
                            "    Risk Difference, Confounder +, Mid-Level:",
                            "               Risk Difference, Confounder -:")
    }
    res <- list(obs.data = tab,
                cfder1.data = tab.cfder1,
                cfder2.data = tab.cfder2,
                nocfder.data = tab.nocfder,
                obs.measures = rmat,
                adj.measures = rmatc, 
                bias.parms = bias_parms)
    class(res) <- c("episensr", "list")
    res
}
