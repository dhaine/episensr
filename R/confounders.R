#' Sensitivity analysis to correct for unknown or unmeasured confounding without
#' effect modification
#'
#' Simple sensitivity analysis to correct for unknown or unmeasured confounding
#' without effect modification. Implementation for ratio measures (relative risk
#' -- RR, or odds ratio -- OR) and difference measures (risk difference -- RD).
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param type Choice of implementation, with no effect measure modification for
#' ratio measures (relative risk -- RR; odds ratio -- OR) or difference measures
#' (risk difference -- RD).
#' @param bias_parms Numeric vector defining the 3 necessary bias parameters. This
#' vector has 3 elements, in the following order:
#' \enumerate{
#' \item the association between the confounder and the outcome among those who were not exposed (RR, OR, or RD according to choice of implementation),
#' \item the prevalence of the confounder among the exposed (between 0 and 1), and
#' \item the prevalence of the confounder among the unexposed (between 0 and 1).
#' }
#' @param alpha Significance level.
#' 
#' @return A list with elements:
#' \item{obs.data}{The analyzed 2 x 2 table from the observed data.}
#' \item{cfder.data}{The same table for Confounder +.}
#' \item{nocfder.data}{The same table for Confounder -.}
#' \item{obs.measures}{A table of relative risk with confidence intervals; for
#' Total, Confounder +, and Confounder -.}
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
#' confounders(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RR",
#' bias_parms = c(.63, .8, .05))
#' 
#' confounders(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "OR",
#' bias_parms = c(.63, .8, .05))
#' 
#' confounders(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RD",
#' bias_parms = c(-.37, .8, .05))
#' @export
#' @importFrom stats qnorm
confounders <- function(case,
                        exposed,
                        type = c("RR", "OR", "RD"),
                        bias_parms = NULL,
                        alpha = 0.05){
    if(length(type) > 1)
        stop('Choose between RR, OR, or RD implementation.')

    if(is.null(bias_parms))
        bias_parms <- c(1, 0, 0)
    else bias_parms <- bias_parms
    if(length(bias_parms) != 3)
        stop('The argument bias_parms should be made of the following components: (1) Association between the confounder and the outcome among those who were not exposed, (2) Prevalence of the confounder among the exposed, and (3) Prevalence of the confounder among the unexposed.')
    if(!all(bias_parms[-1] >= 0 & bias_parms[-1] <= 1))
        stop('Prevalences should be between 0 and 1.')
    if(bias_parms[1] <= 0 & type != "RD")
        stop('Association between the confounder and the outcome among those who were not exposed should be greater than 0.')
    
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

        M1 <- (a + c) * bias_parms[2]
        N1 <- (b + d) * bias_parms[3]
        A1 <- (bias_parms[1] * M1 * a) / (bias_parms[1] * M1 + (a + c) - M1)
        B1 <- (bias_parms[1] * N1 * b) / (bias_parms[1] * N1 + (b + d) - N1)
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
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')

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
        rmat <- rbind(c(crude.rr, lci.crude.rr, uci.crude.rr))
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(SMRrr, RRadj.smr), c(MHrr, RRadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio:",
                             "             Mantel-Haenszel:")
        colnames(rmatc) <- c(" ", "Adjusted RR")
        rmat <- rbind(rmat, c(cfder.rr, NA, NA), c(nocfder.rr, NA, NA))
        rownames(rmat) <- c("        Crude Relative Risk:",
                            "Relative Risk, Confounder +:",
                            "Relative Risk, Confounder -:")
    }

    if (type == "OR"){
        crude.or <- (a/b) / (c/d)
        se.log.crude.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.crude.or <- exp(log(crude.or) - qnorm(1 - alpha/2) * se.log.crude.or)
        uci.crude.or <- exp(log(crude.or) + qnorm(1 - alpha/2) * se.log.crude.or)

        C1 <- c * bias_parms[2]
        D1 <- d * bias_parms[3]
        A1 <- (bias_parms[1] * C1 * a) / (bias_parms[1] * C1 + c - C1)
        B1 <- (bias_parms[1] * D1 * b) / (bias_parms[1] * D1 + d - D1)
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
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')

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
        rmat <- rbind(c(crude.or, lci.crude.or, uci.crude.or))
        colnames(rmat) <- c(" ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(SMRor, ORadj.smr), c(MHor, ORadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio:",
                             "             Mantel-Haenszel:")
        colnames(rmatc) <- c(" ", "Adjusted OR")
        rmat <- rbind(rmat, c(cfder.or, NA, NA), c(nocfder.or, NA, NA))
        rownames(rmat) <- c("        Crude Odds Ratio:",
                            "Odds Ratio, Confounder +:",
                            "Odds Ratio, Confounder -:")
    }

    if (type == "RD"){
        crude.rd <- (a / (a + c)) - (b / (b + d))
        se.log.crude.rd <- sqrt((a * c) / (a + c)^3 + (b * d) / (b + d)^3)
        lci.crude.rd <- crude.rd - qnorm(1 - alpha/2) * se.log.crude.rd
        uci.crude.rd <- crude.rd + qnorm(1 - alpha/2) * se.log.crude.rd

        M1 <- (a + c) * bias_parms[2]
        N1 <- (b + d) * bias_parms[3]
        M0 <- (a + c) - M1
        N0 <- (b + d) - N1
        A1 <- (bias_parms[1] * M1 * M0 + M1 * a) / (a + c)
        B1 <- (bias_parms[1] * N1 * N0 + N1 * b) / (b + d)
        C1 <- M1 - A1
        D1 <- N1 - B1
        A0 <- a - A1
        B0 <- b - B1
        C0 <- c - C1
        D0 <- d - D1

        if(A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')
        
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
        rmat <- rbind(c(crude.rd, lci.crude.rd, uci.crude.rd))
        colnames(rmat) <- c("     ",
                            paste(100 * (alpha/2), "%", sep = ""),
                            paste(100 * (1 - alpha/2), "%", sep = ""))
        rmatc <- rbind(c(MHrd, RDadj.mh))
        rownames(rmatc) <- "Mantel-Haenszel:"
        colnames(rmatc) <- c(" ", "Adjusted RD")
        rmat <- rbind(rmat, c(cfder.rd, NA, NA), c(nocfder.rd, NA, NA))
        rownames(rmat) <- c("        Crude Risk Difference:",
                            "Risk Difference, Confounder +:",
                            "Risk Difference, Confounder -:")
    }
    res <- list(obs.data = tab,
                cfder.data = tab.cfder,
                nocfder.data = tab.nocfder,
                obs.measures = rmat,
                adj.measures = rmatc, 
                bias.parms = bias_parms)
    class(res) <- c("episensr", "list")
    res
}
