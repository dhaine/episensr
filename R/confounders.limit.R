confounders.limit <- function(p = NA,
                              RR = NA,
                              OR = NA,
                              crude.RR = NULL,
                              dec = 4,
                              print = TRUE){
    if(is.null(crude.RR))
        stop('Please provide crude relative risk.')
    if(is.null(p) & is.null(RR) & is.null(OR))
        stop('Not enough information.')

    q <- ifelse(is.null(p), NA, 1 - p)
    lower.bound <- crude.RR /
        min(RR, OR, 1/p, RR/(q+RR*p), OR/(q+OR*p), na.rm = TRUE)
    upper.bound <- crude.RR

    if (print)
        cat("\nLimits on adjusted RR:", round(lower.bound, dec),
            "<= RRadj <=", round(upper.bound, dec), "\n")
    if (print)
        cat("\nInput Bias Parameters:",
            "\n----------------\n\n")
    if (print)
        cat("  p(Confounder+|Exposure-):", p,
            "\n    RR(Confounder-Disease):", RR,
            "\n   OR(Confounder-Exposure):", OR,
            "\nCrude RR(Exposure-Disease):", crude.RR, "\n")
    invisible(list(conf.limits = c(lower.bound, upper.bound),
                   bias.parms = c(p, RR, OR, crude.RR)))
}
