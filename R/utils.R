# Given a contingency table (2-by-2) (as a data frame or matrix)
# it will untabulate it by repeating each row by the number of times
# it was repeated.
#
# @param matrix or data.frame to untable
# @param vector of counts (of same length as \code{df})
# @keywords internal
callback_df <- function(df, num) {
    df[rep(1:nrow(df), num), ]
}

# Logit-log distribution function
#
# @param sesp Distribution parameters as list
# @keywords internal
logitlog.dstr <- function(sesp) {
    u <- runif(sesp[[1]])
    w <- sesp[[2]] + sesp[[3]] * (log(u / (1 - u)))
    p <- sesp[[4]] + (sesp[[5]] - sesp[[4]]) * exp(w) / (1 + exp(w))
    return(p)
}

# Logit-normal distribution function
#
# @param sesp Distribution parameters as list
# @keywords internal
logitnorm.dstr <- function(sesp) {
    u <- runif(sesp[[1]])
    w <- sesp[[2]] + sesp[[3]] * qnorm(u)
    p <- sesp[[4]] + (sesp[[5]] - sesp[[4]]) * exp(w) / (1 + exp(w))
    return(p)
}
