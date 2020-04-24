#' Given a contingency table (2-by-2) (as a data frame or matrix)
#' it will untabulate it by repeating each row by the number of times
#' it was repeated.
#'
#' @param matrix or data.frame to untable
#' @param vector of counts (of same length as \code{df})
#' @keywords internal
callback_df <- function(df, num) {
    df[rep(1:nrow(df), num), ]
}
