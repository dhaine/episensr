#' Pipe bias functions
#'
#' episensr also uses the pipe function, \code{\%>\%} to turn
#' function composition into a series of imperative statements.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs Data or bias function and a function to apply to it
#' @examples
#' # Instead of
#' misclass(matrix(c(118, 832, 103, 884),
#' dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")), nrow = 2, byrow = TRUE),
#' type = "exposure", bias_parms = c(.56, .58, .99, .97))
#' # you can write
#' dat <- matrix(c(118, 832, 103, 884),
#' dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")), nrow = 2, byrow = TRUE)
#' dat %>% misclass(., type = "exposure", bias_parms = c(.56, .58, .99, .97))
NULL
