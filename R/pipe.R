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
#' misclassification(matrix(c(118, 832, 103, 884),
#' dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")), nrow = 2, byrow = TRUE),
#' type = "exposure", bias_parms = c(.56, .58, .99, .97))
#' # you can write
#' dat <- matrix(c(118, 832, 103, 884),
#' dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")), nrow = 2, byrow = TRUE)
#' dat %>% misclassification(., type = "exposure", bias_parms = c(.56, .58, .99, .97))
#' # also for multiple bias:
#' chien %>%
#' misclassification(., type = "exposure", bias_parms = c(.56, .58, .99, .97)) %>%
#' multiple.bias(., bias_function = "selection", bias_parms = c(.73, .61, .82, .76))
NULL
