#' Extract adjusted 2-by-2 table from episensr object
#'
#' Extract the adjusted 2-by-2 table from an episensr function, so that it can be re-used
#' into an other episensr function.
#'
#' @param x An object of class 'episensr'.
#' 
#' @return A matrix.
#'
#' @examples
#' # The data for this example come from:
#' # Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N.
#' # et al.
#' # Population-based incidence estimates of uveal melanoma in Germany. Supplementing
#' # cancer registry data by case-control data.
#' # Eur J Cancer Prev 2006;15:165-70.
#' selection(matrix(c(136, 107, 297, 165),
#' dimnames = list(c("UM+", "UM-"), c("Mobile+", "Mobile-")),
#' nrow = 2, byrow = TRUE),
#' bias_parms = c(.94, .85, .64, .25))
#' @export
adj_tab <- function(x) {
    if(!inherits(x, "episensr"))
        stop('Not an episensr class object.')

    if(!inherits(x, "episensr.probsens")) {
        corr.data <- round(x$corr.data, 1)
        class(corr.data) <- c("episensr", "matrix")
    } else {
        obs.data <- x$obs.data
        obs.measures <- x$obs.measures
        adj.cells <- x$sim.df[, c("A1", "B1", "C1", "D1")]
        reps <- x$reps
        corr.data <- list(obs.data,
                          obs.measures,
                          adj.cells,
                          reps)
        class(corr.data) <- c("episensr", "episensr.probsens", "list")
    }
    corr.data
}
