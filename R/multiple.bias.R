#' Extract adjusted 2-by-2 table from episensr object
#'
#' Extract the adjusted 2-by-2 table from an \code{episensr} function, so that it can
#' be re-used into an other \code{episensr} function when performing multiple (combined)
#' bias analysis.
#' Allowed functions are: \code{selection}, \code{misclassification}, \code{confounders},
#' \code{probsens}, \code{probsens.sel}, and \code{probsens.conf}.
#'
#' For probabilistic bias analyses, median of cells are passed to the next function as
#' starting 2-by-2 table.
#'
#' @param x An object of class 'episensr' or 'episensr.probsens'.
#' @param bias_function Bias function to be called. Choices between 'selection',
#' 'misclassification', 'confounders', 'probsens', 'probsens.sel', 'probsens.conf'.
#' @param ... Additional arguments passed on to methods.
#'
#' @return A list with the elements corresponding to the bias function called.
#'
#' @examples
#' dat <- matrix(c(118, 832, 103, 884),
#' dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")), nrow = 2, byrow = TRUE)
#'
#' dat %>%
#' misclassification(., type = "exposure", bias_parms = c(.56, .58, .99, .97)) %>%
#' multiple.bias(., bias_function = "selection", bias_parms = c(.73, .61, .82, .76))
#' @seealso \code{\link{selection}}, \code{\link{misclassification}},
#' \code{\link{confounders}}, \code{\link{probsens}}, \code{\link{probsens.sel}},
#' \code{\link{probsens.conf}}
#' @export
multiple.bias <- function(x,
                          bias_function = c("selection", "misclassification",
                                            "confounders", "probsens.sel",
                                            "probsens.conf", "probsens"),
                          ...) {
    if (!inherits(x, "episensr"))
        stop("Not an episensr class object.")

    if (missing(bias_function)) stop("Please provide a valid bias function.")

    if (!(bias_function %in% c("selection", "misclassification", "confounders",
                               "probsens.sel", "probsens.conf", "probsens")))
        stop("Please provide a valid bias function.")

    if (!inherits(x, "episensr.probsens")) {
        corr.data <- round(x[["corr.data"]], 1)
        class(corr.data) <- c("episensr", "matrix")
    } else {
        obs.data <- x[["obs.data"]]
        obs.measures <- x[["obs.measures"]]
        adj.cells <- x[["sim.df"]][, c("A1", "B1", "C1", "D1")]
        reps <- x[["reps"]]
        corr.data <- round(matrix(c(median(adj.cells[, "A1"], na.rm = TRUE),
                                    median(adj.cells[, "B1"], na.rm = TRUE),
                                    median(adj.cells[, "C1"], na.rm = TRUE),
                                    median(adj.cells[, "D1"], na.rm = TRUE)),
                                  dimnames = list(rownames(obs.data), colnames(obs.data)),
                                  nrow = 2, byrow = TRUE), 1)
        class(corr.data) <- c("episensr", "matrix")
    }

    arguments <- list(...)

    if (bias_function == "selection") {
        spec <- c("bias_parms", "alpha")
        if (length(arguments) > 2 | length(arguments) < 1)
            stop("Please provide valid arguments to selection bias function.")
        if ((length(arguments) == 1) && (names(arguments) != "bias_parms"))
            stop("Please provide valid arguments to selection bias function.")
        if ((length(arguments) == 2) && (!(names(arguments)[1] %in% spec) |
                                         !(names(arguments)[2] %in% spec)))
            stop("Please provide valid arguments to selection bias function.")
        res <- do.call("selection", c(list(corr.data),
                                      arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "confounders") {
        spec <- c("type", "bias_parms", "alpha")
        if (length(arguments) > 3 | length(arguments) < 2)
            stop("Please provide valid arguments to confounder bias function.")
        if ((length(arguments) == 2) && (!(names(arguments)[1] %in% spec) |
                                         !names(arguments)[2] %in% spec))
            stop("Please provide valid arguments to confounder bias function.")
        if ((length(arguments) == 3) && (!(names(arguments)[1] %in% spec) |
                                         !(names(arguments)[2] %in% spec) |
                                         !(names(arguments)[3] %in% spec)))
            stop("Please provide valid arguments to confounder bias function.")
        res <- do.call("confounders", c(list(corr.data),
                                        arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "misclassification") {
        spec <- c("type", "bias_parms", "alpha")
        if (length(arguments) > 3 | length(arguments) < 2)
            stop("Please provide valid arguments to misclassification bias function.")
        if (length(arguments) == 2 & (!(names(arguments)[1] %in% spec) |
                                      !names(arguments)[2] %in% spec))
            stop("Please provide valid arguments to misclassification bias function.")
        if (length(arguments) == 3 & (!(names(arguments)[1] %in% spec) |
                                      !(names(arguments)[2] %in% spec) |
                                      !(names(arguments)[3] %in% spec)))
            stop("Please provide valid arguments to misclassification bias function.")
        res <- do.call("misclassification", c(list(corr.data),
                                              arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "probsens") {
        spec <- c("type", "reps", "seca.parms", "seexp.parms", "spca.parms", "spexp.parms",
                  "corr.se", "corr.sp", "discard", "alpha")
        if (length(arguments) < 2)
            stop("Please provide valid arguments to probsens function.")
        for (i in length(arguments)) {
            if (!(names(arguments)[i] %in% spec))
                stop("Please provide valid arguments to probsens function.")
        }
        res <- do.call("probsens", c(list(corr.data),
                                     arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    if (bias_function == "probsens.sel") {
        spec <- c("reps", "or.parms", "case.exp", "case.nexp", "ncase.exp",
                  "ncase.nexp", "alpha")
        if (length(arguments) < 1)
            stop("Please provide valid arguments to probsens function.")
        for (i in length(arguments)) {
            if (!(names(arguments)[i] %in% spec))
                stop("Please provide valid arguments to probsens function.")
        }
        res <- do.call("probsens.sel", c(list(corr.data),
                                         arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    if (bias_function == "probsens.conf") {
        spec <- c("reps", "prev.exp", "prev.nexp", "risk", "corr.p", "discard", "alpha")
        if (length(arguments) < 3)
            stop("Please provide valid arguments to probsens function.")
        for (i in length(arguments)) {
            if(!(names(arguments)[i] %in% spec))
                stop("Please provide valid arguments to probsens function.")
        }
        res <- do.call("probsens.conf", c(list(corr.data),
                                          arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    res
}
