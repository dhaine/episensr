#' Extract adjusted 2-by-2 table from episensr object
#'
#' Extract the adjusted 2-by-2 table from an \code{episensr} function, so that it can
#' be re-used into an other \code{episensr} function when performing multiple (combined)
#' bias analyses.
#' Allowed functions are: \code{selection}, \code{misclass}, \code{confounders},
#' \code{probsens}, \code{probsens.sel}, and \code{probsens_conf}.
#'
#' For probabilistic bias analyses, median of cells are passed to the next function as
#' starting 2-by-2 table.
#'
#' @param x An object of class 'episensr' or 'episensr.probsens'.
#' @param bias_function Bias function to be called. Choices between 'selection',
#' 'misclass', 'confounders', 'probsens', 'probsens.sel', 'probsens_conf'.
#' @param ... Additional arguments passed on to methods.
#'
#' @return A list with the elements corresponding to the bias function called.
#'
#' @examples
#' dat <- matrix(c(118, 832, 103, 884),
#' dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")), nrow = 2, byrow = TRUE)
#'
#' dat %>%
#' misclass(., type = "exposure", bias_parms = c(.56, .58, .99, .97)) %>%
#' multiple.bias(., bias_function = "selection", bias_parms = c(.73, .61, .82, .76))
#' @seealso \code{\link{selection}}, \code{\link{misclass}},
#' \code{\link{confounders}}, \code{\link{probsens}}, \code{\link{probsens.sel}},
#' \code{\link{probsens_conf}}
#' @export
multiple.bias <- function(x,
                          bias_function = c("selection", "misclass",
                                            "confounders", "probsens.sel",
                                            "probsens_conf", "probsens"),
                          ...) {
    if (!inherits(x, "episensr"))
        stop(cli::format_error(c("x" = "Not an episensr class object.")))

    if (missing(bias_function)) stop(cli::format_error(c("x" = "Please provide a valid bias function.")))

    if (!(bias_function %in% c("selection", "misclass", "confounders",
                               "probsens.sel", "probsens_conf", "probsens")))
        stop(cli::format_error(c("x" = "Please provide a valid bias function.")))

    if (!inherits(x, "episensr.probsens")) {
        corr_data <- round(x[["corr_data"]], 1)
        class(corr_data) <- c("episensr", "matrix")
    } else {
        obs_data <- x[["obs_data"]]
        obs_measures <- x[["obs_measures"]]
        adj_cells <- x[["sim_df"]][, c("A1", "B1", "C1", "D1")]
        reps <- x[["reps"]]
        corr_data <- round(matrix(c(median(adj_cells[, "A1"], na.rm = TRUE),
                                    median(adj_cells[, "B1"], na.rm = TRUE),
                                    median(adj_cells[, "C1"], na.rm = TRUE),
                                    median(adj_cells[, "D1"], na.rm = TRUE)),
                                  dimnames = list(rownames(obs_data), colnames(obs_data)),
                                  nrow = 2, byrow = TRUE), 1)
        class(corr_data) <- c("episensr", "matrix")
    }

    arguments <- list(...)

    if (bias_function == "selection") {
        spec <- c("bias_parms", "alpha")
        if (length(arguments) > 2 | length(arguments) < 1)
            stop(cli::format_error(c("x" = "Please provide valid arguments to selection bias function.")))
        if ((length(arguments) == 1) && (names(arguments) != "bias_parms"))
            stop(cli::format_error(c("x" = "Please provide valid arguments to selection bias function.")))
        if ((length(arguments) == 2) && (!(names(arguments)[1] %in% spec) |
                                         !(names(arguments)[2] %in% spec)))
            stop(cli::format_error(c("x" = "Please provide valid arguments to selection bias function.")))
        res <- do.call("selection", c(list(corr_data),
                                      arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "confounders") {
        spec <- c("type", "bias_parms", "alpha")
        if (length(arguments) > 3 | length(arguments) < 2)
            stop(cli::format_error(c("x" = "Please provide valid arguments to confounder bias function.")))
        if ((length(arguments) == 2) && (!(names(arguments)[1] %in% spec) |
                                         !names(arguments)[2] %in% spec))
            stop(cli::format_error(c("x" = "Please provide valid arguments to confounder bias function.")))
        if ((length(arguments) == 3) && (!(names(arguments)[1] %in% spec) |
                                         !(names(arguments)[2] %in% spec) |
                                         !(names(arguments)[3] %in% spec)))
            stop(cli::format_error(c("x" = "Please provide valid arguments to confounder bias function.")))
        res <- do.call("confounders", c(list(corr_data),
                                        arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "misclass") {
        spec <- c("type", "bias_parms", "alpha")
        if (length(arguments) > 3 | length(arguments) < 2)
            stop(cli::format_error(c("x" = "Please provide valid arguments to misclassification bias function.")))
        if (length(arguments) == 2 & (!(names(arguments)[1] %in% spec) |
                                      !names(arguments)[2] %in% spec))
            stop(cli::format_error(c("x" = "Please provide valid arguments to misclassification bias function.")))
        if (length(arguments) == 3 & (!(names(arguments)[1] %in% spec) |
                                      !(names(arguments)[2] %in% spec) |
                                      !(names(arguments)[3] %in% spec)))
            stop(cli::format_error(c("x" = "Please provide valid arguments to misclassification bias function.")))
        res <- do.call("misclassification", c(list(corr_data),
                                              arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "probsens") {
        spec <- c("type", "reps", "seca", "seexp", "spca", "spexp",
                  "corr_se", "corr_sp", "alpha")
        if (length(arguments) < 2)
            stop(cli::format_error(c("x" = "Please provide valid arguments to probsens function.")))
        for (i in length(arguments)) {
            if (!(names(arguments)[i] %in% spec))
                stop(cli::format_error(c("x" = "Please provide valid arguments to probsens function.")))
        }
        res <- do.call("probsens", c(list(corr_data),
                                     arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    if (bias_function == "probsens.sel") {
        spec <- c("reps", "case_exp", "case_nexp", "ncase_exp", "ncase_nexp", "alpha")
        if (length(arguments) < 1)
            stop(cli::format_error(c("x" = "Please provide valid arguments to probsens function.")))
        for (i in length(arguments)) {
            if (!(names(arguments)[i] %in% spec))
                stop(cli::format_error(c("x" = "Please provide valid arguments to probsens function.")))
        }
        res <- do.call("probsens.sel", c(list(corr_data),
                                         arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    if (bias_function == "probsens_conf") {
        spec <- c("reps", "prev_exp", "prev_nexp", "risk", "corr_p", "alpha")
        if (length(arguments) < 3)
            stop(cli::format_error(c("x" = "Please provide valid arguments to probsens function.")))
        for (i in length(arguments)) {
            if(!(names(arguments)[i] %in% spec))
                stop(cli::format_error(c("x" = "Please provide valid arguments to probsens function.")))
        }
        res <- do.call("probsens_conf", c(list(corr_data),
                                          arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    res
}
