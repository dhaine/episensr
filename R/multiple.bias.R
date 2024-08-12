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
    } else {
        names_data <- c(rownames(x[["obs_data"]]), colnames(x[["obs_data"]]))
        adj_cells <- x[["sim_df"]][, c("A1", "B1", "C1", "D1")]
        reps <- x[["reps"]]
        corr_data <- round(matrix(c(median(adj_cells[, "A1"], na.rm = TRUE),
                                    median(adj_cells[, "B1"], na.rm = TRUE),
                                    median(adj_cells[, "C1"], na.rm = TRUE),
                                    median(adj_cells[, "D1"], na.rm = TRUE)),
                                  dimnames = list(names_data[1:2], names_data[3:4]),
                                  nrow = 2, byrow = TRUE), 1)
    }



    arguments <- list(...)

    if (bias_function == "selection") {
        spec <- c("bias_parms", "alpha")
        if (length(arguments) > 2 | length(arguments) < 1)
            stop(cli::format_error(c("x" = "Please provide valid arguments to `selection()` bias function.")))
        if ((length(arguments) == 1) && (names(arguments) != "bias_parms"))
            stop(cli::format_error(c("x" = "Please provide valid arguments to `selection()` bias function.")))
        if ((length(arguments) == 2) && (!(names(arguments)[1] %in% spec) |
                                         !(names(arguments)[2] %in% spec)))
            stop(cli::format_error(c("x" = "Please provide valid arguments to `selection()` bias function.")))
        res <- do.call("selection", c(list(corr_data),
                                      arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "confounders") {
        spec <- c("type", "bias_parms", "alpha")
        if (length(arguments) > 3 | length(arguments) < 2)
            stop(cli::format_error(c("x" = "Please provide valid arguments to `confounders()` bias function.")))
        if ((length(arguments) == 2) && (!(names(arguments)[1] %in% spec) |
                                         !names(arguments)[2] %in% spec))
            stop(cli::format_error(c("x" = "Please provide valid arguments to `confounders()` bias function.")))
        if ((length(arguments) == 3) && (!(names(arguments)[1] %in% spec) |
                                         !(names(arguments)[2] %in% spec) |
                                         !(names(arguments)[3] %in% spec)))
            stop(cli::format_error(c("x" = "Please provide valid arguments to `confounders()` bias function.")))
        res <- do.call("confounders", c(list(corr_data),
                                        arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "misclass") {
        spec <- c("type", "bias_parms", "alpha")
        if (length(arguments) > 3 | length(arguments) < 2)
            stop(cli::format_error(c("x" = "Please provide valid arguments to `misclass()` bias function.")))
        if (length(arguments) == 2 & (!(names(arguments)[1] %in% spec) |
                                      !names(arguments)[2] %in% spec))
            stop(cli::format_error(c("x" = "Please provide valid arguments to `misclass()` bias function.")))
        if (length(arguments) == 3 & (!(names(arguments)[1] %in% spec) |
                                      !(names(arguments)[2] %in% spec) |
                                      !(names(arguments)[3] %in% spec)))
            stop(cli::format_error(c("x" = "Please provide valid arguments to `misclass()` bias function.")))
        res <- do.call("misclassification", c(list(corr_data),
                                              arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "list")
    }

    if (bias_function == "probsens") {
        spec <- c("type", "seca", "seexp", "spca", "spexp",
                  "corr_se", "corr_sp", "alpha")
        if (length(arguments) < 2)
            stop(cli::format_error(c("x" = "Please provide valid arguments to `probsens()` function.")))
        for (i in length(arguments)) {
            if (!(names(arguments)[i] %in% spec))
                stop(cli::format_error(c("x" = "Please provide valid arguments to `probsens()` function.")))
        }
#        res <- do.call("probsens", c(list(corr_data),
#                                     arguments[names(arguments) %in% spec]))
#        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    if (bias_function == "probsens.sel") {
        spec <- c("case_exp", "case_nexp", "ncase_exp", "ncase_nexp", "alpha")
        if (length(arguments) < 1)
            stop(cli::format_error(c("x" = "Please provide valid arguments to `probsens.sel()` function.")))
        for (i in length(arguments)) {
            if (!(names(arguments)[i] %in% spec))
                stop(cli::format_error(c("x" = "Please provide valid arguments to `probsens.sel()` function.")))
        }
        dist_list <- c("constant", "uniform", "triangular",
                       "trapezoidal", "normal", "beta")
        call <- match.call()
        s11 <- as.character(c(call[names(call) == "case_exp"]))
        s11_dist <- stringr::str_extract(s11, regex(dist_list))
        s11_dist <- s11_dist[!is.na(s11_dist)]
        s11_parms <- capture.output(cat(sub('.*c((.*))).*', "\\1", s11[[1]])))
        s11_parms <- unlist(regmatches(s11_parms, gregexpr('\\(?[0-9,.]+', s11_parms)))
        s11_parms <- as.numeric(gsub('\\(', '', gsub(',', '', s11_parms)))
        s11 <- list(s11_dist, s11_parms)
        s01 <- as.character(call[names(call) == "case_nexp"])
        s01_dist <- stringr::str_extract(s01, regex(dist_list))
        s01_dist <- s01_dist[!is.na(s01_dist)]
        s01_parms <- capture.output(cat(sub('.*c((.*))).*', "\\1", s01[[1]])))
        s01_parms <- unlist(regmatches(s01_parms, gregexpr('\\(?[0-9,.]+', s01_parms)))
        s01_parms <- as.numeric(gsub('\\(', '', gsub(',', '', s01_parms)))
        s01 <- list(s01_dist, s01_parms)
        s10 <- as.character(call[names(call) == "ncase_exp"])
        s10_dist <- stringr::str_extract(s10, regex(dist_list))
        s10_dist <- s10_dist[!is.na(s10_dist)]
        s10_parms <- capture.output(cat(sub('.*c((.*))).*', "\\1", s10[[1]])))
        s10_parms <- unlist(regmatches(s10_parms, gregexpr('\\(?[0-9,.]+', s10_parms)))
        s10_parms <- as.numeric(gsub('\\(', '', gsub(',', '', s10_parms)))
        s10 <- list(s10_dist, s10_parms)
        s00 <- as.character(call[names(call) == "ncase_nexp"])
        s00_dist <- stringr::str_extract(s00, regex(dist_list))
        s00_dist <- s00_dist[!is.na(s00_dist)]
        s00_parms <- capture.output(cat(sub('.*c((.*))).*', "\\1", s00[[1]])))
        s00_parms <- unlist(regmatches(s00_parms, gregexpr('\\(?[0-9,.]+', s00_parms)))
        s00_parms <- as.numeric(gsub('\\(', '', gsub(',', '', s00_parms)))
        s00 <- list(s00_dist, s00_parms)
        res <- probsens.sel(corr_data, reps = reps,
                            case_exp = list(s11_dist, s11_parms),
                            case_nexp = list(s01_dist, s01_parms),
                            ncase_exp = list(s10_dist, s10_parms),
                            ncase_nexp = list(s00_dist, s00_parms)
                            )
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    if (bias_function == "probsens_conf") {
        spec <- c("reps", "prev_exp", "prev_nexp", "risk", "corr_p", "alpha")
        if (length(arguments) < 3)
            stop(cli::format_error(c("x" = "Please provide valid arguments to `probsens_conf()` function.")))
        for (i in length(arguments)) {
            if(!(names(arguments)[i] %in% spec))
                stop(cli::format_error(c("x" = "Please provide valid arguments to `probsens_conf()` function.")))
        }
        res <- do.call("probsens_conf", c(list(corr_data),
                                          arguments[names(arguments) %in% spec]))
        class(res) <- c("episensr.multiple", "episensr", "episensr.probsens", "list")
    }

    res
}
