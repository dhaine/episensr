#' Bootstrap resampling for selection and misclassification bias models.
#'
#' Generate \code{R} bootstrap replicates of either selection or misclassification bias functions.
#' It then generates a confidence interval of the parameter, by first order normal approximation or the bootstrap percentile interval.
#' Replicates giving negative cell(s) in the adjusted 2-by-2 table are silently ignored.
#'
#' @param bias_model An object of class "episensr.boot", i.e. either selection bias function or misclassification bias function.
#' @param R The number of bootstrap replicates.
#' @param conf Confidence level.
#' @param ci_type A character string giving the type of interval required. Values can be either "norm" or "perc", default to "norm".
#'
#' @return A list with elements:
#' \item{model}{Model ran.}
#' \item{boot_mod}{Bootstrap resampled object, of class \code{boot}.}
#' \item{nrep}{Number of replicates used.}
#' \item{bias_ciRR}{Bootstrap confidence interval object for relative risk.}
#' \item{bias_ciOR}{Bootstrap confidence interval object for odds ratio.}
#' \item{ci}{Confidence intervals for the bias adjusted association measures.}
#' \item{conf}{Confidence interval.}
#'
#' @seealso \code{\link{boot}, \link{selection}, \link{misclass}}
#'
#' @examples
#' misclass_eval <- misclass(matrix(c(215, 1449, 668, 4296),
#' dimnames = list(c("Breast cancer+", "Breast cancer-"),
#' c("Smoker+", "Smoker-")),
#' nrow = 2, byrow = TRUE),
#' type = "exposure",
#' bias_parms = c(.78, .78, .99, .99))
#'
#' set.seed(123)
#' boot.bias(misclass_eval)
#' @export
boot.bias <- function(bias_model,
                      R = 1000,
                      conf = 0.95,
                      ci_type = c("norm", "perc")
                      ) {
    if (!inherits(bias_model, "episensr.boot"))
        stop(cli::format_error(c("x" = 'Not an episensr.boot class object.')))
    if (R < 1)
        stop(cli::format_error(c("i" = 'Please provide a sensible number of replicates to run.')))

    model <- bias_model$model
    if (model == "misclassification") type <- bias_model$type
    obs_table <- bias_model$obs_data
    obs_df <- callback_df(data.frame(x = gl(2, 1), y = gl(2, 2)),
                          c(obs_table[2, 2], obs_table[2, 1],
                            obs_table[1, 2], obs_table[1, 1]))
    bias <- bias_model$bias_parms
    ci_type <- match.arg(ci_type)

    if (model == "misclassification") {
        boot_fun <- function(data, indices) {
            d <- data[indices, ]
            bias_boot <- tryCatch({
                                      misclass(d$y, d$x,
                                               type = type,
                                               bias_parms = bias)$adj_measures
                                  },
                error = function(err) {
                    return(c(NA, NA))
                }
            )
            return(bias_boot)
        }
    }

    if (model == "selection") {
        boot_fun <- function(data, indices) {
            d <- data[indices, ]
            bias_boot <- tryCatch({
                                      selection(d$y, d$x,
                                                bias_parms = bias)$adj_measures
                                  },
                error = function(err) {
                    return(c(NA, NA))
                }
            )
            return(bias_boot)
        }
    }

    boot_mod <- boot::boot(data = obs_df, statistic = boot_fun, R = R)
    nrep <- length(which(!is.na(boot_mod$t[, 1])))

    if (ci_type == "norm") {
        bias_ci1 <- boot::boot.ci(boot_mod,
                                  t0 = log(boot_mod$t0[1]),
                                  t = log(boot_mod$t[, 1]),
                                  conf = conf,
                                  type = "norm",
                                  hinv = exp)
        bias_ci2 <- boot::boot.ci(boot_mod,
                                  t0 = log(boot_mod$t0[2]),
                                  t = log(boot_mod$t[, 2]),
                                  conf = conf,
                                  type = "norm",
                                  hinv = exp)
        rci <- rbind(c(bias_ci1[[4]][2:3]), c(bias_ci2[[4]][2:3]))
    } else if (ci_type == "perc") {
        bias_ci1 <- boot::boot.ci(boot_mod,
                                  t0 = boot_mod$t0[1],
                                  t = boot_mod$t[, 1],
                                  conf = conf,
                                  type = "perc")
        bias_ci2 <- boot::boot.ci(boot_mod,
                                  t0 = boot_mod$t0[2],
                                  t = boot_mod$t[, 2],
                                  conf = conf,
                                  type = "perc")
        rci <- rbind(c(bias_ci1[[4]][4:5]), c(bias_ci2[[4]][4:5]))
        rownames(rci) <- c("RR:", "OR:")
        colnames(rci) <- NULL
    }

    res <- list(model = model,
                boot_mod = boot_mod,
                nrep = nrep,
                bias_ciRR = bias_ci1,
                bias_ciOR = bias_ci2,
                ci = rci,
                conf = conf)
    class(res) <- c("episensr.booted", "list")
    res
}
