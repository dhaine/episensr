#' Boot for bias model.
#'
#' Bootstrap sample.
#'
#' @param model.
#'
#' @return
#'
#' @references
#'
#' @examples
#'
#' @export
#' @importFrom stats
boot.bias <- function(model){
    if(class(model) != episensr.boot)
        stop('Not an episensr.boot class object.')

    model <- t
    model <- t2
    obs_table <- model$obs.data
    obs_df <- untable(data.frame(x = gl(2, 1), y = gl(2, 2)),
                      c(obs_table[2, 2], obs_table[1, 2],
                        obs_table[2, 1], obs_table[1, 1]))
    bias <- model$bias.parms

    f2 <- function(data, indices){
        d <- data[indices, ]
        bias_mod <- tryCatch(
        {
            misclassification(d$y, d$x, type = "exposure",
                              bias_parms = model$bias.parms)$adj.measures[1]
        },
        error = function(err){
            return(NA)
        }
        )
        #bias <- bias_mod$adj.measures
        return(bias_mod)
    }
    f <- function(data, indices){
        d <- data[indices, ]
        bias_mod <- misclassification(d$y, d$x,
                                      type = "exposure",
                                      bias_parms = model$bias.parms)$adj.measures[1]
        return(bias_mod)
    }

    set.seed(123)
    boot_bias <- boot(data = obs_df, statistic = f, R = 1000)
    boot_bias2 <- boot(data = obs_df, statistic = f2, R = 1000)
    plot(boot_bias2)
    ci <- boot.ci(boot_bias2, index = 1, conf=.95,type='perc')
    ci2 <- ci$percent[, c(4, 5)]
    hist(boot_bias$t[, 1], main = 'Bias adjusted relative risk',
         xlab = 'Relative risk', col = 'grey')
    hist(boot_bias2$t[, 1], main = 'Bias adjusted relative risk',
         xlab = 'Relative risk', col = 'grey', prob = T)
    lines(density(boot_bias2$t[, 1], na.rm=T), col = 'blue')
    abline(v = ci2, col = 'red')

    bootstrap rr=r(rrdx_mie), reps(1000): episensi  24 16 60 100, st(cs) dseca(c(.75)) dspca(c(.9)) dsenc(c(.85)) dspnc(c(.95))
    
    t2 <- misclassification(icu$STA, icu$INF, type = "exposure", bias_parms = c(.75, .85, .9, .95))
}
