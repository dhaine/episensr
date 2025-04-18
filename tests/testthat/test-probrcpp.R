context("check probabilistic sensitivity analysis when using Rcpp")

test_that("rbinom procedure is correct", {
              nrow <- 200
              reps <- 10^3
              set.seed(1234)
              p <- matrix(rbeta(nrow, 1, 1), nrow = nrow, ncol = reps)
              mat1 <-matrix(NA, nrow = nrow, ncol = reps)
              mat2 <-matrix(NA, nrow = nrow, ncol = reps)
              set.seed(1234)
              for (i in 1:reps) {
                  mat1[, i] <- suppressWarnings({
                                               rbinom(nrow, 1, p[, i])})
              }
              set.seed(1234)
              for (i in 1:reps) {
                  mat2[, i] <- cpprbinom(nrow, 1, p[, i])
              }
              expect_equal(mat1, mat2)
          })

test_that("matrix loop is correct", {
              a <- 40; b <- 20; c <- 60; d <- 80
              D <- data.frame(e_obs = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
                              d = c(rep(1, a), rep(1, b), rep(0, c), rep(0, d)))

              set.seed(1234)
              mat1 <- test_fun(D, x = e_obs, y = d, type = "exposure", reps = 3,
#                               measure = "RR", type = "exposure",
                               seca = list("beta", c(25, 3)),
                               spca = list("trapezoidal", c(.9, .93, .97, 1)),
                               seexp = list("beta", c(45, 7)),
                               spexp = list("trapezoidal", c(.8, .83, .87, .9)),
                               corr_se = .8, corr_sp = .8)
              draws <- mat1$draws
              obs_mat <- mat1$obs_mat
              obs_mat2 <- mat1$obs_mat[, c(1:2, 14:17)]

#              set.seed(1234)
#              mat1b <- test_fun2(obs_mat, draws, reps = 3)

              set.seed(1234)
              mat2 <- RRexpo_e(3, obs_mat2, draws)

#              set.seed(1234)
#              mat3 <- probrcpp(D, x = e_obs, y = d, reps = 3,
#                               measure = "RR", type = "exposure",
#                               seca = list("beta", c(25, 3)),
#                               spca = list("trapezoidal", c(.9, .93, .97, 1)),
#                               seexp = list("beta", c(45, 7)),
#                               spexp = list("trapezoidal", c(.8, .83, .87, .9)),
#                               corr_se = .8, corr_sp = .8)

#              expect_equal(mat1$ptest, mat2$ptest)
#              expect_equal(mat1b$ptest, mat2$ptest)
#              expect_equal(mat1$ptest, mat3$ptest)
              expect_equal(mat1$etest, mat2)
              expect_equal(mat1$etest, mat1b$etest)
              expect_equal(mat1b$etest, mat2$etest)
              expect_equal(mat1$etest, mat3$etest)
              expect_equal(mat1b$etest, mat3$etest)
              expect_equal(mat2$etest, mat3$etest)
          })
