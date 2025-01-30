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
