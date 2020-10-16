context("check bias due to unmeasured confounders based on external adjustment")

test_that("correct if null bias parameters", {
    expect_output(confounders.ext(RR = 1))
})

test_that("correct number of arguments for bias parameters", {
    expect_error(confounders.ext(RR = 1,
                                 bias_parms = c(0.1, 0.9, 0.1)))
})

test_that("prevalence of confounder between 0 and 1", {
    expect_error(confounders.ext(RR = 1,
                                 bias_parms = c(0.1, 0.9, 2, 0.4)))
})

test_that("prevalence of exposure between 0 and 1", {
    expect_error(confounders.ext(RR = 1,
                                 bias_parms = c(0.1, 0.9, 0.1, 2)))
})

test_that("association between confounder and outcome >= 0", {
    expect_error(confounders.ext(RR = 1,
                                 bias_parms = c(-1, 0.9, 0.1, 0.4)))
})

test_that("association between exposure and confounder >= 0", {
    expect_error(confounders.ext(RR = 1,
                                 bias_parms = c(0.1, -1, 0.1, 0.4)))
})

test_that("true RR > 0", {
    expect_error(confounders.ext(RR = -1,
                                 bias_parms = c(0.1, 0.9, 0.1, 0.4)))
})

test_that("results ok!", {
    model <- confounders.ext(RR = 1,
                             bias_parms = c(0.1, 0.9, 0.1, 0.4),
                             print = FALSE)
    expect_equal(as.numeric(model[6]), 1.0093, tolerance = 1e-4)
    expect_equal(as.numeric(model[7]), 0.9328, tolerance = 1e-4)
})
