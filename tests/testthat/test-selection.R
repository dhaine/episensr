context("check selection bias")

test_that("correct number of arguments for selection probabilities", {
    expect_that(selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          bias_parms = c(.94, .85, .64)),
                throws_error())
})

test_that("selection probabilities between 0 and 1", {
    expect_that(selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          bias_parms = c(.94, -1, .64, 3)),
                throws_error())
})

test_that("RR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                              nrow = 2, byrow = TRUE),
                          bias_parms = c(.94, .85, .64, .25))
    expect_equal(model$obs_measures[1, 1], 0.7984, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 0.6518, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 0.9780, tolerance = 1e-4, scale = 1)
})

test_that("OR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                              nrow = 2, byrow = TRUE),
                          bias_parms = c(.94, .85, .64, .25))
    expect_equal(model$obs_measures[2, 1], 0.7061, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 0.5144, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 0.9693, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted RR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                              nrow = 2, byrow = TRUE),
                          bias_parms = c(.94, .85, .64, .25))
    expect_equal(model$adj_measures[1], 1.4838, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted OR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                              nrow = 2, byrow = TRUE),
                          bias_parms = c(.94, .85, .64, .25))
    expect_equal(model$adj_measures[2], 1.6346, tolerance = 1e-4, scale = 1)
})

test_that("bias-factor works (1)", {
              model <- selection(matrix(c(136, 107, 297, 165),
                                        nrow = 2, byrow = TRUE),
                                 bias_parms = 0.43)
              expect_equal(model$adj_measures[1], 1.8568, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2], 1.6421, tolerance = 1e-4, scale = 1)
          })

test_that("bias-factor works (2)", {
              model <- selection(matrix(c(136, 107, 297, 165),
                                        nrow = 2, byrow = TRUE),
                                 bias_parms = 0.99)
              expect_equal(model$adj_measures[1], 0.8065, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2], 0.7132, tolerance = 1e-4, scale = 1)
          })

test_that("correct number of arguments for M-bias", {
    expect_that(mbias(or = c(2, 5.4, 2.5)),
                throws_error())
})

test_that("M-bias OR positive", {
    expect_that(mbias(or = c(-2, 5.4, 2.5, 1.5, 1)),
                throws_error())
})

test_that("M-bias OR is correct", {
    model <- mbias(or = c(2, 5.4, 2.5, 1.5, 1))
    expect_equal(model$adj_measures, 0.9938, tolerance = 1e-4, scale = 1)
    expect_equal(model$mbias_parms[1], 1.3149, tolerance = 1e-4, scale = 1)
    expect_equal(model$mbias_parms[2], 1.0952, tolerance = 1e-4, scale = 1)
    expect_equal(model$mbias_parms[3], 1.0062, tolerance = 1e-4, scale = 1)
})
