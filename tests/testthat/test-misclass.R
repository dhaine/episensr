context("check misclassification bias")

test_that("correct number of arguments for bias parameters", {
    expect_that(misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  type = "exposure",
                                  bias_parms = c(.99, .99)),
                throws_error())
})

test_that("bias parameters between 0 and 1", {
    expect_that(misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  type = "exposure",
                                  bias_parms = c(-1, .78, .99, 2)),
                throws_error())
})

test_that("Observed measures are correct for exposure misclassification", {
    model <- misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  type = "exposure",
                                  bias_parms = c(.78, .78, .99, .99))
    expect_equal(model$obs.measures[1, 1], 0.9654, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.8524, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 1.0933, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 0.9542, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 0.8093, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 1.1252, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for exposure misclassification", {
    model <- misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  type = "exposure",
                                  bias_parms = c(.78, .78, .99, .99))
    expect_equal(model$adj.measures[1], 0.9614, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2], 0.9491, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.7896, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 1.1408, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for exposure misclassification", {
    model <- misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  type = "exposure",
                                  bias_parms = c(.13, .78, .91, 1))
    expect_equal(model$adj.measures[1], 82.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2], 237.0529, tolerance = 1e-4, scale = 1)
})

test_that("Observed measures are correct for outcome misclassification", {
    model <- misclassification(matrix(c(4558, 3428, 46305, 46085),
                                         nrow = 2, byrow = TRUE),
                                  type = "outcome",
                                  bias_parms = c(.53, .53, .99, .99))
    expect_equal(model$obs.measures[1, 1], 1.2944, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.2404, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 1.3506, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.3233, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2636, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 1.3858, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for outcome misclassification", {
    model <- misclassification(matrix(c(4558, 3428, 46305, 46085),
                                         nrow = 2, byrow = TRUE),
                                  type = "outcome",
                                  bias_parms = c(.53, .53, .99, .99))
    expect_equal(model$adj.measures[1], 1.3440, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2], 1.4062, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for outcome misclassification", {
    model <- misclassification(matrix(c(4558, 3428, 46305, 46085),
                                         nrow = 2, byrow = TRUE),
                                  type = "outcome",
                                  bias_parms = c(.78, .75, .92, .97))
    expect_equal(model$adj.measures[1], 0.2520, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2], 0.2416, tolerance = 1e-4, scale = 1)
})

test_that("CI is correct for OR when exposure misclassification", {
    model <- misclassification(matrix(c(126, 92, 71, 224),
                                      nrow = 2, byrow = TRUE),
                               type = "exposure",
                               bias_parms = c(.94, .94, .97, .97))
    expect_equal(model$adj.measures[2, 2], 3.2825, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 7.6909, tolerance = 1e-4, scale = 1)
})

test_that("correct number of arguments for bias parameters", {
    expect_that(misclassification_cov(array(c(1319, 38054, 5641, 405546, 565,
                                              3583, 781, 21958, 754, 34471,
                                              4860, 383588),
                                            dimnames = list(c("Twins+", "Twins-"),
                                                            c("Folic acid+", "Folic acid-"),
                                                            c("Total", "IVF+", "IVF-")),
                                            dim = c(2, 2, 3)),
                                      bias_parms = c(.6, .6, .95)),
                throws_error())
})

test_that("bias parameters between 0 and 1", {
    expect_that(misclassification_cov(array(c(1319, 38054, 5641, 405546, 565,
                                              3583, 781, 21958, 754, 34471,
                                              4860, 383588),
                                            dimnames = list(c("Twins+", "Twins-"),
                                                            c("Folic acid+", "Folic acid-"),
                                                            c("Total", "IVF+", "IVF-")),
                                            dim = c(2, 2, 3)),
                                      bias_parms = c(-1, .6, .95, 2)),
                throws_error())
})

test_that("Observed measures are correct for covariate misclassification", {
    model <- misclassification_cov(array(c(1319, 38054, 5641, 405546, 565,
                                              3583, 781, 21958, 754, 34471,
                                              4860, 383588),
                                            dimnames = list(c("Twins+", "Twins-"),
                                                            c("Folic acid+", "Folic acid-"),
                                                            c("Total", "IVF+", "IVF-")),
                                            dim = c(2, 2, 3)),
                                      bias_parms = c(.6, .6, .95, .95))
    expect_equal(model$obs.measures[1, 1], 2.4419, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 2.3019, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.5904, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 2.4919, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 2.3448, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.6483, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for covariate misclassification", {
    model <- misclassification_cov(array(c(1319, 38054, 5641, 405546, 565,
                                              3583, 781, 21958, 754, 34471,
                                              4860, 383588),
                                            dimnames = list(c("Twins+", "Twins-"),
                                                            c("Folic acid+", "Folic acid-"),
                                                            c("Total", "IVF+", "IVF-")),
                                            dim = c(2, 2, 3)),
                                      bias_parms = c(.6, .6, .95, .95))
    expect_equal(model$adj.measures[1, 1], 2.2617, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.0002, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.0797, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.4413, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 2.2288, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 1.0002, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 1.0956, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 2.4415, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[5, 1], 2.3379, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[5, 2], 1.0003, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[6, 1], 1.0659, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[6, 2], 2.4911, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[7, 1], 2.2905, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[7, 2], 1.0002, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[8, 1], 1.0879, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[8, 2], 2.4914, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for covariate misclassification", {
    model <- misclassification_cov(array(c(1319, 38054, 5641, 405546, 565,
                                              3583, 781, 21958, 754, 34471,
                                              4860, 383588),
                                            dimnames = list(c("Twins+", "Twins-"),
                                                            c("Folic acid+", "Folic acid-"),
                                                            c("Total", "IVF+", "IVF-")),
                                            dim = c(2, 2, 3)),
                                      bias_parms = c(.78, .78, .96, .96))
    expect_equal(model$adj.measures[1, 2], 1.7755, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.3753, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 1.7259, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 1.4149, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[5, 2], 1.8756, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[6, 2], 1.3286, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[7, 2], 1.7903, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[8, 2], 1.3919, tolerance = 1e-4, scale = 1)
})
