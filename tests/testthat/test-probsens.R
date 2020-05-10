context("check probabilistic sensitivity analysis")

test_that("correct arguments --- list", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca.parms = c("uniform", 1, 2),
                         spca.parms = list("beta", c(153, 6))),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca.parms = list(c("uniform", "beta"), c(908, 16)),
                         spca.parms = list("beta", c(153, 6))),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca.parms = list("beta", c(908, 16)),
                         seexp.parms = list("beta", "trapezoidal", c(156, 56)),
                         spca.parms = list("beta", c(153, 6)),
                         spexp.parms = list("beta", c(205, 18)),
                         corr.se = .8, corr.sp = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca.parms = list("beta", c(-1, 16)),
                         seexp.parms = list("logit-normal", c(156, 56)),
                         spca.parms = list("uniform", c(1, .5)),
                         spexp.parms = list("triangular", c(1, 0.5, 0.2)),
                         corr.se = .8, corr.sp = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca.parms = list("beta", c(908, 16)),
                         seexp.parms = list("logit-normal", c(156, 56)),
                         spca.parms = list("uniform", c(1, .5)),
                         spexp.parms = list("triangular", c(1, 0.5, 0.2)),
                         corr.se = .8, corr.sp = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca.parms = list("beta", c(908, 16)),
                         seexp.parms = list("logit-normal", c(156, 56)),
                         spca.parms = list("uniform", c(.5, 1)),
                         spexp.parms = list("triangular", c(1, 0.5, 0.2)),
                         corr.se = .8, corr.sp = .8),
                throws_error())
})

test_that("correct arguments --- parameters", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca.parms = list("beta", "beta", c(908, 16)),
                         seexp.parms = list("beta", c(156, 56)),
                         spca.parms = list("beta", c(153, 6)),
                         spexp.parms = list("beta", c(205, 18)),
                         corr.se = .8, corr.sp = .8),
                throws_error())
})

test_that("exposure misclassification (ND): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)))
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)))
    expect_equal(model$adj.measures[1, 1], 2.1765, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.7378, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 6.4922, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.4549, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.8723, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 13.8513, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("beta", c(908, 16)),
                      spca.parms = list("beta", c(153, 6)))
    expect_equal(model$adj.measures[1, 1], 1.7479, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.6857, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.8909, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.8878, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.8086, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 2.0740, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      seexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      seexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      spca.parms = list("trapezoidal", c(.75, .85, .95, 1)),
                      spexp.parms = list("trapezoidal", c(.7, .8, .9, .95)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 2.9161, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.6866, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 10.2742, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 3.5341, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.8127, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 54.5525, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("beta", c(908, 16)),
                      seexp.parms = list("beta", c(156, 56)),
                      spca.parms = list("beta", c(153, 6)),
                      spexp.parms = list("beta", c(205, 18)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 1.6084, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.4714, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.9361, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.7140, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.5468, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 2.1316, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      spca.parms = list("uniform", c(.8, 1)))
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      spca.parms = list("uniform", c(.8, 1)))
    expect_equal(model$adj.measures[1, 1], 2.2803, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.6629, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 20.5547, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.4596, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.7932, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 22.1893, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      seexp.parms = list("uniform", c(.7, .95)),
                      spca.parms = list("uniform", c(.8, 1)),
                      spexp.parms = list("uniform", c(.7, .95)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("uniform", c(.8, 1)),
                      seexp.parms = list("uniform", c(.7, .95)),
                      spca.parms = list("uniform", c(.8, 1)),
                      spexp.parms = list("uniform", c(.7, .95)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 4.8303, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 2.1173, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 50.2278, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 5.4494, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.2103, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 57.4525, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, .8)),
                      spca.parms = list("logit-logistic", c(0, .8)))
    expect_equal(model$adj.measures[1, 1], 0.9354, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4488, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 16.8135, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.5596, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.0387, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 23.7620, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, 0.8)),
                      seexp.parms = list("logit-logistic", c(0, .5)),
                      spca.parms = list("logit-logistic", c(0, 0.8)),
                      spexp.parms = list("logit-logistic", c(0, .5)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 0.9916, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.1969, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 37.2372, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.9416, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.0100, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 86.3129, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---beta): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(173, 602, 134, 663), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("beta", c(100, 5)),
                      seexp.parms = list("beta", c(110, 10)),
                      spca.parms = list("beta", c(120, 15)),
                      spexp.parms = list("beta", c(130, 30)),
                      corr.se = .8,
                      corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 1.3736, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.3138, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.5105, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.8097, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.6795, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 2.0774, tolerance = 1e-4, scale = 1)
})


test_that("exposure misclassification (ND---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, .8)),
                      spca.parms = list("logit-logistic", c(0, .8)))
    expect_equal(model$adj.measures[1, 1], 0.5855, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.0771, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 12.3942, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.5447, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.0269, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 31.4069, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D---logit-logistic): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-logistic", c(0, 0.8)),
                      seexp.parms = list("logit-logistic", c(0, .5)),
                      spca.parms = list("logit-logistic", c(0, 0.8)),
                      spexp.parms = list("logit-logistic", c(0, .5)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 0.9410, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.0460, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 22.4035, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.9344, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.0108, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 90.9576, tolerance = 1e-4, scale = 1)
})
####
test_that("outcome misclassification (ND---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, 0.28)),
                      spca.parms = list("logit-normal", c(2.159, .28)))
    expect_equal(model$adj.measures[1, 1], 5.9926, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 2.6443, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 48.1013, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 6.4646, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.8488, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 51.8446, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, 0.28)),
                      seexp.parms = list("logit-normal", c(0, .5)),
                      spca.parms = list("logit-normal", c(2.159, 0.28)),
                      spexp.parms = list("logit-normal", c(0, .5)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 3.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 3.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 3.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 3.6503, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 3.6503, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 3.6503, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, .28)),
                      spca.parms = list("logit-normal", c(2.159, .28)))
    expect_equal(model$adj.measures[1, 1], 2.1230, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.8816, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 3.1177, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.3822, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.0578, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 3.9231, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("logit-normal", c(2.159, 0.28)),
                      seexp.parms = list("logit-normal", c(0, .5)),
                      spca.parms = list("logit-normal", c(2.159, 0.28)),
                      spexp.parms = list("logit-normal", c(0, .5)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 0.2808, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.0349, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 10.9991, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.2389, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.0032, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 55.4015, tolerance = 1e-4, scale = 1)
})

test_that("Selection bias: observed measures are correct", {
    set.seed(123)
    model <- probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      or.parms = list("triangular", c(.35, 1.1, .43)))
    expect_equal(model$obs.measures[1, 1], 0.7061, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.5144, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.9693, tolerance = 1e-4, scale = 1)
})

test_that("Selection bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      or.parms = list("triangular", c(.35, 1.1, .43)))
    expect_equal(model$adj.measures[1, 1], 1.1850, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7160, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.8153, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.1793, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.6456, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 2.0608, tolerance = 1e-4, scale = 1)
})

test_that("Selection bias: adjusted measures are correct (logit)", {
    set.seed(123)
    model <- probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      or.parms = list("log-logistic", c(15, 20)))
    expect_equal(model$adj.measures[1, 1], 14.1010, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 11.0719, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 17.9912, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 14.0991, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 9.5216, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 20.9828, tolerance = 1e-4, scale = 1)
})

test_that("correct arguments --- list", {
    expect_that(probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         prev.exp = c("uniform", 1, 2),
                         prev.nexp = list("beta", c(10, 16)),
                         risk = list("triangular", c(.6, .7, .63)),
                         corr.p = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         prev.exp = list(c("uniform", "beta"), c(10, 16)),
                         prev.nexp = list("beta", c(10, 16)),
                         risk = list("triangular", c(.6, .7, .63)),
                         corr.p = .8),
                throws_error())
})

test_that("Confounding bias: observed measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 20000,
                      prev.exp = list("triangular", c(.7, .9, .8)),
                      prev.nexp = list("trapezoidal", c(.03, .04, .05, .06)),
                      risk = list("triangular", c(.6, .7, .63)),
                      corr.p = .8)
    expect_equal(model$obs.measures[1, 1], 0.3479, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.2757, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.4390, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 0.2180, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 0.1519, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 0.3128, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 20000,
                      prev.exp = list("triangular", c(.7, .9, .8)),
                      prev.nexp = list("trapezoidal", c(.03, .04, .05, .06)),
                      risk = list("triangular", c(.6, .7, .63)),
                      corr.p = .8)
    expect_equal(model$adj.measures[1, 1], 0.4788, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4530, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.5080, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4782, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3768, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.6074, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.3658, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.3361, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.4020, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.3652, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.2524, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.5292, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 20000,
                      prev.exp = list("triangular", c(.7, .9, .8)),
                      prev.nexp = list("triangular", c(.03, .05, .04)),
                      risk = list("triangular", c(.6, .7, .65)),
                      corr.p = .8)
    expect_equal(model$adj.measures[1, 1], 0.4759, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4518, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.5060, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4755, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3751, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.6041, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.3628, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.3347, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.3991, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.3626, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.2507, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.5250, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 10000,
                      prev.exp = list("constant", c(.8)),
                      prev.nexp = list("constant", c(.05)),
                      risk = list("constant", c(.63)))
    expect_equal(model$adj.measures[1, 1], 0.4850, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4850, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.4850, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4858, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3857, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.6127, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.3726, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.3726, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.3726, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.3736, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.2611, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.5356, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 10000,
                      prev.exp = list("beta", c(200, 56)),
                      prev.nexp = list("beta", c(10, 16)),
                      risk = list("triangular", c(.6, .7, .63)),
                      corr.p = .8)
    expect_equal(model$adj.measures[1, 1], 0.4226, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4079, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.4419, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4231, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3346, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.5328, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.2953, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.2811, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.3178, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.2965, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.2059, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.4239, tolerance = 1e-4, scale = 1)
})
