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
                         seexp.parms = list("normal", c(0.4, 0.7, 0.5, 0.01)),
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
                         seexp.parms = list("normal", c(156, 56)),
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
                         seexp.parms = list("normal", c(0.4, 1, 0.8, 0.01)),
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
    expect_equal(model$adj.measures[1, 1], 2.9223, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.8113, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 10.0058, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 3.5498, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.9705, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 47.4763, tolerance = 1e-4, scale = 1)
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
    expect_equal(model$adj.measures[1, 1], 1.5928, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.3478, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 2.0209, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.6946, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.3990, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 2.2433, tolerance = 1e-4, scale = 1)
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
    expect_equal(model$adj.measures[1, 1], 2.2911, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.6630, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 26.7703, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 2.4694, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.7932, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 28.9868, tolerance = 1e-4, scale = 1)
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
    expect_equal(model$adj.measures[1, 1], 4.8406, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 2.3768, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 78.6928, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 5.4640, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 2.5203, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 91.2977, tolerance = 1e-4, scale = 1)
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
    expect_equal(model$adj.measures[1, 1], 1.3576, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.2472, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.5344, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.7747, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.5464, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 2.1224, tolerance = 1e-4, scale = 1)
})


####
test_that("outcome misclassification (ND---normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("normal", c(0.4, 0.7, 0.5, 0.001)),
                      spca.parms = list("normal", c(0.4, 1, 0.9, 0.005)))
    expect_equal(model$adj.measures[1, 1], 40.1652, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 10.1777, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 718.0268, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 46.8348, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 11.8743, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 837.7270, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca.parms = list("normal", c(0.1, 0.9, 0.8, 0.001)),
                      seexp.parms = list("normal", c(0.1, 0.9, 0.82, 0.001)),
                      spca.parms = list("normal", c(0.4, 1, 0.95, 0.001)),
                      spexp.parms = list("normal", c(0.4, 1, 0.99, 0.001)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 1.3287, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.3112, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.3467, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.3787, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.3586, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 1.3994, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND---normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("normal", c(0.4, 0.7, 0.5, 0.001)),
                      spca.parms = list("normal", c(0.4, 1, 0.9, 0.005)))
    expect_equal(model$adj.measures[1, 1], 2.7860, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 2.6963, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 2.8910, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 3.1919, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 3.0671, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 3.3408, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D---normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca.parms = list("normal", c(0.1, 0.9, 0.8, 0.001)),
                      seexp.parms = list("normal", c(0.1, 0.9, 0.82, 0.001)),
                      spca.parms = list("normal", c(0.4, 1, 0.95, 0.001)),
                      spexp.parms = list("normal", c(0.4, 1, 0.99, 0.001)),
                      corr.se = .8, corr.sp = .8)
    expect_equal(model$adj.measures[1, 1], 1.6076, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 1.5988, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.6165, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.7095, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 1.6988, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 1.7204, tolerance = 1e-4, scale = 1)
})

test_that("correct arguments --- list", {
    expect_that(probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         or.parms = list("beta", c(900, 10))),
                throws_error())
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

test_that("Selection bias: adjusted measures are correct (beta)", {
    set.seed(123)
    model <- probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                      reps = 20000,
                      case.exp = list("beta", c(200, 56)),
                      case.nexp = list("beta", c(100, 16)),
                      ncase.exp = list("beta", c(200, 16)),
                      ncase.nexp = list("beta", c(100, 56)))
    expect_equal(model$adj.measures[1, 1], 1.1266, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.9649, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.3281, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 1.1275, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.7966, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 1.6115, tolerance = 1e-4, scale = 1)
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
    expect_equal(model$adj.measures[1, 1], 0.4790, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4532, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.5067, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4789, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3760, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.6094, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.3660, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.3363, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.3999, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.3661, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.2515, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.5321, tolerance = 1e-4, scale = 1)
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
    expect_equal(model$adj.measures[1, 2], 0.4521, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.5050, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4762, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3737, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.6057, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.3630, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.3350, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.3978, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.3634, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.2498, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.5281, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 10000,
                      prev.exp = list("constant", c(.8)),
                      prev.nexp = list("constant", c(.05)),
                      risk = list("constant", c(.63)))
    expect_equal(model$adj.measures[1, 1], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4858, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3857, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.6128, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.3727, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.3727, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.3727, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.3736, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.2611, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.5357, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                      reps = 10000,
                      prev.exp = list("beta", c(200, 56)),
                      prev.nexp = list("beta", c(10, 16)),
                      risk = list("triangular", c(.6, .7, .63)),
                      corr.p = .8)
    expect_equal(model$adj.measures[1, 1], 0.4165, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.3898, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 0.4418, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4164, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.3261, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 0.5277, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 1], 0.2885, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 2], 0.2566, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[3, 3], 0.3179, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 1], 0.2883, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 2], 0.1964, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[4, 3], 0.4176, tolerance = 1e-4, scale = 1)
})

test_that("correct arguments --- list (IRR)", {
    expect_that(probsens.irr(matrix(c(2, 67232, 58, 10539000), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         seca.parms = c("beta", c(900, 10)),
                         spca.parms = list("constant", 1)),
                throws_error())
})

test_that("correct arguments --- list (IRR)", {
    expect_that(probsens.irr(matrix(c(2, 67232, 58, 10539000), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         seca.parms = list("beta", c(-1, 10)),
                         spca.parms = list("constant", 1)),
                throws_error())
})

test_that("correct arguments --- list (IRR)", {
    expect_that(probsens.irr(matrix(c(2, 67232, 58, 10539000), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         seca.parms = list("beta", "uniform", c(900, 10)),
                         spca.parms = list("constant", 1)),
                throws_error())
    })

test_that("IRR exposure misclassification: observed measures are correct", {
    set.seed(123)
    model <- probsens.irr(matrix(c(2, 67232, 58, 10539000), nrow = 2, byrow = TRUE),
                          reps = 10000,
                          seca.parms = list("trapezoidal", c(.4, .45, .55, .6)),
                          spca.parms = list("constant", 1))
    expect_equal(model$obs.measures[1, 1], 5.4054, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.6546, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 19.5268, tolerance = 1e-4, scale = 1)
})

test_that("IRR exposure misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.irr(matrix(c(2, 67232, 58, 10539000), nrow = 2, byrow = TRUE),
                          reps = 10000,
                          seca.parms = list("trapezoidal", c(.4, .45, .55, .6)),
                          spca.parms = list("constant", 1))
    expect_equal(model$adj.measures[1, 1], 5.2587, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 5.2586, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 5.2587, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 5.3569, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.9784, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 28.3129, tolerance = 1e-4, scale = 1)
})

test_that("correct arguments --- list", {
    expect_that(probsens.irr.conf(matrix(c(77, 10000, 87, 10000), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         prev.exp = c("uniform", 1, 2),
                         prev.nexp = list("beta", c(10, 16)),
                         risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
                         corr.p = .8),
                throws_error())
})

test_that("IRR confounding bias: observed measures are correct", {
    set.seed(123)
    model <- probsens.irr.conf(matrix(c(77, 10000, 87, 10000), nrow = 2, byrow = TRUE),
                      reps = 20000,
                      prev.exp = list("trapezoidal", c(.01, .2, .3, .51)),
                      prev.nexp = list("trapezoidal", c(.09, .27, .35, .59)),
                      risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
                      corr.p = .8)
    expect_equal(model$obs.measures[1, 1], 0.8851, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.6979, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 1.1072, tolerance = 1e-4, scale = 1)
})

test_that("IRR Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.irr.conf(matrix(c(77, 10000, 87, 10000), nrow = 2, byrow = TRUE),
                      reps = 20000,
                      prev.exp = list("trapezoidal", c(.01, .2, .3, .51)),
                      prev.nexp = list("trapezoidal", c(.09, .27, .35, .59)),
                      risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
                      corr.p = .8)
    expect_equal(model$adj.measures[1, 1], 0.9653, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7620, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 3], 1.2656, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.9711, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.7063, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 3], 1.3607, tolerance = 1e-4, scale = 1)
})
