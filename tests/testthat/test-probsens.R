context("check probabilistic sensitivity analysis")

test_that("correct arguments --- list", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = c("uniform", 1, 2),
                         spca = list("beta", c(153, 6))),
                throws_error())
})

test_that("correct arguments --- reps", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 0.5,
                         seca = c("uniform", 1, 2),
                         spca = list("beta", c(153, 6))),
                throws_error())
})

test_that("correct arguments --- dist", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = c("uniform", 1, 2)),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = list(c("uniform", "beta"), c(908, 16)),
                         spca = list("beta", c(153, 6))),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = list("beta", c(908, 16)),
                         seexp = list("beta", "trapezoidal", c(156, 56)),
                         spca = list("beta", c(153, 6)),
                         spexp = list("beta", c(205, 18)),
                         corr_se = .8, corr_sp = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = list("beta", c(-1, 16)),
                         seexp = list("normal", c(0.4, 0.7, 0.5, 0.01)),
                         spca = list("uniform", c(1, .5)),
                         spexp = list("triangular", c(1, 0.5, 0.2)),
                         corr_se = .8, corr_sp = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = list("beta", c(908, 16)),
                         seexp = list("normal", c(156, 56)),
                         spca = list("uniform", c(1, .5)),
                         spexp = list("triangular", c(1, 0.5, 0.2)),
                         corr_se = .8, corr_sp = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = list("beta", c(908, 16)),
                         seexp = list("normal", c(0.4, 1, 0.8, 0.01)),
                         spca = list("uniform", c(.5, 1)),
                         spexp = list("triangular", c(1, 0.5, 0.2)),
                         corr_se = .8, corr_sp = .8),
                throws_error())
    })

test_that("correct arguments --- distribution", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure_pv",
                         reps = 1000,
                         seca = list("beta", c(908, 16)),
                         seexp = list("normal", c(0.4, 1, 0.8, 0.01)),
                         spca = list("uniform", c(.5, 1)),
                         spexp = list("triangular", c(1, 0.5, 0.2))),
                throws_error())
})

test_that("correct arguments --- parameters", {
    expect_that(probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                         type = "exposure",
                         reps = 1000,
                         seca = list("beta", "beta", c(908, 16)),
                         seexp = list("beta", c(156, 56)),
                         spca = list("beta", c(153, 6)),
                         spexp = list("beta", c(205, 18)),
                         corr_se = .8, corr_sp = .8),
                throws_error())
})

test_that("exposure misclassification (ND): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("trapezoidal", c(.75, .85, .95, 1)),
                      spca = list("trapezoidal", c(.75, .85, .95, 1)))
    expect_equal(model$obs_measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("trapezoidal", c(.75, .85, .95, 1)),
                      spca = list("trapezoidal", c(.75, .85, .95, 1)))
    expect_equal(model$adj_measures[2, 1], 2.2431, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.3164, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 6.8389, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 2.5388, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.3584, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 14.6663, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("beta", c(908, 16)),
                      spca = list("beta", c(153, 6)))
    expect_equal(model$adj_measures[2, 1], 1.7520, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.2116, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 2.5463, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 1.8929, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.2356, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 2.9299, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("trapezoidal", c(.75, .85, .95, 1)),
                      seexp = list("trapezoidal", c(.7, .8, .9, .95)),
                      spca = list("trapezoidal", c(.75, .85, .95, 1)),
                      spexp = list("trapezoidal", c(.7, .8, .9, .95)),
                      corr_se = .8, corr_sp = .8)
    expect_equal(model$obs_measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("trapezoidal", c(.75, .85, .95, 1)),
                      seexp = list("trapezoidal", c(.7, .8, .9, .95)),
                      spca = list("trapezoidal", c(.75, .85, .95, 1)),
                      spexp = list("trapezoidal", c(.7, .8, .9, .95)),
                      corr_se = .8, corr_sp = .8)
    expect_equal(model$adj_measures[2, 1], 3.0106, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.4977, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 10.2873, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 3.6676, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.5770, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 48.1868, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("beta", c(908, 16)),
                      seexp = list("beta", c(156, 56)),
                      spca = list("beta", c(153, 6)),
                      spexp = list("beta", c(205, 18)),
                      corr_se = .8, corr_sp = .8)
    expect_equal(model$adj_measures[2, 1], 1.6015, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.0478, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 2.4756, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 1.7046, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.0502, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 2.8417, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: Fox 2023 paper is reproduced", {
    set.seed(1234)
    model <- probsens(matrix(c(40, 20, 60, 80), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 10^6,
                      seca = list("beta", c(25, 3)),
                      spca = list("trapezoidal", c(.9, .93, .97, 1)),
                      seexp = list("beta", c(45, 7)),
                      spexp = list("trapezoidal", c(.8, .83, .87, .9)),
                      corr_se = .8,
                      corr_sp = .8)
    expect_equal(model$obs_measures[1, 1], 2.0000, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.2630, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 3.1671, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 1], 2.7536, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 2.3295, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 5.0172, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 2.8415, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.5355, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 7.7767, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 4.1753, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 3.3127, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 8.5107, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 4.3618, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.8721, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 14.8008, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: Ch.8 is reproduced", {
    set.seed(1234)
    model <- probsens(matrix(c(215, 1449, 668, 4296), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 10^5,
                      seca = list("beta", c(50.6, 14.3)),
                      spca = list("beta", c(70, 1)))
    expect_equal(model$obs_measures[2, 1], 0.9542, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 0.8093, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 1.1252, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 1], 0.9614, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.9450, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 0.9641, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 0.9601, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.8276, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 1.1035, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 0.9490, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 0.9278, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 0.9525, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 0.9472, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 0.7796, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 1.1402, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification: Fox 2023 paper is reproduced", {
    set.seed(1234)
    model <- probsens(matrix(c(40, 20, 60, 80), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 10^5,
                      seca = list("beta", c(254, 24)),
                      spca = list("beta", c(450, 67)),
                      seexp = list("trapezoidal", c(.94, .96, .98, 1)),
                      spexp = list("trapezoidal", c(.9, .92, .93, .95)),
                      corr_se = .8,
                      corr_sp = .8)
    expect_equal(model$obs_measures[1, 1], 2.0000, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.2630, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 3.1671, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 2.6667, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 1.4166, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 5.0199, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 1], 2.4695, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 2.2612, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 2.7380, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 2.4914, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.2725, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 5.6527, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 3.2446, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 2.9278, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 3.6446, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 3.2929, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.3890, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 8.7248, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (PPV/NPV): Ch.8 is reproduced", {
    set.seed(1234)
    model <- probsens(matrix(c(599, 4978, 31175, 391851), nrow = 2, byrow = TRUE),
                      type = "exposure_pv",
                      reps = 10^6,
                      seca = list("beta", c(50, 27)),
                      spca = list("beta", c(120, .5)),
                      seexp = list("beta", c(132, 47)),
                      spexp = list("beta", c(115, 2)))
    expect_equal(model$obs_measures[1, 1], 1.5028, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.3817, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 1.6345, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 1.5125, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 1.3885, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 1.6476, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 1], 1.0712, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.6835, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 1.5144, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca = list("uniform", c(.8, 1)),
                      spca = list("uniform", c(.8, 1)))
    expect_equal(model$obs_measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (ND): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca = list("uniform", c(.8, 1)),
                      spca = list("uniform", c(.8, 1)))
    expect_equal(model$adj_measures[1, 1], 2.2754, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 1.6628, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 18.5378, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 2.3114, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.2928, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 25.9579, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 2.4536, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 1.7930, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 19.9649, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 2.5194, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.3279, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 28.2854, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D): observed measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca = list("uniform", c(.8, 1)),
                      seexp = list("uniform", c(.7, .95)),
                      spca = list("uniform", c(.8, 1)),
                      spexp = list("uniform", c(.7, .95)),
                      corr_se = .8, corr_sp = .8)
    expect_equal(model$obs_measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca = list("uniform", c(.8, 1)),
                      seexp = list("uniform", c(.7, .95)),
                      spca = list("uniform", c(.8, 1)),
                      spexp = list("uniform", c(.7, .95)),
                      corr_se = .8, corr_sp = .8)
    expect_equal(model$adj_measures[1, 1], 4.7529, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 2.3768, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 49.0042, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 4.9584, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.8708, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 76.8991, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 5.3563, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 2.5199, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 56.1394, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 5.6472, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.9727, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 89.0901, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---beta): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(173, 602, 134, 663), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca = list("beta", c(100, 5)),
                      seexp = list("beta", c(110, 10)),
                      spca = list("beta", c(120, 15)),
                      spexp = list("beta", c(130, 30)),
                      corr_se = .8,
                      corr_sp = .8)
    expect_equal(model$adj_measures[1, 1], 1.3576, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 1.2472, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 1.5344, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 1.3601, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.1420, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 1.6539, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 1.7747, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 1.5464, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 2.1224, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 1.7810, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.2739, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 2.5257, tolerance = 1e-4, scale = 1)
})


####
test_that("outcome misclassification (ND---normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca = list("normal", c(0.4, 0.7, 0.5, 0.001)),
                      spca = list("normal", c(0.4, 1, 0.9, 0.005)))
    expect_equal(model$adj_measures[1, 1], 32.0221, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 9.7897, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 241.2247, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 31.3154, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 7.0538, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 396.5734, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 37.3653, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 11.4276, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 281.4657, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 36.5450, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 8.0001, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 464.7237, tolerance = 1e-4, scale = 1)
})

test_that("outcome misclassification (D---logit-normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "outcome",
                      reps = 50000,
                      seca = list("normal", c(0.1, 0.9, 0.8, 0.001)),
                      seexp = list("normal", c(0.1, 0.9, 0.82, 0.001)),
                      spca = list("normal", c(0.4, 1, 0.95, 0.001)),
                      spexp = list("normal", c(0.4, 1, 0.99, 0.001)),
                      corr_se = .8, corr_sp = .8)
    expect_equal(model$adj_measures[1, 1], 1.3287, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 1.3112, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 1.3467, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 1.3267, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.8339, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 2.0168, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 1.3787, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 1.3586, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 1.3994, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 1.3754, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 0.8147, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 2.2280, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (ND---normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("normal", c(0.4, 0.7, 0.5, 0.001)),
                      spca = list("normal", c(0.4, 1, 0.9, 0.005)))
    expect_equal(model$adj_measures[1, 1], 2.7860, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 2.6963, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 2.8910, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 2.7976, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.7019, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 4.6640, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 3.1919, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 3.0671, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 3.3408, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 3.2078, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.8205, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 5.7337, tolerance = 1e-4, scale = 1)
})

test_that("exposure misclassification (D---normal): adjusted measures are correct", {
    set.seed(123)
    model <- probsens(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
                      type = "exposure",
                      reps = 50000,
                      seca = list("normal", c(0.1, 0.9, 0.8, 0.001)),
                      seexp = list("normal", c(0.1, 0.9, 0.82, 0.001)),
                      spca = list("normal", c(0.4, 1, 0.95, 0.001)),
                      spexp = list("normal", c(0.4, 1, 0.99, 0.001)),
                      corr_se = .8, corr_sp = .8)
    expect_equal(model$adj_measures[1, 1], 1.6076, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 1.5988, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 1.6165, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 1.6092, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.0705, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 2.3940, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 1.7095, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 1.6988, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 1.7204, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 1.7110, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.0770, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 2.6977, tolerance = 1e-4, scale = 1)
})

test_that("correct arguments --- list", {
    expect_that(probsens.sel(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
                         reps = 1000,
                         case_exp = list("beta", c(900, 10))),
                throws_error())
})

test_that("Selection bias: adjusted measures are correct (beta)", {
    set.seed(1234)
    model <- probsens.sel(matrix(c(139, 114, 369, 377), nrow = 2, byrow = TRUE),
                      reps = 50000,
                      case_exp = list("beta", c(139, 5.1)),
                      case_nexp = list("beta", c(114, 11.9)),
                      ncase_exp = list("beta", c(369, 96.1)),
                      ncase_nexp = list("beta", c(377, 282.9)))
    expect_equal(model$obs_measures[1, 1], 1.1785, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 0.9511, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 1.4602, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 1.2457, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 0.9357, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 1.6586, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 1], 1.4770, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 1.3543, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 1.6062, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 1.3981, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.2575, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 1.5507, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 1.6247, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 1.4614, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 1.8007, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 1.5100, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.3212, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 1.7204, tolerance = 1e-4, scale = 1)
})

test_that("correct arguments --- list", {
    expect_that(probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                             reps = 1000,
                             prev.exp = c("uniform", 1, 2),
                             prev.nexp = list("beta", c(10, 16)),
                             risk = list("triangular", c(.6, .7, .63)),
                             corr.p = .8),
                throws_error())
})

test_that("correct arguments --- distribution", {
    expect_that(probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                             reps = 1000,
                             prev.exp = list(c("uniform", "beta"), c(10, 16)),
                             prev.nexp = list("beta", c(10, 16)),
                             risk = list("triangular", c(.6, .7, .63)),
                             corr.p = .8),
                throws_error())
})

test_that("Confounding bias: observed measures are correct", {
    set.seed(123)
    model <- probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                          reps = 20000,
                          prev_exp = list("triangular", c(.7, .9, .8)),
                          prev_nexp = list("trapezoidal", c(.03, .04, .05, .06)),
                          risk = list("triangular", c(.6, .7, .63)),
                          corr_p = .8)
    expect_equal(model$obs_measures[1, 1], 0.3479, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 0.2757, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 0.4390, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 0.2180, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 0.1519, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 0.3128, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                          reps = 20000,
                          prev_exp = list("triangular", c(.7, .9, .8)),
                          prev_nexp = list("trapezoidal", c(.03, .04, .05, .06)),
                          risk = list("triangular", c(.6, .7, .63)),
                          corr_p = .8)
    expect_equal(model$adj_measures[1, 1], 0.4792, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.4533, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 0.5068, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 0.4797, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.3366, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 0.6827, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                          reps = 20000,
                          prev_exp = list("triangular", c(.7, .9, .8)),
                          prev_nexp = list("triangular", c(.03, .05, .04)),
                          risk = list("triangular", c(.6, .7, .65)),
                          corr_p = .8)
    expect_equal(model$adj_measures[1, 1], 0.4761, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.4523, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 0.5052, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 0.4756, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.3343, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 0.6790, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                          reps = 10000,
                          prev_exp = list("constant", c(.8)),
                          prev_nexp = list("constant", c(.05)),
                          risk = list("constant", c(.63)))
    expect_equal(model$adj_measures[1, 1], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 0.4858, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.3469, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 0.6779, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                          reps = 10000,
                          prev_exp = list("beta", c(200, 56)),
                          prev_nexp = list("beta", c(10, 16)),
                          risk = list("triangular", c(.6, .7, .63)),
                          corr_p = .8)
    expect_equal(model$adj_measures[1, 1], 0.4165, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.3898, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 0.4418, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 0.4154, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.3188, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 0.5493, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: Fox paper 2023", {
    set.seed(1234)
    model <- probsens_conf(matrix(c(40, 20, 60, 80), nrow = 2, byrow = TRUE),
                          reps = 10^5,
                          prev_exp = list("beta", c(10, 20)),
                          prev_nexp = list("beta", c(5, 20)),
                          risk = list("trapezoidal", c(1.5, 1.7, 2.3, 2.5)))
    expect_equal(model$obs_measures[1, 1], 2.0000, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 1.2630, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 3.1671, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 2.6667, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 1.4166, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 5.0199, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 1], 1.8070, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 1.4907, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 2.1576, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 1.8069, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.0902, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 2.9426, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 2.9709, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 2.0709, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 4.8405, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 3.0076, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 1.6338, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 5.7567, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: book is correct", {
    set.seed(1234)
    model <- probsens_conf(matrix(c(40, 20, 60, 80), nrow = 2, byrow = TRUE),
                          reps = 10^5,
                          prev_exp = list("beta", c(10, 20)),
                          prev_nexp = list("beta", c(5, 20)),
                          risk = list("log-normal", c(log(2), .23)))
    expect_equal(model$adj_measures[1, 1], 1.8188, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 1.4263, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 2.1558, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 1.8037, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 1.0773, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 2.9570, tolerance = 1e-4, scale = 1)
})

test_that("Confounding bias: book is correct", {
    set.seed(1234)
    model <- probsens_conf(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                          reps = 10^5,
                          prev_exp = list("trapezoidal", c(.7, .75, .85, .9)),
                          prev_nexp = list("trapezoidal", c(.03, .04, .07, .1)),
                          risk = list("trapezoidal", c(.5, .6, .7, .8)))
    expect_equal(model$obs_measures[1, 1], 0.3479, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 0.2757, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 0.4390, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 1], 0.2180, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 2], 0.1519, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[2, 3], 0.3128, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 1], 0.4718, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.4194, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 0.5487, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 0.4738, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.3278, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 0.6866, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 1], 0.3785, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 2], 0.3093, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[3, 3], 0.4747, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 1], 0.3798, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 2], 0.2542, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[4, 3], 0.5765, tolerance = 1e-4, scale = 1)
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
                          seca = list("trapezoidal", c(.4, .45, .55, .6)),
                          spca = list("constant", 1))
    expect_equal(model$obs_measures[1, 1], 5.4054, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 0.6546, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 19.5268, tolerance = 1e-4, scale = 1)
})

test_that("IRR exposure misclassification: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.irr(matrix(c(2, 67232, 58, 10539000), nrow = 2, byrow = TRUE),
                          reps = 10000,
                          seca = list("trapezoidal", c(.4, .45, .55, .6)),
                          spca = list("constant", 1))
    expect_equal(model$adj_measures[1, 1], 5.2587, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 5.2586, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 5.2587, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 5.3569, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.9784, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 28.3129, tolerance = 1e-4, scale = 1)
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
                      prev_exp = list("trapezoidal", c(.01, .2, .3, .51)),
                      prev_nexp = list("trapezoidal", c(.09, .27, .35, .59)),
                      risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
                      corr_p = .8)
    expect_equal(model$obs_measures[1, 1], 0.8851, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 2], 0.6979, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs_measures[1, 3], 1.1072, tolerance = 1e-4, scale = 1)
})

test_that("IRR Confounding bias: adjusted measures are correct", {
    set.seed(123)
    model <- probsens.irr.conf(matrix(c(77, 10000, 87, 10000), nrow = 2, byrow = TRUE),
                      reps = 20000,
                      prev_exp = list("trapezoidal", c(.01, .2, .3, .51)),
                      prev_nexp = list("trapezoidal", c(.09, .27, .35, .59)),
                      risk = list("trapezoidal", c(2, 2.5, 3.5, 4.5)),
                      corr_p = .8)
    expect_equal(model$adj_measures[1, 1], 0.9701, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 2], 0.8199, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[1, 3], 1.1804, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 1], 0.9718, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 2], 0.7270, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj_measures[2, 3], 1.3065, tolerance = 1e-4, scale = 1)
})

test_that("Probcase -- Fox paper -- misclassification -- exposure -- RR", {
              a <- 40; b <- 20; c <- 60; d <- 80
              D <- data.frame(e_obs = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
                              d = c(rep(1, a), rep(1, b), rep(0, c), rep(0, d)))
              set.seed(1234)
              model <- probcase(D, x = e_obs, y = d, reps = 10^3,
                                measure = "RR", type = "exposure",
                                seca = list("beta", c(25, 3)),
                                spca = list("trapezoidal", c(.9, .93, .97, 1)),
                                seexp = list("beta", c(45, 7)),
                                spexp = list("trapezoidal", c(.8, .83, .87, .9)),
                                corr_se = .8, corr_sp = .8)
              expect_equal(model$obs_measures[1, 1], 2.0000, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[1, 2], 1.2630, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[1, 3], 3.1671, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[2, 1], 2.6667, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[2, 2], 1.4166, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[2, 3], 5.0199, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[1, 1], 2.7776, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[1, 2], 2.3251, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[1, 3], 4.5794, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2, 1], 2.7853, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2, 2], 1.4233, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2, 3], 9.5569, tolerance = 1e-4, scale = 1)
          })

#test_that("Probcase -- Fox paper -- misclass -- exposure -- RR LONG!", {
#              a <- 40; b <- 20; c <- 60; d <- 80
#              D <- data.frame(e_obs = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
#                              d = c(rep(1, a), rep(1, b), rep(0, c), rep(0, d)))
#              set.seed(1234)
#              model <- probcase(D, x = e_obs, y = d, reps = 10^5,
#                                measure = "RR", type = "exposure",
#                                seca = list("beta", c(25, 3)),
#                                spca = list("trapezoidal", c(.9, .93, .97, 1)),
#                                seexp = list("beta", c(45, 7)),
#                                spexp = list("trapezoidal", c(.8, .83, .87, .9)),
#                                corr_se = .8, corr_sp = .8)
#              expect_equal(model$adj_measures[1, 1], 2.7533, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 2], 2.3317, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 3], 5.0717, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 1], 2.8427, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 2], 1.5391, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 3], 8.0743, tolerance = 1e-4, scale = 1)
#          })

#test_that("Probcase -- Ch.9 QBA book -- misclassification -- exposure -- OR", {
#              a <- 215; b <- 1449; c <- 668; d <- 4296
#              D <- data.frame(e_obs = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
#                              d = c(rep(1, a), rep(1, b), rep(0, c), rep(0, d)))
#              set.seed(1234)
#              model <- probcase(D, x = e_obs, y = d, reps = 10^3,
#                                measure = "OR", type = "exposure",
#                                seca = list("beta", c(50.6, 14.3)),
#                                spca = list("beta", c(70, 1)))
#              expect_equal(model$obs_measures[1, 1], .9654, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[1, 2], .8524, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[1, 3], 1.0934, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[2, 1], .9542, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[2, 2], .8093, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[2, 3], 1.1252, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 1], 0.9491, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 2], 0.9298, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 3], 0.9526, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 1], 0.9546, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 2], 0.7867, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 3], 1.1499, tolerance = 1e-4, scale = 1)
#          })

#test_that("Probcase -- Ch.9 QBA book -- misclassification -- exposure -- OR LONG!", {
#              a <- 215; b <- 1449; c <- 668; d <- 4296
#              D <- data.frame(e_obs = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
#                              d = c(rep(1, a), rep(1, b), rep(0, c), rep(0, d)))
#              set.seed(1234)
#              system.time(model <- probcase(D, x = e_obs, y = d, reps = 10^5,
#                                            measure = "OR", type = "exposure",
#                                            seca = list("beta", c(50.6, 14.3)),
#                                            spca = list("beta", c(70, 1))))
#              expect_equal(model$obs_measures[1, 1], .9654, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[1, 2], .8524, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[1, 3], 1.0934, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[2, 1], .9542, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[2, 2], .8093, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_measures[2, 3], 1.1252, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 1], 0.9491, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 2], 0.9298, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[1, 3], 0.9526, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 1], 0.9546, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 2], 0.7867, tolerance = 1e-4, scale = 1)
#              expect_equal(model$adj_measures[2, 3], 1.1499, tolerance = 1e-4, scale = 1)
#          })

test_that("Probcase -- Fox paper -- misclassification -- outcome", {
              a <- 40; b <- 20; c <- 60; d <- 80
              D <- data.frame(e = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
                              d_obs = c(rep(1, a), rep(1, b), rep(0, c), rep(0, d)))
              set.seed(1234)
              model <- probcase(D, x = e, y = d_obs, reps = 10^3,
                                measure = "RR", type = "outcome",
                                seca = list("beta", c(254, 24)),
                                spca = list("trapezoidal", c(.94, .96, .98, 1)),
                                seexp = list("beta", c(450, 67)),
                                spexp = list("trapezoidal", c(.9, .92, .93, .95)),
                                corr_se = .8, corr_sp = .8)
              expect_equal(model$obs_measures[1, 1], 2.0000, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[1, 2], 1.2630, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[1, 3], 3.1671, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[2, 1], 2.6667, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[2, 2], 1.4166, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_measures[2, 3], 5.0199, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[1, 1], 2.6875, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[1, 2], 1.8095, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[1, 3], 4.3333, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2, 1], 2.6146, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2, 2], 1.4877, tolerance = 1e-4, scale = 1)
              expect_equal(model$adj_measures[2, 3], 5.6669, tolerance = 1e-4, scale = 1)
          })
