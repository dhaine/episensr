context("check bias due to unmeasured and unknown confounders")

test_that("correct number of arguments for bias parameters", {
    expect_that(confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                            type = "RR",
                            bias_parms = c(.63, .8, .05, .5)),
                throws_error())
})

test_that("confounder prevalences between 0 and 1", {
    expect_that(confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                            type = "RR",
                            bias_parms = c(.63, -1, 2)),
                throws_error())
})

test_that("RR is correct", {
    model <- confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         type = "RR",
                         bias_parms = c(.63, .8, .05),
                         print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.3479, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.2757, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.4390, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted RR are correct", {
    model <- confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         type = "RR",
                         bias_parms = c(.63, .8, .05),
                         print = FALSE)
    expect_equal(model$adj.measures[1, 1], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7173, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.4851, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.7173, tolerance = 1e-4, scale = 1)
})

test_that("OR is correct", {
    model <- confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         type = "OR",
                         bias_parms = c(.63, .8, .05),
                         print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.2180, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.1519, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.3128, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted OR are correct", {
    model <- confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         type = "OR",
                         bias_parms = c(.63, .8, .05),
                         print = FALSE)
    expect_equal(model$adj.measures[1, 1], 0.3039, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7173, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.3039, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.7173, tolerance = 1e-4, scale = 1)
})

test_that("RD is correct", {
    model <- confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         type = "RD",
                         bias_parms = c(-.37, .8, .05),
                         print = FALSE)
    expect_equal(model$obs.measures[1, 1], -0.3114, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], -0.3903, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], -0.2325, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted RD are correct", {
    model <- confounders(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         type = "RD",
                         bias_parms = c(-.37, .8, .05),
                         print = FALSE)
    expect_equal(model$adj.measures[1, 1], -0.0339, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], -0.2775, tolerance = 1e-4, scale = 1)
})
