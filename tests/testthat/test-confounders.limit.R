context("check limits on confounding")

test_that("non-missing crude RR", {
    expect_error(confounders.limit(OR = 1.65))
})

test_that("enough information", {
    expect_error(confounders.limit(crude.RR = 1.5))
})

test_that("Limits are correct", {
    model <- confounders.limit(OR = 1.65, crude.RR = 1.5, print = FALSE)
    expect_equal(model$conf.limits[1], 0.9091, tolerance = 1e-4, scale = 1)
    expect_equal(model$conf.limits[2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(p = 0.61, RR = 2.01, crude.RR = 1.5, print = FALSE)
    expect_equal(model$conf.limits[1], 1.2060, tolerance = 1e-4, scale = 1)
    expect_equal(model$conf.limits[2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(RR = 2.01, crude.RR = 1.5, print = FALSE)
    expect_equal(model$conf.limits[1], 0.7463, tolerance = 1e-4, scale = 1)
    expect_equal(model$conf.limits[2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(RR = 2.01, OR = 1.65, crude.RR = 1.5, print = FALSE)
    expect_equal(model$conf.limits[1], 0.9091, tolerance = 1e-4, scale = 1)
    expect_equal(model$conf.limits[2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(p = 0.61, crude.RR = 1.5, print = FALSE)
    expect_equal(model$conf.limits[1], 0.915, tolerance = 1e-4, scale = 1)
    expect_equal(model$conf.limits[2], 1.5000, tolerance = 1e-4, scale = 1)
})
