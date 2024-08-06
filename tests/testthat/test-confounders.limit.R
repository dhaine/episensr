context("check limits on confounding")

test_that("non-missing crude RR", {
    expect_error(confounders.limit(OR = 1.65))
})

test_that("enough information", {
    expect_error(confounders.limit(crude_RR = 1.5))
})

test_that("Limits are correct", {
    model <- confounders.limit(OR = 1.65, crude_RR = 1.5)
    expect_equal(model[[3]][1], 0.9091, tolerance = 1e-4, scale = 1)
    expect_equal(model[[3]][2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(p = 0.61, RR = 2.01, crude_RR = 1.5)
    expect_equal(model[[3]][1], 1.2060, tolerance = 1e-4, scale = 1)
    expect_equal(model[[3]][2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(RR = 2.01, crude_RR = 1.5)
    expect_equal(model[[3]][1], 0.7463, tolerance = 1e-4, scale = 1)
    expect_equal(model[[3]][2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(RR = 2.01, OR = 1.65, crude_RR = 1.5)
    expect_equal(model[[3]][1], 0.9091, tolerance = 1e-4, scale = 1)
    expect_equal(model[[3]][2], 1.5000, tolerance = 1e-4, scale = 1)
})

test_that("Limits are correct", {
    model <- confounders.limit(p = 0.61, crude_RR = 1.5)
    expect_equal(model[[3]][1], 0.915, tolerance = 1e-4, scale = 1)
    expect_equal(model[[3]][2], 1.5000, tolerance = 1e-4, scale = 1)
})
