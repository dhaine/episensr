context("check E-value")

test_that("association greater than 0", {
    expect_that(confounders.evalue(est = -1,
                                   type = "RR"),
                throws_error())
})

test_that("Closest CIs, E-value are correct", {
    model <- confounders.evalue(est = 10.73, lower_ci = 8.02, upper_ci = 14.36,
                                type = "RR")
    expect_equal(model[1, 2], 8.02, tolerance = 1e-4, scale = 1)
    expect_equal(model[2, 1], 20.947, tolerance = 1e-3, scale = 1)
    expect_equal(model[2, 2], 15.523, tolerance = 1e-3, scale = 1)
})

test_that("Closest CIs, E-value are correct (non-null H_0)", {
    model <- confounders.evalue(est = 0.8, lower_ci = 0.71, upper_ci = 0.91,
                                type = "RR", true_est = 1.2)
    expect_equal(model[1, 2], 0.91, tolerance = 1e-4, scale = 1)
    expect_equal(model[2, 1], 2.366, tolerance = 1e-3, scale = 1)
    expect_equal(model[2, 2], 1.967, tolerance = 1e-3, scale = 1)
})

test_that("Closest CIs, E-value are correct", {
    model <- confounders.evalue(est = 1.47, lower_ci = 1.12, upper_ci = 1.93,
                                type = "ORc")
    expect_equal(model[1, 1], 1.212, tolerance = 1e-3, scale = 1)
    expect_equal(model[1, 2], 1.058, tolerance = 1e-3, scale = 1)
    expect_equal(model[2, 1], 1.720, tolerance = 1e-3, scale = 1)
    expect_equal(model[2, 2], 1.307, tolerance = 1e-3, scale = 1)
})

test_that("Closest CIs, E-value are correct", {
    model <- confounders.evalue(est = -0.42, sd = 0.14, type = "diff_RR")
    expect_equal(model[1, 1], 0.682, tolerance = 1e-3, scale = 1)
    expect_equal(model[1, 2], 0.875, tolerance = 1e-3, scale = 1)
    expect_equal(model[2, 1], 2.292, tolerance = 1e-3, scale = 1)
    expect_equal(model[2, 2], 1.545, tolerance = 1e-3, scale = 1)
})
