context("check bootstrap bias analysis")

## Test cases
# Misclassification, ci == norm
misclass_eval <- misclassification(matrix(c(215, 1449, 668, 4296),
                                          dimnames = list(c("Breast cancer+",
                                                            "Breast cancer-"),
                                                          c("Smoker+", "Smoker-")),
                                          nrow = 2, byrow = TRUE),
                                   type = "exposure",
                                   bias_parms = c(.78, .78, .99, .99))
set.seed(123)
mod <- boot.bias(misclass_eval)


test_that("RR and OR are correct, misclass + CI norm", {
    expect_equal(mod$ci[1, 1], 0.8376, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[1, 2], 1.1028, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 1], 0.7904, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 2], 1.1370, tolerance = 1e-4, scale = 1)
    })

# Misclassification, ci == perc
set.seed(123)
mod <- boot.bias(misclass_eval, ci_type = "perc")


test_that("RR and OR are correct, misclass + CI perc", {
    expect_equal(mod$ci[1, 1], 0.8390, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[1, 2], 1.1070, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 1], 0.7942, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 2], 1.1469, tolerance = 1e-4, scale = 1)
    })


# Selection, ci == norm
sel <- selection(matrix(c(136, 107, 297, 165),
                        nrow = 2, byrow = TRUE),
                 bias_parms = c(.94, .85, .64, .25))
set.seed(123)
mod <- boot.bias(sel)


test_that("RR and OR are correct, selection + CI norm", {
    expect_equal(mod$ci[1, 1], 1.1356, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[1, 2], 1.9333, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 1], 1.1754, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 2], 2.2681, tolerance = 1e-4, scale = 1)
    })


# Selection, ci == perc
set.seed(123)
mod <- boot.bias(sel, ci_type = "perc")


test_that("RR and OR are correct, selection + CI perc", {
    expect_equal(mod$ci[1, 1], 1.1398, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[1, 2], 1.9491, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 1], 1.1798, tolerance = 1e-4, scale = 1)
    expect_equal(mod$ci[2, 2], 2.2877, tolerance = 1e-4, scale = 1)
    })
