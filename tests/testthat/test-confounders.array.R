context("check bias due to unmeasured confounders based on confounding imbalance")

test_that("correct number of arguments for type", {
    expect_that(confounders.array(crude.risk = 1.5,
                            type = c("binary", "continuous"),
                            bias_parms = c(5.5, 0.5)),
                throws_error())
})

test_that("correct number of arguments for bias parameters", {
    expect_that(confounders.array(crude.risk = 1.5,
                            type = "binary",
                            bias_parms = c(5.5, 0.5)),
                throws_error())
})

test_that("correct if null bias parameters", {
    expect_output(confounders.array(crude.risk = 1.5,
                                    type = "binary"))
})

test_that("confounder prevalence among exposed between 0 and 1", {
    expect_that(confounders.array(crude.risk = 1.5,
                            type = "binary",
                            bias_parms = c(5.5, -1, 0.1)),
                throws_error())
})

test_that("confounder prevalence among unexposed between 0 and 1", {
    expect_that(confounders.array(crude.risk = 1.5,
                            type = "binary",
                            bias_parms = c(5.5, 0.5, 2)),
                throws_error())
})

test_that("association between confounder and outcome >= 0", {
    expect_that(confounders.array(crude.risk = 1.5,
                            type = "continuous",
                            bias_parms = c(-1, 7.8, 7.9)),
                throws_error())
})

test_that("crude risk > 0 (binary)", {
    expect_that(confounders.array(crude.risk = -1,
                            type = "binary",
                            bias_parms = c(5.5, 0.5, 0.1)),
                throws_error())
})

test_that("crude risk > 0 (continuous)", {
    expect_that(confounders.array(crude.risk = -1,
                            type = "continuous",
                            bias_parms = c(1.009, 7.8, 7.9)),
                throws_error())
})

test_that("crude risk between -1 and 1 (RD)", {
    expect_error(confounders.array(crude.risk = -2,
                            type = "RD",
                            bias_parms = c(0.009, 8.5, 8)))
    })

test_that("RR and percent are correct (binary)", {
    model <- confounders.array(crude.risk = 1.5,
                         type = "binary",
                         bias_parms = c(5.5, 0.5, 0.1),
                         print = FALSE)
    expect_equal(as.numeric(model[5]), 0.6692, tolerance = 1e-4)
    expect_equal(as.numeric(model[6]), 124.1379, tolerance = 1e-4)
})

test_that("RR and percent are correct (continuous)", {
    model <- confounders.array(crude.risk = 1.5,
                         type = "continuous",
                         bias_parms = c(1.009, 7.8, 7.9),
                         print = FALSE)
    expect_equal(as.numeric(model[5]), 1.5013, tolerance = 1e-4)
    expect_equal(as.numeric(model[6]), -0.0895, tolerance = 1e-4)
})

test_that("RR and percent are correct (RD)", {
    model <- confounders.array(crude.risk = 0.05,
                         type = "RD",
                         bias_parms = c(0.009, 8.5, 8),
                         print = FALSE)
    expect_equal(as.numeric(model[5]), 0.0455, tolerance = 1e-4)
    expect_equal(as.numeric(model[6]), 9.89011, tolerance = 1e-4)
})
