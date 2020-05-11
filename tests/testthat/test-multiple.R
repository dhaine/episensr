context("check multiple bias analysis")

test_that("Correct functions are selected", {
              dat <- matrix(c(118, 832, 103, 884),
                            dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
                            nrow = 2, byrow = TRUE)
              expect_that(dat %>%
                          misclassification(., type = "exposure",
                                            bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "confounders.ext" ,
                                        RR = 1, bias_parms = c(0.1, 0.9, 0.1, 0.4)),
                          throws_error())
})

test_that("Correct table passed on", {
              dat <- matrix(c(118, 832, 103, 884),
                            dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
                            nrow = 2, byrow = TRUE)
              model <- dat %>%
                  misclassification(., type = "exposure",
                                    bias_parms = c(.56, .58, .99, .97)) %>%
                  multiple.bias(., bias_function = "selection",
                                bias_parms = c(.73, .61, .82, .76))
              expect_equal(model$obs.data[1, 1], 197.3, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs.data[1, 2], 752.7, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs.data[2, 1], 133.4, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs.data[2, 2], 853.6, tolerance = 1e-4, scale = 1)
})

test_that("Correct table passed on (probabilistic)", {
              dat <- matrix(c(118, 832, 103, 884),
                            dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
                            nrow = 2, byrow = TRUE)
              set.seed(123)
              model <- dat %>%
                  probsens(., type = "exposure", reps = 5000,
                           seca.parms = list("trapezoidal", c(.45, .5, .6, .65)),
                           seexp.parms = list("trapezoidal", c(.4, .48, .58, .63)),
                           spca.parms = list("trapezoidal", c(.95, .97, .99, 1)),
                           spexp.parms = list("trapezoidal", c(.96, .98, .99, 1)),
                           corr.se = .8, corr.sp = .8) %>%
                  multiple.bias(., bias_function = "probsens.sel",
                                case.exp = list("logit-normal", c(-1.1, 0, 0, 1)),
                                case.nexp = list("trapezoidal", c(.75, .85, .95, 1)),
                                ncase.exp = list("logit-normal", c(-1.2, 0, 0, 1)),
                                ncase.nexp = list("trapezoidal", c(0.7, 0.8, 0.9, 1)))
              expect_equal(model$obs.data[1, 1], 182.6, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs.data[1, 2], 767.4, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs.data[2, 1], 169.3, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs.data[2, 2], 817.7, tolerance = 1e-4, scale = 1)
})
