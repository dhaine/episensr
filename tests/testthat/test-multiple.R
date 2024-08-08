context("check multiple bias analysis")

dat <- matrix(c(118, 832, 103, 884),
              dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
              nrow = 2, byrow = TRUE)

chien <- matrix(c(118, 832, 103, 884),
                dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
                nrow = 2, byrow = TRUE)

test_that("episensr object is provided", {
              expect_error(table(dat) %>%
                           multiple.bias(., bias_function = "misclass",
                                         type = "exposure",
                                         bias_parms = c(.56, .58, .99, .97)))
          })

test_that("bias function is provided", {
             expect_error(dat %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(.,
                                        RR = 1, bias_parms = c(0.1, 0.9, 0.1, 0.4)))
         })

test_that("Correct functions are selected", {
              dat <- matrix(c(118, 832, 103, 884),
                            dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
                            nrow = 2, byrow = TRUE)
              expect_that(dat %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "confounders.ext",
                                        RR = 1, bias_parms = c(0.1, 0.9, 0.1, 0.4)),
                          throws_error())
})

test_that("Correct table passed on", {
              dat <- matrix(c(118, 832, 103, 884),
                            dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
                            nrow = 2, byrow = TRUE)
              model <- dat %>%
                  misclass(., type = "exposure",
                           bias_parms = c(.56, .58, .99, .97)) %>%
                  multiple.bias(., bias_function = "selection",
                                bias_parms = c(.73, .61, .82, .76))
              expect_equal(model$obs_data[1, 1], 197.3, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_data[1, 2], 752.7, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_data[2, 1], 133.4, tolerance = 1e-4, scale = 1)
              expect_equal(model$obs_data[2, 2], 853.6, tolerance = 1e-4, scale = 1)
})

#test_that("Correct table passed on (probabilistic)", {
#              dat <- matrix(c(118, 832, 103, 884),
#                            dimnames = list(c("BC+", "BC-"), c("AD+", "AD-")),
#                            nrow = 2, byrow = TRUE)
#              set.seed(123)
#              model2 <- dat %>%
#                  probsens(., type = "exposure", reps = 5000,
#                           seca = list("trapezoidal", c(.45, .5, .6, .65)),
#                           seexp = list("trapezoidal", c(.4, .48, .58, .63)),
#                           spca = list("trapezoidal", c(.95, .97, .99, 1)),
#                           spexp = list("trapezoidal", c(.96, .98, .99, 1)),
#                           corr_se = .8, corr_sp = .8) %>%
#                  multiple.bias(., bias_function = "probsens.sel",
#                                case_exp = list("beta", c(114, 12)),
#                                case_nexp = list("trapezoidal", c(.75, .85, .95, 1)),
#                                ncase_exp = list("beta", c(370, 100)),
#                                ncase_nexp = list("trapezoidal", c(0.7, 0.8, 0.9, 1)))
#              expect_equal(model$obs_data[1, 1], 182.9, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_data[1, 2], 767.1, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_data[2, 1], 168.9, tolerance = 1e-4, scale = 1)
#              expect_equal(model$obs_data[2, 2], 818.1, tolerance = 1e-4, scale = 1)
#})

test_that("selection bias has correct arguments", {
             expect_error(chien %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "selection"))
             expect_error(chien %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "selection",
                                        bias = c(0.5, 0.5, 0.5, 0.5)))
             expect_error(chien %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "selection",
                                        bias_parms = c(.73, .61, .82, .76),
                                        beta = .8))
         })

test_that("confounding bias has correct arguments", {
             expect_error(chien %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "selection",
                                        bias_parms = c(.73, .61, .82, .76)) %>%
                          multiple.bias(., bias_function = "confounders",
                                        type = "OR"))
             expect_error(chien %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "selection",
                                        bias_parms = c(.73, .61, .82, .76)) %>%
                          multiple.bias(., bias_function = "confounders",
                                        tipe = "OR", bias_parms = c(.92, .3, .44)))
             expect_error(chien %>%
                          misclass(., type = "exposure",
                                   bias_parms = c(.56, .58, .99, .97)) %>%
                          multiple.bias(., bias_function = "selection",
                                        bias_parms = c(.73, .61, .82, .76)) %>%
                          multiple.bias(., bias_function = "confounders",
                                        type = "OR", bias_parms = c(.92, .3, .44),
                                        beta = 0.8))
         })

test_that("misclassification bias has correct arguments", {
             expect_error(chien %>%
                          selection(., bias_parms = c(.73, .61, .82, .76)) %>%
                          multiple.bias(., bias_function = "misclass",
                                        bias_parms = c(.56, .58, .99, .97)))
             expect_error(chien %>%
                          selection(., bias_parms = c(.73, .61, .82, .76)) %>%
                          multiple.bias(., bias_function = "misclass",
                                        tipo = "exposure",
                                        bias_parms = c(.56, .58, .99, .97)))
         })

test_that("probsens has correct arguments", {
             set.seed(123)
             expect_error(chien %>%
                          probsens.sel(., reps = 1000,
                                       case_exp = list("beta", c(8.08, 24.25)),
                                       case_nexp = list("trapezoidal", c(.75, .85, .95, 1)),
                                       ncase_exp = list("beta", c(12.6, 50.4)),
                                       ncase_nexp = list("trapezoidal",
                                                         c(0.7, 0.8, 0.9, 1))) %>%
                          multiple.bias(., bias_function = "probsens",
                                        seca = list("trapezoidal",
                                                    c(.45, .5, .6, .65))))
             expect_error(chien %>%
                          probsens.sel(., reps = 1000,
                                       case_exp = list("beta", c(8.08, 24.25)),
                                       case_nexp = list("trapezoidal", c(.75, .85, .95, 1)),
                                       ncase_exp = list("beta", c(12.6, 50.4)),
                                       ncase_nexp = list("trapezoidal",
                                                         c(0.7, 0.8, 0.9, 1))) %>%
                          multiple.bias(., bias_function = "probsens",
                                        type = "exposure",
                                        seca = list("trapezoidal",
                                                    c(.45, .5, .6, .65)),
                                        seexp = list("trapezoidal",
                                                     c(.4, .48, .58, .63)),
                                        shibidi = list("trapezoidal",
                                                       c(.95, .97, .99, 1)),
                                        spexp = list("trapezoidal",
                                                     c(.96, .98, .99, 1)),
                                        corr_se = .8, corr_sp = .8))
         })

test_that("probsens.sel has correct arguments", {
             set.seed(123)
             expect_error(chien %>%
                          probsens(., type = "exposure", reps = 1000,
                                   seca = list("trapezoidal", c(.45, .5, .6, .65)),
                                   seexp = list("trapezoidal", c(.4, .48, .58, .63)),
                                   spca = list("trapezoidal", c(.95, .97, .99, 1)),
                                   spexp = list("trapezoidal", c(.96, .98, .99, 1)),
                                   corr_se = .8, corr_sp = .8) %>%
                          multiple.bias(., bias_function = "probsens.sel"))
         })

test_that("probsens.conf has correct arguments", {
             set.seed(123)
             expect_error(chien %>%
                          probsens(., type = "exposure", reps = 1000,
                                   seca = list("trapezoidal", c(.45, .5, .6, .65)),
                                   seexp = list("trapezoidal", c(.4, .48, .58, .63)),
                                   spca = list("trapezoidal", c(.95, .97, .99, 1)),
                                   spexp = list("trapezoidal", c(.96, .98, .99, 1)),
                                   corr_se = .8, corr_sp = .8) %>%
                          multiple.bias(., bias_function = "probsens.sel",
                                        case_exp = list("beta", c(8.08, 24.25)),
                                        case_nexp = list("trapezoidal", c(.75, .85, .95, 1)),
                                        ncase_exp = list("beta", c(12.6, 50.4)),
                                        ncase_nexp = list("trapezoidal",
                                                          c(0.7, 0.8, 0.9, 1))) %>%
                          multiple.bias(., bias_function = "probsens_conf",
                                        prev_nexp = list("beta", c(42.9, 54.6)),
                                        risk = list("trapezoidal", c(.2, .58, 1.01, 1.24))))
             expect_error(chien %>%
                          probsens(., type = "exposure", reps = 1000,
                                   seca = list("trapezoidal", c(.45, .5, .6, .65)),
                                   seexp = list("trapezoidal", c(.4, .48, .58, .63)),
                                   spca = list("trapezoidal", c(.95, .97, .99, 1)),
                                   spexp = list("trapezoidal", c(.96, .98, .99, 1)),
                                   corr_se = .8, corr_sp = .8) %>%
                          multiple.bias(., bias_function = "probsens.sel",
                                        case_exp = list("beta", c(8.08, 24.25)),
                                        case_nexp = list("trapezoidal", c(.75, .85, .95, 1)),
                                        ncase_exp = list("beta", c(12.6, 50.4)),
                                        ncase_nexp = list("trapezoidal",
                                                          c(0.7, 0.8, 0.9, 1))) %>%
                          multiple.bias(., bias_function = "probsens_conf",
                                        whoaa = list("beta", c(24.9, 58.1)),
                                        prev_nexp = list("beta", c(42.9, 54.6)),
                                        risk = list("trapezoidal", c(.2, .58, 1.01, 1.24))))
         })
