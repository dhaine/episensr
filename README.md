
<!-- README.md is generated from README.Rmd. Please edit that file -->

# episensr <img src="man/figures/logo.png" align="right" width=120 />

[![Build
Status](https://travis-ci.org/dhaine/episensr.svg?branch=master)](https://travis-ci.org/dhaine/episensr)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/episensr)](https://cran.r-project.org/package=episensr)
[![Codecov test
coverage](https://codecov.io/gh/dhaine/episensr/branch/master/graph/badge.svg)](https://codecov.io/gh/dhaine/episensr?branch=master)
[\!Active](https://www.repostatus.org/badges/latest/active.svg)
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Total CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/episensr)](https://cran.r-project.org/package=episensr)

The R package **episensr** allows to do basic sensitivity analysis of
epidemiological results as described in **Applying Quantitative Bias
Analysis to Epidemiological Data** by Timothy L. Lash, Matthew P. Fox,
and Aliza K. Fink (ISBN: 978-0-387-87960-4,
[bias.analysis](https://sites.google.com/site/biasanalysis/)). A similar
function is available in Stata
([episens](http://ideas.repec.org/c/boc/bocode/s456792.html)).

## License

This package is free and open source software, licensed under GPL2.

## Example

We will use a case-control study by [Stang et
al.](http://www.ncbi.nlm.nih.gov/pubmed/16523014) on the relation
between mobile phone use and uveal melanoma. The observed odds ratio for
the association between regular mobile phone use vs. no mobile phone use
with uveal melanoma incidence is 0.71 \[95% CI 0.51-0.97\]. But there
was a substantial difference in participation rates between cases and
controls (94% vs 55%, respectively) and so selection bias could have an
impact on the association estimate. The 2X2 table for this study is the
following:

|          | regular use | no use |
| -------- | ----------- | ------ |
| cases    | 136         | 107    |
| controls | 297         | 165    |

We use the function `selection` as shown below.

``` r
library(episensr)
#> Loading required package: ggplot2

selection(matrix(c(136, 107, 297, 165),
                 dimnames = list(c("UM+", "UM-"), c("Mobile+", "Mobile-")),
                 nrow = 2, byrow = TRUE),
          bias_parms = c(.94, .85, .64, .25))
#> --Observed data-- 
#>          Outcome: UM+ 
#>        Comparing: Mobile+ vs. Mobile- 
#> 
#>     Mobile+ Mobile-
#> UM+     136     107
#> UM-     297     165
#> 
#>                                        2.5%     97.5%
#> Observed Relative Risk: 0.7984287 0.6518303 0.9779975
#>    Observed Odds Ratio: 0.7061267 0.5143958 0.9693215
#> ---
#>                                                 
#> Selection Bias Corrected Relative Risk: 1.483780
#>    Selection Bias Corrected Odds Ratio: 1.634608
```

The 2X2 table is provided as a matrix and selection probabilities given
with the argument `bias_parms`, a vector with the 4 probabilities
(guided by the participation rates in cases and controls) in the
following order: among cases exposed, among cases unexposed, among
noncases exposed, and among noncases unexposed. The output shows the
observed 2X2 table, the observed odds ratio (and relative risk) followed
by the corrected ones.

## Installation

You can get the latest release from **CRAN**:

``` r
install.packages('episensr')
```

Or install the development version from **GitHub** with **devtools**
package:

``` r
devtools::install_github('dhaine/episensr', ref = "develop")
```
