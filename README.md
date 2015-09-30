---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# episensr

[![Build Status](https://travis-ci.org/dhaine/episensr.svg?branch=master)](https://travis-ci.org/dhaine/episensr)

The R package **episensr** allows to do basic sensitivity analysis of
epidemiological results as described in **Applying Quantitative Bias Analysis to
Epidemiological Data** by Timothy L. Lash, Matthew P. Fox, and Aliza K. Fink
(ISBN: 978-0-387-87960-4,
[bias.analysis](https://sites.google.com/site/biasanalysis/)). A similar
function is available in Stata
([episens](http://ideas.repec.org/c/boc/bocode/s456792.html)).

## License

This package is free and open source software, licensed under GPL2.

## Example

We will use a case-control study by
[Stang et al.](http://www.ncbi.nlm.nih.gov/pubmed/16523014) on the relation
between mobile phone use and uveal melanoma.
The observed odds ratio for the association between regular mobile phone use vs.
no mobile phone use with uveal melanoma incidence is 0.71 [95% CI 0.51-0.97].
But there was a substantial difference in participation rates between cases and
controls (94% vs 55%, respectively) and so selection bias could have an impact
on the association estimate.
The 2X2 table for this study is the following:

|          | regular use | no use |
|----------|-------------|--------|
| cases    | 136         | 107    |
| controls | 297         | 165    |

We use the function `selection` as shown below.


```r
library(episensr)

selection(matrix(c(136, 107, 297, 165),
                 dimnames = list(c("UM+", "UM-"), c("Mobile+", "Mobile-")),
                 nrow = 2, byrow = TRUE),
          selprob = c(.94, .85, .64, .25))
#> Observed Data: 
#> --------------------------------------------------- 
#> Outcome   : UM+ 
#> Comparing : Mobile+ vs. Mobile- 
#> 
#>     Mobile+ Mobile-
#> UM+     136     107
#> UM-     297     165
#> 
#> Data Corrected for Selected Proportions: 
#> ---------------------------------------------------
#> 
#>      Mobile+  Mobile-
#> UM+ 144.6809 125.8824
#> UM- 464.0625 660.0000
#> 
#>                                95% conf. interval
#> Observed Relative Risk: 0.7984    0.6518   0.9780
#>    Observed Odds Ratio: 0.7061    0.5144   0.9693
#> 
#>                                           [,1]
#> Selection Bias Corrected Relative Risk: 1.4838
#>    Selection Bias Corrected Odds Ratio: 1.6346
#> 
#>                                                 [,1]
#>      Selection probability among cases exposed: 0.94
#>    Selection probability among cases unexposed: 0.85
#>   Selection probability among noncases exposed: 0.64
#> Selection probability among noncases unexposed: 0.25
```

The 2X2 table is provided as a matrix and selection probabilities given with the
argument selprob, a vector with the 4 probabilities (guided by the participation
rates in cases and controls) in the following order: among cases exposed, among
cases unexposed, among noncases exposed, and among noncases unexposed.
The output shows the observed 2X2 table, the same table corrected for the
selection proportions, the observed odds ratio (and relative risk) followed by
the corrected ones, and the input parameters.

Here's an other example to correct for selection bias caused by M bias.


```r
mbias(or = c(2, 5.4, 2.5, 1.5, 1),
      var = c("HIV", "Circumcision", "Muslim", "Low CD4", "Participation"))
#> Correction for selection bias: 
#> ---------------------------------------- 
#> OR observed between the exposure and the outcome: 1 
#>              Maximum bias from conditioning on P: 1.006236 
#>                  OR corrected for selection bias: 0.9938024
```

## Installation

You can get the latest release from **CRAN**:


```r
install.packages('episensr')
```

Or install the development version from **GitHub** with **devtools** package:


```r
devtools::install_github('dhaine/episensr', ref = "develop")
```
