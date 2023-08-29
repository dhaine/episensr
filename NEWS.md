# episensr 1.3.0

# episensr 1.2.0
- Fix CI for OR in `plot.booted` (@codiewood, #2)
- Replaced dot-dot notation by after_stat(density) in `plot.probsens` and `plot.booted`
- `misclassification_cov` replaced by `misclassification.cov`
- Standard output for `confounders.ext`, `confounders.limit`, and `confounders.array`
- Suggest `tidyr` instead of `tidyverse`

# episensr 1.1.0
- Fix confidence interval values in plotting output of `probsens`.
- Fix bug that created integer overflow in rare cases when 2-by-2 table cells
  were very large.
- Update documentation.

# episensr 1.0.0
- Provide default values for variables in `mbias`.
- Update M-bias plot by using `ggdag` (and `dagitty`).
- Add beta distribution to probabilistic bias analyses.
- Plotting functions for probabilistic bias analyses.
- Review criteria to remove negative cell counts in probabilistic bias analyses.
- Add citation and DOI.

# episensr 0.9.6
- Fix formulas in `probsens.conf`.

# episensr 0.9.5
- Fix selection-bias factor in `selection` which was returning a constant value.

# episensr 0.9.4
- Fix `bias_parms` in `selection` to use a single selection-bias factor (was
  skipping it and returning NAs).

# episensr 0.9.3
- Update description of `parms` in `confounders`, `confounders.emm`, and
  `confounders.poly`
- New function `confounders.array`, sensitivity analysis for unmeasured
  confounders based on confounding imbalance among exposed and unexposed
  (Schneeweiss, 2006)
- Deprecated `bias` parameters removed from function `multidimBias`. Please use
  `bias_parms` instead.
- New function `confounders.evalue`, computing E-value to assess bias due to
  unmeasured confounder (VanderWeele and Ding, 2017)
- New function `multiple.bias` allowing to extract 2-by-2 table from an
  `episensr` object to feed another function (multiple bias analysis)

# episensr 0.9.2
- Fix bug for distributions and computations of OR/RR in `probsens.conf`
- Update distributions in `probsens.irr.conf`
- Add example using `probsens.conf` in vignette

# episensr 0.9.1
- Fix bug when using triangular distribution in `probsens.conf` function for
  prevalence of exposure among the non-exposed (as producing NaNs) (#1).

# episensr 0.9.0

- Add `misclassification_cov` for a misclassified covariate (confounder or
  effect modifier).
- As such, this (wrongly) added option into `misclassification` in the previous
  version is now removed.
- Add computation of confidence interval for odds ratio as per Chu et al. for
  exposure misclassification.
- Use `bias_parms` instead of `bias` in `misclassification` function.

# episensr 0.8.0

- Fix bug when building 2-by-2 table.

- Various formatting improvements in output of `confounders`, `confounders.emm`,
  `misclassification` and `selection` functions.
- Standardize use of `bias_parms`.

- Add vignette.

- Selection bias factor now available in output of `selection` function.
- Add bootstrap option

# episensr 0.7.2

- Fix 2-by-2 tables when variables are provided instead of a matrix.

# episensr 0.7.1

- Fix R version dependency (R >= 3.2.0)

# episensr 0.7.0

- Harmonization of arguments across functions.

- New distributions added to `probsens` series of functions: constant,
logit-logistic, logit-normal, log-logistic, and log-normal.

- Probabilistic analysis of person-time data added with `probsens.irr` for
exposure misclassification, and `probsens.irr.conf` for unmeasured confounder.

- Sensitivity analysis to correct for selection bias caused by M bias with
`mbias` function, including DAG plot and print function.

- Fix CI formatting.

- NAMESPACE: add imports to `stats` functions to avoid new R CMD CHECK warnings
