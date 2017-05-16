# episensr 0.8.x

- Add `misclassification_cov` for a misclassified covariate (confounder or
  effect modifier).
- As such, this (wrongly) added option into `misclassification` in the previous
  version is now removed.

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
