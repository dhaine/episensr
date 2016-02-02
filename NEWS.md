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
