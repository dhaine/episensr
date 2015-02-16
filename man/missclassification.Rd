\name{misclassification}
\alias{misclassification}
\title{Simple sensitivity analysis for misclassification}
\description{Simple sensitivity analysis for misclassification.}
\usage{
misclassification(exposed, case, implement = c("exposure", "outcome", "confounder"),
bias = NULL, alpha = 0.05, dec = 4, print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{implement}{Choice of implementation:
    \begin{enumerate}
    \item Exposure: bias analysis of exposure misclassification; corrections using
  sensitivity and specificity: nondifferential and independent errors,
  \item  Outcome: bias analysis of outcome misclassification.
  \end{enumerate}
}
  \item{bias}{Vector defining the bias parameters. This vector has 4
    elements between 0 and 1, in the following order:
    \begin{enumerate}
    \item Sensitivity of exposure classification among those with the
    outcome, 
    \item Sensitivity of exposure classification among those without the
    outcome, 
    \item Specificity of exposure classification among those with the
    outcome,and 
    \item Specificity of exposure classification among those without the outcome.
    \end{enumerate}
  }
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{corr.data}{The expected observed data given the true data
    assuming misclassfication.}
  \item{obs.measures}{A table of observed relative risk and odds ratio
    with confidence intervals.}
  \item{corr.rr}{Corrected relative risk.}
  \item{corr.or}{Corrected odds ratio.}
  \item{bias.params}{Input bias parameters.}
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.79--108, Springer. 
}
\author{Denis Haine \email{denis.haine@gmail.com}
}
\examples{
# The data for this example come from:
# Fink, A.K., Lash,  T.L. A null association between smoking during pregnancy and breast cancer using Massachusetts registry data (United States).
# Cancer Causes Control 2003;14:497-503.
misclassification(matrix(c(215,1449,668,4296),
                   dimnames = list(c("Breast cancer+", "Breast cancer-"),
                       c("Smoker+", "Smoker-")),
                   nrow=2, byrow=TRUE),
            implement = "exposure",
            bias = c(.78, .78, .99, .99))
misclassification(matrix(c(105,85,527,93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow=2, byrow=TRUE),
            implement = "OR",
            prev.cfder = c(.8, .05),
            cfder.dis.OR = .63)
misclassification(matrix(c(105,85,527,93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow=2, byrow=TRUE),
            implement = "RD",
            prev.cfder = c(.8, .05),
            cfder.dis.RD = -.37)
}
