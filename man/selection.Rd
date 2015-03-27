\name{selection}
\alias{selection}
\title{Simple sensitivity analysis to correct for selection bias using
  estimates of the selection proportions}
\description{Simple sensitivity analysis to correct for selection bias using
  estimates of the selection proportions.}
\usage{selection(exposed, case, selprob = NULL, alpha = 0.05, dec = 4,
print = TRUE)}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{selprob}{Vector defining the selection probabilities. This
    vector has 4 elements between 0 and 1, in the following order:
    \enumerate{
    \item Selection probability among cases exposed,
    \item Selection probability among cases unexposed,
    \item Selection probabillity among noncases exposed, and
    \item Selection probability among noncases unexposed.
    }
  }
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{corr.data}{The same table corrected for  selection proportions.}
  \item{obs.measures}{A table of odds ratios and relative risk with
    confidence intervals.}
  \item{adj.measures}{Selection bias corrected measures of
    outcome-exposure relationship.}
  \item{bias.parms}{Input bias parameters: selection probabilities.}
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.43--58, Springer.}
\author{Denis Haine \email{denis.haine@gmail.com}}
\examples{
# The data for this example come from:
# Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N.
# et al.
# Population-based incidence estimates of uveal melanoma in Germany. Supplementing
# cancer registry data by case-control data.
# Eur J Cancer Prev 2006;15:165-70.
selection(matrix(c(136, 107, 297, 165),
                 dimnames = list(c("UM+", "UM-"), c("Mobile+", "Mobile-")),
                 nrow = 2, byrow = TRUE),
          selprob = c(.94, .85, .64, .25))
}
