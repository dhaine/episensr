\name{probsens.sel}
\alias{probsens.sel, probsens, probsens.conf}
\title{Probabilistic sensitivity analysis for selection bias.}
\description{Probabilistic sensitivity analysis to correct for selection
  bias.}
\usage{
probsens.sel(exposed, case, reps = 1000, or.parms = list(dist =
c("uniform", "triangular", "trapezoidal"), parms = NULL), alpha = 0.05, dec = 4, print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{reps}{Number of replications to run.}
  \item{or.parms}{List defining the selection bias odds. The first
    argument provides the probability distribution function (uniform, triangular,
    or trapezoidal) and the second its parameters as a vector:
    \enumerate{
    \item Uniform: min, max,
    \item Triangular: lower limit, upper limit, mode,
    \item Trapezoidal: min, lower mode, upper mode, max.
    }
  }
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{obs.measures}{A table of observed odds ratio with confidence intervals.}
  \item{adj.measures}{A table of corrected odds ratios.}
  \item{sim.df}{Data frame of random parameters and computed values.}
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.}
\author{Denis Haine \email{denis.haine@gmail.com}}
\examples{
# The data for this example come from:
# Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N.
# et al.
# Population-based incidence estimates of uveal melanoma in Germany. Supplementing
# cancer registry data by case-control data.
# Eur J Cancer Prev 2006;15:165-70.
set.seed(123)
probsens.sel(matrix(c(136, 107, 297, 165),
                   dimnames = list(c("Melanoma+", "Melanoma-"), c("Mobile+", "Mobile-")),
                   nrow = 2, byrow = TRUE),
         reps = 20000,
         or.parms = list("triangular", c(.35, 1.1, .43)))
}
