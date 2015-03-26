\name{probsens.conf}
\alias{probsens.conf}
\alias{probsens}
\alias{probsens.sel}
\title{Probabilistic sensitivity analysis for unmeasured confounding.}
\description{Probabilistic sensitivity analysis to correct for unknown
  or unmeasurred confounding and random error simultaneously.}
\usage{
probsens.conf(exposed, case, reps = 1000, prev.exp = list(dist =
c("uniform", "triangular", "trapezoidal"), parms = NULL), prev.nexp =
list(dist = c("uniform", "triangular", "trapezoidal"), parms = NULL),
risk = list(dist = c("uniform", "triangular", "trapezoidal"), parms = NULL), corr.p
= NULL, alpha = 0.05, dec = 4, print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{reps}{Number of replications to run.}
  \item{prev.exp}{List defining the prevalence of exposure among the
    exposed. The first argument provides the probability distribution
    function (uniform, triangular, or trapezoidal) and the second its
    parameters as a vector:
    \begin{enumerate}
    \item Uniform: min, max,
    \item Triangular: lower limit, upper limit, mode,
    \item Trapezoidal: min, lower mode, upper mode, max.
    \end{enumerate}
  }
  \item{prev.nexp}{List defining the prevalence of exposure among the unexposed.}
  \item{risk}{List defining the confounder-disease relative risk or the
    confounder-exposure odds ratio.}
  \item{corr.p}{Correlation between the exposure-specific confounder prevalences.}
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{obs.measures}{A table of observed relative risk and odds ratio
    with confidence intervals.}
  \item{adj.measures}{A table of corrected relative risks and odds ratios.}
  \item{sim.df}{Data frame of random parameters and computed values.}
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.117--150, Springer.}
\author{Denis Haine \email{denis.haine@gmail.com}}
\examples{
# The data for this example come from:
# Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O. et al.
# Increased risk of infection with human immunodeficiency virus type 1 among uncircumcised men presenting with genital ulcer disease in Kenya.
# Clin Infect Dis 1996;23:449-53.
set.seed(123)
probsens.conf(matrix(c(105, 85, 527, 93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow = 2, byrow = TRUE),
         reps = 20000,
         prev.exp = list("trapezoidal", c(.7, .8, .8, .9)),
         prev.nexp = list("trapezoidal", c(.03, .04, .05, .06)),
         risk = list("trapezoidal", c(.6, .63, .63, .7)),
         corr.p = .8)
}
