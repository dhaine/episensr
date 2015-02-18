\name{multidimBias}
\alias{multidimBias}
\title{Multidimensional sensitivity analysis for different sources of
  bias}
\description{Multidimensional sensitivity analysis for different sources of bias}
\usage{multidimBias(exposed, case,
se = NULL, sp = NULL, alpha = 0.05, dec = 4,
print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{se}{Vector of sensitivities.}
  \item{sp}{Vector of specificities.}  
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{obs.measures}{A table of odds ratios and relative risk with
    confidence intervals.}
  \item{rr.mat}{Multidimensional corrected relative risk data.}
  \item{or.mat}{Multidimensional corrected odds ratio data.}  
  \item{sesp}{Sensitivities and specificities.}
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.109--116, Springer. 
}
\author{Denis Haine \email{denis.haine@gmail.com}
}
\examples{
multidimBias(matrix(c(45, 94, 257, 945),
                    dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                    nrow = 2, byrow = TRUE),
          se = c(1, 1, 1, .9, .9, .9, .8, .8, .8),
          sp = c(1, .9, .8, 1, .9, .8, 1, .9, .8))
}
