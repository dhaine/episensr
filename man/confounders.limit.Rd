\name{confounders.limit}
\alias{confounders}
\alias{confounders.emm}
\alias{confounders.poly}
\title{Bounding the bias limits of unmeasured confounding.}
\description{Function to elicit the limits on measures of effect
  corrected for an unmeasured confounder when only some of the bias
  parameters are known.}
\usage{
confounders.limit(p = NA, RR = NA, OR = NA, crude.RR = NULL, dec = 4, print = TRUE)
}
\arguments{
  \item{p}{Proportion with the confounder among the unexposed group.}
  \item{RR}{Relative risk between the confounder and the outcome.}
  \item{OR}{Odds ratio between the confounder and the outcome.}
  \item{crude.RR}{Crude relative risk between the exposure and the outcome.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{conf.limits}{Limits on confounding.}
  \item{bias.params}{Input bias parameters.}
}
\references{
  Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
    Bias Analysis to Epidemiologic Data}, pp.59--78, Springer.
  Flanders, W. Dana, Khoury, Muin J., 1990. Indirect Assessment of Confounding:
  Graphic Description and Limits on Effect of Adjusting for Covariates.
  \emph{Epidemiology} 1(3): 239--246.
}
\author{Denis Haine \email{denis.haine@gmail.com}}
\examples{
confounders.limit(OR = 1.65, crude.RR = 1.5)
}
