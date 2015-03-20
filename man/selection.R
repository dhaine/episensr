\name{selection}
\alias{selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple sensitivity analysis to correct for selection bias using
  estimates of the selection proportions
%%  ~~function to do ... ~~
}
\description{Simple sensitivity analysis to correct for selection bias using
  estimates of the selection proportions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
selection(exposed, case, selprob = NULL, alpha = 0.05, dec = 4,
print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{selprob}{Vector defining the selection probabilities. This
    vector has 4 elements between 0 and 1, in the following order:
    \begin{enumerate}
    \item Selection probability among cases exposed,
    \item Selection probability among cases unexposed,
    \item Selection probabillity among noncases exposed, and
    \item Selection probability among noncases unexposed.
    \end{enumerate}
  }
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{corr.data}{The same table corrected for  selection proportions.}
  \item{obs.measures}{A table of odds ratios and relative risk with
    confidence intervals.}
  \item{corr.measures}{Selection bias corrected measures of
    outcome-exposure relationship.}
  \item{probs}{Input bias parameters.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.43--58, Springer. 
%% ~put references to the literature/web site here ~
}
\author{Denis Haine \email{denis.haine@gmail.com}
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# The data for this example come from:
# Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N. et al.
# Population-based incidence estimates of uveal melanoma in Germany. Supplementing cancer registry data by case-control data.
# Eur J Cancer Prev 2006;15:165-70.
selection(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE), selprob = c(.94, .85, .64, .25))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
