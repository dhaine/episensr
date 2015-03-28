\name{confounders.poly}
\alias{confounders.poly}
\alias{confounders, confounders.emm, confounders.limit}
\title{Simple sensitivity analysis to correct for unknown or unmeasured
  polychotomous confounding without effect modification}
\description{Simple sensitivity analysis to correct for unknown or unmeasured
  polychotomous (3-level) confounding without effect modification.
  Implementation for ratio measures (relative risk -- RR, or odds ratio
  -- OR) and difference measures (risk difference -- RD).}
\usage{
confounders.poly(exposed, case, implement = c("RR", "OR", "RD"),
p = NULL, RR.cd = NULL, OR.cd = NULL, RD.cd = NULL, alpha = 0.05,
dec = 4, print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{implement}{Choice of implementation, with no effect measure
    modification for ratio measures (relative risk -- RR; odds
    ratio -- OR) or difference measures (risk difference -- RD).}
  \item{p}{Vector defining the prevalence of the confounder. This vector
    has 4 elements between 0 and 1, in the following order:
    \enumerate{
    \item Prevalence of the highest level confounder among the exposed,
    \item Prevalence of the highest level confounder among the
    unexposed,
    \item Prevalence of the mid-level confounder among the exposed, and
    \item Prevalence of the mid-level confounder among the unexposed.
    }
  }
  \item{RR.cd}{Vector defining the confounder-disease relative risk.
    This vector has two elements in the following order:
    \enumerate{
    \item Relative risk of having the highest-level confounder in
    diseased, and
    \item Relative risk of having the mid-level confounder in diseased.
    }
  }
  \item{OR.cd}{Vector defining the confounder-disease odds ratio. This
    vector has two elements in the following order:
    \enumerate{
    \item Odds ratio of having the highest-level confounder in diseased,
    and
    \item Odds ratio of having the mid-level confounder in diseased.
    }
  }
  \item{RD.cd}{Vector defining the confounder-disease risk
    difference. This vector has two elements in the following order:
    \enumerate{
    \item Risk difference of having the highest-level confounder in
    diseased, and
    \item Risk difference of having the mid-level confounder in diseased.
    }
  }
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{cfder1.data}{The same table for Mid-level Confounder +.}
  \item{cfder2.data}{The same table for Highest-level Confounder +.}
  \item{nocfder.data}{The same table for Confounder -.}
  \item{obs.measures}{A table of relative risk with confidence
    intervals; Total and by confounders.}
  \item{adj.measures}{A table of Standardized Morbidity Ratio and
    Mantel-Haenszel estimates.}
  \item{bias.parms}{Input bias parameters.}
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.59--78, Springer.}
\author{Denis Haine \email{denis.haine@gmail.com}}
\examples{
# The data for this example come from:
# Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O.
# et al.
# Increased risk of infection with human immunodeficiency virus type 1 among
# uncircumcised men presenting with genital ulcer disease in Kenya.
# Clin Infect Dis 1996;23:449-53.
confounders.poly(matrix(c(105, 85, 527, 93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow = 2, byrow = TRUE),
            implement = "RR",
            p = c(.6, .05, .2, .2),
            RR.cd = c(.4, .8))
confounders.poly(matrix(c(105, 85, 527, 93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow = 2, byrow = TRUE),
            implement = "OR",
            p = c(.6, .05, .2, .2),
            OR.cd = c(.4, .8))
confounders.poly(matrix(c(105, 85, 527, 93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow = 2, byrow = TRUE),
            implement = "RD",
            p = c(.6, .05, .2, .2),
            RD.cd = c(-.4, -.2))
}
