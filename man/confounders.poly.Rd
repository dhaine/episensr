\name{confounders.poly}
\alias{confounders.poly}
\alis{confounders}
\title{Simple sensitivity analysis to correct for unknown or unmeasured
  polychotomous confounding without effect modification}
\description{Simple sensitivity analysis to correct for unknown or unmeasured
  polychotomous confounding without effect modification. Implementation for ratio
  measures (relative risk -- RR, or odds ratio -- OR) and difference
  measures (risk difference -- RD).}
\usage{
confounders.poly(exposed, case, implement = c("RR", "OR", "RD"),
prev.cfder = NULL, cfder.dis.RR = NULL, cfder.dis.OR = NULL,
cfder.dis.RD = NULL, alpha = 0.05, dec = 4, print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{implement}{Choice of implementation, with no effect measure
    modification for ratio measures (relative risk -- RR; odds
    ratio -- OR) or difference measures (risk difference -- RD).}
  \item{prev.cfder}{Vector defining the prevalence of the confounder. This
    vector has 4 elements between 0 and 1, in the following order:
    \begin{enumerate}
    \item Prevalence of the highest level confounder among the exposed,
    \item Prevalence of the highest level confounder among the
    unexposed,
    \item Prevalence of the mid-level confounder among the exposed, and
    \item Prevalence of the mid-level confounder among the unexposed.
    \end{enumerate}
  }
  \item{cfder.dis.RR}{Vector defining the confounder-disease relative
    risk. This vector has two elements in the following order:
    \begin{enumerate}
  \item Relative risk of having the highest-level confounder in
  diseased, and
  \item Relative risk of having the mid-level confounder in disased.
  \end{enumerate}}
  \item{cfder.dis.OR}{Confounder-disease odds ratio.}
  \item{cfder.dis.RD}{Confounder-disease risk difference.}
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{cfder1.data}{The same table for Mid-level Confounder +.}
  \item{cfder2.data}{The same table for Highest-level Confounder +.}
  \item{nocfder.data}{The same table for Confounder -.}
  \item{obs.measures}{A table of relative risk with confidence intervals.}
  \item{SMR}{Standardized Morbidity Ratio.}
  \item{MH}{Mantel-Haenszel.}
  \item{cfder1.rr}{Relative risk, Mid-level Confounder +.}
  \item{cfder2.rr}{Relative risk, Highest-level confounder +.}
  \item{cfder.or}{Odds ratio, Confounder +.}
  \item{cfder.rd}{Risk difference, Confounder +.}
  \item{nocfder.rr}{Relative risk, Confounder -.}
  \item{nocfder.or}{Odds ratio, Confounder -.}
  \item{nocfder.rd}{Risk difference, Confounder -.}
  \item{RRadj.smr}{Relative risk adjusted for the effect of the
    confounder using SMR estimate.}
  \item{ORadj.smr}{Odds ratio adjusted for the effect of the confounder using
    SMR estimate.}
  \item{RRadj.mh}{Relative risk adjusted for the effect of the
    confounder using MH estimate.}
  \item{ORadj.mh}{Odds ratio adjusted for the effect of the confounder
    using MH estimate.}
  \item{RDadj.mh}{Risk difference adjusted for the effect of the
    confounder using MH estimate.}
  \item{bias.params}{Input bias parameters.}
}
\references{Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying
    Quantitative Bias Analysis to Epidemiologic Data}, pp.43--58, Springer. 
}
\author{Denis Haine \email{denis.haine@gmail.com}
}
\examples{
# The data for this example come from:
# Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O. et al.
# Increased risk of infection with human immunodeficiency virus type 1 among uncircumcised men presenting with genital ulcer disease in Kenya.
# Clin Infect Dis 1996;23:449-53.
confounders(matrix(c(105,85,527,93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow=2, byrow=TRUE),
            implement = "RR",
            prev.cfder = c(.6, .05, .2, .2),
            cfder.dis.RR = c(.4, .8))
confounders(matrix(c(105,85,527,93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow=2, byrow=TRUE),
            implement = "OR",
            prev.cfder = c(.8, .05),
            cfder.dis.OR = .63)
confounders(matrix(c(105,85,527,93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow=2, byrow=TRUE),
            implement = "RD",
            prev.cfder = c(.8, .05),
            cfder.dis.RD = -.37)
}
