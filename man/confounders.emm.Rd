\name{confounders.emm}
\alias{confounders.emm}
\alias{confounders}
\alias{confounders.poly}
\title{Simple sensitivity analysis to correct for unknown or unmeasured
  confounding with effect measure modification}
\description{Simple sensitivity analysis to correct for unknown or unmeasured
  confounding with effect measure modification. Implementation for ratio
  measures (relative risk -- RR, or odds ratio -- OR) and difference
  measures (risk difference -- RD).}
\usage{
confounders(exposed, case, implement = c("RR", "OR", "RD"),
p = NULL, RR.cd = NULL, OR.cd = NULL, RD.cd = NULL, alpha = 0.05, dec =
4, print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{implement}{Choice of implementation, with no effect measure
    modification for ratio measures (relative risk -- RR; odds
    ratio -- OR) or difference measures (risk difference -- RD).}
  \item{p}{Vector defining the prevalence of the confounder. This vector
    has 2 elements between 0 and 1, in the following order:
    \begin{enumerate}
    \item Prevalence of the confounder among the exposed, and
    \item Prevalence of the confounder among the unexposed.
    \end{enumerate}
  }
  \item{RR.cd}{Vector defining the confounder-disease relative risk.
    This vector has two elements between 0 and 1, in the following
    order:
    \begin{enumerate}
    \item Among the exposed, and
    \item among the unexposed.
    \end{enumerate}}
  \item{OR.cd}{Vector defining the confounder-disease odds ratio. This
    vector has two elements between 0 and 1, in the folowing order:
    \begin{enumerate}
    \item Among the exposed, and
    \item among the unexposed.
    \end{enumerate}}
  \item{RD.cd}{Vector dfining the confounder-disease risk difference. This
    vector has two elements between 0 and 1, in the following 
    \begin{enumerate}
    \item Among the exposed, and
    \item among the unexposed.
    \end{enumerate}}
  \item{alpha}{Significance level.}
  \item{dec}{Number of decimals in the printout.}
  \item{print}{Should the results be printed?}
}
\value{A list with elements:
  \item{obs.data}{The analysed 2 x 2 table from the observed data.}
  \item{cfder.data}{The same table for Confounder +.}
  \item{nocfder.data}{The same table for Confounder -.}
  \item{obs.measures}{A table of relative risk with confidence intervals.}
  \item{SMR}{Standardized Morbidity Ratio.}
  \item{MH}{Mantel-Haenszel.}
  \item{cfder.rr}{Relative risk, Confounder +.}
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
    Quantitative Bias Analysis to Epidemiologic Data}, pp.59--78, Springer.}
\author{Denis Haine \email{denis.haine@gmail.com}}
\examples{
# The data for this example come from:
# Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O. et al.
# Increased risk of infection with human immunodeficiency virus type 1 among uncircumcised men presenting with genital ulcer disease in Kenya.
# Clin Infect Dis 1996;23:449-53.
confounders.emm(matrix(c(105, 85, 527, 93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow = 2, byrow = TRUE),
            implement = "RR",
            p = c(.8, .05),
            RR.cd = c(.4, .7))
confounders.emm(matrix(c(105, 85, 527, 93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow = 2, byrow = TRUE),
            implement = "OR",
            p = c(.8, .05),
            OR.cd = c(.4, .7))
confounders.emm(matrix(c(105, 85, 527, 93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow = 2, byrow = TRUE),
            implement = "RD",
            p = c(.8, .05),
            RD.cd = c(-.6, -.3))
}
