\name{confounders}
\alias{confounders}
\title{Simple sensitivity analysis to correct for unknown or unmeasured
  confounding (RR)}
\description{Simple sensitivity analysis to correct for unknown or unmeasured
  confounding (RR).}
\usage{
confounders(exposed, case, implement = c("noEM.RR", "noEM.OR", "noEM.RD"),
prev.cfder = NULL, cfder.dis.RR = NULL, cfder.dis.OR = NULL,
cfder.dis.RD = NULL, alpha = 0.05, dec = 4, print = TRUE)
}
\arguments{
  \item{exposed}{Exposure variable. If a variable, this variable is tabulated against.}
  \item{case}{Outcome variable.}
  \item{implement}{Choice of implementation, with no effect measure
    modification (noEM) for ratio measures (relative risk -- RR; odds
    ratio -- OR) or differenmce measures (risk difference -- RD).}
  \item{prev.cfder}{Vector defining the prevalence of the confounder. This
    vector has 2 elements between 0 and 1, in the following order:
    \begin{enumerate}
    \item Prevalence of the confounder among the exposed, and
    \item Prevalence of the confounder among the unexposed.
    \end{enumerate}
  }
  \item{cfder.dis.RR}{Confounder-disease relative risk.}
  \item{cfder.dis.OR}{Confounder-disease odds ratio.}
  \item{cfder.dis.RD}{Confounder-disease risk difference.}
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
  \item{RR0}{Relative risk corrected by SMR.}
  \item{OR0}{Odds ratio corrected by SMR.}
  \item{RRc}{Relative risk corrected by MH.}
  \item{ORc}{Odds ratio corrected by MH.}
  \item{RDc}{Risk difference corrected by MH.}
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
            implement = "noEM.RR",
            prev.cfder = c(.8, .05),
            cfder.dis.RR = .63)
confounders(matrix(c(105,85,527,93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow=2, byrow=TRUE),
            implement = "noEM.OR",
            prev.cfder = c(.8, .05),
            cfder.dis.OR = .63)
confounders(matrix(c(105,85,527,93),
                   dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
                   nrow=2, byrow=TRUE),
            implement = "noEM.RD",
            prev.cfder = c(.8, .05),
            cfder.dis.OR = -.37)
}
