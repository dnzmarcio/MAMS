\name{stepdown.mams}
\alias{stepdown.mams}

\title{
Function to find stopping boundaries for a 2- or 3-stage (step-down) multiple-comparisons-with-control test.
}

\description{
The function determines stopping boundaries for all intersection hypothesis tests in a multi-arm multi-stage study, given the amount of alpha (familywise error rate) to be spent at each analysis.
}

\usage{
stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4),
              alpha.star=c(0.01, 0.025), lb=0,
              selection="all.promising")
}

\arguments{
  \item{nMat}{Matrix containing the cumulative sample sizes in each treatment arm (columns: control, trt 1, \dots, trt K), at each analysis (rows). The number of analyses must be either 2 or 3 (default=matrix(c(10, 20), nrow=2, ncol=4)).}
  \item{alpha.star}{Cumulative familywise error rate to be spent at each analysis (default=c(0.01, 0.025)).}
  \item{lb}{Fixed lower boundary (default=0).}
  \item{selection}{How are treatments selected for the next stage? Using the default "all.promising" method, all treatments with a test statistic exceeding the lower boundary are taken forward to the next stage. If "select.best", only the treatment with the largest statistic may be selected for future stages. (default="all.promising").}
}

\details{
The function implements the methods described in Magirr et al. (2014) to find individual boundaries for all intersection hypotheses.
}

\value{
An object of the class MAMS.stepdown containing the following components: \cr 
  \item{l}{Lower boundaries.} 
  \item{u}{Upper boundaries.}
  \item{nMat}{Cumulative sample sizes on each treatment arm.}
  \item{K}{Number of experimental treatments.}
  \item{J}{Number of stages in the trial.}
  \item{alpha.star}{Cumulative familywise error rate spent at each analysis.}
  \item{selection}{Pre-specified method of treatment selection.}
  \item{zscores}{A list containing the observed test statistics at analyses so far (at the design stage this is NULL).}
  \item{selected.trts}{A list containing the treatments selected for each stage.}

}

\references{
Jaki T, Pallmann P and Magirr D (2019). "The R Package MAMS for Designing Multi-Arm Multi-Stage Clinical Trials."" Journal of Statistical Software, 88(4), pp. 1-25. doi: 10.18637/jss.v088.i04
(URL: http://doi.org/10.18637/jss.v088.i04)

Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection. Biometrika, 99(2), 494-501.\cr

Magirr D, Stallard N, Jaki T (2014) Flexible sequential designs for multi-arm clinical trials. Statistics in Medicine, 33(19), 3269-3279.

Stallard N, Todd S (2003) Sequential designs for phase III clinical trials incorporating treatment selection. Statistics in Medicine, 22(5), 689-703.\cr
}

\author{Dominic Magirr}

\examples{\donttest{
# Note that some of these examples may take a few minutes to run

# 2-stage 3-treatments versus control design, all promising treatments are selected:
stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4), 
              alpha.star=c(0.01, 0.05), lb=0, 
              selection="all.promising")

# select the best treatment after the first stage:
stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4), 
              alpha.star=c(0.01, 0.05), lb=0, 
              selection="select.best")

# 3 stages and unequal randomization:
stepdown.mams(nMat=matrix(c(20, 40, 60, rep(c(10, 20, 30), 3)), nrow=3, ncol=4), 
              alpha.star=c(0.01, 0.025, 0.05), lb=c(0, 0.75), 
              selection="all.promising")
}}

\keyword{ design }
