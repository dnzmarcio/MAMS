\name{stepdown.update}
\alias{stepdown.update}

\title{
Update the stopping boundaries of multi-arm multi-stage study at an interim analysis, allowing for unplanned treatment selection and/or sample-size reassessment. 
}

\description{
Function to update a planned multi-arm multi-stage design to account for unplanned adaptations.
}

\usage{
stepdown.update(current.mams=stepdown.mams(), nobs=NULL, 
                zscores=NULL, selected.trts=NULL, nfuture=NULL)
}

\arguments{
  \item{current.mams}{The planned step-down MAMS design prior to the current interim analysis (default=stepdown.mams()).}
  \item{nobs}{Cumulative sample sizes observed on each treatment arm up to and including the current interim analysis.}
  \item{zscores}{Observed vector of test statistics at the current interim analysis.}
  \item{selected.trts}{The set of experimental treatments to be taken forward to the next stage of testing. This argument should be omitted at the final analysis.}
  \item{nfuture}{A matrix of future cumulative sample sizes. The number of rows must be equal to the originally planned number of stages (2 or 3) minus the number of stages already observed. The number of columns must be equal to the number of treatment arms (default=NULL).}
}

\details{
The function implements the ideas described in Magirr et al. (2014) to update a design according to unplanned design modifications. It takes as input the planned multi-arm multi-stage design prior to the interim analysis, together with the actually observed cumulative sample sizes and test statistics. Treatments to be included in future stages, as well as future sample sizes, can be chosen without following pre-specified rules. The output is a new multi-arm multi-stage design for the remaining stages such that the familywise error remains controlled at the pre-specified level. 
}

\value{
An object of the class MAMS.stepdown containing the following components: \cr 
  \item{l}{Lower boundaries.} 
  \item{u}{Upper boundaries.}
  \item{sample.sizes}{Cumulative sample sizes on each treatment arm.}
  \item{K}{Number of experimental treatments.}
  \item{J}{Number of stages in the trial.}
  \item{alpha.star}{Cumulative familywise error rate spent at each analysis, conditional on results so far.}
  \item{selection}{Pre-specified method of treatment selection.}
  \item{zscores}{A list containing the observed test statistics at analyses so far (at the design stage this is NULL).}
  \item{selected.trts}{A list containing the treatments selected for each stage.}
}

\references{
Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection. Biometrika, 99(2), 494-501.\cr

Stallard N, Todd S (2003) Sequential designs for phase III clinical trials incorporating treatment selection. Statistics in Medicine, 22(5), 689-703.\cr

Magirr D, Stallard N, Jaki T (2014) Flexible sequential designs for multi-arm clinical trials. Statistics in Medicine, 33(19), 3269-3279. 
}

\author{Dominic Magirr}

\examples{\donttest{
# 2-stage 3-treatments versus control design
# all promising treatments are selected:
orig_mams <- stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4), 
                           alpha.star=c(0.01, 0.05), lb=0, 
                           selection="all.promising")

# make adjustment for the observed sample sizes 
# not being exactly as planned:
stepdown.update(orig_mams, nobs=c(9, 8, 13, 11), zscores=c(1.1, -0.5, 0.2), 
                selected.trts=1:3, nfuture=NULL)

# make adjustment for the observed sample sizes 
# not being exactly as planned. In addition, drop treatment 2:
stepdown.update(orig_mams, nobs=c(9, 8, 13, 11), zscores=c(1.1, -0.5, 0.2), 
                selected.trts=c(1, 3), nfuture=NULL)

# make adjustment for the observed sample sizes not being 
# exactly as planned. In addition, drop treatment 2. In addition, 
# double the planed cumulative second stage sample sizes:
updated_mams1 <- stepdown.update(orig_mams, nobs=c(9, 8, 13, 11), 
                                 zscores=c(1.1, -0.5, 0.2), selected.trts=c(1, 3), 
                                 nfuture=matrix(c(40, 40, 13, 40), nrow=1, ncol=4))

# Account for the observed second stage sample sizes:
stepdown.update(updated_mams1, nobs=c(38, 41, 13, 36), zscores=c(1.9, -Inf, 1.2), 
                selected.trts=NULL)

# 'select.best' design. Account for actually observed sample sizes 
# in first stage, and drop treatment 2:
orig_mams <- stepdown.mams(nMat=matrix(c(10, 20), nrow=2, ncol=4), 
                           alpha.star=c(0.01, 0.05), lb=0, selection="select.best")

stepdown.update(orig_mams, nobs=c(9, 8, 13, 11), zscores=c(1.1, -0.5, 0.2), 
                selected.trts=c(1, 3), nfuture=NULL)
}}

\keyword{design}
