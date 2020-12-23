\name{ordinal.mams}
\alias{ordinal.mams}

\title{
Function to design multi-arm multi-stage studies with ordinal or binary endpoints
}

\description{
The function determines (approximately) the boundaries of a multi-arm multi-stage study with ordinal or binary endpoints for a given boundary shape and finds the required number of subjects.
}

\usage{
ordinal.mams(prob=c(0.35, 0.4, 0.25), or=2, or0=1.2, K=4, J=2, alpha=0.05,
     power=0.9, r=1:2, r0=1:2, ushape="obf", lshape="fixed", ufix=NULL,
     lfix=0, nstart=1, nstop=NULL, sample.size=TRUE, N=20)
}

\arguments{
  \item{prob}{Vector of expected probabilities of falling into each category under control conditions. The elements must sum up to one (default=c(0.35, 0.4, 0.25)).}
  \item{or}{Interesting treatment effect on the scale of odds ratios (default=2).}
  \item{or0}{Uninteresting treatment effect on the scale of odds ratios (default=1.2).}
  \item{K}{Number of experimental treatments (default=4).}
  \item{J}{Number of stages (default=2).}
  \item{alpha}{One-sided familywise error rate (default=0.05).}
  \item{power}{Desired power (default=0.9).}
  \item{r}{Vector of allocation ratios (default=1:2).}
  \item{r0}{Vector ratio on control (default=1:2).}
  \item{ushape}{Shape of upper boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}.}
  \item{lshape}{Shape of lower boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"}, \code{"triangular"} and \code{"fixed"} (the default).}
  \item{ufix}{Fixed upper boundary (default=NULL). Only used if \code{shape="fixed"}.}
  \item{lfix}{Fixed lower boundary (default=0). Only used if \code{shape="fixed"}.}
  \item{nstart}{Starting point for finding the sample size (default=1).}
  \item{nstop}{Stopping point for finding the sample size (default=NULL).}
  \item{sample.size}{Logical if sample size should be found as well (default=TRUE).}
  \item{N}{Number of quadrature points per dimension in the outer integral (default=20).}
}

\details{
This function finds the (approximate) boundaries and sample size of a multi-arm multi-stage study with ordinal or binary endpoints with K active treatments plus control in which all promising treatments are continued at interim analyses as described in Magirr et al (2012). It is a wrapper around the basic \code{mams} function to facilitate its use with ordinal and binary endpoints, following ideas of Whitehead & Jaki (2009) and Jaki & Magirr (2013). For a binary endpoint the vector \code{prob} has only two elements (success/failure, yes/no, etc.). See \code{?mams} for further details on the basic methodology.\cr
}

\value{
An object of the class MAMS containing the following components:\cr 
  \item{l}{Lower boundary.}
  \item{u}{Upper boundary.}
  \item{n}{Sample size on control in stage 1.}
  \item{N}{Maximum total sample size.}
  \item{K}{Number of experimental treatments.}
  \item{J}{Number of stages in the trial.}
  \item{alpha}{Familywise error rate.}
  \item{alpha.star}{Cumulative familywise error rate spent by each analysis.}
  \item{power}{Power under least favorable configuration.}
  \item{rMat}{Matrix of allocation ratios. First row corresponds to control while subsequent rows are for the experimental treatments.}

}

\references{
Jaki T, Magirr D (2013) Considerations on covariates and endpoints in multi-arm multi-stage clinical trials selecting all promising treatments. Statistics in Medicine, 32(7), 1150-1163.\cr

Jaki T, Pallmann P and Magirr D (2019). "The R Package MAMS for Designing Multi-Arm Multi-Stage Clinical Trials."" Journal of Statistical Software, 88(4), pp. 1-25. doi: 10.18637/jss.v088.i04
(URL: http://doi.org/10.18637/jss.v088.i04)

Magirr D, Stallard N, Jaki T. (2014) Flexible sequential designs for multi-arm clinical
trials. Statistics in Medicine. 33(19):3269-3279.

Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection. Biometrika, 99(2), 494-501.\cr

Whitehead J, Jaki T (2009) One- and two-stage design proposals for a phase II trial comparing three active treatments with control using an ordered categorical endpoint. Statistics in Medicine, 28(5), 828-847.
}

\author{Philip Pallmann}

\examples{\donttest{
## An example based on the example in Whitehead & Jaki (2009)
# 2-stage design with triangular efficacy and futility boundaries
prob <- c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166)
ordinal.mams(prob=prob, or=3.06, or0=1.32, K=3, J=2, alpha=0.05,
             power=0.9, r=1:2, r0=1:2, ushape="triangular",
             lshape="triangular")
}}

\keyword{design}
