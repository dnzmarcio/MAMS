\name{new.bounds}
\alias{new.bounds}

\title{
Function to update boundaries based on observed sample sizes
}

\description{
The function determines updated boundaries of a multi-arm multi-stage study based on observed number of observations per arm.
}

\usage{
new.bounds(K=3, J=2, alpha=0.05, nMat=matrix(c(10,20),nrow=2,ncol=4), 
           u=NULL, l=NULL, u.shape="obf", l.shape="fixed", lfix=0, 
           ufix=NULL, N=20)
}

\arguments{
  \item{K}{Number of experimental treatments (default=4).}
  \item{J}{Number of stages (default=2).}
  \item{alpha}{One-sided familywise error rate (default=0.05).}
  \item{nMat}{Jx(K+1) dimensional matrix of observed/expected sample sizes. Rows correspond to stages and columns to arms. First column is control (default: 2x4 matrix with 10 subjects per stage and arm). }
  \item{u}{ Vector of previously used upper boundaries (default=NULL).}
  \item{l}{ Vector of previously used upper boundaries (default=NULL).}
  \item{u.shape}{ Shape of upper boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}. See details.}
  \item{l.shape}{ Shape of lower boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}. See details.}
  \item{lfix}{ Fixed lower boundary (default=NULL). Only used if \code{shape="fixed"}.}
  \item{ufix}{ Fixed upper boundary (default=NULL). Only used if \code{shape="fixed"}.}
  \item{N}{ Number of quadrature points per dimension in the outer integral (default=20).}
}

\details{
This function finds the boundaries for a given matrix of sample sizes in multi-arm multi-stage study with K active treatments plus control. The vectors \code{u} and \code{l} are the boundaries used so far while \code{u.shape} and \code{l.shape} specify the shape to the boundaries for the remaining analysis. By specifying \code{u} and \code{l} as \code{NULL}, a design using only the shapes given by \code{u.shape} and \code{l.shape} can be found for any sample sizes per stage and arm.\cr

The shape of the boundaries (\code{u.shape}, \code{l.shape}) are either using the predefined shapes following Pocock (1977), O'Brien & Fleming (1979) or the triangular Test (Whitehead, 1997) using options \code{"pocock"}, \code{"obf"} or \code{"triangular"} respectively, are constant (option \code{"fixed"}) or supplied in as a function. If a function is passed it should require exactly one argument specifying the number of stages and return a vector of the same length. The lower boundary shape is required to be non-decreasing while the upper boundary shape needs to be non-increasing. If a fixed lower boundary is used, \code{lfix} must be smaller than \eqn{\Phi^{-1}(1-\alpha)/2}{Phi(1-alpha)/2} to ensure that it is smaller than the upper boundary.\cr
}

\value{
An object of the class MAMS containing the following components: \cr 
  \item{l}{Lower boundary. } 
  \item{u}{Upper boundary.}
  \item{n}{Sample size on control in stage 1.}
  \item{N}{Maximum total sample size.}
  \item{K}{Number of experimental treatments.}
  \item{J}{Number of stages in the trial. }
  \item{alpha}{Familywise error rate.}
  \item{power}{Power under least favorable configuration.}

  \item{rMat}{Matrix of allocation ratios. First row corresponds to control and second row to experimental treatments.}

}

\references{
Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett Test for Multi-arm Multi-stage Clinical Studies with Treatment Selection. Biometrika. 99(2):494-501.\cr

Pocock SJ (1977) Group sequential methods in the design and analysis of clinical trials. Biometrika 64, 191-199.\cr

O'Brien PC, Fleming TR (1979) A multiple testing procedure for clinical trials. Biometrics 35, 549-556.\cr 

Whitehead J (1997) The Design and Analysis of Sequential Clinical Trials. Wiley: Chichester.\cr 
}

\author{Thomas Jaki and Dominic Magirr}

\examples{
######
#### Note that some of these examples may take a few minutes to run
######
\dontrun{
# Note that some of these examples may take a few minutes to run
# 2-stage design with O'Brien & Fleming efficacy and zero futility boundary with 
# equal sample size per arm and stage. Results are equivalent to using 
# mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, u.shape="obf", l.shape="fixed", 
#      lfix=0, sample.size=FALSE)

new.bounds(K=4, J=2, alpha=0.05, nMat=matrix(c(10,20),nrow=2,ncol=5), u=NULL, l=NULL, 
           u.shape="obf", l.shape="fixed", lfix=0)

# A 2-stage design that was designed to use an O'Brien & Fleming efficacy and 
# zero futility boundary with equal sample size per arm and stage (n=14).
# The observed sample size after stage one are 10,10,18,10,13 for each arm 
# while the original upper bounds used are (3.068,2.169) for stage 1. 
# The updated bounds are (3.068,2.167).

new.bounds(K=4, J=2, alpha=0.05, 
           nMat=matrix(c(10,28,10,28,18,28,10,28,13,28),nrow=2,ncol=5), 
           u=3.068, l=0, u.shape="obf", l.shape="fixed", lfix=0)
}
}

\keyword{ design }
