\name{mams}
\alias{mams}

\title{
Function to design multi-arm multi-stage studies with normal endpoints
}

\description{
The function determines the boundaries of a multi-arm multi-stage study for a given boundary shape and finds the required number of subjects.
}

\usage{
mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.75 , p0=0.5, 
     u.shape="obf", l.shape="fixed", lfix=0, ufix=NULL, nstart=1, 
     sample.size=TRUE, N=20)
}

\arguments{
  \item{K}{Number of experimental treatments (default=4).}
  \item{J}{Number of stages (default=2).}
  \item{alpha}{One-sided familywise error rate (default=0.05).}
  \item{power}{ Desired power (default=0.9).}
  \item{r}{ Vector of allocation ratios (default=1:2).}
  \item{r0}{ Vector ratio on control (default=1:2).}
  \item{p}{ Interesting treatment effect. See Details (default=0.75).}
  \item{p0}{ Uninteresting treatment effect. See Details (default=0.5).}
  \item{u.shape}{ Shape of upper boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}. See details.}
  \item{l.shape}{ Shape of lower boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}. See details.}
  \item{lfix}{ Fixed lower boundary (default=0). Only used if \code{shape="fixed"}.}
  \item{ufix}{ Fixed upper boundary (default=NULL). Only used if \code{shape="fixed"}.}
  \item{nstart}{ Starting point for finding the sample size (default=1).}
  \item{sample.size}{ Logical if sample size should be found as well (default=TRUE). }
  \item{N}{ Number of quadrature points per dimension in the outer integral (default=20).}
}

\details{
This function finds the boundaries and sample size of a multi-arm multi-stage study with K active treatments plus control in which all promising treatments are continued at interim analyses as described in Magirr et al (2012). At each interim analysis the test statistics are compared to the lower (futility) bound and any treatment whose corresponding test statistic falls below that bound is discontinued. Similarly if any test statistic exceeds the upper (efficacy) bound the null hypothesis corresponding to that treatment can be rejected and superiority of that treatment over control claimed. At the same time the study is stopped. If at least one test statistic exceeds the lower bound and none exceeds the upper bound the study is continued and further patients are recruited to all remaining experimental treatments plus control.\cr

The design is found under the least favorable configuration which requires an interesting treatment effect \code{p} that if present we would like to find with high probability and an uninteresting effect \code{p0}. Both \code{p} and \code{p0} are parameterized as \eqn{P(X_k > X_0 ) = p}{P(X_k > X_0 ) = p}, that is the probability of a randomly selected person on treatment k observing a better outcome than a random person on control. For \code{p=0.5} the experimental treatment and control perform equally well. The advantage of this paramterization is that no knowledge about the variance is required. To convert traditional effect sizes, \eqn{\delta}{delta} to this format use \eqn{p=\Phi(\frac{\delta}{\sqrt{2}\sigma})}{Phi(delta/(2^0.5*sigma))}.\cr

The shape of the boundaries (\code{u.shape}, \code{l.shape}) are either using the predefined shapes following Pocock (1977), O'Brien & Fleming (1979) or the triangular Test (Whitehead, 1997) using options \code{"pocock"}, \code{"obf"} or \code{"triangular"} respectively, are constant (option \code{"fixed"}) or supplied in as a function. If a function is passed it should require exactly one argument specifying the number of stages and return a vector of the same length. The lower boundary shape is required to be non-decreasing while the upper boundary shape needs to be non-increasing. If a fixed lower boundary is used, \code{lfix} must be smaller than \eqn{\Phi^{-1}(1-\alpha)/2}{Phi(1-alpha)/2} to ensure that it is smaller than the upper boundary.\cr

Computation of designs with more than four stages are very time consuming and not advised. The parameter \code{sample.size} controls whether the required sample size is computed as well. Setting to \code{FALSE} approximately halves the compuation time.\cr
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
  \item{alpha_star}{Cumulative familywise error rate spent by each analysis.}
  \item{power}{Power under least favorable configuration.}

  \item{rMat}{Matrix of allocation ratios. First row corresponds to control while subsequent rows are for the experimental treatments.}

}

\references{
Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett Test for Multi-arm Multi-stage Clinical Studies with Treatment Selection. Biometrika. 99(2):494-501.\cr

Pocock SJ (1977) Group sequential methods in the design and analysis of clinical trials. Biometrika 64, 191-199.\cr

O'Brien PC, Fleming TR (1979) A multiple testing procedure for clinical trials. Biometrics 35, 549-556.\cr 

Whitehead J (1997) The Design and Analysis of Sequential Clinical Trials. Wiley: Chichester.\cr 
}

\author{Thomas Jaki and Dominic Magirr}

\examples{
## A fixed sample (single stage) design
mams(K=4, J=1, alpha=0.05, power=0.9, r=1, r0=1, p=0.65 , p0=0.55)


## An example in Table 1 of Magirr et al (2012)
# 2-stage design with O'Brien & Fleming efficacy and zero futility boundary
mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.65 , p0=0.55, 
     u.shape="obf", l.shape="fixed", lfix=0, nstart=40)

\dontrun{
# Note that these examples may take a few minutes to run
# 3-stage design with Triangular efficacy and futility boundary
mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.65 , p0=0.55, 
     u.shape="triangular", l.shape="triangular", nstart=30)

## Example of a custom boundary function without sample size evaluation
mams(K=6, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.7 , p0=0.5, 
     u.shape=function(x)return(x:1), l.shape="fixed", lfix=0, 
     sample.size=FALSE)

## Different allocation ratios between control and experimental treatments. 
## Twice as many patients are randomized to control at each stage.
mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=c(2, 4), p=0.65 , p0=0.55, 
     u.shape="obf", l.shape="fixed", lfix=0, nstart=30)
}

}

\keyword{ design }
