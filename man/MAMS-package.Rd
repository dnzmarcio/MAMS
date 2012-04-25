\name{MAMS-package}
\alias{MAMS-package}
\alias{MAMS}
\docType{package}
\title{
Designing Multi-Arm Multi-Stage Studies
}
\description{
A package to design multi-arm multi-stage clinical studies with J stages and K treatments that allows early stopping for efficacy and futility. Normally distributed endpoints with known variance are assumed.
}
\details{
\tabular{ll}{
Package: \tab MAMS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-04-23\cr
License: \tab GPL-2 \cr
}

}
\author{
Thomas Jaki and Dominic Magirr

Maintainer: Thomas Jaki <jaki.thomas@gmail.com>

}
\references{
Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett Test for Multi-arm Multi-stage Clinical Studies with Treatment Selection. Biometrika. Published Online ahead of print.
}
\keyword{ package }
\examples{
# 2-stage design with O'Brien & Fleming efficacy and zero futility boundary
mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.65 , p0=0.55, u.shape="obf", 
     l.shape="fixed", lfix=0, nstart=40)
}
