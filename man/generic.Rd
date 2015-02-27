\name{plot}
\alias{plot.MAMS}
\alias{print.MAMS}
\alias{summary.MAMS}
\alias{print.MAMS.sim}
\alias{summary.MAMS.sim}
\alias{print.MAMS.step_down}
\alias{summary.MAMS.step_down}
\alias{plot.MAMS.step_down}


\title{Different generic functions for class MAMS.}
\description{Generic functions for summarizing an object of class MAMS.}
\usage{
\method{print}{MAMS}(x, digits=max(3, getOption("digits") - 4), \dots)

\method{summary}{MAMS}(object, digits=max(3, getOption("digits") - 4), \dots)

\method{plot}{MAMS}(x, col=NULL, pch=NULL, lty=NULL, main=NULL, xlab="Analysis", 
     ylab="Test statistic", ylim=NULL, type=NULL, \dots)

\method{print}{MAMS.sim}(x, digits=max(3, getOption("digits") - 4), \dots)

\method{summary}{MAMS.sim}(object, digits=max(3, getOption("digits") - 4), \dots)

\method{print}{MAMS.step_down}(x, digits=max(3, getOption("digits") - 4), ...)

\method{summary}{MAMS.step_down}(object, digits=max(3, getOption("digits") - 4), \dots)

\method{plot}{MAMS.step_down}(x, col=NULL, pch=NULL, lty=NULL, main=NULL, xlab="Analysis", 
     ylab="Test statistic", ylim=NULL, type=NULL, \dots)


}

\arguments{
  \item{x}{ An output object of class MAMS.}
  \item{digits}{Number of significant digits to be printed.}
  \item{object}{ An output object of class MAMS.}
  \item{col}{A specification for the default plotting color (default=\code{NULL}). See \code{par} for more details.}
  \item{pch}{Either an integer specifying a symbol or a single character to be used as the default in plotting points (default=\code{NULL}). See \code{par} for more details.}
  \item{lty}{A specification for the default line type to be used between analyses (default=\code{NULL}). Setting to zero supresses ploting of the lines. See \code{par} for more details.}
  \item{main}{An overall title for the plot (default=\code{NULL}).}
  \item{xlab}{A title for the x axis (default=\code{"Analysis"}).}
  \item{ylab}{A title for the y axis (default=\code{"Test statistic"}).}
  \item{ylim}{Numeric vector of length 2, giving the y coordinates range (default=\code{"NULL"}).}
  \item{type}{Type of plot to be used (default=\code{NULL}). See \code{plot} for more details.}
  \item{\dots}{Further (graphical) arguments to be passed to methods.}
}

\details{
  \code{print.MAMS} produces a summary of an object from class MAMS including boundaries and requires sample size if initially requested.\cr

  \code{summary.MAMS} produces same output as \code{print.MAMS}.\cr

  \code{plot.MAMS} produces as plot of the boundaries.

  \code{print.MAMS.sim} produces a summary of an object from class MAMS.sim including type-I-error and expected sample size.\cr

  \code{summary.MAMS.sim} produces same output as \code{print.MAMS.sim}.\cr

  \code{print.MAMS.step_down} produces a summary of an object from class MAMS including boundaries and requires sample size if initially requested.\cr

  \code{summary.MAMS.step_down} produces same output as \code{print.step_down_mams}.\cr

  \code{plot.MAMS.step_down} produces a plot of the boundaries.

}

\value{Screen or graphics output.}

\references{
Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett Test for Multi-arm Multi-stage Clinical Studies with Treatment Selection. Biometrika. 99(2):494-501.\cr


Stallard N, Todd S (2003) Sequential designs for phase III clinical trials incorporating treatment selection. Statistics in Medicine. 22, 689-703. \cr

Magirr D, Stallard N, Jaki T (2014) Flexible sequential designs for multi-arm clinical trials. Statistics in Medicine. Published online. 
}

\author{Thomas Jaki and Dominic Magirr}

\examples{
\donttest{
# 2-stage design with triangular boundaries
res <- mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.65 , p0=0.55, 
            u.shape="triangular", l.shape="triangular", nstart=30)

print(res)
summary(res)
plot(res)
}
res <- mams.sim(nsim=10000,nMat=matrix(c(44,88),nrow=2,ncol=5),u=c(3.068,2.169),
                l=c(0.000,2.169),pv=c(0.65,0.55,0.55,0.55),ptest=c(1:2,4))

print(res)

# 2-stage 3-treatments versus control design, all promising treatments are selected:
res <- step_down_mams(nMat = matrix(c(10, 20), nrow = 2, ncol = 4), 
                      alpha_star = c(0.01, 0.05), lb = 0, 
                      selection = "all_promising")

print(res)
summary(res)
plot(res)

}

\keyword{classes}
