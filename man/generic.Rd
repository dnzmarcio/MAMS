\name{plot}
\alias{plot.MAMS}
\alias{print.MAMS}
\alias{summary.MAMS}
\title{Different generic functions for class MAMS.}
\description{Generic functions for summarizing an object of class PK}
\usage{
\method{print}{MAMS}(x, digits=max(3, getOption("digits") - 4), \dots)

\method{summary}{MAMS}(object, digits=max(3, getOption("digits") - 4), \dots)

\method{plot}{MAMS}(x, col=NULL, pch=NULL, lty=NULL, main=NULL, xlab="Analysis", 
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
  \code{print.MAMS} produces a minimal summary of an estimation object from class PK including point estimate, standard error and confidence interval. The confidence interval is the first of \code{"boott"}, \code{"fieller"}, \code{"t"} or \code{"z"} that was originally requested. \cr\cr

  \code{summary.MAMS} prints a more detailed summary of an estimation object from class PK. Most notably all confidence intervals originally requested are printed.\cr\cr

  \code{plot.MAMS} produces as concentration versus time plot of the data used of an estimation object from class PK.
}

\value{Screen or graphics output.}

\references{
Magirr D, Jaki T, Whitehead J (2012) A generalized Dunnett Test for Multi-arm Multi-stage Clinical Studies with Treatment Selection. Biometrika. Published online ahead of print.\cr
}

\author{Thomas Jaki}

\examples{
\donttest{
# 3-stage design with triangular boundaries
res <- mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.65 , p0=0.55, u.shape="triangular", 
            l.shape="triangular", nstart=30)

print(res)
summary(res)
plot(res)
}
}

\keyword{classes}
