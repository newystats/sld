\name{lmom.sample}
\alias{lmom.sample}
\title{Calculate sample L-Moments}
\description{Calculate sample L-Moments of a dataset}
\usage{
lmom.sample(data,max.mom=3)
}
\arguments{
  \item{data}{A vector containing a dataset}
  \item{max.mom}{The maximum order of L-Moment to estimate}
}
\details{
This function is a wrapper around the \code{\link[lmom]{samlmu}} function from
Hosking's \code{lmom} package, to give different argument defaults.  
It calculates sample L Moments.
}
\value{
  A vector of length \code{max.mom}.  The first two elements are the first two
  L moments.  If \code{max.mom} is greater than two, the following elements are
  the corresponding L moment ratios (the L Moment divided by the 2nd L Moment).
}
\references{
Hosking, J. R. M. (1990). \emph{L-moments: analysis and estimation of distributions using linear combinations of order statistics.} Journal of the Royal Statistical Society, Series B, \bold{52}, 105--124.

Hosking, J. R. M., and Wallis, J. R. (1997). \emph{Regional frequency analysis: an approach based on L-moments.} Cambridge University Press.

\url{http://tolstoy.newcastle.edu.au/rking/SLD/SLD.html}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/} and Paul van Staden}

\seealso{\code{\link{sld}}}
\examples{
generated.data <- rsl(300,c(0,1,.4))
lmom.sample(data=generated.data,max.mom=3)
}
\keyword{estimation}