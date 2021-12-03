\name{plotsl}
\alias{plotsl}
\alias{plotslc}
\alias{plotsld}
\title{Plots of density and distribution function for the quantile based skew logistic distribution}
\description{
Produces plots of density and distribution function for the  quantile based skew 
logistic distribution.  Although you could use \code{curve(dsl(x))} 
to do this, the fact that the density and quantiles of this distribution are
defined in terms of the depth, \eqn{u}, means that a seperate function that
uses the depths to produce the values to plot is more accurate and efficient.
}
\usage{
plotsld(parameters, add=FALSE, granularity = 10000,
  xlab = "x", ylab="density", quant.probs = seq(0,1,.25), ...)
plotslc(parameters, add=FALSE, granularity = 10000,
  xlab = "quantile", ylab="depth", quant.probs = seq(0,1,.25), ...)}
\arguments{
  \item{parameters}{A vector of length 3, giving the parameters of the 
  quantile-based skew logistic distribution.  The 3 elements are 
  \eqn{\alpha}{alpha} (location), \eqn{\beta}{beta} (scale) and
  \eqn{\delta}{delta} (skewing).  \eqn{\alpha}{alpha} can take on any real
  value, \eqn{\beta}{beta} can take on any positive value and 
  \eqn{\delta}{delta} must satisfy \eqn{0 \leq \delta \leq 1}{0 <= delta <= 1}.
  \eqn{\delta = 0.5}{delta = 0.5} gives the logistic distribution,
  \eqn{\delta = 0}{delta = 0} gives the reflected exponential distribution and
  \eqn{\delta = 1}{delta = 1} gives the exponential distribution.
  }
  \item{add}{Should this add to an existing plot (using \code{lines}) or 
  produce a new plot (using \code{plot}).}
  \item{granularity}{Number of points at which quantiles and density will be 
  calculated.--- see \emph{details}}
  \item{xlab}{X axis label}
  \item{ylab}{Y axis label}
  \item{quant.probs}{Quantiles of distribution to return (see \emph{value}
  below).  Set to NULL to suppress this return entirely.}
  \item{...}{arguments that get passed to \code{plot} if this is a new plot}
}
\details{
  The quantile-based skew logistic distribution is defined by its quantile
  function.  The density of the distribution is available explicitly as a
  function of depths, \eqn{u}, but not explicitly available as a function of 
  \eqn{x} (except for the special cases at \eqn{\delta}{delta}=0,0.5 and 1).
  This function calculates quantiles and density as a function of depths to
  produce a density plot \code{plotsld} or cumulative probability plot
  \code{plotslc}.
}
\value{
  A number of quantiles from the distribution, the default being the minimum,
maximum and quartiles.
}
\references{
van Staden, P.J. and King, Robert A.R. (2015) \emph{The quantile-based skew logistic distribution}, Statistics and Probability Letters \bold{96},
109--116. 
\doi{10.1016/j.spl.2014.09.001}

van Staden, Paul J. 2013 \emph{Modeling of generalized families of probability distribution in the quantile statistical universe}. 
PhD thesis, University of Pretoria.
\url{http://hdl.handle.net/2263/40265}

\url{https://github.com/newystats/sld}
}
\author{Robert King, \email{robert.king.newcastle@gmail.com},
\url{https://github.com/newystats/sld}}

\seealso{\code{\link{sld}}}
\examples{
plotsld(c(0,1,1),main="Exponential Distribution")
plotsld(c(0,1,0.5),main="Logistic Distribution")
plotsld(c(0,1,0.7))
plotslc(c(0,1,0.7))
}
\keyword{distribution}
\keyword{hplot}
\keyword{aplot}

