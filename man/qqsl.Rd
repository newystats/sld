\name{qqsl}
\alias{qqsl}
\title{Quantile-Quantile plot against the skew logistic distribution}
\description{\code{qqgl} produces a Quantile-Quantile plot of data against the
quantile-based skew logistic distribution, 
or a Q-Q plot to compare two sets of parameter values
for the quantile-based skew logistic distribution.  This function does for the 
skew logistic distribution what \code{\link{qqnorm}} does for the normal.
}
\usage{
qqsl(y=NULL,parameters1,parameters2=NULL,abline=TRUE,
      granularity.for.2.dists=4000,use.endpoints=FALSE,...)
}
\arguments{
  \item{y}{The data sample}
  \item{parameters1}{A vector of length 3, containg the parameters of the 
  skew logistic distribution, \eqn{\alpha}{alpha}, \eqn{\beta}{beta} 
  and \eqn{\delta}{delta}.
  }
  \item{parameters2}{Second set of parameters of the skew logistic distribution. A vector of length 3, as described above for \code{parameters1}.}
  \item{abline}{A logical value, TRUE adds a line through the origian
with a slope of 1 to the plot}
  \item{granularity.for.2.dists}{Number of quantiles to use in a Q-Q plot comparing two sets of parameter values}
  \item{use.endpoints}{logical. When comparing two sets of parameter values, should Q(0) and Q(1) be used?  TRUE will give QQ plots including the theoretical
  minimum and maximum of the distribution, which is arguably not equivalent to
  what would be seen in QQ plots based on data.  FALSE will give QQ plots based on 
  ideal depths (type 8 quantiles, see \code{\link[stats]{quantile}} 
  documentation), where \code{n} is the \code{granularity.for.2.dists}.}
  \item{\dots}{graphical parameters, passed to \code{\link{qqplot}}}
}
\details{
See \code{\link{sld}} for more details on the Skew Logistic
Distribution.  A Q-Q plot provides a way to visually assess the
correspondence between a dataset and a particular distribution, or between two 
distributions.
}
\value{
  A list of the same form as that returned by \code{\link{qqline}}
  \item{x}{The x coordinates of the points that were/would be plotted,
corresponding to a skew logistic distibution with parameters 
\eqn{\alpha}{alpha}, \eqn{\beta}{beta} and \eqn{\delta}{delta}.}
  \item{y}{The original \code{y} vector, i.e., the corresponding y
    coordinates, or a corresponding set of quantiles from a skew logistic
distribution with the second set of parameters}
}
\references{
van Staden, P.J. and King, R.A.R. (2015) \emph{The quantile-based skew logistic distribution}, Statistics and Probability Letters \bold{96} 109--116.
\url{http://dx.doi.org/10.1016/j.spl.2014.09.001}

van Staden, Paul J. 2013 \emph{Modeling of generalized families of probability distribution in the quantile statistical universe}. 
PhD thesis, University of Pretoria.
\url{http://hdl.handle.net/2263/40265}

\url{http://tolstoy.newcastle.edu.au/rking/SLD/SLD.html}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/} and Paul van Staden}

\seealso{\code{\link{sld}}}
\examples{
qqsl(y=rsl(100,c(0,1,0.7)),parameters1=c(0,1,0.7))
qqsl(parameters1=c(0,1,0.7),parameters2=c(0,0.9,0.5),col="blue")
}
\keyword{distribution}
\keyword{hplot}
\keyword{aplot}
