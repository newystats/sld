\name{fit.sld.lmom}
\alias{fit.sld.lmom.given}
\alias{fit.sld.lmom}
\title{Fit the skew logistic distribution using L-Moments}
\description{Fits the quantile-based Skew Logistic Distribution using L-Moments.
\code{fit.sld.lmom} calculates the sample L Moments of a dataset and uses the
method of L Moments to estimate the parameters of the skew logistic distribution.
\code{fit.sld.lmom.given} fits the skew logistic using user-supplied values
of the first three L Moments.
}
\usage{
fit.sld.lmom.given(lmoms)
fit.sld.lmom(data)
}
\arguments{
  \item{lmoms}{A vector of length 3, containing the first and second (sample)
  L Moments and the 3rd (sample) L Moment ratio (\eqn{\tau_3}{tau3 })}
  \item{data}{A vector containing a dataset}
}
\details{
  The method of L-Moments estimates of the parameters of the quantile-based
  skew logistic distribution are:
  \deqn{\hat\alpha=L_1 - 6 L_3}{alpha = L1 - 6L3}
  \deqn{\hat\beta = 2 L_2}{beta = 2 L2}
  \deqn{\hat\delta= \frac 12 \left( 1 + 3\tau_3\right)}{delta = 0.5*(1+3*tau3)}
  Note that L3 in the \eqn{\hat\alpha}{alpha} estimate is the 3rd L-Moment, not the 3rd L-Moment 
  ratio (\eqn{\tau_3 = L3/L2}{tau3 = L3/L2}).
  
  \code{fit.sld.lmom} uses the \code{\link[lmom]{samlmu}} function to
  calculate the sample L moments, then \code{fit.sld.lmom.given} to calculate
  the estimates.
}
\value{
  A vector of length 3, with the estimated parameters, 
  \eqn{\alpha}{alpha}, \eqn{\beta}{beta} and \eqn{\delta}{delta}.
}
\references{
van Staden, P.J. and King, R.A.R. (2014) The quantile-based skew logistic distribution. \emph{Under Review}

van Staden, Paul J. 2013 \emph{Modeling of generalized families of probability distribution in the quantile statistical universe}. 
PhD thesis, University of Pretoria.
\url{http://hdl.handle.net/2263/40265}

\url{http://tolstoy.newcastle.edu.au/rking/SLD/SLD.html}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/} and Paul van Staden}

\seealso{\code{\link{sld}}}
\examples{
generated.data <- rsl(300,c(0,1,.4))
estimate1 <- fit.sld.lmom(generated.data)
estimate2 <- fit.sld.lmom.given(c(0,1,.3))
}
\keyword{estimation}

