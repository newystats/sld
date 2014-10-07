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
fit.sld.lmom.given(lmoms,n=NULL)
fit.sld.lmom(data)
}
\arguments{
  \item{lmoms}{A vector of length 3, containing the first and second (sample)
  L Moments and the 3rd (sample) L Moment ratio (\eqn{\tau_3}{tau3 })}
  \item{n}{The sample size}
  \item{data}{A vector containing a dataset}
}
\details{
  The method of L-Moments estimates of the parameters of the quantile-based
  skew logistic distribution are:
  \deqn{\hat\alpha=L_1 - 6 L_3}{alpha = L1 - 6L3}
  \deqn{\hat\beta = 2 L_2}{beta = 2 L2}
  \deqn{\hat\delta= \frac 12 \left( 1 + 3\tau_3\right)}{delta = 0.5*(1+3*tau3)}
  Note that \eqn{L_3}{L3} in the \eqn{\hat\alpha}{alpha} estimate is the 3rd L-Moment, not the 3rd L-Moment 
  ratio (\eqn{\tau_3 = L_3/L_2}{tau3 = L3/L2}).
  
  \code{fit.sld.lmom} uses the \code{\link[lmom]{samlmu}} function (from 
  the \code{lmom} package) to calculate the sample L moments, then 
  \code{fit.sld.lmom.given} to calculate the estimates.
}
\value{
  If the sample size is unknown (via using \code{\link{fit.sld.lmom.given}} and not specifying the sample size), a vector of length 3, with the estimated parameters, 
  \eqn{\hat\alpha}{alpha}, \eqn{\hat\beta}{beta} and \eqn{\hat\delta}{delta}.
  
  If the sample size is known, a 3 by 2 matrix.  The first column 
contains the estimated parameters, 
  \eqn{\hat\alpha}{alpha}, \eqn{\hat\beta}{beta} and \eqn{\hat\delta}{delta},
  and the second column provides asymptotic standard errors for these.
  
Note that if \eqn{|\tau_3| > \frac 13}{abs(tau3) > 1/3}, 
\eqn{\hat\delta}{delta hat} is beyond its allowed value of [0,1]
and the function returns an error.  
Values of \eqn{|\tau_3|}{abs(tau3)}, beyond 
\eqn{\frac 13}{1/3} correspond to distributions with 
greater skew than the exponential / reflected exponential, 
which form the limiting cases of the skew logistic distribution.
}
\references{
van Staden, P.J. and King, Robert A.R. (2015) \emph{The quantile-based skew logistic distribution}, Statistics and Probability Letters \bold{96} 109--116.

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
estimate2 <- fit.sld.lmom.given(c(0,1,.3),n=300)
}
\keyword{estimation}

