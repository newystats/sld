\name{SkewLogisticDistribution}
\alias{sld}
\alias{dqsl}
\alias{dsl}
\alias{psl}
\alias{qsl}
\alias{rsl}
\title{The quantile-based Skew Logistic Distribution}
\description{Density, density quantile, distribution and quantile 
functions and random generation for the quantile-based skew logistic 
distribution.}
\usage{
dsl(x,parameters,inverse.eps=.Machine$double.eps,max.iterations=500)
dqsl(p,parameters)
psl(q,parameters,inverse.eps=.Machine$double.eps,max.iterations=500)
qsl(p,parameters)
rsl(n,parameters)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
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
  \item{inverse.eps}{Accuracy of calculation for the numerical determination 
  of \eqn{F(x)}, defaults to \code{.Machine$double.eps}}
  \item{max.iterations}{Maximum number of iterations in the numerical 
	determination of \eqn{F(x)}, defaults to 500}
}
\details{
  The quantile-based skew logistic distribution is a generalisation of the 
  logistic distribution, defined by its quantile funtion, Q(u), the inverse 
  of the distribution function.
\deqn{Q(u)= \alpha + \beta \left[ (1- \delta) \log (u)  - \delta (\log (1-u)) \right] )}{Q(u) = alpha + beta ( (1 - delta)*(log(u)) - delta * (log(1-u)) ) }
for \eqn{\beta > 0}{beta >0} and \eqn{0 \leq \delta \leq 1}{0 <= delta <= 1}.

  The distribution was first used by Gilchrist (2000) in the book 
\emph{Statistical Modelling with Quantile Functions}.  Full details of 
the properties of the distributions, including
moments, L-moments and estimation via L-Moments are given in 
van Staden and King (2015).
  
The distribution is defined by its quantile function and its distribution and
density functions do not exist in closed form (except for some special cases).
Accordingly, the results from \code{psl} and \code{dsl} are the result of
numerical solutions to the quantile function, using the Newton-Raphson method.
Since the density quantile function, \eqn{f(Q(u))}, does exist, an 
additional function, \code{dqsl}, computes this.

  The distribution has closed form method of L-Moment estimates 
(see \code{\link{fit.sld.lmom}} for details).  The 4th L-Moment ratio of the
the distribution is constant \eqn{\tau_4 = \frac{1}{6}}{tau4 = 1/6} for 
all values of \eqn{\delta}{delta}.  
The 3rd L-Moment ratio of the distribution is restricted to
\eqn{\frac{-1}{3} \leq \tau_3 \leq \frac{1}{3}}{-1/3 <= tau3 <= 1/3}, being 
the the 3rd L-moment ratio values of the reflected exponential and 
the exponential distributions respectively.
}
\value{
  \code{dsl} gives the density (based on the quantile density and a 
  numerical solution to \eqn{Q(u)=x}{Q(u)=x}), 

  \code{dqsl} gives the density quantile,
  
  \code{psl} gives the distribution function (based on a numerical
	solution to Q(u)=x and \code{\link{dqsl}}

  \code{qsl} gives the quantile function, and

  \code{rsl} generates random deviates.
}
\references{
Gilchrist, W.G. (2000)
\emph{Statistical Modelling with Quantile Functions}
Chapman & Hall, print 978-1-58488-174-2, e-book 978-1-4200-3591-9.

van Staden, P.J. and King, Robert A.R. (2015) \emph{The quantile-based skew logistic distribution}, Statistics and Probability Letters \bold{96} 109--116.
\doi{10.1016/j.spl.2014.09.001}

van Staden, Paul J. 2013 \emph{Modeling of generalized families of probability distribution in the quantile statistical universe}. 
PhD thesis, University of Pretoria.
\url{http://hdl.handle.net/2263/40265}

\url{https://github.com/newystats/sld}
}
\author{Robert King, \email{robert.king.newcastle@gmail.com}, \url{https://github.com/newystats/}}

\examples{
qsl(seq(0,1,0.02),c(0,1,0.123))
psl(seq(-2,2,0.2),c(0,1,.1),inverse.eps=1e-10)
rsl(21,c(3,2,0.3))}
\keyword{distribution}
