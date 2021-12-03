# sld
R package for the Quantile based skew logistic distribution.  This is available on CRAN as https://cran.r-project.org/package=sld

The quantile-based skew logistic distribution is a skewed distribution, produced by a skewing function of Gilchrist.  It is defined by its quantile function;

Q(u) = alpha + beta ( (1 - delta)*(log(u)) - delta * (log(1-u)) )

for beta >0 and 0 <= delta <= 1.

The distribution was first used by Gilchrist (2000) in the book Statistical Modelling with Quantile Functions. 
Full details of the properties of the distributions, including moments, L-moments and estimation via L-Moments are given in van Staden and King (2015).

The distribution is defined by its quantile function and its distribution and density functions do not exist in closed form (except for some special cases). 
Accordingly, the results from psl and dsl are the result of numerical solutions to the quantile function, using the Newton-Raphson method. 
Since the density quantile function, f(Q(u)), does exist, an additional function, dqsl, computes this.

The distribution has closed form method of L-Moment estimates (see fit.sld.lmom for details). 
The 4th L-Moment ratio of the the distribution is constant tau4 = 1/6 for all values of delta. 
The 3rd L-Moment ratio of the distribution is restricted to -1/3 <= tau3 <= 1/3, being the the 3rd L-moment ratio values of 
the reflected exponential and the exponential distributions respectively.

The package provides random number generation, probability density function, cumulative distribution function, plotting of these functions, estimation using L-Moments, 
standard errors.

## References

Gilchrist, W.G. (2000) Statistical Modelling with Quantile Functions Chapman & Hall, print 978-1-58488-174-2, e-book 978-1-4200-3591-9.

van Staden, P.J. and King, Robert A.R. (2015) The quantile-based skew logistic distribution, Statistics and Probability Letters 96 109–116. http://dx.doi.org/10.1016/j.spl.2014.09.001

van Staden, Paul J. 2013 Modeling of generalized families of probability distribution in the quantile statistical universe. PhD thesis, University of Pretoria. http://hdl.handle.net/2263/40265
