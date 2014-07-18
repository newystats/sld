# Incomplete
qqsl <- function (y=NULL,parameters1,parameters2=NULL,abline=TRUE,
      granularity.for.2.dists=4000,use.endpoints=FALSE,...)
{
  # Check the parameter values are OK
  if(!sl.check.pars(parameters1)) {
    stop(paste("The parameter values", paste(parameters1,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Use the names for the parameters
  alpha1 <- parameters1[1]
  p.beta1 <- parameters1[2]
  delta1 <- parameters1[3]
if(is.null(parameters2)) { # No second set of parameters - compare to a dataset
	n <- length(y)
	}
else { # A second set of parameters - compare two distributions
	if(!is.null(y)){stop("Can\'t produce a QQ plot comparing 2 sets of parameter values and also a dataset")}
	if(!sl.check.pars(parameters2)) {
	  stop(paste("The (2nd set of) parameter values", paste(parameters2,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n")) 
	}
	n <- granularity.for.2.dists
}
if (use.end.points){
  u <- seq(from=0,to=1,length.out=n)
} else {
  # uses ideal depths - this is quantile type 8 - perhaps I should use the quantile function
  u <- seq(from = (2/3)/(n + 1/3), by = 1/(n+ 1/3), length = n)
}

## Done up to here
if(is.null(lambda.pars2)) { # No second set of parameters - compare to a dataset
	Theoretical.Quantiles <- qgl(u, lambda1=lambdas, param=param)
	Data <- y
	ret <- qqplot(Theoretical.Quantiles,Data,...)
	}
else	{ # A second set of parameters - compare two distributions
	Quantiles1 <- qgl(u, lambda1=lambdas, param=param)
	Quantiles2 <- qgl(u, lambda1=lambdas2, param=param2)
	ret <- qqplot(Quantiles1,Quantiles2,...)
	}
if(abline) { 
	abline(0,1)
	}
invisible(ret)
}
