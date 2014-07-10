qsl <- function(p,parameters)
{
# Check the parameter values are OK
if(!sl.check.pars(parameters)) {
        stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
	}
# Use the names for the parameters
alpha <- parameters[1]
p.beta <- parameters[2]
delta <- parameters[3]
# Do something sensible with stupid ps
outside.range <- !as.logical((p <= 1) * (p >= 0))
u <- p[!outside.range]
# Special cases at delta=0,1 require 
if (delta == 0){ # These special cases are here in case u=1 when delta is 0 and lambda is negative see delta zero question in Robert Kings gld package notes
  # reflected exponential
  quants <- alpha + p.beta * ( log(u) )
  } else {
    if (delta ==1) { # exponential
      quants <- alpha - (p.beta * log(1-u) ) # beta * -1 * log(1-u)
      } else { # skew logistic
        quants <- alpha + p.beta * ( (1-delta)*log(u) - delta*log(1-u))
      }
  }
result <- p * NaN
result[!outside.range] <- quants
result
}

dqgl <- function(p,parameters){
  # This is the density quantile function of the skew logistic distribution
  # Check the parameter values are OK
  if(!sl.check.pars(parameters)) {
    stop(paste("The parameter values", paste(parameters,collapse=" "),"\ndo not produce a proper skew logistic distrbution.\nNote that beta must be positive and delta needs to be in the range [0,1]\n"))
  }
  # Use the names for the parameters
  alpha <- parameters[1]
  p.beta <- parameters[2]
  delta <- parameters[3]
  # Do something sensible with stupid ps  
  outside.range <- !as.logical((p<=1)*(p>=0))
  # u gets only the probabilities in [0,1]
  dens <- p*0
### Fix this with somtehign like the qsl function
  u <- p[!outside.range]	
dens <-  lambdas[2]/(lambdas[3] * (p^(lambdas[3] -1)) + lambdas[4] * ((1 - p)^(lambdas[4] -1)))
dens
}

.dqgl.fmkl <- function(p,lambdas)
{
# Check the values are OK)
if(!gl.check.lambda(lambdas,param="fkml",vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas,collapse=" "),
"\ndo not produce a proper distribution with the FMKL",
"parameterisation - see \ndocumentation for gl.check.lambda"))
	}
outside.range <- !as.logical((p<=1)*(p>=0))
# u gets only the probabilities in [0,1]
u <- p[!outside.range]
# The density is given by 1/Q'(u)
dens <- lambdas[2]/(p^(lambdas[3] - 1) + (1 - p)^(lambdas[4] - 1))
dens
}

.dqgl.fm5 <- function(p,lambdas)
{
# Check the values are OK)
if(!gl.check.lambda(lambdas,param="fm5",vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas,collapse=" "),
"\ndo not produce a proper distribution with the FM5",
"parameterisation - see \ndocumentation for gl.check.lambda"))
	}
outside.range <- !as.logical((p<=1)*(p>=0))
# u gets only the probabilities in [0,1]
u <- p[!outside.range]
# The density is given by 1/Q'(u)
dens <- lambdas[2]/((1-lambdas[5])*(u^(lambdas[3] - 1)) + (1+lambdas[5])*((1 - u)^(lambdas[4] - 1)) )
dens
}

.dqgl.vsk <- function(p,lambdas)
{
	# lambdas is a parameter containing (alpha,beta,lambda,delta)
	alpha <- lambdas[1]
	p.beta <- lambdas[2]
	delta <- lambdas[3]
	lambda <- lambdas[4]
  outside.range <- !as.logical((p <= 1) * (p >= 0))
	u <- p[!outside.range]
	if (lambda == 0){
		dens <- u*(1-u)/(p.beta* (delta*u + (1-delta)*(1-u)))
	} else {
		dens <- 1/(p.beta * ( (1-delta)*(u^(lambda -1)) + delta*( (1-u)^(lambda -1))) )
	}
	result <- p * 0
	result[!outside.range] <- dens
	result
}
