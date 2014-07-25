fit.sld.lmom.given <- function(lmoms){
  lmomest <- c(lmoms[1] - 6*(lmoms[3]*lmoms[2]),2*lmoms[2],0.5*(1+3*lmoms[3]))
  names(lmomest) <- c("alpha","beta","delta")
  lmomest
}

fit.sld.lmom <- function(data){
  fit.sld.lmom.given(lmom.sample(data,3))
}

lmom.sample <- function(data,max.mom=3){
  LM <- samlmu(data,max.mom)
  LM
}