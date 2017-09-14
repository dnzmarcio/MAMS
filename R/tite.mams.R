tite.mams <- function(hr=1.5, hr0=1.1, K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2,
                      ushape="obf", lshape="fixed", ufix=NULL, lfix=0, nstart=1,
                      nstop=NULL, sample.size=TRUE, N=20){
 
  if (hr0 < 1) {
    stop("The uninteresting effect must be 1 or larger.")
  }
  
  if (hr <= hr0) {
    stop("The interesting effect must be larger than the uninteresting one.")
  }
  
  p <- pnorm(log(hr)/sqrt(2))
  p0 <- pnorm(log(hr0)/sqrt(2))
  
  mams(K=K, J=J, alpha=alpha, power=power, r=r, r0=r0, p=p, p0=p0, delta=NULL, delta0=NULL, sd=NULL,
       ushape=ushape, lshape=lshape, ufix=ufix, lfix=lfix, nstart=nstart, nstop=nstop, sample.size=sample.size,
       N=N, type="tite")
  
}