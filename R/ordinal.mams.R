ordinal.mams <- function(prob=c(0.35, 0.4, 0.25), or=2, or0=1.2, K=4, J=2, alpha=0.05, power=0.9,
                         r=1:2, r0=1:2, ushape="obf", lshape="fixed", ufix=NULL, lfix=0, nstart=1,
                         nstop=NULL, sample.size=TRUE, N=20, 
                         parallel = FALSE, ncores = NULL){
  
  if(sum(prob)!=1){
    stop("The elements of prob must sum to one.")
  }
  if(or0 < 1){
    stop("The uninteresting effect must be 1 or larger.")
  }
  if(or <= or0){
    stop("The interesting effect must be larger than the uninteresting one.")
  }
  
  q <- (1 - sum(prob^3))/3
  sigma <- 1/sqrt(q)
  
  p <- pnorm(log(or)/sqrt(2 * sigma^2))
  p0 <- pnorm(log(or0)/sqrt(2 * sigma^2))
  
  mams(K=K, J=J, alpha=alpha, power=power, r=r, r0=r0, p=p, p0=p0, delta=NULL, delta0=NULL, sd=NULL,
       ushape=ushape, lshape=lshape, ufix=ufix, lfix=lfix, nstart=nstart, nstop=nstop, sample.size=sample.size,
       N=N, type="ordinal", parallel = parallel, ncores = ncores)
  
}