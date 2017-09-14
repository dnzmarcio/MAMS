new.bounds <- function(K=3, J=2, alpha=0.05, nMat=matrix(c(10,20), nrow=2, ncol=4), u=NULL, l=NULL, ushape="obf", lshape="fixed", ufix=NULL, lfix=0, N=20){


  #require(mvtnorm) ## the function pmvnorm is required to evaluate multivariate normal probabilities

  ############################################################################################
  ## 'mesh' creates the points and respective weights to use in the outer quadrature integrals
  ## you provide one-dimensional set of points (x) and weights (w) (e.g. midpoint rule, gaussian quadrature, etc)
  ## and the number of dimensions (d) and it gives d-dimensional collection of points and weights
  ###############################################################################################

  mesh<-function(x,d,w=1/length(x)+x*0){
    n<-length(x)
    W<-X<-matrix(0,n^d,d)
    for (i in 1:d){
      X[,i]<-x
      W[,i]<-w
      x<-rep(x,rep(n,length(x)))
      w<-rep(w,rep(n,length(w)))
    }
    w<-exp(rowSums(log(W)))
    list(X=X,w=w)
  }

  prodsum<-function(x,l,u,R,r0,r0diff,J,K,Sigma){
    ########################################################################################################
    ## x is vector of dummy variables ( the t_j in gdt paper ), l and u are boundary vectors
    ########################################################################################################
    int<-prod(sapply(x,dnorm))
    for (k in 1:K){
      L<-sqrt(R[1,k])/r0[1]*sqrt(r0diff[1])*x[1]+l[1]*sqrt(1+R[1,k]/r0[1])
      U<-sqrt(R[1,k])/r0[1]*sqrt(r0diff[1])*x[1]+u[1]*sqrt(1+R[1,k]/r0[1])
      insum<-pnorm(L)
      if(J>1){
        for (j in 2:J){
          L[j]<-sqrt(R[j,k])/r0[j]*(sqrt(r0diff[1:j])%*%x[1:j])+l[j]*sqrt(1+R[j,k]/r0[j])
          U[j]<-sqrt(R[j,k])/r0[j]*(sqrt(r0diff[1:j])%*%x[1:j])+u[j]*sqrt(1+R[j,k]/r0[j])
          insum<-insum+pmvnorm(lower=c(L[1:j-1],-Inf),upper=c(U[1:j-1],L[j]),sigma=Sigma[1:j,k,1:j])[1]
        }
      }
      int<-int*insum
    }
    return(int)
  }

  #############################################################################################################
  ##  'typeI' performs the outer quadrature integral in type I error equation using 'mesh' and 'prodsum'
  ##  and calculates the difference with the nominal alpha.
  ##  The accuracy of the quadrature integral will depend on the choice of points and weights.
  ##  Here, the number of points in each dimension is an input, N.
  ##  The midpoint rule is used with a range of -6 to 6 in each dimension. 
  #############################################################################################################

  typeI<-function(C,alpha,N,R,r0,r0diff,J,K,Sigma,u,l,ushape,lshape,lfix=NULL,ufix=NULL){

    ## number of stages already done
    if(!is.null(l)){
      j <- length(l)
    }else{
      j <- 0
    }
    ########################################################################
    ## the form of the boundary constraints are determined as functions of C. 
    ######################################################################## 
    Jleft <- J-j
    if(Jleft==1){
      ind <- J
    }else{
      ind <- (j+1):J
    }
    if(!is.function(ushape)){
      if (ushape=='obf'){
        ub<-c(u,C*sqrt(J/(1:J))[ind])
      }
      else if (ushape=='pocock'){
        ub<-c(u,rep(C,Jleft))
      }
      else if (ushape=='fixed'){
        ub<-c(u,rep(ufix,Jleft-1),C)
      } 
      else if (ushape=='triangular') {
        #ub<-c(u,C*(1+(1:J)/J)/sqrt(1:J)[(j+1):J])
        ub<-c(u,(C*(1+(1:J)/J)/sqrt(1:J))[ind])
      }
    }else{
      ub <- c(u,C*ushape(J)[(j+1):J])
    }

    if(!is.function(lshape)){
      if (lshape=='obf'){
        lb<- c(l,-C*sqrt(J/(1:(J-1)))[(j+1):(J-1)],ub[J])
      }
      else if (lshape=='pocock'){
        lb<-c(l,rep(-C,Jleft-1),ub[J])
      }
      else if (lshape=='fixed'){
        lb<-c(l,rep(lfix,Jleft-1),ub[J])
      } else if (lshape=='triangular') {
        #lb<-c(l,-C*(1-3*(1:J)/J)/sqrt(1:J)[(j+1):J])
        lb<-c(l,(-C*(1-3*(1:J)/J)/sqrt(1:J))[ind])
      }
    }else{
      lb <- c(l,C*lshape(J)[(j+1):(J-1)],ub[J])
    }

    mmp<-mesh((1:N-.5)/N*12-6,J,rep(12/N,N))
    evs<-apply(mmp$X,1,prodsum,l=lb,u=ub,R=R,r0=r0,r0diff=r0diff,J=J,K=K,Sigma=Sigma)
    truealpha<-1-mmp$w%*%evs

    return(truealpha-alpha)
  }

  # check input

  if(K%%1>0 | J%%1>0){stop("K and J need to be integers.")}
  if(K < 1 | J < 1){stop("The number of stages and treatments must be at least 1.")}
  if(N<=3){stop("Number of points for integration by quadrature to small or negative.")}
  if(N>3 & N<=10){warning("Number of points for integration by quadrature is small which may result in inaccurate solutions.")}
  if(alpha<0 | alpha>1){stop("Error rate not between 0 and 1.")}

  if(!is.function(ushape)){
    if(!ushape%in%c("pocock","obf","triangular","fixed")){stop("Upper boundary does not match the available options")}
    if(ushape=="fixed" & is.null(ufix)){stop("ufix required when using a fixed upper boundary shape.")}
  }else{
    b <- ushape(J)
    if(!all(sort(b,decreasing=TRUE)==b)){stop("Upper boundary shape is increasing")}
  }
  if(!is.function(lshape)){
   if(!lshape%in%c("pocock","obf","triangular","fixed")){stop("Lower boundary does not match the available options")}
   if(lshape=="fixed" & is.null(lfix)){stop("lfix required when using a fixed lower boundary shape.")}
  }else{
    b <- lshape(J)
    if(!all(sort(b,decreasing=FALSE)==b)){stop("Lower boundary shape is decreasing")}
  }
  
  
  if(!all(diff(nMat)>=0)){stop("Total sample size per arm can not decrease between stages")}
  if(ncol(nMat)!=K+1){stop("Number of columns in nMat not equal to K+1")}
  if(nrow(nMat)!=J){stop("Number of rows in nMat not equal to J")}
  if(length(l)!=length(u)){stop("Length of l must be the same as length of u")}
  if(length(l)>(J-1)){stop("Maximum number of stages, J, greater or equal to length of boundary vector")}


  r0 <- nMat[,1]/nMat[1,1]
  R <-  nMat[,-1]/nMat[1,1]

  ####################################################################
  ## Create the variance covariance matrix from allocation proportions:
  ####################################################################

  bottom<-array(R,dim=c(J,K,J))
  top<-array(rep(t(R),rep(J,K*J)),dim=c(J,K,J))
  for (i in 1:K){
    top[,i,][upper.tri(top[,i,])]<-t(top[,i,])[upper.tri(top[,i,])]
    bottom[,i,][upper.tri(bottom[,i,])]<-t(bottom[,i,])[upper.tri(bottom[,i,])] 
  }
  Sigma<-sqrt(top/bottom)

  ###############################################################################
  ## Create r0diff: the proportion of patients allocated to each particular stage
  ###############################################################################

  r0lag1<-c(0,r0[1:J-1])
  r0diff<-r0-r0lag1 


  ################################
  ## Find boundaries using 'typeI'
  ################################

  uJ <- NULL
  try(uJ<-uniroot(typeI,c(0,5),alpha=alpha,N=N,R=R,r0=r0,r0diff=r0diff,J=J,K=K,Sigma=Sigma,u=u,l=l,ushape=ushape,lshape=lshape,lfix=lfix,ufix=ufix,tol=0.001)$root,silent=TRUE) 
 
  if(is.null(uJ)) {
    stop("No solution found")
  }


  ## number of stages already done
  if(!is.null(l)){
    j <- length(l)
  }else{
    j <- 0
  }
  ## number of stages left
  Jleft <- J-j
  if(Jleft==1){
    ind <- J
  }else{
    ind <- (j+1):J
  }

  ########################################################################
  ## the form of the boundary constraints are determined as functions of C. 
  ######################################################################## 

  if(!is.function(ushape)){
    if (ushape=='obf'){
      ub<-c(u,uJ*sqrt(J/(1:J))[ind])
    }
    else if (ushape=='pocock'){
      ub<-c(u,rep(uJ,Jleft))
    }
    else if (ushape=='fixed'){
      ub<-c(u,rep(ufix,Jleft),uJ)
    } 
    else if (ushape=='triangular') {
      ub<-c(u,(uJ*(1+(1:J)/J)/sqrt(1:J))[ind])
    }
  }else{
    ub <- c(u,uJ*ushape(J)[(j+1):J])
  }

  if(!is.function(lshape)){
    if (lshape=='obf'){
      lb<- c(l,-uJ*sqrt(J/(1:(J-1)))[(j+1):(J-1)],ub[J])
    }
    else if (lshape=='pocock'){
      lb<-c(l,rep(-uJ,Jleft-1),ub[J])
    }
    else if (lshape=='fixed'){
      lb<-c(l,rep(lfix,Jleft-1),ub[J])
    } else if (lshape=='triangular') {
      lb<-c(l,(-uJ*(1-3*(1:J)/J)/sqrt(1:J))[ind])
    }
  }else{
    lb <- c(l,uJ*lshape(J)[(j+1):(J-1)],ub[J])
  }



  #########################################################
  ## Find alpha.star
  #########################################################
  alpha.star <- numeric(J)
  alpha.star[1] <- typeI(ub[1], alpha = 0, N = N, R = t(as.matrix(R[1,])), r0 = r0[1], r0diff = r0diff[1], J = 1, K = K, Sigma = Sigma, u = NULL, l = NULL, ushape = "fixed", lshape = "fixed", lfix = NULL, ufix = NULL)
  if (J > 1){
      for (j in 2:J){
          alpha.star[j] <- typeI(ub[j], alpha = 0, N = N, R = R[1:j,], r0 = r0[1:j], r0diff = r0diff[1:j], J = j, K = K, Sigma = Sigma, u = NULL, l = NULL, ushape = "fixed", lshape = "fixed", lfix = lb[1:(j - 1)], ufix = ub[1:(j - 1)])
      }
  }

  

  res <- NULL
  res$l <- lb  
  res$u <- ub
  res$n <- nMat[1,1]
  res$rMat <- rbind(r0,t(R)) ## allocation ratios
  res$N <- sum(res$rMat[,J]*res$n) ## maximum total sample size

  res$K <- K
  res$J <- J
  res$alpha <- alpha
  res$alpha.star <- alpha.star
  res$power <- NA

  
  class(res)<-"MAMS"

  return(res)

}

