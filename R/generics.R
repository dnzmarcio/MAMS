MAMSNews <- function() file.show(system.file("NEWS", package="MAMS"))

print.MAMS <- function (x, digits=max(3, getOption("digits") - 4), ...) {

  cat(paste("Design parameters for a ", x$J, " stage trial with ", x$K, " treatments\n\n",sep=""))

  if(!is.na(x$power)){
    res <- matrix(NA,nrow=2,ncol=x$J)
    colnames(res)<-paste("Stage",1:x$J)
    rownames(res) <- c("Cumulative sample size per stage (control):", "Cumulative sample size per stage (active):")

    res[1,] <- x$n*x$rMat[1,]
    res[2,] <- x$n*x$rMat[2,]

    print(res)
  
    cat(paste("\nMaximum total sample size: ", x$N,"\n\n"))

  }

  res <- matrix(NA,nrow=2,ncol=x$J)
  colnames(res)<-paste("Stage",1:x$J)
  rownames(res) <- c("Upper bound:", "Lower bound:")
  res[1,] <- round(x$u,digits)
  res[2,] <- round(x$l,digits)
  
  print(res)

}


summary.MAMS<-function(object, digits=max(3, getOption("digits") - 4), ...){

  cat(paste("Design parameters for a ", object$J, " stage trial with ", object$K, " treatments\n\n",sep=""))

  if(!is.null(object$n)){
    res <- matrix(NA,nrow=2,ncol=object$J)
    colnames(res)<-paste("Stage",1:object$J)
    rownames(res) <- c("Cumulative sample size per stage (control):", "Cumulative sample size per stage (active):")

    res[1,] <- object$n*object$rMat[1,]
    res[2,] <- object$n*object$rMat[2,]

    print(res)
  
    cat(paste("\nMaximum total sample size: ", object$N,"\n\n"))
  }

  res <- matrix(NA,nrow=2,ncol=object$J)
  colnames(res)<-paste("Stage",1:object$J)
  rownames(res) <- c("Upper bound:", "Lower bound:")
  res[1,] <- round(object$u,digits)
  res[2,] <- round(object$l,digits)
  
  print(res)
}

plot.MAMS <- function (x, col=NULL, pch=NULL, lty=NULL, main=NULL, xlab="Analysis", ylab="Test statistic", ylim=NULL, type=NULL, ...) {

  if(is.null(type))type<-"p"
  if(is.null(pch))pch<-1
  if(is.null(col))col<-1
  if(is.null(lty))lty<-2
  if(is.null(ylim)){
    r<-range(x$l,x$u)
    ylim <- c(r[1]-diff(r)/6,r[2]+diff(r)/6)
  }

  matplot(1:x$J, cbind(x$l,x$u), type=type, pch=pch, col=col, ylab=ylab, xlab=xlab, ylim=ylim, axes=FALSE, ...)
  mtext(1:x$J,side=1,at=1:x$J)
#  axis(side=2)
  axis(side=2,at=seq(-10,10,1))
  lines(x$u,lty=lty)
  lines(x$l[1:(x$J)],lty=lty)

}


print.MAMS.sim <- function (x, digits=max(3, getOption("digits") - 4), ...) {

  cat(paste("Simulated error rates based on ", x$nsim," simulations\n\n",sep=""))

  res <- matrix(NA,nrow=4,ncol=1)

  res[1,1] <- round(x$typeI,digits)
  res[2,1] <- round(x$power,digits)
  res[3,1] <- round(x$prop.rej,digits)
  res[4,1] <- round(x$exss,digits)

  if(length(x$ptest)==1){
  rownames(res) <- c("Prop. rejecting at least 1 hypothesis:", "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
                      paste("Prop. rejecting hypothesis ",x$ptest,":",sep=""),"Expected sample size:")
  }else{
  rownames(res) <- c("Prop. rejecting at least 1 hypothesis:", "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
                      paste("Prop. rejecting hypotheses ",paste(as.character(x$ptest),collapse=" or "),":",sep=""),"Expected sample size:")
  }
  colnames(res)<-""
  
  print(res)

}

summary.MAMS.sim<-function(object, digits=max(3, getOption("digits") - 4), ...){

  print(object)

}
