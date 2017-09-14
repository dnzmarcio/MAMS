stepdown.mams <- function(nMat = matrix(c(10, 20), nrow=2, ncol=4), alpha.star = c(0.01, 0.025), lb = 0, selection = "all.promising"){

    # checking input parameters
    if (!all(diff(nMat) >= 0)) {stop("total sample size per arm cannot decrease between stages.")}
    J <- dim(nMat)[1]
    K <- dim(nMat)[2] - 1
    if ((J != 2) && (J != 3)) {stop("number of stages must be 2 or 3")}
    if (K < 2) {stop("must have at least two experimental treatments")}
    if (length(alpha.star) != J) {stop("length of error spending vector must be same as number of stages")}
    if (!all(diff(alpha.star) >= 0)) {stop("cumulative familywise error must increase.")}
    if (length(lb) != J - 1) {stop("lower boundary must be specified at all analysis points except the last")}
    match.arg(selection,c("all.promising","select.best"))
    
  get.hyp <- function(n){ # find the nth intersection hypothesis (positions of 1s in binary n)
    indlength = ceiling(log(n)/log(2)+.0000001)
    ind = rep(0,indlength)
    newn=n
        
    for (h in seq(1,indlength)){
      ind[h] = (newn/(2^(h-1))) %% 2
      newn = newn - ind[h]*2^(h-1)
    }
    seq(1,indlength)[ind==1]
  }

  create.block <- function(control.ratios = 1:2, active.ratios = matrix(1:2, 2, 3)){ # for argument c(i,j) this gives covariance between statistics in stage i with statistics in stage j
    K <- dim(active.ratios)[2]
    block <- matrix(NA, K, K)
    for(i in 1:K){
      block[i, i] <- sqrt(active.ratios[1, i] * control.ratios[1] * (active.ratios[2, i] + control.ratios[2]) / (active.ratios[1, i] + control.ratios[1]) / active.ratios[2, i] / control.ratios[2])
    }
    for (i in 2:K){
      for (j in 1:(i - 1)){
        block[i, j] <- sqrt(active.ratios[1, i] * control.ratios[1] * active.ratios[2, j] / (active.ratios[1, i] + control.ratios[1]) / (active.ratios[2, j] + control.ratios[2]) / control.ratios[2])
        block[j, i] <- sqrt(active.ratios[1, j] * control.ratios[1] * active.ratios[2, i] / (active.ratios[1, j] + control.ratios[1]) / (active.ratios[2, i] + control.ratios[2]) / control.ratios[2])
      }
    }
    block
  }

  create.cov.matrix <- function(control.ratios = 1:2, active.ratios = matrix(1:2, 2, 3)){ # create variance-covariance matrix of the test statistics

    J <- dim(active.ratios)[1]
    K <- dim(active.ratios)[2]

    cov.matrix <- matrix(NA, J * K, J * K)
    for (i in 1:J){
      for (j in i:J){
        cov.matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)] <- create.block(control.ratios[c(i, j)], active.ratios[c(i, j), ])
        cov.matrix[((j - 1) * K + 1):(j * K), ((i - 1) * K + 1):(i * K)] <- t(cov.matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)])
      }
    }
    cov.matrix
  }

  get.path.prob <- function(surviving.subset1, surviving.subset2 = NULL, cut.off, treatments, cov.matrix, lb, upper.boundary, K, stage){ # find the probability that no test statistic crosses the upper boundary + only treatments in surviving_subsetj reach the jth stage
    treatments2 <- treatments[surviving.subset1]
    if (stage == 2){
      lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)))
      lower[surviving.subset1] <- lb[1]

      upper <- c(rep(lb[1], length(treatments)), rep(cut.off, length(treatments2)))
      upper[surviving.subset1] <- upper.boundary[1]

      return(pmvnorm(lower = lower, upper = upper, sigma = cov.matrix[c(treatments, K + treatments2), c(treatments, K + treatments2)])[1])
    }
    treatments3 <- treatments2[surviving.subset2]

    lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)), rep(-Inf, length(treatments3)))
    lower[surviving.subset1] <- lb[1]
    lower[length(treatments) + surviving.subset2] <- lb[2]

    upper <- c(rep(lb[1], length(treatments)), rep(lb[2], length(treatments2)), rep(cut.off, length(treatments3)))
    upper[surviving.subset1] <- upper.boundary[1]
    upper[length(treatments) + surviving.subset2] <- upper.boundary[2]

    pmvnorm(lower = lower, upper = upper, sigma = cov.matrix[c(treatments, K + treatments2, 2 * K + treatments3), c(treatments, K + treatments2, 2 * K + treatments3)])[1]
  }


  rejection.paths <- function(selected.treatment, cut.off, treatments, cov.matrix, lb, upper.boundary, K, stage){ # for the "select.best" method, find the probability that "select.treatment" is selected and subsequently crosses the upper boundary

    contrast <- diag(-1, K + stage - 1)
    contrast[1:K, selected.treatment] <- 1
    for (i in 1:(stage - 1)) contrast[K + i, K + i] <- 1
  
    bar.cov.matrix <- contrast %*% cov.matrix[c(1:K, 1:(stage - 1) * K + selected.treatment), c(1:K, 1:(stage - 1) * K + selected.treatment)] %*% t(contrast)

    lower <- c(rep(0, length(treatments)), cut.off)
    if (stage > 2) lower <- c(rep(0, length(treatments)), lb[2:(stage - 1)], cut.off)
    lower[which(treatments == selected.treatment)] <- lb[1]

    upper <- c(rep(Inf, length(treatments)), Inf)
    if (stage > 2) upper <- c(rep(Inf, length(treatments)), upper.boundary[2:(stage - 1)], Inf)
    upper[which(treatments == selected.treatment)] <- upper.boundary[1]

    pmvnorm(lower = lower, upper = upper, sigma = bar.cov.matrix[c(treatments, K + 1:(stage - 1)), c(treatments, K + 1:(stage - 1))])[1]
    
  }

  excess.alpha <- function(cut.off, alpha.star, treatments, cov.matrix, lb, upper.boundary, selection, K, stage){ # for "all.promising" rule, this gives the cumulative typeI error for 'stage' stages

      # for "select.best" rule, this gives the Type I error spent at the 'stage'th stage

      if (stage == 1) return(1 - alpha.star[1] - pmvnorm(lower = rep(-Inf, length(treatments)), upper = rep(cut.off, length(treatments)), sigma = cov.matrix[treatments, treatments])[1])
      if (selection == "select.best") return(alpha.star[stage] - alpha.star[stage - 1] - sum(unlist(lapply(treatments, rejection.paths, cut.off = cut.off, treatments = treatments, cov.matrix = cov.matrix, lb = lb, upper.boundary = upper.boundary, K = K, stage = stage)))) # any of 'treatments' could be selected, so we add all these probabilities
      if (stage == 2){
          surviving.subsets <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp)) # list all possible subsets of surviving treatments after the first stage
          return(1 - alpha.star[2] - sum(unlist(lapply(surviving.subsets, get.path.prob, cut.off = cut.off, treatments = treatments, cov.matrix = cov.matrix, lb = lb, upper.boundary = upper.boundary, K = K, stage = stage))))
      }
      surviving.subsets1 <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp)) # all possible subsets of surviving treatments after the first stage
      surviving.subsets2 <- c(list(list(numeric(0))), lapply(surviving.subsets1[-1], function(x) c(list(numeric(0)), lapply(as.list(1:(2 ^ length(x) - 1)), get.hyp)))) # for each possible subset of survivng subsets after stage 1, list the possible subsets still surviving after stage 2
      1 - alpha.star[3] - sum(unlist(Map(function(x, y) sum(unlist(lapply(y, get.path.prob, surviving.subset1 = x, cut.off = cut.off, treatments = treatments, cov.matrix = cov.matrix, lb = lb, upper.boundary = upper.boundary, K = K, stage = stage))), surviving.subsets1, surviving.subsets2)))
  }


  # get sample size ratios
  R <- nMat[, -1] / nMat[1, 1]
  r0 <- nMat[, 1] / nMat[1, 1]

  cov.matrix <- create.cov.matrix(r0, R)
  
  l <- u <- as.list(1:(2 ^ K - 1))

  alpha.star <- rep(list(alpha.star), 2 ^ K - 1)
        
  for (i in 1:(2 ^ K - 1)){

    names(u)[i] <- paste("U_{",paste(get.hyp(i), collapse = " "),"}",sep="")
    names(l)[i] <- paste("L_{",paste(get.hyp(i), collapse = " "),"}",sep="")
    names(alpha.star)[i] <- paste("alpha.star.{",paste(get.hyp(i), collapse = " "),"}",sep="")

    for (j in 1:J){

      try(new.u <- uniroot(excess.alpha, c(0, 10), alpha.star = alpha.star[[i]], treatments = get.hyp(i), cov.matrix = cov.matrix, lb = lb, upper.boundary = u[[i]], selection = selection, K = K, stage = j)$root, silent = TRUE)
      if (is.null(new.u)) {stop("upper boundary not between 0 and 10")}
      u[[i]][j] <- round(new.u, 2)
                
    }

    l[[i]] <- c(lb, u[[i]][J])
  }

  res <- NULL
  res$l <- l
  res$u <- u
  res$sample.sizes <- nMat
  res$K <- K
  res$J <- J
  res$alpha.star <- alpha.star
  res$selection <- selection
  res$zscores <- NULL
  res$selected.trts <- list(1:K)
  class(res) <- "MAMS.stepdown"

  return(res)
}
