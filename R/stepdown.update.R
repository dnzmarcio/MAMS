stepdown.update <- function(current.mams = stepdown.mams(), nobs = NULL, zscores = NULL, selected.trts = NULL, nfuture = NULL) {

    zscores <- c(current.mams$zscores, list(zscores))

    if (!is.null(selected.trts)) selected.trts <- c(current.mams$selected.trts, list(selected.trts))
    
    # checking input parameters
    if (!is(current.mams, "MAMS.stepdown")) {stop("current.mams must be a 'MAMS.stepdown' object")}
    if (length(nobs) != current.mams$K + 1) {stop("must provide observed cumulative sample size for each treatment")}

    completed.stages <- length(zscores)
    for (i in 1:completed.stages) {
        if (length(zscores[[i]]) != current.mams$K) {stop("vector of statistics is wrong length")}
    }

    if (is.null(selected.trts)){
        if (current.mams$J > completed.stages) {stop("must specify treatments selected for next stage")}
    }
    
    for (i in seq_along(selected.trts)){
        if (length(setdiff(selected.trts[[i]], 1:current.mams$K) > 0)) {stop("inappropriate treatment selection")}
    }
    if (is.matrix(nfuture)){
        if (dim(nfuture)[1] != current.mams$J - completed.stages) {stop("must provide future sample sizes for all remaining stages")}
        if (dim(nfuture)[2] != current.mams$K + 1) {stop("must provide future sample sizes for all treatment arms")}
    }

    # load all necessary functions
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
    create.block <- function(control.ratios = 1:2, active.ratios = matrix(1:2, 2, 3)){  # for argument c(i,j) this gives covariance between statistics in stage i with statistics in stage j

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
  create.cond.cov.matrix <- function(cov.matrix, K, completed.stages){ # find the conditional covariance of future test statistics given data so far

    sigma_1_1 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages * K), ((completed.stages - 1) * K + 1):(completed.stages * K)]
    sigma_1_2 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages * K), -(1:(completed.stages * K))]
    sigma_2_1 <- t(sigma_1_2)
    sigma_2_2 <- cov.matrix[-(1:(completed.stages * K)), -(1:(completed.stages * K))]
    sigma_2_2 - sigma_2_1 %*% solve(sigma_1_1) %*% sigma_1_2

  }
  create.cond.mean <- function(cov.matrix, K, completed.stages, zscores){ # find the conditional mean of future test statistics given data so far

    sigma_1_1 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages * K), ((completed.stages - 1) * K + 1):(completed.stages * K)]
    sigma_1_2 <- cov.matrix[((completed.stages - 1) * K + 1):(completed.stages * K), -(1:(completed.stages * K))]
    sigma_2_1 <- t(sigma_1_2)
    sigma_2_1 %*% solve(sigma_1_1) %*% zscores

  } 
  get.path.prob <- function(surviving.subset1, surviving.subset2 = NULL, cut.off, treatments, cov.matrix, lower.boundary, upper.boundary, K, stage, z.means){ # find the probability that no test statistic crosses the upper boundary + only treatments in surviving_subsetj reach the jth stage

    treatments2 <- treatments[surviving.subset1]
    if (stage == 2){
        lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)))
        lower[surviving.subset1] <- lower.boundary[1]

        upper <- c(rep(lower.boundary[1], length(treatments)), rep(cut.off, length(treatments2)))
        upper[surviving.subset1] <- upper.boundary[1]

        return(pmvnorm(lower = lower, upper = upper, mean = z.means[c(treatments, K + treatments2)], sigma = cov.matrix[c(treatments, K + treatments2), c(treatments, K + treatments2)])[1])
    }
    treatments3 <- treatments2[surviving.subset2]

    lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)), rep(-Inf, length(treatments3)))
    lower[surviving.subset1] <- lower.boundary[1]
    lower[length(treatments) + surviving.subset2] <- lower.boundary[2]

    upper <- c(rep(lower.boundary[1], length(treatments)), rep(lower.boundary[2], length(treatments2)), rep(cut.off, length(treatments3)))
    upper[surviving.subset1] <- upper.boundary[1]
    upper[length(treatments) + surviving.subset2] <- upper.boundary[2]

    pmvnorm(lower = lower, upper = upper, mean = z.means[c(treatments, K + treatments2, 2 * K + treatments3)], sigma = cov.matrix[c(treatments, K + treatments2, 2 * K + treatments3), c(treatments, K + treatments2, 2 * K + treatments3)])[1]
  }
  rejection.paths <- function(selected.treatment, cut.off, treatments, cov.matrix, lower.boundary, upper.boundary, K, stage, z.means){  # for the "select.best" method, find the probability that "select.treatment" is selected and subsequently crosses the upper boundary

    contrast <- diag(-1, K + stage - 1)
    contrast[1:K, selected.treatment] <- 1
    for (i in 1:(stage - 1)) contrast[K + i, K + i] <- 1

    bar.mean <- contrast %*% z.means[c(1:K, 1:(stage - 1) * K + selected.treatment)]
  
    bar.cov.matrix <- contrast %*% cov.matrix[c(1:K, 1:(stage - 1) * K + selected.treatment), c(1:K, 1:(stage - 1) * K + selected.treatment)] %*% t(contrast)
    
    lower <- c(rep(0, length(treatments)), cut.off)
    if (stage > 2) lower <- c(rep(0, length(treatments)), lower.boundary[2:(stage - 1)], cut.off)
    lower[which(treatments == selected.treatment)] <- lower.boundary[1]

    upper <- c(rep(Inf, length(treatments)), Inf)
    if (stage > 2) upper <- c(rep(Inf, length(treatments)), upper.boundary[2:(stage - 1)], Inf)
    upper[which(treatments == selected.treatment)] <- upper.boundary[1]
    
    pmvnorm(lower = lower, upper = upper, mean = bar.mean[c(treatments, K + 1:(stage - 1))], sigma = bar.cov.matrix[c(treatments, K + 1:(stage - 1)), c(treatments, K + 1:(stage - 1))])[1]
  }
  excess.alpha <- function(cut.off, alpha.star, treatments, cov.matrix, lower.boundary, upper.boundary, selection.method, K, stage, z.means){# for "all.promising" rule, this gives the cumulative typeI error for 'stage' stages

      # for "select.best" rule, this gives the Type I error spent at the 'stage'th stage
      if (stage == 1) return(1 - alpha.star[1] - pmvnorm(lower = rep(-Inf, length(treatments)), upper = rep(cut.off, length(treatments)), mean = z.means[treatments], sigma = cov.matrix[treatments, treatments])[1])
      if (selection.method == "select.best") return(sum(unlist(lapply(treatments, rejection.paths, cut.off = cut.off, treatments = treatments, cov.matrix = cov.matrix, lower.boundary = lower.boundary, upper.boundary = upper.boundary, K = K, stage = stage, z.means = z.means))) - (alpha.star[stage] - alpha.star[stage - 1])) # any of 'treatments' could be selected, so we add all these probabilities
      if (stage == 2){
          surviving.subsets <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp))  # list all possible subsets of surviving treatments after the first stage
          return(1 - alpha.star[2] - sum(unlist(lapply(surviving.subsets, get.path.prob, cut.off = cut.off, treatments = treatments, cov.matrix = cov.matrix, lower.boundary = lower.boundary, upper.boundary = upper.boundary, K = K, stage = stage, z.means = z.means))))
      }
      surviving.subsets1 <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp)) # all possible subsets of surviving treatments after the first stage
      surviving.subsets2 <- c(list(list(numeric(0))), lapply(surviving.subsets1[-1], function(x) c(list(numeric(0)), lapply(as.list(1:(2 ^ length(x) - 1)), get.hyp)))) # for each possible subset of survivng subsets after stage 1, list the possible subsets still surviving after stage 2
      1 - alpha.star[3] - sum(unlist(Map(function(x, y) sum(unlist(lapply(y, get.path.prob, surviving.subset1 = x, cut.off = cut.off, treatments = treatments, cov.matrix = cov.matrix, lower.boundary = lower.boundary, upper.boundary = upper.boundary, K = K, stage = stage, z.means = z.means))), surviving.subsets1, surviving.subsets2)))
  }

    # give everything the correct name
    alpha.star <- current.mams$alpha.star
    l <- current.mams$l
    u <- current.mams$u
    selection.method <- current.mams$selection
    sample.sizes <- current.mams$sample.sizes
    sample.sizes[completed.stages, ] <- nobs  # Update given the sample sizes actually observed
    if (!all(diff(sample.sizes) >= 0)) {stop("total sample size per arm cannot decrease between stages.")}
    J <- dim(sample.sizes)[1]
    K <- dim(sample.sizes)[2] - 1
    R <- sample.sizes[, -1] / sample.sizes[1, 1]
    r0 <- sample.sizes[, 1] / sample.sizes[1, 1]

    # get conditional distributions BEFORE seeing the new z scores
    cond.cov.matrix <- cov.matrix <- create.cov.matrix(r0, R)
    cond.mean <- rep(0, K * J)
    if (completed.stages > 1){
        cond.cov.matrix <- create.cond.cov.matrix(cov.matrix, K, completed.stages - 1)
        cond.mean <- create.cond.mean(cov.matrix, K, completed.stages - 1, zscores = zscores[[completed.stages - 1]])
    }
    

    # adjust upper boundaries in light of observed sample sizes:
    for (i in 1:(2 ^ K - 1)){ 
        treatments <- intersect(selected.trts[[completed.stages]], get.hyp(i))
        if ((length(treatments > 0)) && (alpha.star[[i]][J] > 0) && (alpha.star[[i]][J] < 1)){
            for (j in completed.stages:J){

                try(new.u <- uniroot(excess.alpha, c(0, 10), alpha.star = alpha.star[[i]][completed.stages:J], treatments = treatments, cov.matrix = cond.cov.matrix, lower.boundary = l[[i]][completed.stages:J], upper.boundary = u[[i]][completed.stages:J], selection.method = selection.method, K = K, stage = j - completed.stages + 1, z.means = cond.mean)$root, silent = TRUE)
                if (is.null(new.u)) {stop("upper boundary not between 0 and 10")}
                u[[i]][j] <- round(new.u, 2)

            }
            l[[i]][J] <- u[[i]][J]
        }
    }
    if (J > completed.stages) {
        cond.cov.matrix <- create.cond.cov.matrix(cov.matrix, K, completed.stages)
        cond.mean <- create.cond.mean(cov.matrix, K, completed.stages, zscores[[completed.stages]])
    }
    for (i in 1:(2 ^ K - 1)) { # get conditional errors
        treatments <- intersect(selected.trts[[completed.stages]], get.hyp(i))
        if ((length(treatments > 0)) && (alpha.star[[i]][J] > 0) && (alpha.star[[i]][J] < 1)){
            max.z <- max(zscores[[completed.stages]][treatments])
            best.treatment <- treatments[which.max(zscores[[completed.stages]][treatments])]
            if (max.z <= u[[i]][completed.stages]) alpha.star[[i]][completed.stages] <- 0
            if (max.z > u[[i]][completed.stages]) {
                alpha.star[[i]][completed.stages:J] <- 1
                if (J > completed.stages) {
                    l[[i]][(completed.stages + 1):J] <- u[[i]][(completed.stages + 1):J] <- -Inf
                }
            }
            else if (max.z <= l[[i]][completed.stages]){
                alpha.star[[i]][completed.stages:J] <- 0
                if (J > completed.stages) {
                    l[[i]][(completed.stages + 1):J] <- u[[i]][(completed.stages + 1):J] <- Inf
                }
            }
            else if (selection.method == "select.best") {
                for (j in (completed.stages + 1):J){
                    alpha.star[[i]][j] <- excess.alpha(cut.off = u[[i]][j], alpha.star = rep(0, J - completed.stages), treatments = best.treatment, cov.matrix = cond.cov.matrix, lower.boundary = l[[i]][(completed.stages + 1):J], upper.boundary = u[[i]][(completed.stages + 1):J], selection.method = selection.method, K = K, stage = j - completed.stages, z.means = cond.mean) + alpha.star[[i]][j - 1]
                }
            }
            else {
                for (j in (completed.stages + 1):J){
                    alpha.star[[i]][j] <- excess.alpha(cut.off = u[[i]][j], alpha.star = rep(0, J - completed.stages), treatments = treatments, cov.matrix = cond.cov.matrix, lower.boundary = l[[i]][(completed.stages + 1):J], upper.boundary = u[[i]][(completed.stages + 1):J], selection.method = selection.method, K = K, stage = j - completed.stages, z.means = cond.mean) 
                }
            }
        }
    }
    if (is.matrix(nfuture)){
        sample.sizes[(completed.stages + 1):J, ] <- nfuture
        if (!all(diff(sample.sizes) >= 0)) {stop("total sample size per arm cannot decrease between stages.")}
        R <- sample.sizes[, -1] / sample.sizes[1, 1]
        r0 <- sample.sizes[, 1] / sample.sizes[1, 1]
        cov.matrix <- create.cov.matrix(r0, R)
        cond.cov.matrix <- create.cond.cov.matrix(cov.matrix, K, completed.stages)
        cond.mean <- create.cond.mean(cov.matrix, K, completed.stages, zscores = zscores[[completed.stages]])
    }
    if (J > completed.stages){
        for (i in 1:(2 ^ K - 1)){ 
            treatments <- intersect(selected.trts[[completed.stages + 1]], get.hyp(i))
            if ((length(treatments > 0)) && (alpha.star[[i]][J] > 0) && (alpha.star[[i]][J] < 1)){
                for (j in (completed.stages + 1):J){
                    try(new.u <- uniroot(excess.alpha, c(0, 10), alpha.star = alpha.star[[i]][(completed.stages + 1):J], treatments = treatments, cov.matrix = cond.cov.matrix, lower.boundary = l[[i]][(completed.stages + 1):J], upper.boundary = u[[i]][(completed.stages + 1):J], selection.method = selection.method, K = K, stage = j - completed.stages, z.means = cond.mean)$root, silent = TRUE)
                    if (is.null(new.u)) {stop("upper boundary not between 0 and 10")}
                    u[[i]][j] <- round(new.u, 2)
                    
                }
                l[[i]][J] <- u[[i]][J]
            }
        }
    }


    res <- NULL
    res$l <- l
    res$u <- u
    res$sample.sizes <- sample.sizes
    res$K <- K
    res$J <- J
    res$alpha.star <- alpha.star
    res$selection <- selection.method
    res$zscores <- zscores
    res$selected.trts <- selected.trts
    class(res) <- "MAMS.stepdown"

    return(res)
      
}
