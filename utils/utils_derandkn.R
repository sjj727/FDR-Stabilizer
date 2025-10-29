#####################################
## The eBH procedure
#####################################
### Input: 
###   E: e-values
###   alpha: target FDR level
### Output:
###   Variables selected by the e-BH procedure
ebh <- function(E, alpha){
  
  p <- length(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (p / alpha / (1:p))
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(list(rej = rej))
}

#######################################
## Computing the knockoff statistics ##
#######################################
### Input: 
###   X: a n-by-p covariate matrix of covariates 
###   Y: a vector of responses of length n 
###   M: the number of knokcoff realizations 
###   mu: mean under PX
###   Sigma: covariance matrix under PX 
###   diags: the diagonal elements for knockoffs construction 
### Output: 
###   an M-by-p matrix of knockoff statistics from multiple runs
kn_stat <- function(X, Y, M, mu, Sigma, 
                    diags, family = "gaussian", 
                    stat_method = stat.glmnet_coefdiff){

  n <- dim(X)[1]
  p <- dim(X)[2]
  stat_mat <- matrix(0, M, p)
  
  for(m in 1:M){
    Xk <- create.gaussian(X, mu, Sigma, diag_s = diags)
    W <- stat_method(X, Xk, Y, family = family)
    stat_mat[m,] <- W
  }
  
  return(stat_mat)
}

#######################################################
## Compute stopping time w/ diff alpha_kn and offset ##
#######################################################
### Input:
###   W: a length p vector of knockoff feature importance statistics
###   fdr: the target FDR level
###   offset: 0 or 1 
### Output: 
###   the knockoff selection threshold
alphakn_threshold <- function(W, fdr, offset) {
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
 

##################################################
### Compute knockoff evals with knockoff stats ###
##################################################
### Input:
###   stat_mat: an M-by-p matrix of multiple knockoff feature importance statistics 
###   gamma: alpha_kn
###   offset: value between 0 and 1
### Output: 
###   E: the list of knockoff e-values 
kn_evals <- function(stat_mat, gamma, offset){

  M <- dim(stat_mat)[1]
  p <- dim(stat_mat)[2]
  E_list <- c()
  
  for(m in 1:M){
    W <- stat_mat[m,]
    tau <- stop_early(W, gamma, offset) 
    ## knockoff e-values for the m-th run
    E <- (W >= tau) / (1 + sum(W <= -tau))
    E_list <- rbind(E_list, E)

  }
  
  E <- p * colMeans(E_list)
  return(E)
}


########################
## Knockoffs e-values ## 
########################
### Input:
###   X: a n-by-p covariate matrix of covariates 
###   Y: a vector of responses of length n 
###   M: the number of knokcoff realizations 
###   gamma: alpha_kn
###   mu: mean under PX
###   Sigma: covariance matrix under PX 
###   diags: the diagonal elements for knockoffs construction 
### Output: 
###   E: the list of knockoff e-values 
ekn <- function(X, Y, M, gamma,
                family = "gaussian", offset = 1){
  is_posdef <- function (A, tol = 1e-09) 
  {
    p = nrow(matrix(A))
    if (p < 500) {
      lambda_min = min(eigen(A)$values)
    }
    else {
      oldw <- getOption("warn")
      options(warn = -1)
      lambda_min = RSpectra::eigs(A, 1, which = "SM", opts = list(retvec = FALSE, 
                                                                  maxitr = 100, tol))$values
      options(warn = oldw)
      if (length(lambda_min) == 0) {
        lambda_min = min(eigen(A)$values)
      }
    }
    return(lambda_min > tol * 10)
  }
  n <- dim(X)[1]
  p <- dim(X)[2]

  E <- matrix(0, M, p)
  
  shrink = F
  mu = colMeans(X)
  if (!shrink) {
    Sigma = cov(X)
    if (!is_posdef(Sigma)) {
      shrink = TRUE
    }
  }
  if (shrink) {
    if (!requireNamespace("corpcor", quietly = T))
      stop("corpcor is not installed", call. = F)
    Sigma <- tryCatch({
      suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X,
                                                             verbose = F)), nrow = ncol(X)))
    }, warning = function(w) {
    }, error = function(e) {
      stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
    }, finally = {
    }) 
  }
  diag_s = NULL
  if (is.null(diag_s)) {
    diag_s <- diag(create.solve_equi(Sigma))
  }
  if (is.null(dim(diag_s))) {
    diag_s = diag(diag_s, length(diag_s))
  }
  if (all(diag_s == 0)) {
    warning("The conditional knockoff covariance matrix is not positive definite. Knockoffs will have no power.")
    return(X)
  }
  #SigmaInv_s = solve(Sigma, diag_s)
  SigmaInv_s <- solve_linear_equation(Sigma, diag_s)
  mu_k <- X - sweep(X, 2, mu, "-") %*% SigmaInv_s
  Sigma_k <- 2 * diag_s - diag_s %*% SigmaInv_s
  
  for(m in 1:M){
    #Xk <- create.second_order(X, method = "equi")
    Xk <- mu_k + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*%
      cholcpp(Sigma_k)
    
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family, cores=1, nlambda=100,)
    tau <- stop_early(W, gamma, offset)
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))

  }
  E <- p * colMeans(E)
  return(list(E = E))

}

#######################################
## Computing the early stopping time ##
#######################################
### Input:
###   W: vector of knockoff feature importance statistics 
###   gamma: alpha_kn 
###   offset: value between 0 and 1
### Output: 
###   The modified knockoff stopping time defined in (14)
stop_early <- function(W, gamma, offset){
  
  tau <- alphakn_threshold(W, fdr =  gamma, offset = offset) 
  ord_W <- order(abs(W), decreasing = TRUE)
  sorted_W <- W[ord_W]
  
  if(sum(W>0) >= 1 / gamma){
    pos_ind <- which(sorted_W > 0)
    tau1 <- sorted_W[pos_ind[ceiling(1/gamma)-1]]
  }else{
    tau1 <- 0
  }
  tau <- min(tau,tau1) 

  return(tau)
}


weighted_ekn <- function(X, Y, M, alpha, gamma, weight, 
                         mu, Sigma, diags, offset = 1, 
                         family = "gaussian"){

  n <- dim(X)[1]
  p <- dim(X)[2]

  E <- matrix(0, M, p)

  for(m in 1:M){
    Xk <- create.gaussian(X, mu, Sigma, diag_s = diags)
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family)        
    tau <- stop_early(W, gamma, offset)
    E[m,] <- (W >= tau) * weight / (weight + sum((W <= -tau) * weight))
  }

  E <- p * colMeans(E)
  rej <- ebh(E, alpha)$rej

  return(list(rej = rej))
}


adaptive_ekn <- function(X, Y, M, alpha, gamma, U, 
                         Sigma, diags, family = "gaussian",
                         offset = 1){

  n <- dim(X)[1]
  p <- dim(X)[2]

  E <- matrix(0, M, p)

  for(m in 1:M){
    Xk <- create.gaussian(X, rep(0,p), Sigma, diag_s = diags)
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family)    

    ## Test different threshold
    akn_res <- filter_gam(W, U, gamma)
    if_remain <- (1:p) %in% akn_res$unrevealed_id
    E[m,] <- if_remain * (W > 0) / (1 + sum((W < 0) * if_remain))
  }

  E <- colMeans(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (1 / alpha / (1:p))
  id <- max(which(comp > 0))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }

  return(list(rej = rej, E = E))

}

pfer_kn <- function(X, Y, v0, M, tau, 
                    mu, Sigma, diags, family = "gaussian"){

  n <- nrow(X)
  p <- ncol(X)
  
  freq <- rep(0, p)
  for (m in 1:M){
    v <- floor(v0)+rbinom(1,1,v0-floor(v0))
    Xk <- create.gaussian(X,mu,Sigma,diag_s = diags)
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family)
    order_w <- order(abs(W),decreasing = TRUE)
    sorted_w <- W[order_w]
    negid <- which(sorted_w<0)
    if(v>0){
      if(sum(W<0)<v){
        S <- which(W>0)
        freq[S] <- freq[S]+1
      }else{
        TT <- negid[v]
        S <- which(sorted_w[1:TT]>0)
        S <- order_w[S]
        freq[S] <- freq[S]+1
        }
      }
  }
  freq <- freq/M
  rej <- which(freq>=tau)

return(list(rej = rej))
}


multienv_ekn <- function(Xlist, Ylist, nenv, 
                         M, alpha, gamma, mu, Sigma, diags, 
                         family = "gaussian", offset = 1){

  n <- dim(Xlist[[1]])[1]
  p <- dim(Xlist[[1]])[2]

  E <- matrix(0, M, p)

  for(m in 1:M){
    
    Xklist <- list()
    for(i in 1:nenv){
      Xklist[[i]] <- create.gaussian(Xlist[[i]], mu, Sigma, diag_s = diags)
    }
    W_matrix <- compute_stats_with_prior(Xlist, Xklist, Ylist)
    W_sign <- 2*((apply(W_matrix, 1, min) > 0) - 0.5)
    W_prod <- abs(apply(W_matrix, 1, prod))
    W <- W_sign * W_prod
    tau <- stop+early(W, gamma, offset)
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))

  }
  E <- p*colMeans(E)
  rej <- ebh(E,alpha)$rej

  return(list(rej = rej))

}


derand_kn <- function(X,Y,M,alpha){
  res <- ekn(X, Y, M, alpha / 2, "gaussian", offset = 1)
  E <- res$E
  rej <- ebh(E, alpha)$rej
  r <- fdp_power(rej,signal_index)
  fdp <- r$fdp
  power <- r$power
  num <- length(rej)
  return(list(fdp=fdp,power=power,num=num,e_values=E,rej=rej))
}

derand_kn_res <- function(X,Y,M,alpha){
  res <- ekn(X, Y, M, alpha / 2, "gaussian", offset = 1)
  E <- res$E
  rej <- ebh(E, alpha)$rej
  
  return(list(rej=rej))
}
