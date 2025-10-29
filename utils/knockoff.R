### model-X knockoff
create.solve_equi <- function(Sigma) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  p = nrow(Sigma)
  tol = 1e-10
  
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  
  # Check that the input matrix is positive-definite
  # if (length(which(eigen(G)$values<0))>0) {
  #   stop('The covariance matrix is not positive-definite: cannot solve SDP', immediate. = T)
  # }
  
  if (p > 2) {
    converged = FALSE
    maxitr = 10000
    while (!converged) {
      lambda_min = RSpectra::eigs(G, 1, which = "SR", opts = list(retvec = FALSE, maxitr = 100000, tol = 1e-8))$values
      if (length(lambda_min) == 1) {
        converged = TRUE
      } else {
        if (maxitr > 1e8) {
          warning('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                  covariance matrix. RSpectra::eigs did not converge. Giving up and computing full SVD with built-in R function.',immediate. = T)
          lambda_min = eigen(G, symmetric = T, only.values = T)$values[p]
          converged=TRUE
        } else {
          warning('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                  covariance matrix. RSpectra::eigs did not converge. Trying again with increased number of iterations.',immediate. = T)
          maxitr = maxitr*10
        }
      }
      }
  } else {
      lambda_min = eigen(G, symmetric = T, only.values = T)$values[p]
  }
  
  if (lambda_min < 0) {
    stop('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
         covariance matrix. The covariance matrix is not positive-definite.')
  }
  
  s = rep(1, nrow(Sigma)) * min(lambda_min, 1)
  
  # Compensate for numerical errors (feasibility)
  psd = 1;
  s_eps = 1e-8;
  while (psd > 0) {
    psd = length(which(eigen(2*G - diag(s*(1 - s_eps)), symmetric = T, only.values = T)$values<0))
    if (psd!=0) {
      s_eps = s_eps*10
    }
  }
  s = s*(1 - s_eps)
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}

knockoff <- function(X, y, q, M=1){
  # equi without shrinkage
  #Xk <- create.gaussian(X, mu = rep(0, p), Sigma = covariance, method = 'equi')
  #M  <- stat.glmnet_coefdiff(X, Xk, y)
  #selected_index <- analys(M, abs(M), q)
  #result <- fdp_power(selected_index, signal_index)
  #fdp <- result$fdp
  #power <- result$power
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
  # equi with shrinkage
  e <- matrix(0, nrow = M, ncol = ncol(X))
  num <- 0
  
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
   
  for (j in 1:M) {
    #Xk <- create.second_order(X, method = 'equi')
    
    Xk <- mu_k + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*%
      cholcpp(Sigma_k)

    
    stat  <- stat.glmnet_coefdiff(X, Xk, y, cores=1, nlambda=100,)
 
    selected_index <- analys(stat, abs(stat), q)
    e[j,selected_index] <- dim(X)[2]/(q*length(selected_index))
    num <- num + length(selected_index)
  }
  num <- num/M
  result <- fdp_power(selected_index, signal_index)
  fdp <- result$fdp
  power <- result$power
  
  e_avg <- apply(e, 2, mean)
  #prevent too much ties
  if(length(which(e_avg == sort(e_avg,decreasing = T)[ceiling(num)]))-1>=num/5){
    e_avg <- e_avg+runif(length(e_avg),0,1e-6)
  }
  sel <- which(e_avg >= sort(e_avg,decreasing = T)[ceiling(num)])
  result_stab <- fdp_power(sel, signal_index)
  kn_fdp_stab <- result_stab$fdp
  kn_power_stab <- result_stab$power
  
  

  return(list(fdp = fdp, power = power,
              stab_fdp = kn_fdp_stab, stab_power = kn_power_stab,
              kn_num=length(selected_index),
              kn_stab_num=ceiling(num)))
}



knockoff_res <- function(X, y, q, M=1, statistic="coef_diff"){

  # equi with shrinkage
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
  
  e <- matrix(0, nrow = M, ncol = ncol(X))
  e_ren <- matrix(0, nrow = M, ncol = ncol(X))
  num <- 0
  
  shrink = F
  mu = colMeans(X)
  if (!shrink) {
    # Sigma = cov(X)
    # if (length(which(eigen(Sigma, symmetric = T, only.values = T)$values<0))>0) {
    #   shrink = TRUE
    # }
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
  
  for (j in 1:M) {
    #Xk <- create.second_order(X, method = 'equi')
    
    Xk <- mu_k + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*%
      cholcpp(Sigma_k)
    
    if(statistic=="lambdasmax"){
      stat  <- stat.glmnet_lambdasmax(X, Xk, y)
    }
    if(statistic=="coef_diff"){
      stat  <- stat.glmnet_coefdiff(X, Xk, y, cores=1, nlambda=100,)
    }
     
    selected_index <- analys(stat, abs(stat), q)
    cutoff_set <- max(abs(stat))
    for(t in (abs(stat)-0.00000001)){
      ps <- length(stat[stat > t])
      ng <- length(na.omit(stat[stat < -t]))
      rto <- (1+ng)/max(ps, 1)
      if(rto <= q){
        cutoff_set <- c(cutoff_set, t)
      }
    }
    e_ren[j,selected_index] <- dim(X)[2]/(1+ng)
    e[j,selected_index] <- dim(X)[2]/(q*length(selected_index))
    num <- num + length(selected_index)
  }
  num <- num/M
  
  e_avg <- apply(e, 2, mean)
  #prevent too much ties
  if(length(which(e_avg == sort(e_avg,decreasing = T)[ceiling(num)]))-1>=num/5){
    e_avg <- e_avg+runif(length(e_avg),0,1e-6)
  }
  sel <- which(e_avg >= sort(e_avg,decreasing = T)[ceiling(num)])
  return(list(stab_sel = sel, e_avg = e_avg, kn_sel = selected_index,e_value = e))
}

