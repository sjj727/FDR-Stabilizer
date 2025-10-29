### High dimension linear model
### base procedure: data splitting
### compare different g
rm(list = ls())
library(snowfall)
library(here)
library(MASS)
library(ExtMallows)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)
library(TAF)
library(plyr)
library(stabs)
library(lars)
library(EnvStats)
library(VGAM)


### read files
setwd(here::here())
folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  source(file_path)
}



derand_ds <- function(X, y, M,target_q = 0.1, q = 0.05) {
  n <- nrow(X); p <- ncol(X)
  
  e_values_mat <- matrix(0, nrow = p, ncol = M)
  
  for (m in 1:M) {
    fit <- DS_single(X, y, q = q)
    e_vals <- rep(0, p)
    
    denom <- q*length(fit$sel)
    e_vals[fit$sel] <- 1
    e_vals <- (p * e_vals) / max(denom,1)
    
    e_values_mat[, m] <- e_vals
  }
  
  # average e-values
  e_values_avg <- rowMeans(e_values_mat, na.rm = TRUE)
  
  # eBH
  ord <- order(e_values_avg, decreasing = TRUE)
  e_sorted <- e_values_avg[ord]
  k_hat <- 0
  if (any(e_sorted > 0)) {
    thr_seq <- (p / (target_q * (1:p)))
    meet <- which(e_sorted >= thr_seq)
    if (length(meet) > 0) k_hat <- max(meet)
  }
  
  if (k_hat == 0) {
    sel_idx <- integer(0)
  } else {
    threshold <- p / (target_q * k_hat)
    sel_idx <- which(e_values_avg >= threshold)
    sel_idx <- sort(sel_idx)
  }
  
  return(sel_idx)
}


aggr_MDS <- function(rej_mat, q){
  shat <- rowSums(rej_mat)
  inclusion_rate <- rej_mat/pmax(shat,1)
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  
  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    
    ### backtracking 
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    
    selected_index <- setdiff(feature_rank, null_feature)
  }
  else{
    return(numeric())
  }
}

stab_eBH <- function(rej_mat,
                     aggr_method = "eavg",
                     statistic = NA,
                     trim_prop = 0.1,
                     tie.method = "random") {
  
  stopifnot(is.matrix(rej_mat))
  M <- nrow(rej_mat); p <- ncol(rej_mat)
  
  shat <- rowSums(rej_mat)
  sbar <- ceiling(mean(shat))
  
  # average e-values
  if (aggr_method == "eavg") {
    
    s_safe <- ifelse(shat > 0, shat, 1)
    e_mat <- sweep(rej_mat, 1, s_safe, FUN = "/")
    e_mat[shat == 0, ] <- 0
    pi_tilde <- rank(-colMeans(e_mat), ties.method = tie.method)
    return(which(pi_tilde <= sbar))
  }
  # selection probability
  if (aggr_method == "sel_prob"){
    prob <- colMeans(rej_mat)
    pi_tilde <- rank(-prob, ties.method = tie.method)
    return(which(pi_tilde <= sbar))
  }
  
  if (anyNA(statistic)) {
    stop("For rank-based aggregation methods (other than 'eavg' and 'sel_prob'), please provide 'statistic' as an M x p matrix.")
  }
  stopifnot(nrow(statistic) == M && ncol(statistic) == p)
  
  rank_stat <- apply(statistic, 1, function(z) rank(-z, ties.method = tie.method))
  if (!is.matrix(rank_stat)) rank_stat <- matrix(rank_stat, nrow = p)
  rank_stat <- t(rank_stat)  # M x p
  
  score <- switch(aggr_method,
    "mean"  = -colMeans(statistic),                                
    "median"= -apply(statistic, 2, median),                        
    "rank_mean" = {                                               
      colMeans(rank_stat)
    },
    "rank_min"   = apply(rank_stat, 2, min),                           
    "rank_geom_mean" = exp(colMeans(log(rank_stat))),                  
    "rank_harm_mean" = {                                               
      1 / colMeans(1 / rank_stat)
    },
    "trimmed_mean" = {                                            
      a <- max(0, min(0.5, trim_prop))
      -apply(statistic, 2, function(v) mean(v, trim = a))
    },
    "rra" = {
      if (requireNamespace("RobustRankAggreg", quietly = TRUE)) {
        perm_mat <- t(apply(rank_stat, 1, order))
        L <- lapply(seq_len(nrow(perm_mat)), function(i) as.character(perm_mat[i, ]))
        out <- RobustRankAggreg::aggregateRanks(L)
        
        ord <- as.integer(out$Name) 
        sc <- integer(p); sc[ord] <- seq_len(p)
        sc
      } else {
        warning("Package 'RobustRankAggreg' not found.")
      }
    },
    "huber" = {
      huber <- function(x) {
        mu <- mean(x)
        s  <- mad(x, constant = 1)
        if (s <= 0) return(mu)
        delta <- 1.345 * s
        r <- x - mu
        mu + mean(pmin(pmax(r, -delta), delta))
      }
      -apply(statistic, 2, huber)
    },
    "quantile" = {
      qq <- max(0, min(1, 0.25))
      apply(-statistic, 2, function(v) quantile(v, probs = 0.25, type = 7, names = FALSE))
    },
    "MM" = {
      # Mallows Model
      if (requireNamespace("ExtMallows", quietly = TRUE)) {
        order_stat <- apply(statistic, 1, function(z) order(-z+runif(length(z), -1e-12, 1e-12)))
        
        res=MM(rankings = order_stat[1:(sbar+10),],initial.method = "mean",it.max = 50)
        return(as.numeric(res$op.pi0[1:sbar]))
      } else{
        warning("Package 'ExtMallows' not found.")
      }
    },
    "EMM" = {
      # Extended Mallows Model
      if (requireNamespace("ExtMallows", quietly = TRUE)) {
        order_stat <- apply(statistic, 1, function(z) order(-z+runif(length(z), -1e-12, 1e-12)))
        
        res=EMM(rankings = order_stat[1:(sbar+10),],initial.method = "mean",it.max = 50)
        return(as.numeric(res$op.pi0[1:sbar]))
      } else{
        warning("Package 'ExtMallows' not found.")
      }
    },
    stop("Unknown aggr_method: ", aggr_method)
  )
  
  pi_tilde <- rank(score, ties.method = tie.method)
  which(pi_tilde <= sbar)
}










DS_single <- function(X, y, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  t <- 0
  while(t<=20){
    t <- t+1
    sample_index1 <- sample(x = c(1:n), size = trunc(0.5 * n), replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    cvfit = cv.glmnet(X[sample_index1, ], y[sample_index1], nfolds = 5, nlambda = 30, lambda.min.ratio = 1e-2)
    beta1 = as.vector(coef(cvfit$glmnet.fit, cvfit$lambda.1se))[-1]
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0)
      break
  }
  if (t>=21){
    return(list(sel=c(1)[0],stat=rep(0,p)))
  }
  ### randomly split the data
  
  if(length(nonzero_index)!=0){
    ### run OLS on the second half of the data, restricted on the selected features
    beta2 <- rep(0, p)
    beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
    
    ### calculate the mirror statistics
    stat <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
    #stat <- abs(beta1 + beta2) - abs(beta1 - beta2)
    stat[is.na(stat)] = 0
    selected_index <- analys(stat, abs(stat), q)
  }
  
  
  return(list(sel=selected_index,stat=stat))
}



run_all_stab_eBH <- function(rej_mat, stat_mat, signal_index, q) {
  methods_vec <- c("eavg","sel_prob","mean","median","rank_mean",
                   "rank_min","rank_geom_mean","rank_harm_mean","quantile",
                   "huber","rra","MM","EMM")
  out <- lapply(methods_vec, function(meth) {
    sel <- stab_eBH(rej_mat, aggr_method = meth, stat_mat)
    fdp_power(sel, signal_index)
  })
  data.frame(method = methods_vec,
             fdp    = sapply(out, `[[`, "fdp"),
             power  = sapply(out, `[[`, "power"))
}




### correlation_Structure = "toeplitz","Compound_Symmetry"
### X_Structure = "Normal","Mix_Normal","t"
### epsilon_Structure = "Normal","Mix_Normal","Heteroscedasticity"

# fix delta, change rho
simulation_corr <- function(n = 800, p = 2000, p0 = 50, q = 0.1, delta = 5,
                            correlation_Structure = "toeplitz", 
                            X_Structure = "Normal",
                            epsilon_Structure = "Normal"){
  res <- data.frame()
  ### algorithmic settings
  for (rho in (1:9)/10) {
    ### correlation
    if (correlation_Structure == "toeplitz"){
      sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
      covariance <- bdiag(rep(list(sig1), 10))
      diag(covariance) <- rep(1, p)
    }
    if (correlation_Structure == "Compound_Symmetry"){
      covariance <- matrix(rep(0,p^2),nrow = p)
      for(i in 1 : p){
        for(j in 1 : p){
          covariance[i,j] <- rho ^ (i!= j)
        }
      }
    }
    ### generate X
    if (X_Structure == "Normal"){
      X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
    }
    if (X_Structure == "Mix_Normal"){
      X1 <- mvrnorm(n/2, mu = rep(0.5, p),  Sigma = covariance)
      X2 <- mvrnorm(n/2, mu = rep(-0.5, p), Sigma = covariance)
      X  <- rbind(X1, X2)
    }
    if (X_Structure == "t"){
      X <- rmvt(n, sigma = as.matrix(covariance), df = 3)
    }
    ### randomly generate the true beta
    beta_star <- rep(0, p)
    signal_index <<- sample(c(1:p), size = p0, replace = F)
    beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
    
    ### generate y
    if (epsilon_Structure == "Normal"){
      y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
    }
    if (epsilon_Structure == "Heteroscedasticity"){
      y <- X%*%beta_star + rowSums(abs(X))/p * rexp(n,2)
    }
    if (epsilon_Structure == "Mix_Normal"){
      y <- X%*%beta_star + rnormMix(n, mean1 = -1, sd1 = 2.5, mean2 = 1, sd2 = 0.5, p.mix = 0.5)
    }
    
    ### number of splits
    M <- 20
    
    ### DS
    rej_mat <- stat_mat <- matrix(0, M, p)
    for (m in 1:M) {
      tmp <- DS_single(X, y, q)
      rej_mat[m, tmp$sel] <- 1
      stat_mat[m, ]      <- tmp$stat
    }
    stab_ds_res <- run_all_stab_eBH(rej_mat, stat_mat, signal_index, q)
    derand_sel_idx <- derand_ds(X, y, M = M, target_q = q, q = q/2)
    derand_res <- fdp_power(derand_sel_idx, signal_index)
    mds_idx <- aggr_MDS(rej_mat, q)
    mds_res <- fdp_power(mds_idx, signal_index)
    
    
    ds_fdp  <- fdp_power(which(rej_mat[1,]==1),signal_index)$fdp
    ds_power<- fdp_power(which(rej_mat[1,]==1),signal_index)$power
    
    ds_row <- data.frame(method = c("DS","derand_DS","MDS"),
                         fdp    = c(ds_fdp,derand_res$fdp,mds_res$fdp),
                         power  = c(ds_power,derand_res$power,mds_res$power))
    
    data_save <- rbind(stab_ds_res, ds_row)
    m   <- data_save$method
    fdp <- data_save$fdp
    pwr <- data_save$power
    
    vals  <- c(rho, fdp[1], pwr[1], fdp[2], pwr[2], fdp[3], pwr[3],
               fdp[4], pwr[4], fdp[5], pwr[5], fdp[6], pwr[6],
               fdp[7], pwr[7], fdp[8], pwr[8], fdp[9], pwr[9],
               fdp[10], pwr[10], fdp[11], pwr[11],
               fdp[12], pwr[12], fdp[13], pwr[13], fdp[14], pwr[14],
               fdp[15], pwr[15], fdp[16], pwr[16])
    
    cols <- c("rho",
              paste(rep(m, each = 2), c("fdp","power"), sep = "_"))
    
    data_save <- as.data.frame(matrix(vals, nrow = 1))
    colnames(data_save) <- cols
    res <- rbind(res,data_save)
  }
  return(res)
}



# fix rho, change delta
simulation_signal <- function(n = 800, p = 2000, p0 = 50, q = 0.1,d, 
                              correlation_Structure = "toeplitz", 
                              X_Structure = "Normal",
                              epsilon_Structure = "Normal"){
  res <- data.frame()
  for (delta in 3:7) {
    
    ### correlation
    if (correlation_Structure == "toeplitz"){
      sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
      covariance <- bdiag(rep(list(sig1), 10))
      diag(covariance) <- rep(1, p)
    }
    if (correlation_Structure == "Compound_Symmetry"){
      covariance <- matrix(rep(0,p^2),nrow = p)
      for(i in 1 : p){
        for(j in 1 : p){
          covariance[i,j] <- rho ^ (i!= j)
        }
      }
    }
    ### generate X
    if (X_Structure == "Normal"){
      X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
    }
    if (X_Structure == "Mix_Normal"){
      X1 <- mvrnorm(n/2, mu = rep(0.5, p),  Sigma = covariance)
      X2 <- mvrnorm(n/2, mu = rep(-0.5, p), Sigma = covariance)
      X  <- rbind(X1, X2)
    }
    if (X_Structure == "t"){
      X <- rmvt(n, sigma = as.matrix(covariance), df = 3)
    }
    ### randomly generate the true beta
    beta_star <- rep(0, p)
    signal_index <<- sample(c(1:p), size = p0, replace = F)
    beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
    
    ### generate y
    if (epsilon_Structure == "Normal"){
      y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
    }
    if (epsilon_Structure == "Heteroscedasticity"){
      y <- X%*%beta_star + rowSums(abs(X))/p * rexp(n,2)
    }
    if (epsilon_Structure == "Mix_Normal"){
      y <- X%*%beta_star + rnormMix(n, mean1 = -1, sd1 = 2.5, mean2 = 1, sd2 = 0.5, p.mix = 0.5)
    }
    
    ### number of splits
    M <- 20
    
    
    
    
    ### DS_single + stab_eBH
    rej_mat <- stat_mat <- matrix(0, M, p)
    for (m in 1:M) {
      tmp <- DS_single(X, y, q)
      rej_mat[m, tmp$sel] <- 1
      stat_mat[m, ]      <- tmp$stat
    }
    stab_ds_res <- run_all_stab_eBH(rej_mat, stat_mat, signal_index, q)
    
    derand_sel_idx <- derand_ds(X, y, M = M, target_q = q, q = q/2)
    derand_res <- fdp_power(derand_sel_idx, signal_index)
    
    mds_idx <- aggr_MDS(rej_mat, q)
    mds_res <- fdp_power(mds_idx, signal_index)
    
    
    ds_fdp  <- fdp_power(which(rej_mat[1,]==1),signal_index)$fdp
    ds_power<- fdp_power(which(rej_mat[1,]==1),signal_index)$power
    
    ds_row <- data.frame(method = c("DS","derand_DS","MDS"),
                         fdp    = c(ds_fdp,derand_res$fdp,mds_res$fdp),
                         power  = c(ds_power,derand_res$power,mds_res$power))
    
    data_save <- rbind(stab_ds_res, ds_row)
    m   <- data_save$method
    fdp <- data_save$fdp
    pwr <- data_save$power
    
    vals  <- c(delta, fdp[1], pwr[1], fdp[2], pwr[2], fdp[3], pwr[3],
               fdp[4], pwr[4], fdp[5], pwr[5], fdp[6], pwr[6],
               fdp[7], pwr[7], fdp[8], pwr[8], fdp[9], pwr[9],
               fdp[10], pwr[10], fdp[11], pwr[11],
               fdp[12], pwr[12], fdp[13], pwr[13], fdp[14], pwr[14],
               fdp[15], pwr[15], fdp[16], pwr[16])
    
    cols <- c("delta",
              paste(rep(m, each = 2), c("fdp","power"), sep = "_"))
    
    data_save <- as.data.frame(matrix(vals, nrow = 1))
    colnames(data_save) <- cols
    res <- rbind(res,data_save)
    
  }
  return(res)
}

Fun <- function(num, type = "signal", n, p, p0, q, arg = 0.5,
                correlation_Structure = "toeplitz", 
                X_Structure = "Normal", 
                epsilon_Structure = "Normal"){
  set.seed(num)
  if (type == "signal"){
    ret <- simulation_signal(n = n, p = p, p0 = p0, q = q, rho = arg, 
                             correlation_Structure = correlation_Structure, 
                             X_Structure = X_Structure,
                             epsilon_Structure = epsilon_Structure)
  }
  if (type == "corr"){
    ret <- simulation_corr(n = n, p = p, p0 = p0, q = q, delta = arg, 
                           correlation_Structure = correlation_Structure, 
                           X_Structure = X_Structure,
                           epsilon_Structure = epsilon_Structure)
  }
  return(ret)
  
}

simulation_par <- function(rep_num, n, p, p0, q, arg = 0.5, type = "signal", 
                           correlation_Structure = "toeplitz", 
                           X_Structure = "Normal",
                           epsilon_Structure = "Normal"){
  
  num <- 1:rep_num
  res <- sfLapply(num, fun = Fun, type, 
                  n, p, p0, q, arg, 
                  correlation_Structure, 
                  X_Structure, 
                  epsilon_Structure)
  
  data_save <- Reduce("+",res)/max(num)
  if(type == "corr"){
    save(data_save, file = paste("results/results_figure2/delta", arg, "_n", n, 
                                 "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                                 "_eps_", epsilon_Structure, ".RData", sep = ""))
  }
  if(type == "signal"){
    save(data_save, file = paste("results/results_figure2/rho", arg, "_n", n, 
                                 "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                                 "_eps_", epsilon_Structure, ".RData", sep = ""))
  }
  
}





avg_jaccard <- function(list_sets) {
  if (length(list_sets) < 2) return(0)
  mean(combn(length(list_sets), 2,
             function(i) {
               length(intersect(list_sets[[i[1]]], list_sets[[i[2]]])) / max(length(union(list_sets[[i[1]]], list_sets[[i[2]]])), 1)
             }, simplify = TRUE))
}

run_all_jaccard <- function(X, y, q, M, B, signal_index) {
  stab_methods <- c("eavg","sel_prob","mean","median","rank_mean",
                    "rank_min","rank_geom_mean","rank_harm_mean","quantile",
                    "huber","rra","MM","EMM")
  
  extra_methods <- c("DS", "MDS", "derand_DS")
  
  all_methods <- c(stab_methods, extra_methods)
  jac_list    <- setNames(vector("list", length(all_methods)), all_methods)
  
  for (b in 1:B) {
    rej_mat <- stat_mat <- matrix(0, M, ncol(X))
    
    for (m in 1:M) {
      tmp <- DS_single(X, y, q)
      rej_mat[m, tmp$sel] <- 1
      stat_mat[m, ]      <- tmp$stat
    }
    
    for (meth in stab_methods) {
      sel <- stab_eBH(rej_mat, aggr_method = meth, stat_mat)
      jac_list[[meth]][[b]] <- sel
    }
    
    jac_list[["DS"]][[b]] <- which(rej_mat[1, ] == 1)
    
    mds_idx <- aggr_MDS(rej_mat, q)
    jac_list[["MDS"]][[b]] <- mds_idx
    
    derand_idx <- derand_ds(X, y, M, target_q = q, q = q/2)
    jac_list[["derand_DS"]][[b]] <- derand_idx
  }
  
  sapply(jac_list, function(lst) avg_jaccard(lst))
}




simulation_corr_jaccard <- function(B=50, n = 800, p = 2000, p0 = 50, q = 0.1, delta = 5,
                                    correlation_Structure = "toeplitz", 
                                    X_Structure = "Normal",
                                    epsilon_Structure = "Normal"){
  res <- data.frame()
  ### algorithmic settings
  for (rho in (1:9)/10) {
    ### correlation
    if (correlation_Structure == "toeplitz"){
      sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
      covariance <- bdiag(rep(list(sig1), 10))
      diag(covariance) <- rep(1, p)
    }
    if (correlation_Structure == "Compound_Symmetry"){
      covariance <- matrix(rep(0,p^2),nrow = p)
      for(i in 1 : p){
        for(j in 1 : p){
          covariance[i,j] <- rho ^ (i!= j)
        }
      }
    }
    ### generate X
    if (X_Structure == "Normal"){
      X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
    }
    if (X_Structure == "Mix_Normal"){
      X1 <- mvrnorm(n/2, mu = rep(0.5, p),  Sigma = covariance)
      X2 <- mvrnorm(n/2, mu = rep(-0.5, p), Sigma = covariance)
      X  <- rbind(X1, X2)
    }
    if (X_Structure == "t"){
      X <- rmvt(n, sigma = as.matrix(covariance), df = 3)
    }
    ### randomly generate the true beta
    beta_star <- rep(0, p)
    signal_index <<- sample(c(1:p), size = p0, replace = F)
    beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
    
    ### generate y
    if (epsilon_Structure == "Normal"){
      y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
    }
    if (epsilon_Structure == "Heteroscedasticity"){
      y <- X%*%beta_star + rowSums(abs(X))/p * rexp(n,2)
    }
    if (epsilon_Structure == "Mix_Normal"){
      y <- X%*%beta_star + rnormMix(n, mean1 = -1, sd1 = 2.5, mean2 = 1, sd2 = 0.5, p.mix = 0.5)
    }
    
    jac_vec <- run_all_jaccard(X, y, q, 20, B = B)
    cols  <- c("rho", paste0(names(jac_vec), "_jaccard"))
    vals  <- c(rho, jac_vec)
    data_save <- as.data.frame(matrix(vals, nrow = 1))
    colnames(data_save) <- cols
    res <- rbind(res, data_save)
  }
  return(res)
}

# fix rho, change delta
simulation_signal_jaccard <- function(B=50, n = 800, p = 2000, p0 = 50, q = 0.1, rho = 0.5, 
                                      correlation_Structure = "toeplitz", 
                                      X_Structure = "Normal",
                                      epsilon_Structure = "Normal"){
  res <- data.frame()
  for (delta in 3:7) {
    
    ### correlation
    if (correlation_Structure == "toeplitz"){
      sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
      covariance <- bdiag(rep(list(sig1), 10))
      diag(covariance) <- rep(1, p)
    }
    if (correlation_Structure == "Compound_Symmetry"){
      covariance <- matrix(rep(0,p^2),nrow = p)
      for(i in 1 : p){
        for(j in 1 : p){
          covariance[i,j] <- rho ^ (i!= j)
        }
      }
    }
    ### generate X
    if (X_Structure == "Normal"){
      X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
    }
    if (X_Structure == "Mix_Normal"){
      X1 <- mvrnorm(n/2, mu = rep(0.5, p),  Sigma = covariance)
      X2 <- mvrnorm(n/2, mu = rep(-0.5, p), Sigma = covariance)
      X  <- rbind(X1, X2)
    }
    if (X_Structure == "t"){
      X <- rmvt(n, sigma = as.matrix(covariance), df = 3)
    }
    ### randomly generate the true beta
    beta_star <- rep(0, p)
    signal_index <<- sample(c(1:p), size = p0, replace = F)
    beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
    
    ### generate y
    if (epsilon_Structure == "Normal"){
      y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
    }
    if (epsilon_Structure == "Heteroscedasticity"){
      y <- X%*%beta_star + rowSums(abs(X))/p * rexp(n,2)
    }
    if (epsilon_Structure == "Mix_Normal"){
      y <- X%*%beta_star + rnormMix(n, mean1 = -1, sd1 = 2.5, mean2 = 1, sd2 = 0.5, p.mix = 0.5)
    }
    jac_vec <- run_all_jaccard(X, y, q, 20, B = B)
    cols  <- c("delta", paste0(names(jac_vec), "_jaccard"))
    vals  <- c(delta, jac_vec)
    data_save <- as.data.frame(matrix(vals, nrow = 1))
    colnames(data_save) <- cols
    res <- rbind(res, data_save)
    
  }
  return(res)
}

Fun_jaccard <- function(num, type = "signal", B, n, p, p0, q, arg = 0.5,
                        correlation_Structure = "toeplitz", 
                        X_Structure = "Normal", 
                        epsilon_Structure = "Normal"){
  set.seed(num)
  if (type == "signal"){
    ret <- simulation_signal_jaccard(B = B, n = n, p = p, p0 = p0, q = q, rho = arg, 
                                     correlation_Structure = correlation_Structure, 
                                     X_Structure = X_Structure,
                                     epsilon_Structure = epsilon_Structure)
  }
  if (type == "corr"){
    ret <- simulation_corr_jaccard(B = B, n = n, p = p, p0 = p0, q = q, delta = arg, 
                                   correlation_Structure = correlation_Structure, 
                                   X_Structure = X_Structure,
                                   epsilon_Structure = epsilon_Structure)
  }
  return(ret)
  
}

simulation_par_jaccard <- function(rep_num, B, n, p, p0, q, arg = 0.5, type = "signal", 
                                   correlation_Structure = "toeplitz", 
                                   X_Structure = "Normal",
                                   epsilon_Structure = "Normal"){
  
  num <- 1:rep_num
  res <- sfLapply(num, fun = Fun_jaccard, type, B, 
                  n, p, p0, q, arg, 
                  correlation_Structure, 
                  X_Structure, 
                  epsilon_Structure)
  
  data_save <- Reduce("+",res)/max(num)
  if(type == "corr"){
    save(data_save, file = paste("results/results_figure2/jaccard_delta", arg, "_n", n, 
                                 "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                                 "_eps_", epsilon_Structure, ".RData", sep = ""))
  }
  if(type == "signal"){
    save(data_save, file = paste("results/results_figure2/jaccard_rho", arg, "_n", n, 
                                 "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                                 "_eps_", epsilon_Structure, ".RData", sep = ""))
  }
  
}


sfInit(parallel = TRUE, cpus = 50)

sfLibrary(snowfall)
sfLibrary(MASS)
sfLibrary(glmnet)
sfLibrary(knockoff)
sfLibrary(ExtMallows)
sfLibrary(TAF)
sfLibrary(mvtnorm)
sfLibrary(hdi)
sfLibrary(plyr)
sfLibrary(stabs)
sfLibrary(lars)
sfLibrary(EnvStats)
sfLibrary(VGAM)


folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  sfSource(file_path)
}
sfExport("simulation_signal")
sfExport("simulation_corr")
sfExport("DS_single")
sfExport("stab_eBH")
sfExport("run_all_stab_eBH")
sfExport("aggr_MDS")
sfExport("derand_ds")
sfExport("simulation_signal_jaccard")
sfExport("simulation_corr_jaccard")
sfExport("run_all_jaccard")
sfExport("avg_jaccard")

#rep_num = 1000 in paper, 50 for fast
rep_num <- 50


n <- 800
p <- 2000
p0 <- 50
q <- 0.1

# Figure 2 FDR and Power in paper
system.time(simulation_par(rep_num, n, p, p0, q, arg = 2, type = "corr",
                           correlation_Structure = "toeplitz",X_Structure = "Normal",
                           epsilon_Structure = "Normal"))
system.time(simulation_par(rep_num, n, p, p0, q, arg = 8, type = "corr",
                           correlation_Structure = "toeplitz",X_Structure = "Normal",
                           epsilon_Structure = "Normal"))

# B: repeat times to calculate Jaccard Index, B=50 in paper
B = 5
# Figure 2 Jaccard Index in paper
system.time(simulation_par_jaccard(rep_num, B, n, p, p0, q, arg = 2, type = "corr",
                           correlation_Structure = "toeplitz",X_Structure = "Normal",
                           epsilon_Structure = "Normal"))
system.time(simulation_par_jaccard(rep_num, B, n, p, p0, q, arg = 8, type = "corr",
                           correlation_Structure = "toeplitz",X_Structure = "Normal",
                           epsilon_Structure = "Normal"))



############# plot ##############

library(ggplot2)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(scales)

plot_result <- function(n, p, arg){
  load(paste("results/results_figure2/delta", arg, "_n", n, 
             "_p", p, "_corr_X_Normal_Xcor_toeplitz_eps_Normal.RData", sep = ""))
  
  df_jaccard <- data_save[1:8,c(1,15:17,4:6,3,2,13,14)]
  colnames(df_jaccard)=c("rho","DS","MDS", "derand_DS","stab_mean","stab_median","stab_rank_mean","stab_sel_prob","stab_e_avg","stab_MM","stab_EMM")
  df_jaccard <- melt(df_jaccard,id.vars="rho",variable.name="method",value.name = "Jaccard")
  
  load(paste("Z:/User/JiajunSun/Stability e-BH/simulations/result_fdp_power_g/delta", arg, "_n", n, 
             "_p", p, "_corr_X_Normal_Xcor_toeplitz_eps_Normal.RData", sep = ""))
  
  df_fdp <- data_save[1:8,c(1,2*c(15,17,16,4:6,3,2,13,14)-2)]
  colnames(df_fdp)=c("rho","DS","MDS", "derand_DS","stab_mean","stab_median","stab_rank_mean","stab_sel_prob","stab_e_avg","stab_MM","stab_EMM")
  
  df_power <- data_save[1:8,c(1,2*c(15,17,16,4:6,3,2,13,14)-1)]
  colnames(df_power)=c("rho","DS","MDS", "derand_DS","stab_mean","stab_median","stab_rank_mean","stab_sel_prob","stab_e_avg","stab_MM","stab_EMM")
  
  df_fdp <- melt(df_fdp,id.vars="rho",variable.name="method",value.name = "FDR")
  df_power <- melt(df_power,id.vars="rho",variable.name="method",value.name = "Power")
  
  n_method <- length(unique(df_fdp$method))
  shapes <- rep(15:19, length.out = n_method)
  linetypes <- rep(c("solid", "dashed", "dotted", "dotdash",
                     "longdash", "twodash"), length.out = n_method)
  
  plot_jaccard <- ggplot(df_jaccard, aes(rho,Jaccard, colour = method, linetype = method, shape = method)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) + scale_linetype_manual(values = linetypes) +scale_shape_manual(values = shapes) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    theme_bw(base_size = 14) +
    theme(legend.position = "right",legend.key.size = unit(1, "cm"))
  
  
  
  plot_fdp <-ggplot(df_fdp, aes(rho, FDR, colour = method, linetype = method, shape = method)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.1, linetype = "dotted", colour = "red") +
    scale_linetype_manual(values = linetypes) +scale_shape_manual(values = shapes) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    theme_bw(base_size = 14) +
    theme(legend.position = "right",legend.key.size = unit(1, "cm"))
  
  plot_power <- ggplot(df_power, aes(rho,Power, colour = method, linetype = method, shape = method)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) + scale_linetype_manual(values = linetypes) +scale_shape_manual(values = shapes) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    theme_bw(base_size = 14) +
    theme(legend.position = "right",legend.key.size = unit(1, "cm"))
  
  g <- ggarrange(plot_fdp,plot_power,plot_jaccard,common.legend = T,legend = "right",ncol = 3, nrow = 1)
  
  ggsave(file=paste("figures/figure2/delta", arg, "_n", n, 
                    "_p", p, "_corr_X_Normal_Xcor_toeplitz_eps_Normal.pdf", sep = ""), g, width = 14, height = 5)
  
  g
}





n=800
p=2000
arg=2
plot_result(n=n, p=p, arg=arg)

arg=8
plot_result(n=n, p=p, arg=arg)



# FDR and Power results in Supplementary Material
# n <- 800
# p <- 1000
# p0 <- 50
# q <- 0.1
# system.time(simulation_par(rep_num, n, p, p0, q, arg = 2, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# system.time(simulation_par(rep_num, n, p, p0, q, arg = 8, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par_jaccard(rep_num, n, p, p0, q, arg = 2, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# system.time(simulation_par_jaccard(rep_num, n, p, p0, q, arg = 8, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# n <- 500
# p <- 500
# p0 <- 50
# q <- 0.1
# system.time(simulation_par(rep_num, n, p, p0, q, arg = 2, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# system.time(simulation_par(rep_num, n, p, p0, q, arg = 8, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# system.time(simulation_par_jaccard(rep_num, n, p, p0, q, arg = 2, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# system.time(simulation_par_jaccard(rep_num, n, p, p0, q, arg = 8, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))

sfStop()



