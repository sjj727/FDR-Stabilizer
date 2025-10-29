### data-splitting methods (DS and MDS)
DS <- function(X, y, num_split, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  MM <- matrix(rep(0,num_split*p),nrow=num_split,ncol=p)
  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    
    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.min
    
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      
      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      M[which(is.na(M)==1)] <- 0
      #M <- abs(beta1 + beta2) - abs(beta1 - beta2)
      selected_index <- analys(M, abs(M), q)
      MM[iter,] <- M
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
        
        ### calculate fdp and power
        result <- fdp_power(selected_index, signal_index)
        fdp[iter] <- result$fdp
        power[iter] <- result$power
      }
    }
  }
  
  ### single data-splitting (DS) result
  DS_fdp <- fdp[1]
  DS_power <- power[1]
  

  selectnum <- rep(0,num_split)
  for (k in 1:(num_split)) {
   tau <- 100
   for (t in sort(abs(MM[k,which(abs(MM[k,])!=0)])+1e-8)){
     fdp_hat <- (1+length(which(MM[k,]<= -t)))/max(length(which(MM[k,]>=t)),1)
     if (fdp_hat<=q){
       tau <- t
       break
     }

   }
   selnum <- which(MM[k,]>=tau)
   selectnum[k] <- length(selnum)
  }
  num <- sum(selectnum)/num_split

  ### multiple data-splitting (MDS) result
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
    #prevent too much ties
    if(length(which(inclusion_rate == sort(inclusion_rate,decreasing = T)[ceiling(num)]))-1>=num/5){
      inclusion_rate <- inclusion_rate+runif(length(inclusion_rate),0,1e-6)
    }
    selected_index_stab <- which(inclusion_rate>=sort(inclusion_rate,decreasing = T)[ceiling(num)])
    ### calculate fdp and power
    result <- fdp_power(selected_index, signal_index)
    MDS_fdp <- result$fdp
    MDS_power <- result$power
    result1 <- fdp_power(selected_index_stab, signal_index)
    DS_stab_fdp <- result1$fdp
    DS_stab_power <- result1$power
  }
  else{
    selected_index <- NULL
    selected_index_stab <- NULL
    MDS_fdp <- 0
    MDS_power <- 0
    DS_stab_fdp <- 0
    DS_stab_power <- 0
  }
  

  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power,
              stab_DS_fdp = DS_stab_fdp, stab_DS_power = DS_stab_power, 
              DS_num=length(analys(MM[1,], abs(MM[1,]), q)),
              MDS_num=length(selected_index),
              DS_stab_num=length(selected_index_stab)))
}

DS_different_aggr <- function(X, y, num_split, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  MM <- matrix(rep(0,num_split*p),nrow=num_split,ncol=p)
  sel_prob <- matrix(0, nrow = num_split, ncol = p)
  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    
    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.min
    
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      
      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      M[which(is.na(M)==1)] <- 0
      #M <- abs(beta1 + beta2) - abs(beta1 - beta2)
      selected_index <- analys(M, abs(M), q)
      MM[iter,] <- M
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
        sel_prob[iter, selected_index] <- 1
        ### calculate fdp and power
        result <- fdp_power(selected_index, signal_index)
        fdp[iter] <- result$fdp
        power[iter] <- result$power
      }
    }
  }
  
  ### single data-splitting (DS) result
  DS_fdp <- fdp[1]
  DS_power <- power[1]
  
  
  selectnum <- rep(0,num_split)
  for (k in 1:(num_split)) {
    tau <- 100
    for (t in sort(abs(MM[k,which(abs(MM[k,])!=0)])+1e-8)){
      fdp_hat <- (1+length(which(MM[k,]<= -t)))/max(length(which(MM[k,]>=t)),1)
      if (fdp_hat<=q){
        tau <- t
        break
      }
      
    }
    selnum <- which(MM[k,]>=tau)
    selectnum[k] <- length(selnum)
  }
  num <- sum(selectnum)/num_split
  
  sel_prob <- apply(sel_prob, 2, mean)
  ### multiple data-splitting (MDS) result
  inclusion_rate_2 <- apply(inclusion_rate^2, 2, mean)
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  
  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    
    ### backtracking 
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        #l <- feature_index
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    
    selected_index <- setdiff(feature_rank, null_feature)
    
    selected_index_eavg <- which(inclusion_rate>=sort(inclusion_rate,decreasing = T)[ceiling(num)])
    
    selected_index_mean <- which(colMeans(MM)>=sort(colMeans(MM),decreasing = T)[ceiling(num)])
    
    selected_index_median <- which(apply(MM,2,median)>=sort(apply(MM,2,median),decreasing = T)[ceiling(num)])
    
    selected_index_rankmean <- which(apply(apply(t(-MM), 2, rank),1,mean)<=sort(apply(apply(t(-MM), 2, rank),1,mean),decreasing = F)[ceiling(num)])
    
    selected_index_prob <- which(sel_prob>=sort(sel_prob,decreasing = T)[ceiling(num)])
    
    selected_index_eavg2 <- which(inclusion_rate_2>=sort(inclusion_rate_2,decreasing = T)[ceiling(num)])
    ### calculate fdp and power
    result <- fdp_power(selected_index, signal_index)
    MDS_fdp <- result$fdp
    MDS_power <- result$power
    result1 <- fdp_power(selected_index_eavg, signal_index)
    result2 <- fdp_power(selected_index_mean, signal_index)
    result3 <- fdp_power(selected_index_prob, signal_index)
    result4 <- fdp_power(selected_index_median, signal_index)
    result5 <- fdp_power(selected_index_rankmean, signal_index)
    result6 <- fdp_power(selected_index_eavg2, signal_index)
    DS_stab_eavg_fdp <- result1$fdp
    DS_stab_eavg_power <- result1$power
    
    DS_stab_mean_fdp <- result2$fdp
    DS_stab_mean_power <- result2$power
    
    DS_stab_prob_fdp <- result3$fdp
    DS_stab_prob_power <- result3$power
    
    DS_stab_median_fdp <- result4$fdp
    DS_stab_median_power <- result4$power
    
    DS_stab_rankmean_fdp <- result5$fdp
    DS_stab_rankmean_power <- result5$power
    
    DS_stab_eavg2_fdp <- result6$fdp
    DS_stab_eavg2_power <- result6$power
  }
  else{
    selected_index_eavg <- NULL
    selected_index_eavg2 <- NULL
    selected_index_mean <- NULL
    selected_index_prob <- NULL
    selected_index_median <- NULL
    selected_index_rankmean <- NULL
    MDS_fdp <- 0
    MDS_power <- 0
    DS_stab_eavg_fdp <- 0
    DS_stab_eavg_power <- 0
    DS_stab_mean_fdp <- 0
    DS_stab_mean_power <- 0
    DS_stab_prob_fdp <- 0
    DS_stab_prob_power <- 0
    DS_stab_median_fdp <- 0
    DS_stab_median_power <- 0
    DS_stab_rankmean_fdp <- 0
    DS_stab_rankmean_power <- 0
    DS_stab_eavg2_fdp <- 0
    DS_stab_eavg2_power <- 0
  }
  
  

  return(list(DS_fdp = DS_fdp, DS_power = DS_power, MDS_fdp = MDS_fdp, MDS_power = MDS_power,
              stab_DS_eavg_fdp = DS_stab_eavg_fdp, stab_DS_eavg_power = DS_stab_eavg_power,
              stab_DS_mean_fdp = DS_stab_mean_fdp, stab_DS_mean_power = DS_stab_mean_power,
              stab_DS_prob_fdp = DS_stab_prob_fdp, stab_DS_prob_power = DS_stab_prob_power,
              stab_DS_median_fdp = DS_stab_median_fdp, stab_DS_median_power = DS_stab_median_power,
              stab_DS_rankmean_fdp = DS_stab_rankmean_fdp, stab_DS_rankmean_power = DS_stab_rankmean_power,
              stab_DS_eavg2_fdp = DS_stab_eavg2_fdp, stab_DS_eavg2_power = DS_stab_eavg2_power))
}



### data-splitting methods (DS and MDS)
DS_res <- function(X, y, num_split, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  num_select <- rep(0, num_split)
  MM <- matrix(rep(0,num_split*p),nrow=num_split,ncol=p)
  for(iter in 1:num_split){
    while(TRUE){
      sample_index1 <- sample(x = c(1:n), size = trunc(0.5 * n), replace = F)
      sample_index2 <- setdiff(c(1:n), sample_index1)
      cvfit = cv.glmnet(X[sample_index1, ], y[sample_index1], nfolds = 10)
      beta1 = as.vector(coef(cvfit$glmnet.fit, cvfit$lambda.1se))[-1]
      nonzero_index <- which(beta1 != 0)
      if(length(nonzero_index)!=0)
        break
    }
    ### randomly split the data
    
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      
      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      #M <- abs(beta1 + beta2) - abs(beta1 - beta2)
      M[is.na(M)] = 0
      selected_index <- analys(M, abs(M), q)
      MM[iter,] <- M
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
      }
    }
  }
  selectnum <- rep(0,num_split)
  for (k in 1:(num_split)) {
    tau <- 100
    for (t in sort(abs(MM[k,which(abs(MM[k,])!=0)])+1e-8)){
      fdp_hat <- (1+length(which(MM[k,]< -t)))/max(length(which(MM[k,]>t)),1)
      if (fdp_hat<=q){
        tau <- t
        break
      }
      
    }
    selnum <- which(MM[k,]>tau)
    selectnum[k] <- length(selnum)
  }
  num <- sum(selectnum)/num_split

  ### multiple data-splitting (MDS) result
  e_values <- inclusion_rate*p/q
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
    #prevent too much ties
    if(length(which(inclusion_rate == sort(inclusion_rate,decreasing = T)[ceiling(num)]))-1>=num/5){
      inclusion_rate <- inclusion_rate+runif(length(inclusion_rate),0,1e-6)
    }
    selected_index1 <- which(inclusion_rate>=sort(inclusion_rate,decreasing = T)[ceiling(num)])
    
  }
  e_avg <- inclusion_rate*p/q
  
  return(list(stab_sel = selected_index1, e_avg = e_avg, MDS = selected_index, e_values = e_values))
}
