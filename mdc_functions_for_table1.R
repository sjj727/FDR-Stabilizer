library(randomForest)
library(e1071)
library(EDMeasure)
library(glmnet)
library(xgboost)

mdc_test = function(X, Y, B)
{
  bstp = rep(0, B)
  nl = nrow(X)
  for(ii in 1:B)
  {
    ind1 = sample(1:nl, nl, replace = F)
    bstp[ii] = mdc(X[ind1,], Y, center = "U")
  }
  mean(mdc(X, Y, center = "U")<bstp)
}


getX_SVG = function(image_data, positions){
  
  n = nrow(positions)
  Xmat = cbind(positions, matrix(0, n, 3))
  for(ii in 1:n)
  {
    x = positions[ii, 2]
    y = positions[ii, 3]
    x = round(as.numeric(x))
    y = round(as.numeric(y))
    ind1 = (x-25):(x+25)
    ind2 = (y-25):(y+25)
    Xmat[ii, 4:6] = c(mean(image_data[ind1,ind2, 1]), mean(image_data[ind1,ind2, 2]), 
                      mean(image_data[ind1,ind2, 3]))
  }
  rownames(Xmat) = NULL
  colnames(Xmat) = paste("V", 1:6, sep = "")
  data.frame(Xmat)
}




SVG_test_combine_ab = function(Xa, Xb, Y, testMethod = "two_sample", split.times = 20)
{
  n = nrow(Xa); n1 = n/2
  pval.all <- mse <- mse2 <- pea <-  rep(0, split.times)
  for(ii in 1:split.times){
    ind1 = sample(1:n, size = n1, replace = F)
    Xa1 = Xa[ind1,]
    Xa2 = Xa[-ind1,]
    Xb1 = Xb[ind1,]
    Xb2 = Xb[-ind1,]
    Y1 = Y[ind1]
    Y2 = Y[-ind1]
    
    dtrain <- xgb.DMatrix(data = Xa1, label = Y1)
    dtest <- xgb.DMatrix(data = Xa2)
    bstDMatrix <- xgboost(data = dtrain, max.depth = 2, eta = 1, verbose = F,
                          nthread = 2, nrounds = 50, objective = "reg:squarederror")
    yhat1 = predict(bstDMatrix, newdata = dtest)
    
    dtrain2 <- xgb.DMatrix(data = Xb1, label = Y1)
    dtest2 <- xgb.DMatrix(data = Xb2)
    bstDMatrix2 <- xgboost(data = dtrain2, max.depth = 2, eta = 1, verbose = F,
                           nthread = 2, nrounds = 50, objective = "reg:squarederror")
    yhat2 = predict(bstDMatrix2, newdata = dtest2)
    
    ep1sq = (Y2 - yhat1)^2
    ep2sq = (Y2 - yhat2)^2
    
    if(testMethod=="two_sample"){
      ts = mean(ep1sq - ep2sq)
      pvalue = getPvalueCCK(ep1sq, ep2sq, B = 1000)
    }
    if(testMethod=="rank_sum"){
      ts = rank_sum_test(ep1sq, ep2sq)
      pvalue = getPvalueRS_Boot(ep1sq, ep2sq, B = 1000)
    }
    
    pval.all[ii] = pvalue
    mse[ii] = mean(ep1sq)
    mse2[ii] = mean(ep2sq)
    pea[ii] = cor(Y2, yhat1)
  }
  
  ### Use Cauchy combination test
  t0 = mean(tan((0.5 - pval.all)*pi))
  pval = 0.5 - atan(t0)/pi
  
  mse_mean = mean(mse)
  mse_med = median(mse)
  list(pval = pval, mse_mean = mse_mean, mse_med = mse_med, pearson_mean = mean(pea), pearson_med = median(pea))
}



SVG_test_combine = function(Xa, Xb, Y, testMethod = "two_sample", samp, B = 1000)
{
  n = nrow(Xa); n1 = n/2; times = ncol(samp)
  pval.all = rep(0, times)
  for(ii in 1:times){
    #ind1 = sample(1:n, size = n1, replace = F)
    ind1 = samp[,ii]
    Xa1 = Xa[ind1,]
    Xa2 = Xa[-ind1,]
    Xb1 = Xb[ind1,]
    Xb2 = Xb[-ind1,]
    Y1 = Y[ind1]
    Y2 = Y[-ind1]
    
    dtrain <- xgb.DMatrix(data = Xa1, label = Y1)
    dtest <- xgb.DMatrix(data = Xa2)
    bstDMatrix <- xgboost(data = dtrain, max.depth = 2, eta = 1, verbose = F,
                          nthread = 2, nrounds = 50, objective = "reg:squarederror")
    yhat1 = predict(bstDMatrix, newdata = dtest)
    
    dtrain2 <- xgb.DMatrix(data = Xb1, label = Y1)
    dtest2 <- xgb.DMatrix(data = Xb2)
    bstDMatrix2 <- xgboost(data = dtrain2, max.depth = 2, eta = 1, verbose = F,
                           nthread = 2, nrounds = 50, objective = "reg:squarederror")
    yhat2 = predict(bstDMatrix2, newdata = dtest2)
    
    ep1sq = (Y2 - yhat1)^2
    ep2sq = (Y2 - yhat2)^2
    
    if(testMethod=="two_sample"){
      ts = mean(ep1sq - ep2sq)
      pvalue = getPvalueCCK(ep1sq, ep2sq, B = B)
    }
    if(testMethod=="rank_sum"){
      ts = rank_sum_test(ep1sq, ep2sq)
      pvalue = getPvalueRS_Boot(ep1sq, ep2sq, B = B)
    }
    
    pval.all[ii] = pvalue
  }
  
  ### Use Cauchy combination test
  t0 = mean(tan((0.5 - pval.all)*pi))
  CCT <- 0.5 - atan(t0)/pi
  return(list(CCT = CCT, pval.all = pval.all))
}



SVG_test = function(X, Y, testMethod = "two_sample")
{
  n = nrow(X); n1 = n/2
  ind1 = sample(1:n, size = n1, replace = F)
 
  Xa1 = Xa[ind1,]
  Xa2 = Xa[-ind1,]
  Xb1 = Xb[ind1,]
  Xb2 = Xb[-ind1,]
  Y1 = Y[ind1]
  Y2 = Y[-ind1]
  
  dtrain <- xgb.DMatrix(data = Xa1, label = Y1)
  dtest <- xgb.DMatrix(data = Xa2)
  bstDMatrix <- xgboost(data = dtrain, max.depth = 2, eta = 1, verbose = F,
                        nthread = 2, nrounds = 50, objective = "reg:squarederror")
  yhat1 = predict(bstDMatrix, newdata = dtest)
  
  dtrain2 <- xgb.DMatrix(data = Xb1, label = Y1)
  dtest2 <- xgb.DMatrix(data = Xb2)
  bstDMatrix2 <- xgboost(data = dtrain2, max.depth = 2, eta = 1, verbose = F,
                         nthread = 2, nrounds = 50, objective = "reg:squarederror")
  yhat2 = predict(bstDMatrix2, newdata = dtest2)
  
  ep1sq = (Y2 - yhat1)^2
  ep2sq = (Y2 - yhat2)^2
  
  if(testMethod=="two_sample"){
    ts = mean(ep1sq - ep2sq)
    pvalue = getPvalueCCK(ep1sq, ep2sq, B = 1000)
  }
  if(testMethod=="rank_sum"){
    ts = rank_sum_test(ep1sq, ep2sq)
    pvalue = getPvalueRS_Boot(ep1sq, ep2sq, B = 1000)
  }
  
  list(test_stat = ts, pvalue = pvalue)
}


getPvalueRS_Boot = function(s1, s2, B = 1000)
{
  ts = rank_sum_test(s1, s2)
  nn = length(s1)
  s = c(s1, s2)
  ts_boot = rep(0, B)
  for (ii in 1:B) {
    ind1 = sample(1:(2*nn), nn, replace = F)
    s1b = s[ind1]
    s2b = s[-ind1]
    ts_boot[ii] = rank_sum_test(s1b, s2b)
  }
  1 - pnorm(ts, mean = 0.5, sd = sd(ts_boot))
}

getPvalueCCK = function(s1, s2, B = 1000)
{
  nn = length(s1)
  muhat = mean(s1 - s2)
  epmu = s1-s2-muhat
  epmat = matrix(rep(epmu, B), nn, B)*matrix(rnorm(nn*B), nn, B)
  tsb = colMeans(epmat)
  #mean(tsb<=muhat)
  pnorm(muhat, mean = 0, sd = sd(tsb))
}

getPvaluePermute = function(X, Y, muhat, B = 1000)
{
  n = nrow(X)
  n1 = n/2
  tsb = rep(0, B)
  ind1 = sample(1:n, size = n1, replace = F)
  for(bb in 1:B){
    
    Y = Y[sample(1:n, n, replace = F)]
    X1 = X[ind1,]
    X2 = X[-ind1,]
    Y1 = Y[ind1]
    Y2 = Y[-ind1]
    train = data.frame(X=X1, Y = Y1)
    l1 = randomForest(Y~., data = train, ntree = 100)
    yhat = predict(l1, newdata = data.frame(X=X2))
    r2 = 1 - sum((Y2 - yhat)^2)/sum((Y2 - mean(Y2))^2)
    ep1sq = (Y2 - yhat)^2
    ep2sq = (Y2 - mean(Y2))^2
    tsb[bb] = rank_sum_test(ep1sq, ep2sq)
  }
  hist(tsb)
  mean(tsb<=muhat)
}

rank_sum_test = function(s1, s2)
{
  m = length(s1)
  ss = c(s1, s2)
  s = rank(ss)[1:m]
  sp = sum(s) - m*(m+1)/2
  1 - sp/(m^2)
}

rank_sum_test_slow = function(s1, s2)
{
  m = length(s1)
  s = 0
  for (ii in 1:m) {
    for (jj in 1:m) {
      s = s+(s1[ii]<s2[jj])
    }
  }
  s/(m^2)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

create_sample <- function(split.time, n){
  n1 <- ceiling(n/2)
  result_matrix <- matrix(nrow = n1, ncol = split.time)
  for (i in 1:split.time) {
    result <- sample(1:n, size = n1, replace = F)  
    result_matrix[,i] <- result
  }
  return(result_matrix)
}

get_result <- function(res, q, split_times){
  CCT_all <- list.map(list.map(res,result_all),CCT)
  CCT_marker <- list.map(list.map(res,result_marker),CCT)
  res_CCT_all <- which(p.adjust(CCT_all, method = "BY")<q)
  res_CCT_marker <- which(p.adjust(CCT_marker, method = "BY")<q)
  
  stab_all <- matrix(unlist(list.map(list.map(res,result_all),pval.all)),nrow=split_times,ncol=226)
  stab_marker <- matrix(unlist(list.map(list.map(res,result_marker),pval.all)),nrow=split_times,ncol=226)
  
  e_all <- matrix(0,nrow=split_times,ncol=226)
  e_marker <- matrix(0,nrow=split_times,ncol=226)
  s_all <- 0
  s_marker <- 0
  for (b in 1:split_times){
    select_all <- which(p.adjust(stab_all[b,], method = "BY")<q)
    select_marker <- which(p.adjust(stab_marker[b,], method = "BY")<q)
    s_all <- s_all + length(select_all)
    s_marker <- s_marker + length(select_marker)
    for (protein in 1:226) {
      if(protein %in% select_all){
        e_all[b,protein] <- 226/(length(select_all)*q)
      }
      if(protein %in% select_marker){
        e_marker[b,protein] <- 226/(length(select_marker)*q)
      }
    }
  }
  s_all <- s_all/split_times
  s_marker <- s_marker/split_times
  res_stab_all <- which(colMeans(e_all) >= sort(colMeans(e_all),decreasing = T)[ceiling(s_all)])
  res_stab_marker <- which(colMeans(e_marker) >= sort(colMeans(e_marker),decreasing = T)[ceiling(s_marker)])
  
  return(list(res_stab_all=res_stab_all,res_stab_marker=res_stab_marker,
              res_CCT_all=res_CCT_all,res_CCT_marker=res_CCT_marker))
}
