rm(list = ls())
library(snowfall)
library(here)
library(MASS)
library(glmnet)
library(mvtnorm)
library(knockoff)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(Matrix)
setwd(here::here())
folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  source(file_path)
}
sourceCpp("utils/solve_linear_equation.cpp")
sourceCpp("utils/cholcpp.cpp")




simulation_ds <- function(X, y, q){
  ds_result <- DS_res(X, y, 50, q)
  n <- nrow(X)
  p <- ncol(X)
  ### randomly split the data
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  ### get the penalty lambda for Lasso
  cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
  lambda <- cvfit$lambda.min
  ### run Lasso on the first half of the data
  beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
  nonzero_index <- which(beta1 != 0)
  selected_index <- NULL
  if(length(nonzero_index)!=0){
    ### run OLS on the second half of the data, restricted on the selected features
    beta2 <- rep(0, p)
    beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
    ### calculate the mirror statistics
    M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
    M[which(is.na(M)==1)] <- 0
    selected_index <- analys(M, abs(M), q)
  }
  ds_res <- 0
  if(length(selected_index)>0){
    ds_res <- length(selected_index)
  }
  ds_stab_res <- length(ds_result$stab_sel)
  mds_res <- length(ds_result$MDS)
  ds_power <- fdp_power(selected_index, signal_index)$power
  ds_stab_power <- fdp_power(ds_result$stab_sel, signal_index)$power
  mds_power <- fdp_power(ds_result$MDS, signal_index)$power
  ds_fdp <- fdp_power(selected_index, signal_index)$fdp
  ds_stab_fdp <- fdp_power(ds_result$stab_sel, signal_index)$fdp
  mds_fdp <- fdp_power(ds_result$MDS, signal_index)$fdp
  result <- as.data.frame(list(ds_res=ds_res,mds_res=mds_res,ds_stab_res=ds_stab_res,ds_power=ds_power,mds_power=mds_power,ds_stab_power=ds_stab_power,
                               ds_fdp=ds_fdp,mds_fdp=mds_fdp,ds_stab_fdp=ds_stab_fdp))
  return(result)
}

simulation_kn <- function(X, y, q){
  dkn <- derand_kn_res(X,y,50,q)
  stab_kn <- knockoff_res(X,y,q,50)
  res <- length(stab_kn$kn_sel)
  stab_res <- length(stab_kn$stab_sel)
  dkn_res <- length(dkn$rej)
  kn_power <- fdp_power(stab_kn$kn_sel, signal_index)$power
  stab_power <- fdp_power(stab_kn$stab_sel, signal_index)$power
  dkn_power <- fdp_power(dkn$rej, signal_index)$power
  kn_fdp <- fdp_power(stab_kn$kn_sel, signal_index)$fdp
  stab_fdp <- fdp_power(stab_kn$stab_sel, signal_index)$fdp
  dkn_fdp <- fdp_power(dkn$rej, signal_index)$fdp
  result <- as.data.frame(list(kn_res=res,dkn_res=dkn_res,kn_stab_res=stab_res,kn_power=kn_power,dkn_power=dkn_power,kn_stab_power=stab_power,
                               kn_fdp=kn_fdp,dkn_fdp=dkn_fdp,kn_stab_fdp=stab_fdp))
  return(result)
}

Fun <- function(num, X, y, q){
  set.seed(num)
  
  ret_ds <- simulation_ds(X = X, y = y, q = q)
  ret_kn <- simulation_kn(X = X, y = y, q = q)
  ret <- as.data.frame(list(ds=ret_ds,kn=ret_kn))
  return(ret)
  
}

simulation_par <- function(rep_num, n, p, p0, q, arg = 0.5){
  num <- 1:rep_num
  set.seed(1000)
  phi <- arg
  covariance <- matrix(rep(0,p^2),nrow = p)
  for(i in 1 : p){
    for(j in 1 : p){
      covariance[i,j] <- phi ^ (i!= j)
    }
  }
  X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
  beta_star <- rep(0, p)
  signal_index <- sample(c(1:p), size = p0, replace = F)
  beta_star[signal_index] <- runif(p0, 0.1, 1.5)*sample(c(-1, 1), p0, replace = TRUE)
  y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
  
  sfExport("signal_index")
  res <- sfLapply(num, fun = Fun, X, y, q)
  
  
  save(res, file = paste("results/results_figure1/seed1000_unif1.5_CS_rho", arg, ".RData", sep = ""))
  
}




sfInit(parallel = TRUE, cpus = 50)
sfLibrary(snowfall)
sfLibrary(MASS)
sfLibrary(glmnet)
sfLibrary(mvtnorm)
sfLibrary(knockoff)
sfLibrary(ggplot2)
sfLibrary(Rcpp)
sfLibrary(RcppArmadillo)
sfLibrary(Matrix)


folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  sfSource(file_path)
}

sfExport("simulation_ds")
sfExport("simulation_kn")
sfClusterEval(sourceCpp("utils/solve_linear_equation.cpp"))
sfClusterEval(sourceCpp("utils/cholcpp.cpp"))
#rep_num = 1000 in paper, 50 for fast
rep_num = 50

n <- 500
p <- 200
p0 <- 30
q <- 0.1
system.time(simulation_par(rep_num=rep_num, n=n, p=p, p0=p0, q=q, arg = 0.5))

sfStop()



################################# plot Figure 1#############################################################

load("results/results_figure1/seed1000_unif1.5_CS_rho0.5.RData")

res <- do.call(rbind,res)
kn_r <- res[,10:12]
kn_power <- res[,13:15]
colnames(kn_r) <- c("knockoff","derand_knockoff","stab (knockoff)")
colnames(kn_power) <- c("knockoff","derand_knockoff","stab (knockoff)")
kn_power <- colMeans(kn_power)
ds_r <- res[,1:3]
ds_power <- res[,4:6]
colnames(ds_r) <- c("DS","MDS","stab (DS)")
colnames(ds_power) <- c("DS","MDS","stab (DS)")
ds_power <- colMeans(ds_power)
pdf("figures/figure1/Figure_1_kn_motivation_example.pdf",width = 8, height = 5)

par(mar=c(0,0,0,0)+2.2)
boxplot(kn_r, col = "lightblue",ylim=c(18,43))
for (i in 1:3) {
  text(x = i, y = 33, labels = paste("Power:", round(kn_power[i], 5)), pos = 1)
}
dev.off()

pdf("figures/figure1/Figure_1_ds_motivation_example.pdf",width = 8, height = 5)

par(mar=c(0,0,0,0)+2.2)
boxplot(ds_r, col = "lightblue",ylim=c(18,43))
for (i in 1:3) {
  text(x = i, y = 33, labels = paste("Power:", round(ds_power[i], 5)), pos = 1)
}
dev.off()

