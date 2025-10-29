# multiple t test fdr and power
rm(list = ls())
library(MASS)
library(here)
library(Matrix)
library(snowfall)
setwd(here::here())
# aggr_method = c("eavg","mean_rank","median_rank","rank_mean_rank")
stab_eBH <- function(rej_mat, aggr_method="eavg", statistic=NA){
  shat <- rowSums(rej_mat)
  sbar <- ceiling(mean(shat))
  if (aggr_method=="eavg"){
    p <- ncol(rej_mat)
    e <- p*rej_mat/shat
    e[is.nan(e)] <- 0
    pi_tilde <- rank(-colMeans(e))
  }
  if (aggr_method=="mean_rank" & anyNA(statistic)==F){
    pi_tilde <- rank(-colMeans(statistic))
  }
  if (aggr_method=="median_rank" & anyNA(statistic)==F){
    pi_tilde <- rank(-apply(statistic,2,median))
  }
  if (aggr_method=="rank_mean_rank" & anyNA(statistic)==F){
    rank_stat <- matrix(0, nrow = nrow(statistic), ncol = ncol(statistic))
    for (i in 1:nrow(statistic)) {
      rank_stat[i,] <- rank(-statistic[i,])
    }
    pi_tilde <- rank(colMeans(rank_stat))
  }
  return(which(pi_tilde<=sbar))
}


aggr_MDS <- function(rej_mat, q){
  shat <- rowSums(rej_mat)
  inclusion_rate <- rej_mat/shat
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


CCT_BH <- function(pval,q){
  p_CCT <- rep(0, ncol(pval))
  for (i in 1:ncol(pval)) {
    t0 = mean(tan((0.5 - pval[,i])*pi))
    p_CCT[i] <- 0.5 - atan(t0)/pi
  }
  return(which(p.adjust(p_CCT,method = "BH")<=q))
}



gen_data <- function(n, p, s, mu, rho, Xtype="AR1"){
  true_effects <- rep(0, p)
  effect_indices <- sample(1:p, s, replace =  F)
  true_effects[effect_indices] <- rnorm(s, mean = mu, sd = 1)
  if (Xtype=="AR1"){
    cov_matrix <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        cov_matrix[i, j] <- rho^abs(i - j)  # AR(1)
      }
    }
  }
  
  if (Xtype=="toeplitz"){
    sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
    cov_matrix <- bdiag(rep(list(sig1), 10))
    diag(cov_matrix) <- rep(1, p)
  }
  
  if (Xtype=="Compound_Symmetry"){
    cov_matrix <- matrix(rep(0,p^2),nrow = p)
    for(i in 1 : p){
      for(j in 1 : p){
        cov_matrix[i,j] <- rho ^ (i!= j)
      }
    }
  }
  
  data <- mvrnorm(n, mu = true_effects, Sigma = cov_matrix)
  return(list(data=data,signal_index=effect_indices))
}


fdp_power <- function(selected_index, signal_index){
  num_selected <- length(selected_index)
  tp <- length(intersect(selected_index, signal_index))
  fp <- num_selected - tp
  fdp <- fp / max(num_selected, 1)
  power <- tp / length(signal_index)
  return(list(fdp = fdp, power = power))
}


t_test_simulation_corr <- function(n = 800, p = 2000, p0 = 50, q = 0.1, mu = 0.5, M=50, Xtype="AR1"){
  res <- data.frame()
  for (rho in c(1:9)/10) {
    # generate data
    data_lst <- gen_data(n,p,p0,mu = mu,rho = rho,Xtype=Xtype)
    data <- data_lst$data
    signal_index <- data_lst$signal_index
    # multiple subsampling and aggregation
    rej_mat_BH <- matrix(0,M,p)
    pval <- matrix(0,M,p)
    for (m in 1:M) {
      sample_index <- sample(1:n,ceiling(n/2),F)
      # t test
      p_values <- rep(0,p)
      for (i in 1:p) {
        p_values[i] <- t.test(data[sample_index, i])$p.value
      }
      pval[m,] <- p_values
      # BH
      adjusted_p_values <- p.adjust(p_values, method = "BH")
      adjusted_p_values_BY <- p.adjust(p_values, method = "BY")
      rej_mat_BH[m, adjusted_p_values<=q] <- 1
    }
    
    # stabilized e-BH using average e-values
    rej_stab_ebh_eavg <- stab_eBH(rej_mat_BH, aggr_method = "eavg")
    fdp_eBH_eavg <- fdp_power(rej_stab_ebh_eavg, signal_index)$fdp
    power_eBH_eavg <- fdp_power(rej_stab_ebh_eavg, signal_index)$power
    
    # stabilized e-BH using rank of mean rank of statistics
    rej_stab_ebh_rmr <- stab_eBH(rej_mat_BH, aggr_method = "rank_mean_rank", statistic = -pval)
    fdp_eBH_rmr <- fdp_power(rej_stab_ebh_rmr, signal_index)$fdp
    power_eBH_rmr <- fdp_power(rej_stab_ebh_rmr, signal_index)$power
    
    
    rej_MDS <- aggr_MDS(rej_mat_BH, q)
    fdp_MDS <- fdp_power(rej_MDS, signal_index)$fdp
    power_MDS <- fdp_power(rej_MDS, signal_index)$power
    
    
    # fdp and power of BH
    fdp_BH <- fdp_power(which(adjusted_p_values<=q),signal_index)$fdp
    power_BH <- fdp_power(which(adjusted_p_values<=q),signal_index)$power
    
    # fdp and power of BY
    fdp_BY <- fdp_power(which(adjusted_p_values_BY<=q),signal_index)$fdp
    power_BY <- fdp_power(which(adjusted_p_values_BY<=q),signal_index)$power
    
    # Aggregate p values using Cauchy Combination Test and then use BH to control FDR
    rej_CCT <- CCT_BH(pval,q)
    fdp_CCT <- fdp_power(rej_CCT,signal_index)$fdp
    power_CCT <- fdp_power(rej_CCT,signal_index)$power
    
    
    data_save <- as.data.frame(list(rho = rho, 
                                    BH_fdp = fdp_BH, BH_power = power_BH,
                                    BY_fdp = fdp_BY, BY_power = power_BY,
                                    MDS_fdp = fdp_MDS, MDS_power = power_MDS,
                                    Stab_eBH_eavg_fdp = fdp_eBH_eavg, Stab_eBH_eavg_power = power_eBH_eavg,
                                    Stab_eBH_rmr_fdp = fdp_eBH_rmr, Stab_eBH_rmr_power = power_eBH_rmr,
                                    CCT_fdp = fdp_CCT, CCT_power = power_CCT
    ))
    
    res <- rbind(res,data_save)
  }
  return(res)
}


t_test_simulation_signal <- function(n = 800, p = 2000, p0 = 50, q = 0.1, rho = 0.5, M=50, Xtype="AR1"){
  res <- data.frame()
  for (mu in c(3:10)/10) {
    # generate data
    data_lst <- gen_data(n,p,p0,mu = mu,rho = rho,Xtype=Xtype)
    data <- data_lst$data
    signal_index <- data_lst$signal_index
    
    # multiple subsampling and aggregation
    rej_mat_BH <- matrix(0,M,p)
    pval <- matrix(0,M,p)
    for (m in 1:M) {
      sample_index <- sample(1:n,ceiling(n/2),F)
      # t test
      p_values <- rep(0,p)
      for (i in 1:p) {
        p_values[i] <- t.test(data[sample_index, i])$p.value
      }
      pval[m,] <- p_values
      # BH
      adjusted_p_values <- p.adjust(p_values, method = "BH")
      adjusted_p_values_BY <- p.adjust(p_values, method = "BY")
      rej_mat_BH[m, adjusted_p_values<=q] <- 1
    }
    
    # stabilized e-BH using average e-values
    rej_stab_ebh_eavg <- stab_eBH(rej_mat_BH, aggr_method = "eavg")
    fdp_eBH_eavg <- fdp_power(rej_stab_ebh_eavg, signal_index)$fdp
    power_eBH_eavg <- fdp_power(rej_stab_ebh_eavg, signal_index)$power
    
    # stabilized e-BH using rank of mean rank of statistics
    rej_stab_ebh_rmr <- stab_eBH(rej_mat_BH, aggr_method = "rank_mean_rank", statistic = -pval)
    fdp_eBH_rmr <- fdp_power(rej_stab_ebh_rmr, signal_index)$fdp
    power_eBH_rmr <- fdp_power(rej_stab_ebh_rmr, signal_index)$power
    
    
    rej_MDS <- aggr_MDS(rej_mat_BH, q)
    fdp_MDS <- fdp_power(rej_MDS, signal_index)$fdp
    power_MDS <- fdp_power(rej_MDS, signal_index)$power
    
    
    # fdp and power of BH
    fdp_BH <- fdp_power(which(adjusted_p_values<=q),signal_index)$fdp
    power_BH <- fdp_power(which(adjusted_p_values<=q),signal_index)$power
    
    # fdp and power of BY
    fdp_BY <- fdp_power(which(adjusted_p_values_BY<=q),signal_index)$fdp
    power_BY <- fdp_power(which(adjusted_p_values_BY<=q),signal_index)$power
    
    # Aggregate p values using Cauchy Combination Test and then use BH to control FDR
    rej_CCT <- CCT_BH(pval,q)
    fdp_CCT <- fdp_power(rej_CCT,signal_index)$fdp
    power_CCT <- fdp_power(rej_CCT,signal_index)$power
    
    
    data_save <- as.data.frame(list(mu = mu, 
                                    BH_fdp = fdp_BH, BH_power = power_BH,
                                    BY_fdp = fdp_BY, BY_power = power_BY,
                                    MDS_fdp = fdp_MDS, MDS_power = power_MDS,
                                    Stab_eBH_eavg_fdp = fdp_eBH_eavg, Stab_eBH_eavg_power = power_eBH_eavg,
                                    Stab_eBH_rmr_fdp = fdp_eBH_rmr, Stab_eBH_rmr_power = power_eBH_rmr,
                                    CCT_fdp = fdp_CCT, CCT_power = power_CCT
    ))
    
    res <- rbind(res,data_save)
    
  }
  return(res)
}

Fun <- function(num, type = "signal", n, p, p0, q, arg = 0.5, M=50, Xtype="AR1"){
  set.seed(num)
  if (type == "signal"){
    ret <- t_test_simulation_signal(n = n, p = p, p0 = p0, q = q, rho = arg, M=M, Xtype=Xtype)
  }
  if (type == "corr"){
    ret <- t_test_simulation_corr(n = n, p = p, p0 = p0, q = q, mu = arg, M=M, Xtype=Xtype)
  }
  return(ret)
  
}

simulation_par <- function(rep_num, n, p, p0, q, arg = 0.5, type = "signal", M = 50, Xtype="AR1"){
  num <- 1:rep_num
  res <- sfLapply(num, fun = Fun, type, 
                  n, p, p0, q, arg, M, Xtype)
  
  data_save <- Reduce("+",res)/max(num)
  if(type == "corr"){
    save(data_save, file = paste("results/results_Sec_S3.3/mu", arg, "_n", n, 
                                 "_p", p, "_s", p0, "_q",q,"_", Xtype, ".RData", sep = ""))
  }
  if(type == "signal"){
    save(data_save, file = paste("results/results_Sec_S3.3/rho", arg, "_n", n, 
                                 "_p", p, "_s", p0, "_q",q,"_", Xtype, ".RData", sep = ""))
  }
  return(data_save)
}





sfInit(parallel = TRUE, cpus = 50)
sfLibrary(MASS)
sfLibrary(Matrix)

sfExport("stab_eBH")
sfExport("gen_data")
sfExport("CCT_BH")
sfExport("aggr_MDS")
sfExport("fdp_power")
sfExport("t_test_simulation_signal")
sfExport("t_test_simulation_corr")


# rep_num = 1000 in paper, 50 for fast
rep_num = 50
# M=50 in paper, 20 for fast
M = 20

n = 100
p = 1000
p0 = 50
q = 0.1
Xtype = "AR1"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))

Xtype = "toeplitz"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))

n = 100
p = 2000
q = 0.1
Xtype = "AR1"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))

Xtype = "toeplitz"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))

q = 0.2
Xtype = "AR1"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))

Xtype = "toeplitz"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))

q = 0.05
Xtype = "AR1"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))

Xtype = "toeplitz"
system.time(simulation_par(rep_num = rep_num, n = n, p = p, p0 = p0, q = q, arg = 0.5, type = "corr", M = M, Xtype = Xtype))





# Jaccard Index simulation

jaccard_index <- function(set1, set2) {
  intersect_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  return(intersect_size / max(union_size,1))
}

t_test_simulation_stability <- function(n = 800, p = 2000, p0 = 50, q = 0.1, mu = 0.5, rho=0.5, M=50, B=50, Xtype="AR1"){

  eBH_eavg_sel_ind <- eBH_rmr_sel_ind <- CCT_sel_ind <- BH_sel_ind <- BY_sel_ind <- MDS_sel_ind <- list()
  # generate data
  data_lst <- gen_data(n,p,p0,mu = mu,rho = rho, Xtype="AR1")
  data <- data_lst$data
  signal_index <- data_lst$signal_index
  for (b in 1:B) {
    # multiple subsampling and aggregation
    rej_mat_BH <- rej_mat_BH2 <- matrix(0,M,p)
    pval <- matrix(0,M,p)
    for (m in 1:M) {
      sample_index <- sample(1:n,ceiling(n/2),F)
      # t test
      p_values <- rep(0,p)
      for (i in 1:p) {
        p_values[i] <- t.test(data[sample_index, i])$p.value
      }
      pval[m,] <- p_values
      # BH
      adjusted_p_values <- p.adjust(p_values, method = "BH")
      adjusted_p_values_BY <- p.adjust(p_values, method = "BY")
      rej_mat_BH[m, adjusted_p_values<=q] <- 1
    }
    rej_BH <- which(adjusted_p_values<=q)
    BH_sel_ind[[b]] <- rej_BH
    
    rej_BY <- which(adjusted_p_values_BY<=q)
    BY_sel_ind[[b]] <- rej_BY
    
    rej_MDS <- aggr_MDS(rej_mat_BH, q)
    MDS_sel_ind[[b]] <- rej_MDS
    
    rej_stab_ebh_eavg <- stab_eBH(rej_mat_BH, aggr_method = "eavg")
    eBH_eavg_sel_ind[[b]] <- rej_stab_ebh_eavg
    
    rej_stab_ebh_rmr <- stab_eBH(rej_mat_BH, aggr_method = "rank_mean_rank", statistic = -pval)
    eBH_rmr_sel_ind[[b]] <- rej_stab_ebh_rmr
    
    rej_CCT <- CCT_BH(pval,q)
    CCT_sel_ind[[b]] <- rej_CCT
    
    
  }
  all_jaccard_BH <- combn(BH_sel_ind,2,function(x) jaccard_index(x[[1]], x[[2]]),simplify = F)
  jaccard_BH <- Reduce("+",all_jaccard_BH)/length(all_jaccard_BH)
  
  all_jaccard_BY <- combn(BY_sel_ind,2,function(x) jaccard_index(x[[1]], x[[2]]),simplify = F)
  jaccard_BY <- Reduce("+",all_jaccard_BY)/length(all_jaccard_BY)
  
  all_jaccard_ebh_eavg <- combn(eBH_eavg_sel_ind,2,function(x) jaccard_index(x[[1]], x[[2]]),simplify = F)
  jaccard_ebh_eavg <- Reduce("+",all_jaccard_ebh_eavg)/length(all_jaccard_ebh_eavg)
  
  all_jaccard_ebh_rmr <- combn(eBH_rmr_sel_ind,2,function(x) jaccard_index(x[[1]], x[[2]]),simplify = F)
  jaccard_ebh_rmr <- Reduce("+",all_jaccard_ebh_rmr)/length(all_jaccard_ebh_rmr)
  
  all_jaccard_CCT <- combn(CCT_sel_ind,2,function(x) jaccard_index(x[[1]], x[[2]]),simplify = F)
  jaccard_CCT <- Reduce("+",all_jaccard_CCT)/length(all_jaccard_CCT)
  
  all_jaccard_MDS <- combn(MDS_sel_ind,2,function(x) jaccard_index(x[[1]], x[[2]]),simplify = F)
  jaccard_MDS <- Reduce("+",all_jaccard_MDS)/length(all_jaccard_MDS)
  
  return(as.data.frame(list(jaccard_BH=jaccard_BH, jaccard_BY=jaccard_BY, jaccard_ebh_eavg=jaccard_ebh_eavg, jaccard_ebh_rmr=jaccard_ebh_rmr, 
                            jaccard_CCT=jaccard_CCT, jaccard_ebh_MDS=jaccard_MDS)))
}



Fun_jaccard <- function(num, n, p, p0, q, mu = 0.5, rho=0.5, M=50, B=50, Xtype="AR1"){
  set.seed(num)
  
  ret <- t_test_simulation_stability(n = n, p = p, p0 = p0, q = q, mu = mu, rho=rho, M=M, B=B, Xtype=Xtype)
  
  return(ret)
  
}

simulation_par_jaccard <- function(rep_num, n, p, p0, q, mu = 0.5, rho=0.5, M = 50, B=50, Xtype=Xtype){
  num <- 1:rep_num
  res <- sfLapply(num, fun = Fun_jaccard, n, p, p0, q, mu, rho, M, B, Xtype=Xtype)
  
  data_save <- Reduce("+",res)/max(num)
  
  return(data_save)
}
sfExport("jaccard_index")
sfExport("t_test_simulation_stability")


# rep_num = 1000 in paper, 50 for fast
rep_num = 50
# M=50 in paper, 20 for fast
M = 20
# B represents the number of repeated runs under the same data set, used to calculate the Jaccard index.
# B=50 in paper, 5 for fast 
B = 5
mu = 0.5

Xtype = "toeplitz"

n = 100
p = 1000
p0 = 50
q = 0.1


res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))

n = 100
p = 2000
q = 0.1
res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))

q = 0.2
res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))

q = 0.05
res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))



Xtype = "AR1"
n = 100
p = 1000
p0 = 50
q = 0.1

mu = 0.5
res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))


n = 100
p = 2000
q = 0.1
res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))

q = 0.2
res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))

q = 0.05
res <- data.frame()
for (rho in c(1:9)/10) {
  res <- rbind(res, simulation_par_jaccard(rep_num=rep_num, n = n, p = p, p0 = p0, q = q, mu = mu, rho = rho, M = M, B = B, Xtype=Xtype))
}
data_save <- cbind(rho=c(1:9)/10, res)
save(data_save, file = paste("results/results_Sec_S3.3/jaccard_mu", mu, "_n", n, 
                             "_p", p, "_s", p0, "_q",q,"_",Xtype, ".RData", sep = ""))




# plot 



library(ggplot2)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(scales)

plot_result <- function(n, p, p0, arg, q, correlation_Structure){
  load(paste("results/results_Sec_S3.3/jaccard_mu", arg, "_n", n, 
             "_p", p, "_s", p0, "_q", q, "_", correlation_Structure, ".RData", sep = ""))
  
  df_jaccard <- data_save[,c(1,2,3,6,7,4,5)]
  colnames(df_jaccard)=c("rho","BH","BY", "CCT","MDS","stab e_avg", "stab mean_rank")
  df_jaccard <- melt(df_jaccard,id.vars="rho",variable.name="method",value.name = "jaccard")
  
  load(paste("results/results_Sec_S3.3/mu", arg, "_n", n, 
             "_p", p, "_s", p0, "_q", q, "_", correlation_Structure, ".RData", sep = ""))
  
  df_fdp <- data_save[,c("rho","BH_fdp","BY_fdp","CCT_fdp","MDS_fdp","Stab_eBH_eavg_fdp","Stab_eBH_rmr_fdp")]
  colnames(df_fdp)=c("rho","BH","BY", "CCT","MDS","stab_e_avg", "stab_mean_rank")
  
  df_power <- data_save[,c("rho","BH_power","BY_power","CCT_power","MDS_power","Stab_eBH_eavg_power","Stab_eBH_rmr_power")]
  colnames(df_power)=c("rho","BH","BY", "CCT","MDS","stab_e_avg", "stab_mean_rank")
  
  df_fdp <- melt(df_fdp,id.vars="rho",variable.name="method",value.name = "FDR")
  df_power <- melt(df_power,id.vars="rho",variable.name="method",value.name = "Power")
  
  plot_jaccard <- ggplot(df_jaccard, aes(rho,jaccard))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    scale_y_continuous(labels = label_number(accuracy = 0.01))
  
  plot_fdp <- ggplot(df_fdp, aes(rho,FDR))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    geom_hline(yintercept=q, linetype='dotted', col = 'red')+ 
    scale_y_continuous(labels = label_number(accuracy = 0.01))
  
  plot_power <- ggplot(df_power, aes(rho,Power)) +
    geom_point(aes(color = method, shape = method), cex=2) +
    geom_line(aes(color = method)) + theme_bw()+ 
    scale_y_continuous(labels = label_number(accuracy = 0.01))
  
  
  g <- ggarrange(plot_fdp,plot_power,plot_jaccard,common.legend = T,legend = "top",ncol = 3, nrow = 1)
  
  ggsave(file=paste("figures/figures_Sec_S3.3/mu", arg, "_n", n, 
                    "_p", p, "_s", p0, "_q", q, "_", correlation_Structure, ".pdf", sep = ""), g, width = 10, height = 3)
  
  g
}





n=100
p=2000
p0=50
arg=0.5
q=0.05
correlation_Structure="AR1"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)



n=100
p=2000
p0=50
arg=0.5
q=0.1
correlation_Structure="AR1"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)


n=100
p=2000
p0=50
arg=0.5
q=0.2
correlation_Structure="AR1"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)


n=100
p=2000
p0=50
arg=0.5
q=0.05
correlation_Structure="toeplitz"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)


n=100
p=2000
p0=50
arg=0.5
q=0.1
correlation_Structure="toeplitz"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)


n=100
p=2000
p0=50
arg=0.5
q=0.2
correlation_Structure="toeplitz"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)








n=100
p=1000
p0=50
arg=0.5
q=0.1
correlation_Structure="AR1"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)


n=100
p=1000
p0=50
arg=0.5
q=0.1
correlation_Structure="toeplitz"
plot_result(n=n, p=p, p0=p0, arg=arg, q=q, correlation_Structure=correlation_Structure)




