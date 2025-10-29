rm(list = ls())
library(snowfall)
library(hdi)
library(here)
library(corpcor)
library(glmnet)
library(knockoff)
library(Rcpp)
library(RcppArmadillo)
### read files
setwd(here::here())
folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  source(file_path)
}

sourceCpp("utils/solve_linear_equation.cpp")
sourceCpp("utils/cholcpp.cpp")





Fun_SNP <- function(num, M, X, y, signal_index, p, p0, q){
  set.seed(num)
  signal_index <<- signal_index
  
  DS_result <- DS(X,y, M, q)
  
  derand_kn_result <- derand_kn(X, y, M, q)
  
  knockoff_result <- knockoff(X, y, q, M=M)
  
  BH_result <- MBHq(X, y, q, M)
  
  data_save <- as.data.frame(list(DS_fdp  = DS_result$DS_fdp,  DS_power  = DS_result$DS_power,
                                  MDS_fdp = DS_result$MDS_fdp, MDS_power = DS_result$MDS_power, 
                                  stab_DS_fdp = DS_result$stab_DS_fdp, stab_DS_power = DS_result$stab_DS_power,
                                  knockoff_fdp = knockoff_result$fdp, knockoff_power = knockoff_result$power,
                                  derand_kn_fdp = derand_kn_result$fdp, derand_kn_power = derand_kn_result$power,
                                  stab_kn_fdp = knockoff_result$stab_fdp, stab_kn_power = knockoff_result$stab_power,
                                  BH_fdp = BH_result$BHfdp, BH_power = BH_result$BHpower,
                                  MBH_fdp = BH_result$fdp, MBH_power = BH_result$power,
                                  BH_stab_fdp = BH_result$stab_fdp, BH_stab_power = BH_result$stab_power
  ))
  return(data_save)
}

simulation_SNP_par <- function(rep_num, M, c, p, p0, q){
  set.seed(1234)
  num <- 1:rep_num
  load("GWAS data/tomato.RData")
  X <- tomato[,sample(1:ncol(tomato),p,replace = F)]
  X <- as.matrix(X)
  n <- nrow(X)
  beta_star <- rep(0, p)
  signal_index <- sample(c(1:p), size = p0, replace = F)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = c^2/n)
  
  y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
  
  data_save <- sfLapply(num, fun = Fun_SNP, M, X, y, signal_index, p, p0, q)
  
  save(data_save, file = paste("results/results_figure5/p",p,"_s",p0,"_c",c,".RData", sep = ""))
  return(data_save)
}


sfInit(parallel = TRUE, cpus = 50)

sfLibrary(snowfall)
sfLibrary(hdi)
sfLibrary(corpcor)
sfLibrary(glmnet)
sfLibrary(knockoff)
sfLibrary(Rcpp)
sfLibrary(RcppArmadillo)

folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  sfSource(file_path)
}
sfClusterEval(sourceCpp("utils/solve_linear_equation.cpp"))
sfClusterEval(sourceCpp("utils/cholcpp.cpp"))

sfExport("Fun_SNP")

# rep_num = 100 in paper, 50 for fast
rep_num = 50
# M=50 in paper, 20 for fast
M = 20

system.time(simulation_SNP_par(rep_num, M, 50,1000,60,0.1))
system.time(simulation_SNP_par(rep_num, M, 50,1000,80,0.1))

sfStop()






# plot Figure 5

library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)



plot_box <- function(c, p, p0){
  load(paste("results/results_figure5/p",p,"_s",p0,"_c",c,".RData", sep = ""))
  res_bind <- do.call(rbind, data_save)
  res_bind_power <- res_bind[,c(14,16,18,8,10,12,2,4,6)]
  res_bind_fdp <- res_bind[,c(13,15,17,7,9,11,1,3,5)]
  res_fdp <- melt(res_bind_fdp,variable.name="method",value.name = "FDP")
  plt_fdp = ggplot(res_fdp, aes(x=method, y = FDP)) + 
    geom_boxplot()+    
    stat_boxplot(geom = "errorbar",width=0.25)+
    scale_x_discrete(labels=c("BH", "MBH", "stab (BH)", "knockoff", "derand_kn", "stab (kn)", "DS", "MDS", "stab (DS)"))+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5))+labs(title = paste("s = ", p0, sep=""))+
    geom_hline(yintercept=0.1, linetype='dotted', col = 'red')
  res_power <- melt(res_bind_power,variable.name="method",value.name = "Power")
  plt_power = ggplot(res_power, aes(x=method, y = Power)) + 
    geom_boxplot()+    
    stat_boxplot(geom = "errorbar",width=0.25)+
    scale_x_discrete(labels=c("BH", "MBH", "stab (BH)", "knockoff", "derand_kn", "stab (kn)", "DS", "MDS", "stab (DS)"))+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5))+labs(title = paste("s = ", p0, sep=""))
  ggsave(file=paste("figures/figure5/SNP_p",p,"_s",p0,"_c",c, 
                    "_fdp.pdf", sep = ""), plt_fdp, width = 10, height = 3)
  ggsave(file=paste("figures/figure5/SNP_p",p,"_s",p0,"_c",c, 
                    "_power.pdf", sep = ""), plt_power, width = 10, height = 3)
}

plot_box(50,1000,60)

plot_box(50,1000,80)





