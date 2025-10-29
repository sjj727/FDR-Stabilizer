rm(list = ls())
library(snowfall)
library(MASS)
library(here)
library(glmnet)
library(knockoff)
library(mvtnorm)
library(hdi)
library(TAF)
library(corpcor)
library(doParallel)
library(plyr)
library(stabs)
library(lars)
library(EnvStats)
library(VGAM)
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



FUN_ds <- function(i, rep_num, X, y){
  set.seed(i)
  lst <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
  MDS_fdp <- MDS_power <- DS_stab_fdp <- DS_stab_power <- DS_num <- MDS_num <- DS_stab_num <- rep(0,rep_num)
  q <- 0.1
  num_split <- lst[i]
  for (t in 1:rep_num){
    ds <- DS(X,y, num_split, q)
    MDS_fdp[t] <- ds$MDS_fdp
    MDS_power[t] <- ds$MDS_power
    DS_stab_fdp[t] <- ds$stab_DS_fdp
    DS_stab_power[t] <- ds$stab_DS_power
    DS_num[t] <- ds$DS_num
    MDS_num[t] <- ds$MDS_num
    DS_stab_num[t] <- ds$DS_stab_num
  }
  E <- c(mean(MDS_fdp),mean(MDS_power),mean(DS_stab_fdp),mean(DS_stab_power),mean(DS_num),mean(MDS_num),mean(DS_stab_num))
  V <- c(var(MDS_fdp),var(MDS_power),var(DS_stab_fdp),var(DS_stab_power),var(DS_num),var(MDS_num),var(DS_stab_num))
  return(list(E=E,V=V))
}



FUN_kn <- function(i, rep_num, X, y){
  set.seed(i)
  lst <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
  deran_kn_fdp <- deran_kn_power <- kn_stab_fdp <- kn_stab_power <- kn_num <- deran_kn_num <- kn_stab_num <- rep(0,rep_num)
  q <- 0.1
  M <- lst[i]
  for (t in 1:rep_num){
    derand_kn <- derand_kn(X,y, M, q)
    kn <- knockoff(X, y, q, M = M)
    deran_kn_fdp[t] <- derand_kn$fdp
    deran_kn_power[t] <- derand_kn$power
    kn_stab_fdp[t] <- kn$stab_fdp
    kn_stab_power[t] <- kn$stab_power
    kn_num[t] <- kn$kn_num
    deran_kn_num[t] <- derand_kn$num
    kn_stab_num[t] <- kn$kn_stab_num
  }
  E <- c(mean(deran_kn_fdp),mean(deran_kn_power),mean(kn_stab_fdp),mean(kn_stab_power),mean(kn_num),mean(deran_kn_num),mean(kn_stab_num))
  V <- c(var(deran_kn_fdp),var(deran_kn_power),var(kn_stab_fdp),var(kn_stab_power),var(kn_num),var(deran_kn_num),var(kn_stab_num))
  return(list(E=E,V=V))
}


FUN_bh <- function(i, rep_num, X, y){
  set.seed(i)
  lst <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
  MBH_fdp <- MBH_power <- BH_stab_fdp <- BH_stab_power <- BH_num <- MBH_num <- BH_stab_num <- rep(0,rep_num)
  q <- 0.1
  num_split <- lst[i]
  for (t in 1:rep_num){
    bh <- MBHq(X, y, q, num_split)
    MBH_fdp[t] <- bh$fdp
    MBH_power[t] <- bh$power
    BH_stab_fdp[t] <- bh$stab_fdp
    BH_stab_power[t] <- bh$stab_power
    BH_num[t] <- bh$num_bh
    MBH_num[t] <- bh$num_mbhq
    BH_stab_num[t] <- bh$num_stab
  }
  E <- c(mean(MBH_fdp),mean(MBH_power),mean(BH_stab_fdp),mean(BH_stab_power),mean(BH_num),mean(MBH_num),mean(BH_stab_num))
  V <- c(var(MBH_fdp),var(MBH_power),var(BH_stab_fdp),var(BH_stab_power),var(BH_num),var(MBH_num),var(BH_stab_num))
  return(list(E=E,V=V))
}


sfInit(parallel = TRUE, cpus = 20)

sfLibrary(snowfall)
sfLibrary(MASS)
sfLibrary(glmnet)
sfLibrary(knockoff)
sfLibrary(mvtnorm)
sfLibrary(TAF)
sfLibrary(hdi)
sfLibrary(doParallel)
sfLibrary(corpcor)
sfLibrary(plyr)
sfLibrary(stabs)
sfLibrary(lars)
sfLibrary(EnvStats)
sfLibrary(VGAM)
sfLibrary(Rcpp)
sfLibrary(RcppArmadillo)
### read files
folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  sfSource(file_path)
}
sfClusterEval(sourceCpp("utils/solve_linear_equation.cpp"))
sfClusterEval(sourceCpp("utils/cholcpp.cpp"))

# rep_num = 100 in paper, 10 for fast
rep_num = 10

#case2 n=500,p=500,p0=50
set.seed(100)
n <- 500
p <- 500
p0 <- 50
q <- 0.1
rho <- 0.4
delta <- 7

covariance <- matrix(rep(0,p^2),nrow = p)
for(i in 1 : p){
  for(j in 1 : p){
    covariance[i,j] <- rho ^ (i!= j)
  }
}

X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
beta_star <- rep(0, p)
signal_index <- sample(c(1:p), size = p0, replace = F)
beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)

sfExport("signal_index")


i <- 1:20
system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
save(res, file = paste("results/results_figure4/case2_bh_n", n, 
                       "_p", p, ".RData", sep = ""))

system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
save(res, file = paste("results/results_figure4/case2_kn_n", n, 
                       "_p", p, ".RData", sep = ""))

system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
save(res, file = paste("results/results_figure4/case2_ds_n", n, 
                       "_p", p, ".RData", sep = ""))


# plot


library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
library(ggpubr)

plot_DS_result <- function(n, p, case){
  load(paste("results/results_figure4/case",case,"_ds_n", n, 
             "_p", p, ".RData", sep = ""))
  
  E <- matrix(0, nrow = 20, ncol = 7)
  V <- matrix(0, nrow = 20, ncol = 7)
  for (i in 1:20) {
    E[i,] <- res[[i]]$E
    V[i,] <- res[[i]]$V
  }
  E <- as.data.frame(E)
  E <- cbind(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),E)
  colnames(E) <- c("M","MDS_fdr","MDS_power","DS_stab_fdr","DS_stab_power","DS_num","MDS_num","DS_stab_num")
  
  V <- as.data.frame(V)
  V <- cbind(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),V)
  colnames(V) <- c("M","MDS_fdr","MDS_power","DS_stab_fdr","DS_stab_power","DS_num","MDS_num","DS_stab_num")
  
  
  E_fdp <- E[,c("M","MDS_fdr","DS_stab_fdr")]
  colnames(E_fdp)=c("M","MDS","stab (DS)")
  
  E_power <- E[,c("M","MDS_power","DS_stab_power")]
  colnames(E_power)=c("M","MDS","stab (DS)")
  
  E_fdp <- melt(E_fdp,id.vars="M",variable.name="method",value.name = "FDR")
  E_power <- melt(E_power,id.vars="M",variable.name="method",value.name = "Power")
  
  plot_E_fdp <- ggplot(E_fdp, aes(M,FDR))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    geom_hline(yintercept=0.1, linetype='dotted', col = 'red')
  
  plot_E_power <- ggplot(E_power, aes(M,Power)) +
    geom_point(aes(color = method, shape = method), cex=2)+
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    geom_line(aes(color = method)) + theme_bw()
  
  
  
  V_fdp <- V[,c("M","MDS_fdr","DS_stab_fdr")]
  colnames(V_fdp)=c("M","MDS","stab (DS)")
  
  V_power <- V[,c("M","MDS_power","DS_stab_power")]
  colnames(V_power)=c("M","MDS","stab (DS)")
  
  V_fdp <- melt(V_fdp,id.vars="M",variable.name="method",value.name = "var_FDR")
  V_power <- melt(V_power,id.vars="M",variable.name="method",value.name = "var_Power")
  
  plot_V_fdp <- ggplot(V_fdp, aes(M,var_FDR))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    scale_y_continuous(labels = label_number(scale = 1e3,accuracy = 0.1))+
    ylab(expression(var_FDR%*%10^-3))
  
  plot_V_power <- ggplot(V_power, aes(M,var_Power)) +
    geom_point(aes(color = method, shape = method), cex=2) +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    geom_line(aes(color = method)) + theme_bw()+
    scale_y_continuous(labels = label_number(scale = 1e3,accuracy = 0.1))+
    ylab(expression(var_Power%*%10^-3))
  
  
  
  
  E_num <- E[,c("M","DS_num","MDS_num","DS_stab_num")]
  colnames(E_num)=c("M","DS","MDS","stab (DS)")
  
  V_num <- V[,c("M","DS_num","MDS_num","DS_stab_num")]
  colnames(V_num)=c("M","DS","MDS","stab (DS)")
  
  E_num <- melt(E_num,id.vars="M",variable.name="method",value.name = "num")
  V_num <- melt(V_num,id.vars="M",variable.name="method",value.name = "var_num")
  
  plot_E_num <- ggplot(E_num, aes(M,num))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw()
  
  plot_V_num <- ggplot(V_num, aes(M,var_num)) +
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw()
  
  g1 <- ggarrange(plot_E_fdp,plot_V_fdp,nrow = 2,common.legend = T,legend = "none")
  g2 <- ggarrange(plot_E_power,plot_V_power,nrow = 2,common.legend = T,legend = "none")
  g3 <- ggarrange(plot_E_num,plot_V_num,nrow = 2,common.legend = T,legend = "none")
  legend_A <- get_legend(plot_E_num, position = "top")
  g4 <- ggarrange(g1,g2,g3,ncol = 3,common.legend = T,legend = "none")
  g4 <- ggarrange(legend_A, g4, nrow = 2, heights = c(0.05, 1))
  
  ggsave(file=paste("figures/figure4/case",case,"_DS_n", n, 
                    "_p", p, ".pdf", sep = ""), g4, width = 7, height = 3)
  g4
}


plot_kn_result <- function(n, p, case){
  load(paste("results/results_figure4/case",case,"_kn_n", n, 
             "_p", p, ".RData", sep = ""))
  
  E <- matrix(0, nrow = 20, ncol = 7)
  V <- matrix(0, nrow = 20, ncol = 7)
  for (i in 1:20) {
    E[i,] <- res[[i]]$E
    V[i,] <- res[[i]]$V
  }
  E <- as.data.frame(E)
  E <- cbind(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),E)
  colnames(E) <- c("M","derand_kn_fdr","derand_kn_power","kn_stab_fdr","kn_stab_power","kn_num","derand_kn_num","kn_stab_num")
  
  V <- as.data.frame(V)
  V <- cbind(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),V)
  colnames(V) <- c("M","derand_kn_fdr","derand_kn_power","kn_stab_fdr","kn_stab_power","kn_num","derand_kn_num","kn_stab_num")
  
  
  E_fdp <- E[,c("M","derand_kn_fdr","kn_stab_fdr")]
  colnames(E_fdp)=c("M","derand_kn","stab (kn)")
  
  E_power <- E[,c("M","derand_kn_power","kn_stab_power")]
  colnames(E_power)=c("M","derand_kn","stab (kn)")
  
  E_fdp <- melt(E_fdp,id.vars="M",variable.name="method",value.name = "FDR")
  E_power <- melt(E_power,id.vars="M",variable.name="method",value.name = "Power")
  
  plot_E_fdp <- ggplot(E_fdp, aes(M,FDR))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+ 
    scale_y_continuous(labels = label_number(accuracy = 0.01))
  
  plot_E_power <- ggplot(E_power, aes(M,Power)) +
    geom_point(aes(color = method, shape = method), cex=2) +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    geom_line(aes(color = method)) + theme_bw()+ 
    scale_y_continuous(labels = label_number(accuracy = 0.001))
  
  
  V_fdp <- V[,c("M","derand_kn_fdr","kn_stab_fdr")]
  colnames(V_fdp)=c("M","derand_kn","stab (kn)")
  
  V_power <- V[,c("M","derand_kn_power","kn_stab_power")]
  colnames(V_power)=c("M","derand_kn","stab (kn)")
  
  V_fdp <- melt(V_fdp,id.vars="M",variable.name="method",value.name = "var_FDR")
  V_power <- melt(V_power,id.vars="M",variable.name="method",value.name = "var_Power")
  
  plot_V_fdp <- ggplot(V_fdp, aes(M,var_FDR))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    theme(plot.title = element_text(hjust = 0.5))+
    #scale_y_continuous(labels = label_number(accuracy = 0.01))+
    #ylim(0,0.1)
    scale_y_continuous(labels = label_number(scale = 1e3,accuracy = 0.1))+
    ylab(expression(var_FDR%*%10^-3))
  
  plot_V_power <- ggplot(V_power, aes(M,var_Power)) +
    geom_point(aes(color = method, shape = method), cex=2) +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    geom_line(aes(color = method)) + theme_bw()+
    scale_y_continuous(labels = label_number(scale = 1e3,accuracy = 0.01))+
    ylab(expression(var_Power%*%10^-3))
  #scale_y_continuous(labels = label_number(scale = 1e2,accuracy = 0.1))+
  #ylab(expression(var_Power%*%10^-2))
  
  E_num <- E[,c("M","kn_num","derand_kn_num","kn_stab_num")]
  colnames(E_num)=c("M","knockoff","derand_kn","stab (kn)")
  
  V_num <- V[,c("M","kn_num","derand_kn_num","kn_stab_num")]
  colnames(V_num)=c("M","knockoff","derand_kn","stab (kn)")
  
  E_num <- melt(E_num,id.vars="M",variable.name="method",value.name = "num")
  V_num <- melt(V_num,id.vars="M",variable.name="method",value.name = "var_num")
  
  plot_E_num <- ggplot(E_num, aes(M,num))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw()
  
  plot_V_num <- ggplot(V_num, aes(M,var_num)) +
    geom_point(aes(color = method, shape = method), cex=2) +
    geom_line(aes(color = method)) + theme_bw()+ 
    scale_y_continuous(labels = label_number(accuracy = 0.1))
  
  
  g1 <- ggarrange(plot_E_fdp,plot_V_fdp,nrow = 2,common.legend = T,legend = "none")
  
  g2 <- ggarrange(plot_E_power,plot_V_power,nrow = 2,common.legend = T,legend = "none")
  
  g3 <- ggarrange(plot_E_num,plot_V_num,nrow = 2,common.legend = T,legend = "none")
  
  legend_A <- get_legend(plot_E_num, position = "top")
  g4 <- ggarrange(g1,g2,g3,ncol = 3,common.legend = T,legend = "none")
  g4 <- ggarrange(legend_A, g4, nrow = 2, heights = c(0.05, 1))
  
  ggsave(file=paste("figures/figure4/case",case,"_kn_n", n, 
                    "_p", p, ".pdf", sep = ""), g4, width = 7, height = 3)
  g4
}


plot_bh_result <- function(n, p, case){
  load(paste("results/results_figure4/case",case,"_bh_n", n, 
             "_p", p, ".RData", sep = ""))
  
  E <- matrix(0, nrow = 20, ncol = 7)
  V <- matrix(0, nrow = 20, ncol = 7)
  for (i in 1:20) {
    E[i,] <- res[[i]]$E
    V[i,] <- res[[i]]$V
  }
  E <- as.data.frame(E)
  E <- cbind(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),E)
  colnames(E) <- c("M","MBHq_fdr","MBHq_power","BH_stab_fdr","BH_stab_power","BH_num","MBHq_num","BH_stab_num")
  
  V <- as.data.frame(V)
  V <- cbind(c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),V)
  colnames(V) <- c("M","MBHq_fdr","MBHq_power","BH_stab_fdr","BH_stab_power","BH_num","MBHq_num","BH_stab_num")
  
  
  E_fdp <- E[,c("M","MBHq_fdr","BH_stab_fdr")]
  colnames(E_fdp)=c("M","MBH","stab (BH)")
  
  E_power <- E[,c("M","MBHq_power","BH_stab_power")]
  colnames(E_power)=c("M","MBH","stab (BH)")
  
  E_fdp <- melt(E_fdp,id.vars="M",variable.name="method",value.name = "FDR")
  E_power <- melt(E_power,id.vars="M",variable.name="method",value.name = "Power")
  
  plot_E_fdp <- ggplot(E_fdp, aes(M,FDR))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))
  
  plot_E_power <- ggplot(E_power, aes(M,Power)) +
    geom_point(aes(color = method, shape = method), cex=2)+
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    geom_line(aes(color = method)) + theme_bw()
  
  
  
  V_fdp <- V[,c("M","MBHq_fdr","BH_stab_fdr")]
  colnames(V_fdp)=c("M","MBH","stab (BH)")
  
  V_power <- V[,c("M","MBHq_power","BH_stab_power")]
  colnames(V_power)=c("M","MBH","stab (BH)")
  
  V_fdp <- melt(V_fdp,id.vars="M",variable.name="method",value.name = "var_FDR")
  V_power <- melt(V_power,id.vars="M",variable.name="method",value.name = "var_Power")
  
  plot_V_fdp <- ggplot(V_fdp, aes(M,var_FDR))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw() +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    scale_y_continuous(labels = label_number(scale = 1e3,accuracy = 0.01))+
    ylab(expression(var_FDR%*%10^-3))
  
  plot_V_power <- ggplot(V_power, aes(M,var_Power)) +
    geom_point(aes(color = method, shape = method), cex=2) +
    scale_shape_manual(values = c(17,15))+
    scale_color_manual(values = c("#00BA38","#619CFF"))+
    geom_line(aes(color = method)) + theme_bw()+
    scale_y_continuous(labels = label_number(scale = 1e3,accuracy = 0.1))+
    ylab(expression(var_Power%*%10^-3))
  
  E_num <- E[,c("M","BH_num","MBHq_num","BH_stab_num")]
  colnames(E_num)=c("M","BH","MBH","stab (BH)")
  
  V_num <- V[,c("M","BH_num","MBHq_num","BH_stab_num")]
  colnames(V_num)=c("M","BH","MBH","stab (BH)")
  
  E_num <- melt(E_num,id.vars="M",variable.name="method",value.name = "num")
  V_num <- melt(V_num,id.vars="M",variable.name="method",value.name = "var_num")
  
  plot_E_num <- ggplot(E_num, aes(M,num))+
    geom_point(aes(color = method, shape = method), cex = 2) +
    geom_line(aes(color = method)) + theme_bw()
  
  plot_V_num <- ggplot(V_num, aes(M,var_num)) +
    geom_point(aes(color = method, shape = method), cex=2) +
    geom_line(aes(color = method)) + theme_bw()
  
  g1 <- ggarrange(plot_E_fdp,plot_V_fdp,nrow = 2,common.legend = T,legend = "none")
  g2 <- ggarrange(plot_E_power,plot_V_power,nrow = 2,common.legend = T,legend = "none")
  g3 <- ggarrange(plot_E_num,plot_V_num,nrow = 2,common.legend = T,legend = "none")
  legend_A <- get_legend(plot_E_num, position = "top")
  g4 <- ggarrange(g1,g2,g3,ncol = 3,common.legend = T,legend = "none")
  g4 <- ggarrange(legend_A, g4, nrow = 2, heights = c(0.05, 1))
  
  ggsave(file=paste("figures/figure4/case",case,"_bh_n", n, 
                    "_p", p, ".pdf", sep = ""), g4, width = 7, height = 3)
  g4
}

n=500
p=500

plot_DS_result(n, p, 2)
plot_kn_result(n, p, 2)
plot_bh_result(n, p, 2)




# Results in Supplementary Material
# 
# #case1 n=800,p=2000,p0=50
# set.seed(100)
# n <- 800
# p <- 2000
# p0 <- 50
# q <- 0.1
# rho <- 0.4
# delta <- 7
# 
# sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
# covariance <- bdiag(rep(list(sig1), 10))
# diag(covariance) <- rep(1, p)
# 
# X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
# beta_star <- rep(0, p)
# signal_index <- sample(c(1:p), size = p0, replace = F)
# beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
# y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
# 
# 
# sfExport("signal_index")
# 
# i <- 1:20
# system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_bh_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_kn_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_ds_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# 
# #case2 n=800,p=2000,p0=50
# set.seed(100)
# n <- 800
# p <- 2000
# p0 <- 50
# q <- 0.1
# rho <- 0.4
# delta <- 7
# 
# covariance <- matrix(rep(0,p^2),nrow = p)
# for(i in 1 : p){
#   for(j in 1 : p){
#     covariance[i,j] <- rho ^ (i!= j)
#   }
# }
# 
# X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
# beta_star <- rep(0, p)
# signal_index <- sample(c(1:p), size = p0, replace = F)
# beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
# y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
# 
# 
# sfExport("signal_index")
# 
# i <- 1:20
# system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_bh_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_kn_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_ds_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# #case1 n=800,p=1000,p0=50
# set.seed(100)
# n <- 800
# p <- 1000
# p0 <- 50
# q <- 0.1
# rho <- 0.4
# delta <- 7
# 
# sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
# covariance <- bdiag(rep(list(sig1), 10))
# diag(covariance) <- rep(1, p)
# 
# X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
# beta_star <- rep(0, p)
# signal_index <- sample(c(1:p), size = p0, replace = F)
# beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
# y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
# 
# 
# sfExport("signal_index")
# 
# i <- 1:20
# system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_bh_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_kn_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_ds_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# 
# #case2 n=800,p=1000,p0=50
# set.seed(100)
# n <- 800
# p <- 1000
# p0 <- 50
# q <- 0.1
# rho <- 0.4
# delta <- 7
# 
# covariance <- matrix(rep(0,p^2),nrow = p)
# for(i in 1 : p){
#   for(j in 1 : p){
#     covariance[i,j] <- rho ^ (i!= j)
#   }
# }
# 
# X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
# beta_star <- rep(0, p)
# signal_index <- sample(c(1:p), size = p0, replace = F)
# beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
# y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
# 
# 
# sfExport("signal_index")
# 
# i <- 1:20
# system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_bh_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_kn_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_ds_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# 
# 
# #case1 n=500,p=500,p0=50
# set.seed(100)
# n <- 500
# p <- 500
# p0 <- 50
# q <- 0.1
# rho <- 0.4
# delta <- 7
# 
# sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
# covariance <- bdiag(rep(list(sig1), 10))
# diag(covariance) <- rep(1, p)
# 
# X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
# beta_star <- rep(0, p)
# signal_index <- sample(c(1:p), size = p0, replace = F)
# beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
# y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
# 
# 
# sfExport("signal_index")
# 
# i <- 1:20
# system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_bh_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_kn_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_ds_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# 
# 
# 
# #case1 n=2000,p=800,p0=50
# set.seed(100)
# n <- 2000
# p <- 800
# p0 <- 50
# q <- 0.1
# rho <- 0.4
# delta <- 7
# 
# sig1 <- toeplitz(seq(rho, 0, length.out = p/10))
# covariance <- bdiag(rep(list(sig1), 10))
# diag(covariance) <- rep(1, p)
# 
# X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
# beta_star <- rep(0, p)
# signal_index <- sample(c(1:p), size = p0, replace = F)
# beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
# y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
# 
# 
# sfExport("signal_index")
# 
# i <- 1:20
# system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_bh_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_kn_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case1_ds_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# 
# 
# #case2 n=2000,p=800,p0=50
# set.seed(100)
# n <- 2000
# p <- 800
# p0 <- 50
# q <- 0.1
# rho <- 0.4
# delta <- 7
# 
# covariance <- matrix(rep(0,p^2),nrow = p)
# for(i in 1 : p){
#   for(j in 1 : p){
#     covariance[i,j] <- rho ^ (i!= j)
#   }
# }
# 
# X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
# beta_star <- rep(0, p)
# signal_index <- sample(c(1:p), size = p0, replace = F)
# beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))
# y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)
# 
# 
# sfExport("signal_index")
# 
# i <- 1:20
# system.time(res <- sfLapply(i,fun=FUN_bh,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_bh_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_kn,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_kn_n", n, 
#                        "_p", p, ".RData", sep = ""))
# 
# system.time(res <- sfLapply(i,fun=FUN_ds,rep_num,X,y))
# save(res, file = paste("results/results_figure4/case2_ds_n", n, 
#                        "_p", p, ".RData", sep = ""))

sfStop()
