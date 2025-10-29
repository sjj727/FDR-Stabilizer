### High dimension linear model
rm(list = ls())
library(snowfall)
library(here)
library(MASS)
library(glmnet)
library(knockoff)
library(hdi)
library(TAF)
library(mvtnorm)
library(plyr)
library(doParallel)
library(stabs)
library(lars)
library(EnvStats)
library(corpcor)
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

#sourceCpp("utils/solve_linear_equation.cpp")
#sourceCpp("utils/cholcpp.cpp")


### correlation_Structure = "toeplitz","Compound_Symmetry"
### X_Structure = "Normal","Mix_Normal","t"
### epsilon_Structure = "Normal","Mix_Normal","Heteroscedasticity"

# fix delta, change rho
simulation_corr <- function(M = 50, n = 800, p = 2000, p0 = 50, q = 0.1, delta = 5,
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

    
    ###
    
    DS_result <- DS(X,y, M, q)
    
    derand_kn_result <- derand_kn(X, y, M, q)
    
    knockoff_result <- knockoff(X, y, q, M=M)
    
    BH_result <- MBHq(X, y, q, M)
    
    data_save <- as.data.frame(list(rho = rho, DS_fdp  = DS_result$DS_fdp,  DS_power  = DS_result$DS_power,
                                    MDS_fdp = DS_result$MDS_fdp, MDS_power = DS_result$MDS_power, 
                                    stab_DS_fdp = DS_result$stab_DS_fdp, stab_DS_power = DS_result$stab_DS_power,
                                    knockoff_fdp = knockoff_result$fdp, knockoff_power = knockoff_result$power,
                                    derand_kn_fdp = derand_kn_result$fdp, derand_kn_power = derand_kn_result$power,
                                    stab_kn_fdp = knockoff_result$stab_fdp, stab_kn_power = knockoff_result$stab_power,
                                    BH_fdp = BH_result$BHfdp, BH_power = BH_result$BHpower,
                                    MBH_fdp = BH_result$fdp, MBH_power = BH_result$power,
                                    BH_stab_fdp = BH_result$stab_fdp, BH_stab_power = BH_result$stab_power))
    
    res <- rbind(res,data_save)
  }
  return(res)
}

# fix rho, change delta
simulation_signal <- function(M = 50, n = 800, p = 2000, p0 = 50, q = 0.1, rho = 0.5, 
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
    
    
    DS_result <- DS(X,y, M, q)
    
    derand_kn_result <- derand_kn(X, y, M, q)
    
    knockoff_result <- knockoff(X, y, q, M = M)
    
    BH_result <- MBHq(X, y, q, M)
    
    data_save <- as.data.frame(list(delta = delta, DS_fdp  = DS_result$DS_fdp,  DS_power  = DS_result$DS_power,
                                    MDS_fdp = DS_result$MDS_fdp, MDS_power = DS_result$MDS_power, 
                                    stab_DS_fdp = DS_result$stab_DS_fdp, stab_DS_power = DS_result$stab_DS_power,
                                    knockoff_fdp = knockoff_result$fdp, knockoff_power = knockoff_result$power,
                                    derand_kn_fdp = derand_kn_result$fdp, derand_kn_power = derand_kn_result$power,
                                    stab_kn_fdp = knockoff_result$stab_fdp, stab_kn_power = knockoff_result$stab_power,
                                    BH_fdp = BH_result$BHfdp, BH_power = BH_result$BHpower,
                                    MBH_fdp = BH_result$fdp, MBH_power = BH_result$power,
                                    BH_stab_fdp = BH_result$stab_fdp, BH_stab_power = BH_result$stab_power))
    
    res <- rbind(res,data_save)
    
  }
  return(res)
}

Fun <- function(num, type = "signal", M, n, p, p0, q, arg = 0.5,
                correlation_Structure = "toeplitz", 
                X_Structure = "Normal", 
                epsilon_Structure = "Normal"){
  set.seed(num)
  if (type == "signal"){
    ret <- simulation_signal(M = M, n = n, p = p, p0 = p0, q = q, rho = arg, 
                             correlation_Structure = correlation_Structure, 
                             X_Structure = X_Structure,
                             epsilon_Structure = epsilon_Structure)
  }
  if (type == "corr"){
    ret <- simulation_corr(M = M, n = n, p = p, p0 = p0, q = q, delta = arg, 
                           correlation_Structure = correlation_Structure, 
                           X_Structure = X_Structure,
                           epsilon_Structure = epsilon_Structure)
  }
  return(ret)
  
}

simulation_par <- function(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal", 
                           correlation_Structure = "toeplitz", 
                           X_Structure = "Normal",
                           epsilon_Structure = "Normal"){
  
  num <- 1:rep_num
  
  res <- sfLapply(num, fun = Fun, type, M,
                  n, p, p0, q, arg, 
                  correlation_Structure, 
                  X_Structure, 
                  epsilon_Structure)
  
  data_save <- Reduce("+",res)/max(num)
  if(type == "corr"){
    save(data_save, file = paste("results/results_figure3/delta", arg, "_n", n, 
                                 "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                                 "_eps_", epsilon_Structure, ".RData", sep = ""))
  }
  if(type == "signal"){
    save(data_save, file = paste("results/results_figure3/rho", arg, "_n", n, 
                                 "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                                 "_eps_", epsilon_Structure, ".RData", sep = ""))
  }
  
}



sfInit(parallel = TRUE, cpus = 50)

sfLibrary(snowfall)
sfLibrary(MASS)
sfLibrary(glmnet)
sfLibrary(knockoff)
sfLibrary(mvtnorm)
sfLibrary(TAF)
sfLibrary(hdi)
sfLibrary(plyr)
sfLibrary(stabs)
sfLibrary(lars)
sfLibrary(doParallel)
sfLibrary(EnvStats)
sfLibrary(VGAM)
sfLibrary(corpcor)
sfLibrary(Rcpp)
sfLibrary(RcppArmadillo)

folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  sfSource(file_path)
}
sfClusterEval(sourceCpp("utils/solve_linear_equation.cpp"))
sfClusterEval(sourceCpp("utils/cholcpp.cpp"))
sfExport("simulation_signal")
sfExport("simulation_corr")

# rep_num = 1000 in paper, 50 for fast
rep_num = 50
# M = 50 in paper, 10 for fast
M = 10

n <- 500
p <- 500
p0 <- 50
q <- 0.1

system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
                           correlation_Structure = "toeplitz",X_Structure = "Normal",
                           epsilon_Structure = "Normal"))

n <- 3000
p <- 500
p0 <- 50
q <- 0.1
system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
                           correlation_Structure = "toeplitz",X_Structure = "Normal",
                           epsilon_Structure = "Normal"))



library(ggplot2)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(scales)

plot_DS_result <- function(n, p, type, arg, correlation_Structure, 
                           X_Structure, epsilon_Structure){
  if (type == "signal"){
    load(paste("results/results_figure3/rho", arg, "_n", n, 
               "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
               "_eps_", epsilon_Structure, ".RData", sep = ""))
    
    df_fdp <- data_save[,c("delta","DS_fdp","MDS_fdp","stab_DS_fdp")]
    colnames(df_fdp)=c("delta","DS","MDS","stab (DS)")
    
    df_power <- data_save[,c("delta","DS_power","MDS_power","stab_DS_power")]
    colnames(df_power)=c("delta","DS","MDS","stab (DS)")
    
    df_fdp <- melt(df_fdp,id.vars="delta",variable.name="method",value.name = "FDR")
    df_power <- melt(df_power,id.vars="delta",variable.name="method",value.name = "Power")
    
    plot_fdp <- ggplot(df_fdp, aes(delta,FDR))+
      geom_point(aes(color = method, shape = method), cex = 2) +
      geom_line(aes(color = method)) + theme_bw() +
      geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
    
    plot_power <- ggplot(df_power, aes(delta,Power)) +
      geom_point(aes(color = method, shape = method), cex=2) +
      geom_line(aes(color = method)) + theme_bw()+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
  }
  if (type == "corr"){
    load(paste("results/results_figure3/delta", arg, "_n", n, 
               "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
               "_eps_", epsilon_Structure, ".RData", sep = ""))
    
    df_fdp <- data_save[,c("rho","DS_fdp","MDS_fdp","stab_DS_fdp")]
    colnames(df_fdp)=c("rho","DS","MDS","stab (DS)")
    
    df_power <- data_save[,c("rho","DS_power","MDS_power","stab_DS_power")]
    colnames(df_power)=c("rho","DS","MDS","stab (DS)")
    
    df_fdp <- melt(df_fdp,id.vars="rho",variable.name="method",value.name = "FDR")
    df_power <- melt(df_power,id.vars="rho",variable.name="method",value.name = "Power")
    
    plot_fdp <- ggplot(df_fdp, aes(rho,FDR))+
      geom_point(aes(color = method, shape = method), cex = 2) +
      geom_line(aes(color = method)) + theme_bw() +
      geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
    
    plot_power <- ggplot(df_power, aes(rho,Power)) +
      geom_point(aes(color = method, shape = method), cex=2) +
      geom_line(aes(color = method)) + theme_bw()+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
  }
  
  g <- ggarrange(plot_fdp,plot_power,common.legend = T,legend = "top")
  if (type == "corr"){
    ggsave(file=paste("figures/figure3/DS_delta",arg,"_n", n, 
                      "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                      "_eps_", epsilon_Structure, ".pdf", sep = ""), g, width = 7, height = 3)
  }
  if (type == "signal"){
    ggsave(file=paste("figures/figure3/DS_rho",arg,"_n", n, 
                      "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                      "_eps_", epsilon_Structure, ".pdf", sep = ""), g, width = 7, height = 3)
  }
  g
}


plot_kn_result <- function(n, p, type, arg, correlation_Structure, 
                           X_Structure, epsilon_Structure){
  
  if (type == "signal"){
    load(paste("results/results_figure3/rho", arg, "_n", n, 
               "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
               "_eps_", epsilon_Structure, ".RData", sep = ""))
    
    df_fdp <- data_save[,c("delta","knockoff_fdp","derand_kn_fdp","stab_kn_fdp")]
    colnames(df_fdp)=c("delta","knockoff","derand_kn","stab (kn)")
    
    df_power <- data_save[,c("delta","knockoff_power","derand_kn_power","stab_kn_power")]
    colnames(df_power)=c("delta","knockoff","derand_kn","stab (kn)")
    
    df_fdp <- melt(df_fdp,id.vars="delta",variable.name="method",value.name = "FDR")
    df_power <- melt(df_power,id.vars="delta",variable.name="method",value.name = "Power")
    
    plot_fdp <- ggplot(df_fdp, aes(delta,FDR))+
      geom_point(aes(color = method, shape = method), cex = 2) +
      geom_line(aes(color = method)) + theme_bw() +
      geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
    
    plot_power <- ggplot(df_power, aes(delta,Power)) +
      geom_point(aes(color = method, shape = method), cex=2) +
      geom_line(aes(color = method)) + theme_bw()+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
  }
  if (type == "corr"){
    load(paste("results/results_figure3/delta", arg, "_n", n, 
               "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
               "_eps_", epsilon_Structure, ".RData", sep = ""))
    
    df_fdp <- data_save[,c("rho","knockoff_fdp","derand_kn_fdp","stab_kn_fdp")]
    colnames(df_fdp)=c("rho","knockoff","derand_kn","stab (kn)")
    
    df_power <- data_save[,c("rho","knockoff_power","derand_kn_power","stab_kn_power")]
    colnames(df_power)=c("rho","knockoff","derand_kn","stab (kn)")
    
    df_fdp <- melt(df_fdp,id.vars="rho",variable.name="method",value.name = "FDR")
    df_power <- melt(df_power,id.vars="rho",variable.name="method",value.name = "Power")
    
    plot_fdp <- ggplot(df_fdp, aes(rho,FDR))+
      geom_point(aes(color = method, shape = method), cex = 2) +
      geom_line(aes(color = method)) + theme_bw() +
      geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
    
    plot_power <- ggplot(df_power, aes(rho,Power)) +
      geom_point(aes(color = method, shape = method), cex=2) +
      geom_line(aes(color = method)) + theme_bw()+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
  }
  
  g <- ggarrange(plot_fdp,plot_power,common.legend = T,legend = "top")
  if (type == "corr"){
    ggsave(file=paste("figures/figure3/knockoff_delta",arg,"_n", n, 
                      "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                      "_eps_", epsilon_Structure, ".pdf", sep = ""), g, width = 7, height = 3)
  }
  if (type == "signal"){
    ggsave(file=paste("figures/figure3/knockoff_rho",arg,"_n", n, 
                      "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                      "_eps_", epsilon_Structure, ".pdf", sep = ""), g, width = 7, height = 3)
  }
  g
}

plot_bh_result <- function(n, p, type, arg, correlation_Structure, 
                           X_Structure, epsilon_Structure){
  if (type == "signal"){
    load(paste("results/results_figure3/rho", arg, "_n", n, 
               "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
               "_eps_", epsilon_Structure, ".RData", sep = ""))
    
    df_fdp <- data_save[,c("delta", "BH_fdp", "MBH_fdp", "BH_stab_fdp")]
    colnames(df_fdp)=c("delta", "BH", "MBH", "stab (BH)")
    
    df_power <- data_save[,c("delta", "BH_power", "MBH_power", "BH_stab_power")]
    colnames(df_power)=c("delta", "BH", "MBH", "stab (BH)")
    
    df_fdp <- melt(df_fdp,id.vars="delta",variable.name="method",value.name = "FDR")
    df_power <- melt(df_power,id.vars="delta",variable.name="method",value.name = "Power")
    
    plot_fdp <- ggplot(df_fdp, aes(delta,FDR))+
      geom_point(aes(color = method, shape = method), cex = 2) +
      #scale_shape_manual(values = c(17,15))+
      #scale_color_manual(values = c("#00BA38","#619CFF"))+
      geom_line(aes(color = method)) + theme_bw() +
      geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
    
    plot_power <- ggplot(df_power, aes(delta,Power)) +
      geom_point(aes(color = method, shape = method), cex=2) +
      #scale_shape_manual(values = c(17,15))+
      #scale_color_manual(values = c("#00BA38","#619CFF"))+
      geom_line(aes(color = method)) + theme_bw()+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
  }
  if (type == "corr"){
    load(paste("results/results_figure3/delta", arg, "_n", n, 
               "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
               "_eps_", epsilon_Structure, ".RData", sep = ""))
    
    df_fdp <- data_save[,c("rho", "BH_fdp", "MBH_fdp", "BH_stab_fdp")]
    colnames(df_fdp)=c("rho", "BH", "MBH", "stab (BH)")
    
    df_power <- data_save[,c("rho", "BH_power", "MBH_power", "BH_stab_power")]
    colnames(df_power)=c("rho", "BH", "MBH", "stab (BH)")
    
    df_fdp <- melt(df_fdp,id.vars="rho",variable.name="method",value.name = "FDR")
    df_power <- melt(df_power,id.vars="rho",variable.name="method",value.name = "Power")
    
    plot_fdp <- ggplot(df_fdp, aes(rho,FDR))+
      geom_point(aes(color = method, shape = method), cex = 2) +
      #scale_shape_manual(values = c(17,15))+
      #scale_color_manual(values = c("#00BA38","#619CFF"))+
      geom_line(aes(color = method)) + theme_bw() +
      geom_hline(yintercept=0.1, linetype='dotted', col = 'red')+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
    
    plot_power <- ggplot(df_power, aes(rho,Power)) +
      geom_point(aes(color = method, shape = method), cex=2) +
      #scale_shape_manual(values = c(17,15))+
      #scale_color_manual(values = c("#00BA38","#619CFF"))+
      geom_line(aes(color = method)) + theme_bw()+ 
      scale_y_continuous(labels = label_number(accuracy = 0.01))
  }
  
  g <- ggarrange(plot_fdp,plot_power,common.legend = T,legend = "top")
  if (type == "corr"){
    ggsave(file=paste("figures/figure3/MBHq_delta",arg,"_n", n, 
                      "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                      "_eps_", epsilon_Structure, ".pdf", sep = ""), g, width = 7, height = 3)
  }
  if (type == "signal"){
    ggsave(file=paste("figures/figure3/MBHq_rho",arg,"_n", n, 
                      "_p", p, "_", type, "_X_", X_Structure, "_Xcor_", correlation_Structure, 
                      "_eps_", epsilon_Structure, ".pdf", sep = ""), g, width = 7, height = 3)
  }
  
  g
}


n=500
p=500

plot_DS_result(n, p, type = "corr", arg = 5, 
               correlation_Structure = "toeplitz",X_Structure = "Normal",
               epsilon_Structure = "Normal")
plot_kn_result(n, p, type = "corr", arg = 5, 
               correlation_Structure = "toeplitz",X_Structure = "Normal",
               epsilon_Structure = "Normal")
plot_bh_result(n, p, type = "corr", arg = 5, 
               correlation_Structure = "toeplitz",X_Structure = "Normal",
               epsilon_Structure = "Normal")
n=3000
p=500

plot_DS_result(n, p, type = "corr", arg = 5, 
               correlation_Structure = "toeplitz",X_Structure = "Normal",
               epsilon_Structure = "Normal")
plot_kn_result(n, p, type = "corr", arg = 5, 
               correlation_Structure = "toeplitz",X_Structure = "Normal",
               epsilon_Structure = "Normal")
plot_bh_result(n, p, type = "corr", arg = 5, 
               correlation_Structure = "toeplitz",X_Structure = "Normal",
               epsilon_Structure = "Normal")

# Results in Supplementary Material
# n <- 800
# p <- 2000
# p0 <- 50
# q <- 0.1
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5,type = "signal",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# n <- 800
# p <- 1000
# p0 <- 50
# q <- 0.1
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# 
# n <- 500
# p <- 500
# p0 <- 50
# q <- 0.1
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# n <- 3000
# p <- 500
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# 
# n <- 2000
# p <- 800
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "toeplitz",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 0.5, type = "signal",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
# 
# system.time(simulation_par(rep_num, M, n, p, p0, q, arg = 5, type = "corr",
#                            correlation_Structure = "Compound_Symmetry",X_Structure = "Normal",
#                            epsilon_Structure = "Normal"))
sfStop()


