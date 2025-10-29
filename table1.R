# renv::install("SeuratObject")
# 
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")
# 
# remotes::install_github("satijalab/seurat", ref = "develop")
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github("cran/EDMeasure")

rm(list = ls())
library(Seurat)
library(SeuratDisk)
library(here)
library(xgboost)
library(snowfall)
library(rlist)
setwd(here::here())
source("mdc_functions_for_table1.R")




################################################################################################################################################################
# Data preprocessing steps
# please download the files yourself using the provided addresses: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
# this file is large and have not been included in the folder here. 
# Therefore, we have saved the final data obtained from the preprocessing pipeline to the folder "real data/single cell data".
################################################################################################################################################################

# pbmc <- LoadH5Seurat("real data/single cell data/pbmc_multimodal.h5seurat")
# 
# Idents(pbmc) <- pbmc$celltype.l1
# markers_Mono <- FindMarkers(pbmc, ident.1 = "Mono",min.pct = 0.5,logfc.threshold = log(1.5))
# markers_CD4_T <- FindMarkers(pbmc, ident.1 = "CD4 T",min.pct = 0.5,logfc.threshold = log(1.5))
# markers_CD8_T <- FindMarkers(pbmc, ident.1 = "CD8 T",min.pct = 0.5,logfc.threshold = log(1.5))
# markers_NK <- FindMarkers(pbmc, ident.1 = "NK",min.pct = 0.5,logfc.threshold = log(1.5))
# markers_B <- FindMarkers(pbmc, ident.1 = "B",min.pct = 0.5,logfc.threshold = log(1.5))
# markers_other_T <- FindMarkers(pbmc, ident.1 = "other T",min.pct = 0.5,logfc.threshold = log(1.5))
# markers_DC <- FindMarkers(pbmc, ident.1 = "DC",min.pct = 0.5,logfc.threshold = log(1.5))
# markers_other <- FindMarkers(pbmc, ident.1 = "other",min.pct = 0.5,logfc.threshold = log(1.5))
# 
# 
# 
# celltype<- "DC"
# markers <- markers_DC
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_DC <- list(X=X,Y=Y,X_marker)
# save(data_DC,file="real data/single cell data/DC.RData")
# 
# celltype<- "Mono"
# markers <- markers_Mono
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_Mono <- list(X=X,Y=Y,X_marker)
# save(data_Mono,file="real data/single cell data/Mono.RData")
# 
# celltype<- "CD4_T"
# markers <- markers_CD4_T
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_CD4_T <- list(X=X,Y=Y,X_marker)
# save(data_CD4_T,file="real data/single cell data/CD4_T.RData")
# 
# celltype<- "CD8_T"
# markers <- markers_CD8_T
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_CD8_T <- list(X=X,Y=Y,X_marker)
# save(data_CD8_T,file="real data/single cell data/CD8_T.RData")
# 
# celltype<- "NK"
# markers <- markers_NK
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_NK <- list(X=X,Y=Y,X_marker)
# save(data_NK,file="real data/single cell data/NK.RData")
# 
# 
# celltype<- "B"
# markers <- markers_B
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_B <- list(X=X,Y=Y,X_marker)
# save(data_B,file="real data/single cell data/B.RData")
# 
# celltype<- "other"
# markers <- markers_other
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_other <- list(X=X,Y=Y,X_marker)
# save(data_other,file="real data/single cell data/other.RData")
# 
# celltype<- "other_T"
# markers <- markers_other_T
# pbmc_type <- subset(pbmc,cells = names(which(pbmc$celltype.l1==celltype)))
# X <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% pbmc_type@assays$SCT@var.features),])),'dgCMatrix')
# Y <- as(t(as.matrix(pbmc_type@assays$ADT@data)),'dgCMatrix')
# Y <- Y[,which(!(colnames(Y) %in% c("CD26-1","TSLPR")))]
# X_marker <- as(t(as.matrix(pbmc_type@assays$SCT@data[which(rownames(pbmc_type@assays$SCT@meta.features) %in% rownames(markers)),])),'dgCMatrix')
# data_other_T <- list(X=X,Y=Y,X_marker)
# save(data_other_T,file="real data/single cell data/other_T.RData")
# 


Fun <- function(protein, X, Y, X_marker, samp, B){
  set.seed(protein)
  result_all <- SVG_test_combine(X, matrix(0,ncol=2,nrow=length(Y[,protein])), Y[,protein], testMethod = "rank_sum", samp, B = B)
  result_marker <- SVG_test_combine(X_marker, matrix(0,ncol=2,nrow=length(Y[,protein])), Y[,protein], testMethod = "rank_sum", samp, B = B)
  return(list(result_all=result_all,result_marker=result_marker))
}

run_simu <- function(X, Y, X_marker, split_times, B){
  
  samp <- create_sample(split_times,nrow(X))
  
  protein <- 1:ncol(Y)
  res <- sfLapply(protein, fun = Fun, X, Y, X_marker, samp, B)
  result <- get_result(res, q = 0.05, split_times=split_times)
  return(list(detail = res, result = result))
}



sfInit(parallel = TRUE, cpus = 50)

sfLibrary(snowfall)
sfLibrary(xgboost)
sfLibrary(rlist)
sfSource("mdc_functions_for_table1.R")



# The code takes a long time to run due to the large dataset, and even repeating the process 20 times is time-consuming.
# rep_num = 20 in paper, 5 for fast
rep_num = 5
# split_times = 50 in paper, 10 for fast
split_times=10
############################################################
load("real data/single cell data/Mono.RData")
X <- data_Mono$X
X_marker <- data_Mono[[3]]
Y <- data_Mono$Y

Mono_res <- list()
system.time(
  for (i in 1:rep_num) {
    Mono_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
  }
)


save(Mono_res,file="results/results_table1/Mono_res.RData")
############################################################
load("real data/single cell data/other_T.RData")
X <- data_other_T$X
X_marker <- data_other_T[[3]]
Y <- data_other_T$Y

other_T_res <- list()
for (i in 1:rep_num) {
  other_T_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
}

save(other_T_res,file="results/results_table1/other_T_res.RData")
############################################################
load("real data/single cell data/other.RData")
X <- data_other$X
X_marker <- data_other[[3]]
Y <- data_other$Y

other_res <- list()
for (i in 1:rep_num) {
  other_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
}

save(other_res,file="results/results_table1/other_res.RData")
############################################################
load("real data/single cell data/NK.RData")
X <- data_NK$X
X_marker <- data_NK[[3]]
Y <- data_NK$Y

NK_res <- list()
for (i in 1:rep_num) {
  NK_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
}

save(NK_res,file="results/results_table1/NK_res.RData")
############################################################
load("real data/single cell data/DC.RData")
X <- data_DC$X
X_marker <- data_DC[[3]]
Y <- data_DC$Y

DC_res <- list()
for (i in 1:rep_num) {
  DC_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
}

save(DC_res,file="results/results_table1/DC_res.RData")
############################################################
load("real data/single cell data/B.RData")
X <- data_B$X
X_marker <- data_B[[3]]
Y <- data_B$Y

B_res <- list()
for (i in 1:rep_num) {
  B_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
}

save(B_res,file="results/results_table1/B_res.RData")
############################################################
load("real data/single cell data/CD4_T.RData")
X <- data_CD4_T$X
X_marker <- data_CD4_T[[3]]
Y <- data_CD4_T$Y

CD4_T_res <- list()
for (i in 1:rep_num) {
  CD4_T_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
}

save(CD4_T_res,file="results/results_table1/CD4_T_res.RData")
############################################################
load("real data/single cell data/CD8_T.RData")
X <- data_CD8_T$X
X_marker <- data_CD8_T[[3]]
Y <- data_CD8_T$Y

CD8_T_res <- list()
for (i in 1:rep_num) {
  CD8_T_res[[i]] <- run_simu(X, Y, X_marker,split_times=split_times, B=B) 
}

save(CD8_T_res,file="results/results_table1/CD8_T_res.RData")
############################################################
sfStop()


results <- list(Mono_result=Mono_res,B_result=B_res,CD4_T_result=CD4_T_res,CD8_T_result=CD8_T_res,
                NK_result=NK_res,DC_result=DC_res,other_res=other_res,other_T_res=other_T_res)
save(results,file="results/results_table1/results.RData")
############################################################

# write table
load("results/results_table1/results.RData")
q <- 0.1
# split_times = 50 ,rep_num = 20 in paper
split_times=50
rep_num = 20

# q = 0.05; 0.2; 0.3  in Supplementary Material Table S1
var_result <- matrix(0,nrow = 8,ncol = 4)
mean_result <- matrix(0,nrow = 8,ncol = 4)
for (c_type in 1:8) {
  stab_all <- rep(0,rep_num)
  stab_marker <- rep(0,rep_num)
  CCT_all <- rep(0,rep_num)
  CCT_marker <- rep(0,rep_num)
  for (i in 1:rep_num) {
    stab_all[i] <- length(get_result(results[[c_type]][[i]]$detail, q = q, split_times=split_times)$res_stab_all)
    stab_marker[i] <- length(get_result(results[[c_type]][[i]]$detail, q = q, split_times=split_times)$res_stab_marker)
    CCT_all[i] <- length(get_result(results[[c_type]][[i]]$detail, q = q, split_times=split_times)$res_CCT_all)
    CCT_marker[i] <- length(get_result(results[[c_type]][[i]]$detail, q = q, split_times=split_times)$res_CCT_marker)
  }
  var_result[c_type,1] <- var(CCT_all)
  var_result[c_type,2] <- var(CCT_marker)
  var_result[c_type,3] <- var(stab_all)
  var_result[c_type,4] <- var(stab_marker)
  mean_result[c_type,1] <- mean(CCT_all)
  mean_result[c_type,2] <- mean(CCT_marker)
  mean_result[c_type,3] <- mean(stab_all)
  mean_result[c_type,4] <- mean(stab_marker)
}
colnames(mean_result) <- c("CCT_all","CCT_marker","stab_all","stab_marker")
colnames(var_result) <- c("CCT_all","CCT_marker","stab_all","stab_marker")
rownames(mean_result) <- c("Mono","B","CD4_T","CD8_T","NK","DC","other","other_T")
rownames(var_result) <- c("Mono","B","CD4_T","CD8_T","NK","DC","other","other_T")

#write.csv(var_result,file = "figures/table1/tableS1_q0.3.csv")
write.csv(var_result,file = "figures/table1/table1.csv")




########################################################################################################################
library(Thermimage)
library(ggplot2)
library(ggpubr)
load("results/results_table1/results.RData")

result_CCT_all <- data.frame(matrix(0,nrow=226,ncol=8))
colnames(result_CCT_all) <- c("Mono_CCT_all","B_CCT_all","CD4_T_CCT_all","CD8_T_CCT_all","NK_CCT_all","DC_CCT_all","other_CCT_all","other_T_CCT_all")
result_CCT_marker <- data.frame(matrix(0,nrow=226,ncol=8))
colnames(result_CCT_marker) <- c("Mono_CCT_marker","B_CCT_marker","CD4_T_CCT_marker","CD8_T_CCT_marker","NK_CCT_marker","DC_CCT_marker","other_CCT_marker","other_T_CCT_marker")
result_stab_all <- data.frame(matrix(0,nrow=226,ncol=8))
colnames(result_stab_all) <- c("Mono_stab_all","B_stab_all","CD4_T_stab_all","CD8_T_stab_all","NK_stab_all","DC_stab_all","other_stab_all","other_T_stab_all")
result_stab_marker <- data.frame(matrix(0,nrow=226,ncol=8))
colnames(result_stab_marker) <- c("Mono_stab_marker","B_stab_marker","CD4_T_stab_marker","CD8_T_stab_marker","NK_stab_marker","DC_stab_marker","other_stab_marker","other_T_stab_marker")

for (c_type in 1:8) {
  result_CCT_all[results[[c_type]][[1]]$result$res_CCT_all,c_type] <- 1
  result_CCT_marker[results[[c_type]][[1]]$result$res_CCT_marker,c_type] <- 1
  result_stab_all[results[[c_type]][[1]]$result$res_stab_all,c_type] <- 1
  result_stab_marker[results[[c_type]][[1]]$result$res_stab_marker,c_type] <- 1
  
}
CCT_table <- cbind(result_CCT_all,result_CCT_marker)
stab_table <- cbind(result_stab_all,result_stab_marker)
table_bind <- cbind(CCT_table,stab_table)

table_bind <- table_bind[,c(1,9,17,25,(c(1,9,17,25)+1),(c(1,9,17,25)+2),(c(1,9,17,25)+3),(c(1,9,17,25)+4),(c(1,9,17,25)+5),(c(1,9,17,25)+6),(c(1,9,17,25)+7))]

write.csv(table_bind,file = "figures/table1/result_table.csv")







# Figure S44 and S45 in Supplementary Material

#########################################################################
type <- c("Mono","B","CD4 T","CD8 T","NK","DC","Other","Other T")

pdf("figures/figures_Sec_S4.2/figureS45.pdf", 10, 8)
par(mfrow=c(8,1),mar=c(0.5,6,0.5,2),font=2)
for (i in 1:8) {
  temp <- cbind(stab_table[,i],-stab_table[,i+8])
  image(as.matrix(temp), col = c("orange","yellow", "blue"),
        xaxt='n', yaxt='n')
  ytick<-seq(0, 1, length=1)
  text(par("usr")[1], ytick, 
       labels = type[i],  pos = 2, xpd = TRUE, cex=1.8)
}
dev.off()

pdf("figures/figures_Sec_S4.2/figureS44.pdf", 10, 8)
par(mfrow=c(8,1),mar=c(0.5,6,0.5,2),font=2)
for (i in 1:8) {
  temp <- cbind(CCT_table[,i],-CCT_table[,i+8])
  image(as.matrix(temp), col = c("orange","yellow", "blue"),
        xaxt='n', yaxt='n')
  ytick<-seq(0, 1, length=1)
  text(par("usr")[1], ytick, 
       labels = type[i], pos = 2, xpd = TRUE, cex=1.8)
}
dev.off()

#########################################################################