rm(list = ls())
library(snowfall)
library(here)
library(knockoff)
library(glmnet)
library(hdi)
library(corpcor)
library(Rcpp)
library(RcppArmadillo)
library(rlist)
### read files
setwd(here::here())
folder_path <- "utils"
file_list <- list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
for (file_path in file_list) {
  source(file_path)
}

sourceCpp("utils/solve_linear_equation.cpp")
sourceCpp("utils/cholcpp.cpp")



drug_class = 'PI' # Possible drug types are 'PI', 'NRTI', and 'NNRTI'. In the paper we restrict to PI and NRTI.

### Cleaning and Fetching the data

gene_data = "real data/HIV_data/PI_DATA.txt"
tsm_data = "real data/HIV_data/PI.txt"


gene_df = read.delim(gene_data, na.string = c('NA', ''), stringsAsFactors = FALSE)
tsm_df = read.delim(tsm_data, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) = c('Position', 'Mutations')

### Removing the rows with error flags or nonstandard mutation codes.
grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}

pos_start = which(names(gene_df) == 'P1')
pos_cols = seq.int(pos_start, ncol(gene_df))
valid_rows = grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[,pos_cols])
gene_df = gene_df[valid_rows,]

### Prepare the design matrix and response variable.
# Flatten a matrix to a vector with names from concatenating row/column names.
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}

# Construct preliminary design matrix.
muts = c(LETTERS, 'i', 'd')
X = outer(muts, as.matrix(gene_df[,pos_cols]), Vectorize(grepl))
X = aperm(X, c(2,3,1))
dimnames(X)[[3]] <- muts
X = t(apply(X, 1, flatten_matrix))
mode(X) <- 'numeric'

# Remove any mutation/position pairs that never appear in the data.
X = X[,colSums(X) != 0]

# Extract response matrix.
Y = gene_df[,4:(pos_start-1)]

analys = function(mm, ww, q = 0.1){
  t_set = max(ww)
  for(t in ww){
    ps = length(mm[mm>=t])
    ng = length(na.omit(mm[mm<=-t]))
    rto = (ng+1)/max(ps, 1)
    if(rto<=q){
      t_set = c(t_set, t)
      #      print(t)
    }
  }
  thre = min(t_set)
  nz_est = which(mm>=thre)
  nz_est
}
Split = function(X, y, q){
  n = nrow(X); p = ncol(X)
  num_split <- 50
  selected_index_multiple <- matrix(0, nrow = num_split, ncol = p)
  fdr_multiple <- rep(0, num_split)
  power_multiple <- rep(0, num_split)
  num_select <- rep(0, num_split)
  MM <- matrix(rep(0,num_split*p),nrow=num_split,ncol=p)
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  for(iter in 1:num_split){
    ### randomly split the data and estimate coefficient
    while(TRUE){
      sample_index1 <- sample(x = c(1:n), size = trunc(0.5 * n), replace = F)
      sample_index2 <- setdiff(c(1:n), sample_index1)
      cvfit = cv.glmnet(X[sample_index1, ], y[sample_index1], nfolds = 10)
      beta1 = as.vector(coef(cvfit$glmnet.fit, cvfit$lambda.1se))[-1]
      nonzero_index <- which(beta1 != 0)
      if(length(nonzero_index)!=0)
        break
    }
    beta2 <- rep(0, p)
    beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
    ### calculate the test statistics
    M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
    M[is.na(M)] = 0
    MM[iter,] <- M
    ### find the threshold
    result = analys(M, abs(M), q)
    current_selected_index <- result
    ### number of selected variables
    num_select[iter] <- length(current_selected_index)
    
    if(num_select[iter]!=0)
      selected_index_multiple[iter, current_selected_index] <- selected_index_multiple[iter, current_selected_index] + 1/num_select[iter]
    inclusion_rate[iter, current_selected_index] <- 1/num_select[iter]
  }
  
  ### single splitting result
  single_split_selected = current_selected_index
  
  ### multiple splitting result
  feature_rank <- order(apply(selected_index_multiple, 2, mean))
  feature_rank <- setdiff(feature_rank, which(apply(selected_index_multiple, 2, mean) == 0))
  null_variable <- numeric()
  fdr_replicate <- rep(0, num_split)
  for(feature_index in 1:length(feature_rank)){
    for(split_index in 1:num_split){
      if(selected_index_multiple[split_index, feature_rank[feature_index]]){
        fdr_replicate[split_index] <- fdr_replicate[split_index] + 1/num_select[split_index]
      }
    }
    if(mean(fdr_replicate) > q){
      break
    }else{
      null_variable <- c(null_variable, feature_rank[feature_index])
    }
  }
  ### conservative one
  multiple_selected_index <- setdiff(feature_rank, null_variable)
  
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
    selected_index1 <- which(inclusion_rate>=sort(inclusion_rate,decreasing = T)[ceiling(mean(num_select))])
  }
  list(SDS = single_split_selected, MDS = multiple_selected_index, DS_stab = selected_index1)
}


DS_knockoff_and_bhq <- function (M, X, y, q) {
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 3]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  # Run the knockoff filter.
  # knock.gen = function(x) create.second_order(x,  method='equi')
  # result = knockoff.filter(X, y, fdr=fdr, knockoffs=knock.gen, statistic=stat.glmnet_lambdasmax)
  # knockoff_selected = names(result$selected)
  result = knockoff_res(X, y, fdr,  M)
  knockoff_selected = colnames(X)[result$kn_sel]
  knockoff_stab_selected = colnames(X)[result$stab_sel]
  
  result = derand_kn_res(X, y, M, fdr)
  derand_kn_sel = colnames(X)[result$derand_kn]
  
  result = MBHq_res(X, y, fdr, M)
  BH_stab_selected = colnames(X)[result$stab_sel]
  MBHq_selected = colnames(X)[result$MBHq_sel]
  bhq_selected = colnames(X)[result$BHq_sel]
  
  result = Split(X, y, fdr)
  SDS_selected = colnames(X)[result$SDS]
  MDS_selected = colnames(X)[result$MDS]
  DS_stab <- colnames(X)[result$DS_stab]
  list(SDS = SDS_selected, MDS = MDS_selected, DS_stab = DS_stab, 
       knockoff = knockoff_selected, knockoff_stab = knockoff_stab_selected, derand_kn = derand_kn_sel, 
       BHq = bhq_selected, MBHq = MBHq_selected, BHq_stab = BH_stab_selected)
}




get_selected_name = function(X,y, selected){
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 3]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  colnames(X)[selected]
}




get_position <- function(x){
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)
}

fun <- function(num, M, X, Y, fdr){
  set.seed(num)
  ret <- lapply(Y, function(y) DS_knockoff_and_bhq(M, X, y, fdr))
  return(ret)
  
}

simulation_real_data <- function(rep_num, M, X, Y, fdr){
  num <- 1:rep_num
  data_save <- sfLapply(num, fun = fun, M, X, Y, fdr)
  
  save(data_save, file = "results/results_figure6/results_PI1.RData")
  return(data_save)
}

sfInit(parallel = TRUE, cpus = 50)

sfLibrary(snowfall)
sfLibrary(glmnet)
sfLibrary(knockoff)
sfLibrary(hdi)
sfLibrary(corpcor)
sfLibrary(doParallel)
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

sfExport("DS_knockoff_and_bhq")
sfExport("Split")
sfExport("analys")

fdr=0.1
sfExport("fdr")

# rep_num = 100 in paper, 50 for fast
rep_num = 50
# M=50 in paper, 20 for fast
M = 20
system.time(res <- simulation_real_data(rep_num, M, X, Y, fdr))

sfStop()





###plot


load("results/results_figure6/results_PI1.RData")
results <- data_save
comparisons <- lapply(results, function(r) {lapply(r, function(drug) {
  lapply(drug, function(selected) {
    positions = unique(get_position(selected)) # remove possible duplicates
    discoveries = length(positions)
    false_discoveries = length(setdiff(positions, tsm_df$Position))
    true_discoveries = discoveries - false_discoveries
    fdr = false_discoveries/discoveries
    power = true_discoveries/length(tsm_df$Position)
    list(true_discoveries = true_discoveries,
         false_discoveries = false_discoveries,
         fdr = fdr, power = power, discoveries = discoveries)
  })
})})


for (drug in names(comparisons[[1]])) {
  a <- paste("list.map(comparisons,",drug,")")
  SDS <- unlist(list.map(list.map(eval(parse(text = a)),SDS),fdr))
  SDS[is.na(SDS)] <- 0
  MDS <- unlist(list.map(list.map(eval(parse(text = a)),MDS),fdr))
  DS_stab <- unlist(list.map(list.map(eval(parse(text = a)),DS_stab),fdr))
  knockoff <- unlist(list.map(list.map(eval(parse(text = a)),knockoff),fdr))
  knockoff[is.na(knockoff)] <- 0
  knockoff_stab <- unlist(list.map(list.map(eval(parse(text = a)),knockoff_stab),fdr))
  derand_kn <- unlist(list.map(list.map(eval(parse(text = a)),derand_kn),fdr))
  derand_kn[is.na(derand_kn)] <- 0
  BHq <- unlist(list.map(list.map(eval(parse(text = a)),BHq),fdr))
  MBHq <- unlist(list.map(list.map(eval(parse(text = a)),MBHq),fdr))
  BHq_stab <- unlist(list.map(list.map(eval(parse(text = a)),BHq_stab),fdr))
  RES <- list(DS=SDS,MDS=MDS,DS_stab=DS_stab,knockoff=knockoff, knockoff_stab=knockoff_stab,derand_kn=derand_kn,
              BHq=BHq,MBHq=MBHq,BHq_stab=BHq_stab)
  
  RES <- melt(RES,id.vars="method",variable.name="method",value.name = "FDP")
  names(RES) <- c("FDP","Method")
  plt = ggplot(RES, aes(x=Method, y = FDP)) + 
    geom_boxplot()+    
    stat_boxplot(geom = "errorbar",width=0.25)+
    scale_x_discrete(limits=c("BHq", "MBHq", "BHq_stab", "knockoff", "derand_kn", "knockoff_stab", "DS", "MDS", "DS_stab"),
                     labels=c("BH", "MBH", "stab (BH)", "knockoff", "derand_kn", "stab (kn)", "DS", "MDS", "stab (DS)"))+
    theme_bw()+labs(title = drug)+theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept=0.1, linetype='dotted', col = 'red')
  
  ggsave(file=paste("figures/figure6/PI_", drug, 
                    "_fdp.pdf", sep = ""), plt, width = 10, height = 3)
}

for (drug in names(comparisons[[1]])) {
  a <- paste("list.map(comparisons,",drug,")")
  SDS <- unlist(list.map(list.map(eval(parse(text = a)),SDS),power))
  MDS <- unlist(list.map(list.map(eval(parse(text = a)),MDS),power))
  DS_stab <- unlist(list.map(list.map(eval(parse(text = a)),DS_stab),power))
  knockoff <- unlist(list.map(list.map(eval(parse(text = a)),knockoff),power))
  knockoff_stab <- unlist(list.map(list.map(eval(parse(text = a)),knockoff_stab),power))
  derand_kn <- unlist(list.map(list.map(eval(parse(text = a)),derand_kn),power))
  BHq <- unlist(list.map(list.map(eval(parse(text = a)),BHq),power))
  MBHq <- unlist(list.map(list.map(eval(parse(text = a)),MBHq),power))
  BHq_stab <- unlist(list.map(list.map(eval(parse(text = a)),BHq_stab),power))
  RES <- list(DS=SDS,MDS=MDS,DS_stab=DS_stab,knockoff=knockoff, knockoff_stab=knockoff_stab,derand_kn=derand_kn,
              BHq=BHq,MBHq=MBHq,BHq_stab=BHq_stab)
  #boxplot(RES,main=drug,ylab="number of power")
  RES <- melt(RES,id.vars="method",variable.name="method",value.name = "Power")
  names(RES) <- c("Power","Method")
  plt = ggplot(RES, aes(x=Method, y = Power)) + 
    geom_boxplot()+    
    stat_boxplot(geom = "errorbar",width=0.25)+
    scale_x_discrete(limits=c("BHq", "MBHq", "BHq_stab", "knockoff", "derand_kn", "knockoff_stab", "DS", "MDS", "DS_stab"),
                     labels=c("BH", "MBH", "stab (BH)", "knockoff", "derand_kn", "stab (kn)", "DS", "MDS", "stab (DS)"))+
    theme_bw()+labs(title = drug)+theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file=paste("figures/figure6/PI_", drug, 
                    "_power.pdf", sep = ""), plt, width = 10, height = 3)
}


for (drug in names(comparisons[[1]])) {
  a <- paste("list.map(comparisons,",drug,")")
  SDS <- unlist(list.map(list.map(eval(parse(text = a)),SDS),discoveries))
  MDS <- unlist(list.map(list.map(eval(parse(text = a)),MDS),discoveries))
  DS_stab <- unlist(list.map(list.map(eval(parse(text = a)),DS_stab),discoveries))
  knockoff <- unlist(list.map(list.map(eval(parse(text = a)),knockoff),discoveries))
  knockoff_stab <- unlist(list.map(list.map(eval(parse(text = a)),knockoff_stab),discoveries))
  derand_kn <- unlist(list.map(list.map(eval(parse(text = a)),derand_kn),discoveries))
  BHq <- unlist(list.map(list.map(eval(parse(text = a)),BHq),discoveries))
  MBHq <- unlist(list.map(list.map(eval(parse(text = a)),MBHq),discoveries))
  BHq_stab <- unlist(list.map(list.map(eval(parse(text = a)),BHq_stab),discoveries))
  RES <- list(DS=SDS,MDS=MDS,DS_stab=DS_stab,knockoff=knockoff, knockoff_stab=knockoff_stab,derand_kn=derand_kn,
              BHq=BHq,MBHq=MBHq,BHq_stab=BHq_stab)
  
  RES <- melt(RES,id.vars="method",variable.name="method",value.name = "discoveries")
  names(RES) <- c("discoveries","Method")
  plt = ggplot(RES, aes(x=Method, y = discoveries)) + 
    geom_boxplot()+    
    stat_boxplot(geom = "errorbar",width=0.25)+
    scale_x_discrete(limits=c("BHq", "MBHq", "BHq_stab", "knockoff", "derand_kn", "knockoff_stab", "DS", "MDS", "DS_stab"),
                     labels=c("BH", "MBH", "stab (BH)", "knockoff", "derand_kn", "stab (kn)", "DS", "MDS", "stab (DS)"))+
    theme_bw()+labs(title = drug)+theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file=paste("figures/figure6/PI_", drug, 
                    "_num_discoveries.pdf", sep = ""), plt, width = 10, height = 3)
}
