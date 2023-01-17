rm(list = ls())

### source code
#setwd('~/realdata')
source('knockoff.R')
library(glmnet)
library(knockoff)

replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(replicate)

data <- t(read.table("~/data.txt"))
y <- c(rep(0, 400), rep(1, 2000))

### dimension reduction
freq <- colSums(data != 0)/nrow(data)
remove_index <- which(freq < 0.1)
X <- data[, -remove_index]
std <- apply(X, 2, sd)
gene_index <- order(std, decreasing = T)
p <- 500
n <- nrow(X)
X <- X[, gene_index[1:p]]

### FDR control level 
q <- 0.1

k = as.numeric(Sys.getenv("att"))

analys_pwer = function(w, k){
  sorted_w = sort(w)
  cutoff = -sorted_w[k]
  selected_index = which(w>cutoff)
  selected_index
}

B = 50
selected_matrix = matrix(0, nrow = B, ncol = p)

for(i in 1:B){
  w = KN(X, y)
  index = analys_pwer(w, k)
  selected_matrix[i, index] = selected_matrix[i, index]+1
}

freq = colMeans(selected_matrix)
set = which(freq>=0.5)

data_save = set

save(data_save, file = paste("~/deran_KN/", "k_", k, ".RData", sep = ""))




