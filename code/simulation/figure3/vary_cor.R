rm(list = ls())
library(MASS)
library(mvnfast)
library(knockoff)
library(glmhd)

### Set working directory to 'figure3' folder, e.g., 
# setwd("~/code/simulation/figure3")

### source code
source('utlis.R')
source('DS.R')
source('MDS.R')
source('BHq.R')
source('GM.R')
source('ABHq.R')

### replicate index
#replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#set.seed(replicate)

f <- function(x){
  exp(x)/(1 + exp(x))
}

p = 500
n = 3000
s = 50
q = 0.1
### Change the correlation
#rho = as.numeric(Sys.getenv("att"))
signal_strength = 11
Sigma = matrix(0, nrow = p, ncol = p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)*1/sqrt(n)
signal_index = sample(1:p, size = s, replace = F)
beta = numeric(p)
beta[signal_index] = sample(c(-1, 1), s, replace = T)*signal_strength
mu = X %*% beta
y = rbinom(n, size = 1, prob = f(mu)) 

DS_time = system.time(DS_result <- DS(X, y))[3]
DS_result <- fdp_power(DS_result$select_index)
DS_fdp <- DS_result$fdp
DS_power <- DS_result$power

MDS_time = system.time(MDS_result <- MDS(X, y, num_split = 50))[3]
MDS_result <- fdp_power(MDS_result$select_index)
MDS_fdp <- MDS_result$fdp
MDS_power <- MDS_result$power

BHq_time = system.time(BHq_result <- BHq(X, y))[3]
BHq_result <- fdp_power(BHq_result$select_index)
BHq_fdp <- BHq_result$fdp
BHq_power <- BHq_result$power

GM_time = system.time(GM_result <- GM(X, y))[3]
GM_result <- fdp_power(GM_result$select_index)
GM_fdp <- GM_result$fdp
GM_power <- GM_result$power

ABHq_time = system.time(ABHq_result <- ABHq(X, y))[3]
ABHq_result <- fdp_power(ABHq_result$select_index)
ABHq_fdp <- ABHq_result$fdp
ABHq_power <- ABHq_result$power

### save data
data_save <- list(DS_fdp   = DS_fdp,   DS_power = DS_power, DS_time = DS_time,
                  MDS_fdp  = MDS_fdp,  MDS_power = MDS_power, MDS_time = MDS_time,
                  BHq_fdp  = BHq_fdp,  BHq_power = BHq_power, BHq_time = BHq_time,
                  GM_fdp   = GM_fdp,   GM_power = GM_power, GM_time = GM_time,
                  ABHq_fdp = ABHq_fdp, ABHq_power = ABHq_power, ABHq_time = ABHq_time)















