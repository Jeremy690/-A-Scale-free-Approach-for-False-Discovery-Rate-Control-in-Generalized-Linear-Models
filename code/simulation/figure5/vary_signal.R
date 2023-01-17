rm(list = ls())
library(glmnet)
library(MASS)
library(knockoff)
library(mvnfast)

### Set working directory to 'figure5' folder, e.g., 
# setwd("~/code/simulation/figure5")

### source code
source('utils.R')
source('DS_MDS.R')
source('BHq.R')
source('Knockoff.R')

# algorithmic setting
n = 800
p = 2000

### replicate index
#replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#set.seed(replicate)

sig1 = toeplitz(seq(0.3,0,length.out = p/10))
Sigma = bdiag(rep(list(sig1),10))
diag(Sigma) = 1

s = 60
### Vary the signal strength
#signal = as.numeric(Sys.getenv("att"))

q = 0.1
beta_q = sample(c(-1, 1), s, replace = T)*signal
signal_index = sample(1:p, size = s, replace = F)
beta = numeric(p)
beta[signal_index]=beta_q
X = rmvn(n, mu = rep(0,p), sigma = 0.1*Sigma)
mu = X %*% beta
y = rbinom(n, size = 1, prob = f(mu))

result1 = DS(X, y)
result2 = BHq(X, y)
result3 = KN(X, y)
result4 = MDS(X, y, 100)

data_save = list(DS.fdp = result1$fdp, DS.power = result1$power,
                 BHq.fdp = result2$fdp, BHq.power = result2$power,
                 KN.fdp = result3$fdp, KN.power = result3$power,
                 MDS.fdp = result4$fdp, MDS.power = result4$power)
