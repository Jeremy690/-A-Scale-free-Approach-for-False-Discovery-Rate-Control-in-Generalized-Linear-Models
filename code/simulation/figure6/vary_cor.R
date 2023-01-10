rm(list = ls())
library(glmnet)
library(MASS)
library(knockoff)
library(mvnfast)
source('utlis.R')
source('knockoff.R')
source('DS.R')

replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(replicate)

n = 800
p = 1600
q = 0.1

### covariance matrix
rho = as.numeric(Sys.getenv("att"))
sig1 = toeplitz(seq(rho,rho,length.out = p/8))
Sigma = bdiag(rep(list(sig1),8))
diag(Sigma) = 1

s = 50
signal_index = sample(1:p, s)

beta = numeric(p)
beta[signal_index] = sample(c(-1, 1), s, replace = T)*2
X = rmvn(n, mu = rep(0, p), sigma = 0.01*Sigma)
while(length(which(abs(X%*% beta)>4))>0){
  X[which(abs(X%*%beta)>4),]=mvrnorm(length(which(abs(X%*% beta)>4)), mu=rep(0,p), Sigma=0.01*Sigma)
}
mu = X %*% beta
mu = exp(mu)
y = rpois(n, lambda = mu)

result1 = MDS(X, y, 50)
result2 = KN(X, y)
result3 = MKN(X, y)



data_save = list(DS.fdp = result1$DS.fdp, DS.power = result1$DS.power,
                 MDS.fdp = result1$MDS.fdp, MDS.power = result1$MDS.power,
                 KN.fdp = result2$fdp, KN.power = result2$power,
                 MKN.fdp = result3$fdp, MKN.power = result3$power)

save(data_save, file = paste("~/result_left/rho_", rho, "_replicate_", replicate, ".RData", sep = ""))



