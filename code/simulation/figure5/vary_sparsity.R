rm(list = ls())
library(glmnet)
library(MASS)
library(knockoff)
library(mvnfast)

f = function(x){
  exp(x)/(1+exp(x))
}

fdp_power <- function(selected_index){
  fdp <- (length(selected_index) - length(intersect(selected_index, signal_index)))/max(length(selected_index), 1)
  power <- length(intersect(selected_index, signal_index))/length(signal_index)
  return(list(fdp = fdp, power = power))
}

analys = function(mm, ww, q){
  t_set = max(ww)
  for(t in ww){
    ps = length(mm[mm>=t])
    ng = length(na.omit(mm[mm<=-t]))
    rto = (ng+1)/max(ps, 1)
    if(rto<=q){
      t_set = c(t_set, t)
    }
  }
  thre = min(t_set)
  nz_est = which(mm>thre)
  nz_est
}

get_M = function(data){
  n - dim(data)[1]
  p = dim(data)[2]
  C = matrix(0, nrow = p, ncol = p)
  tau_square = numeric(p)
  for(j in 1:p){
    y = data[,j]
    X <- data[, -j]
    if(j == 1){
      cvfit <- cv.glmnet(X, y, nfolds = 10, nlambda = 200, intercept = F, standardize = F)
      lambda <- cvfit$lambda.min
    }
    beta1 <- as.vector(glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
    C[j, -j] = -beta1
    tau_square[j] = mean((y-X%*%beta1)*y)
  }
  diag(C) = 1
  T_mat = diag(tau_square)
  M = solve(T_mat)%*%C
  M
}

DS = function(x, y){
  n = dim(x)[1]; p = dim(x)[2]
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  X1 = x[sample_index1, ]
  y1 = y[sample_index1]
  X2 = x[sample_index2, ]
  y2 = y[sample_index2]
  
  cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "mse", nfolds = 10, nlambda = 200,  intercept = F, standardize = F)
  lambda <- cvfit$lambda.min/sqrt(2)
  beta1 <- as.vector(glmnet(X1, y1, family = "binomial", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
  beta2 <- as.vector(glmnet(X2, y2, family = "binomial", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
  
  w1 = as.vector(sqrt(f(X1%*%beta1)*(1-f(X1%*%beta1))))
  w2 = as.vector(sqrt(f(X2%*%beta2)*(1-f(X2%*%beta2))))
  X1_beta = diag(w1)%*%X1
  X2_beta = diag(w2)%*%X2
  M1 = get_M(X1_beta)
  M2 = get_M(X2_beta)
  beta1_d = as.numeric(beta1 + 2/n*M1%*%t(x[sample_index1, ])%*%(y[sample_index1] - f(x[sample_index1,] %*% beta1)))
  beta2_d = as.numeric(beta2 + 2/n*M2%*%t(x[sample_index2, ])%*%(y[sample_index2] - f(x[sample_index2,] %*% beta2)))
  sigma1 = sqrt(diag(M1%*%(t(X1)%*%diag(diag((y1-f(X1%*%beta1))%*%t(y1-f(X1%*%beta1))))%*%X1/(n/2))%*%t(M1)))
  sigma2 = sqrt(diag(M2%*%(t(X2)%*%diag(diag((y2-f(X2%*%beta2))%*%t(y2-f(X2%*%beta2))))%*%X2/(n/2))%*%t(M2)))
  M = beta1_d*beta2_d/(sigma1*sigma2)
  selected_index = analys(M, abs(M), q)
  result = fdp_power(selected_index)
  result
}

MDS = function(x, y, num_split){
  n = dim(x)[1]; p = dim(x)[2]
  inclusion_rate_multiple <- matrix(0, nrow = num_split, ncol = p)
  fdr_multiple <- rep(0, num_split)
  power_multiple <- rep(0, num_split)
  num_select <- rep(0, num_split)
  for(iter in 1:num_split){
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    X1 = x[sample_index1, ]
    y1 = y[sample_index1]
    X2 = x[sample_index2, ]
    y2 = y[sample_index2]
    
    cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "mse", nfolds = 10, nlambda = 200, intercept = F, standardize = F)
    lambda <- cvfit$lambda.min/sqrt(2)
    beta1 <- as.vector(glmnet(X1, y1, family = "binomial", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
    beta2 <- as.vector(glmnet(X2, y2, family = "binomial", alpha = 1, lambda = lambda, intercept = F, standardize = F)$beta)
    
    w1 = as.vector(sqrt(f(X1%*%beta1)*(1-f(X1%*%beta1))))
    w2 = as.vector(sqrt(f(X2%*%beta2)*(1-f(X2%*%beta2))))
    X1_beta = diag(w1)%*%X1
    X2_beta = diag(w2)%*%X2
    M1 = get_M(X1_beta)
    M2 = get_M(X2_beta)
    beta1_d = as.numeric(beta1 + 2/n*M1%*%t(x[sample_index1, ])%*%(y[sample_index1] - f(x[sample_index1,] %*% beta1)))
    beta2_d = as.numeric(beta2 + 2/n*M2%*%t(x[sample_index2, ])%*%(y[sample_index2] - f(x[sample_index2,] %*% beta2)))
    sigma1 = sqrt(diag(M1%*%(t(X1)%*%diag(diag((y1-f(X1%*%beta1))%*%t(y1-f(X1%*%beta1))))%*%X1/(n/2))%*%t(M1)))
    sigma2 = sqrt(diag(M2%*%(t(X2)%*%diag(diag((y2-f(X2%*%beta2))%*%t(y2-f(X2%*%beta2))))%*%X2/(n/2))%*%t(M2)))
    M = beta1_d*beta2_d/(sigma1*sigma2)
    current_selected_index = analys(M, abs(M), q)
    num_select[iter] <- length(current_selected_index)
    inclusion_rate_multiple[iter, current_selected_index] <- 1/num_select[iter]
    result = fdp_power(current_selected_index)
    fdr_multiple[iter] = result$fdp
    power_multiple[iter] = result$power
  }
  ### single splitting result
  single_split_fdr <- fdr_multiple[1]
  single_split_power <- power_multiple[1]
  ### multiple splitting result
  inclusion_rate <- apply(inclusion_rate_multiple, 2, mean)
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  null_variable <- numeric()
  for(feature_index in 1:length(feature_rank)){
    if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
      break
    }else{
      null_variable <- c(null_variable, feature_rank[feature_index])
    }
  }
  selected_index <- setdiff(feature_rank, null_variable)
  result = fdp_power(selected_index)
  result
}

BHq = function(X, y){
  tune.1=1.5
  tune.2=1
  n = dim(X)[1]; p = dim(X)[2]
  #cvfit <- cv.glmnet(X, y, family = "binomial", type.measure = "mse", nfolds = 10, nlambda = 200, intercept = F, standardize = F)
  #lambda = cvfit$lambda.1se
  my.logistic.fit = glmnet(x = X, y = y, family = "binomial", alpha = 1,  intercept=F,
                           lambda = 0.1*sqrt(log(p)/n), standardize=F)
  b.hat = coef(my.logistic.fit)
  W.n1 = c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1)
  zeta.try = matrix(ncol = p,nrow = 5)
  tau.try = matrix(ncol = p,nrow = 5)
  V = matrix(ncol=n, nrow = p)
  tau = c()
  for(i in 1:p){
    nodewise.try = glmnet(x= X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5, standardize = F)
    for(lambda.i in 1:5){
      V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.i]
      zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%X[,-i])/sqrt(sum((V[i,])^2*W.n1))))
      tau.try[lambda.i,i] = sqrt(sum((V[i,])^2*W.n1))/(V[i,]%*% X[,i])
    }
    zeta0 = sqrt(2*log(p))
    if(min(zeta.try[,i])>sqrt(2*log(p))) zeta0 = tune.1*min(zeta.try[,i])
    lambda.chosen = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
    tau[i] = tau.try[lambda.chosen,i]
    lambda.chosen = order(nodewise.try$lambda[tau.try[,i]<=tune.2*tau[i]],decreasing = F)[1]
    tau[i] = tau.try[lambda.chosen,i]
    V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.chosen]
  }
  V2 = t((t(V)*W.n1))
  #debaised estimator
  b.check = c()
  for(j in 1:p){
    b.check[j] = b.hat[j+1]+(V2[j,]%*%(y-f(X %*% as.matrix(b.hat)[-1])))/(V[j,] %*% X[,j])
  }
  
  M = b.check/tau
  cutoff = function(x, alpha) p*(1-1*pchisq(x^2, df=1))/max(1,sum(abs(M) > x)) - alpha
  t = sqrt(2*log(p))
  x= seq(0, sqrt(2*log(p)-2*log(log(p))), 0.001)
  for(k in 1:length(x)){
    if(cutoff(x[k], 0.1)<0) t = min(x[k],t)
  }
  selected_index = which(abs(M)>t)
  result = fdp_power(selected_index)
  result
}



KN = function(X, y){
  k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, nlambda=200, family = "binomial")
  Xk = create.second_order(X, method = 'equi')
  M = k_stat(X, Xk, y)
  selected_index <- analys(M, abs(M), q)
  result = fdp_power(selected_index)
  result
}


### Multiple Knockoff
MKN = function(X, y, B = 50, gamma = 0.3){
  k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, nlambda=200, family = "binomial")
  pvalues = matrix(0, B, p)
  for(iter in 1:B){
    pvals_0 = rep(0, p)
    pvals_1 = rep(0, p)
    Xk = create.second_order(X, shrink = T, method = 'equi')
    M = k_stat(X, Xk, y)
    test_score = M
    test_score_inv = -test_score
    for(i in 1:p){
      if(test_score[i]<=0){
        pvals_0[i] = 1
        pvals_1[i] = 1
      }else{
        pvals_0[i] = sum(test_score_inv>=test_score[i])/p
        pvals_1[i] = (1+sum(test_score_inv>=test_score[i]))/p
      }
    }
    pvalues[iter,] = pvals_0
  }
  converted_pvalue = (1 / gamma) * (apply(pvalues,2, FUN = function(x) quantile(x, probs = gamma)))
  sorted_pvalue = sort(converted_pvalue, index.return = T)
  BH_index = max(which(sorted_pvalue$x<=(1:p)*q/p))
  selected_index = sorted_pvalue$ix[1:BH_index]
  result = fdp_power(selected_index)
  result
}


# algorithmic setting
n = 800
p = 2000

replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(replicate)

sig1 = toeplitz(seq(0.3, 0, length.out = p/10))
Sigma = bdiag(rep(list(sig1),10))
diag(Sigma) = 1

s = as.numeric(Sys.getenv("att"))
signal = 4

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
result5 = MKN(X, y)

data_save = list(DS.fdp = result1$fdp, DS.power = result1$power,
                 BHq.fdp = result2$fdp, BHq.power = result2$power,
                 KN.fdp = result3$fdp, KN.power = result3$power,
                 MDS.fdp = result4$fdp, MDS.power = result4$power,
                 MKN.fdp = result5$fdp, MKN.power = result5$power)

### Set the file directory
save(data_save, file = paste("~/result_left/p1_", s, "_replicate_", replicate, ".RData", sep = ""))




