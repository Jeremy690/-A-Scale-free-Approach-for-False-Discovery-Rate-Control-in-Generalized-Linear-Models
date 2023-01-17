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
  result = fdp_power(selected_index, signal_index)
  result
}