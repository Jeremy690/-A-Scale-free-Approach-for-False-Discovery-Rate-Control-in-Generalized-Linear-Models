KN = function(X, y){
  k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, nlambda=200, family = "binomial")
  result = knockoff.filter(X, y, statistic=k_stat, fdr=q)
  result$statistic
}