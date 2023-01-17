KN = function(X, y){
  k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, nlambda=200, family = "binomial")
  Xk = create.second_order(X, method = 'equi')
  M = k_stat(X, Xk, y)
  selected_index <- analys(M, abs(M), q)
  result = fdp_power(selected_index, signal_index)
  result
}