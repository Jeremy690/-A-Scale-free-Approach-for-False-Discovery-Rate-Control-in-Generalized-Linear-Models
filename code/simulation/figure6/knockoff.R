KN = function(X, y){
  k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, nlambda=200, family = "poisson")
  Xk = create.second_order(X, shrink = T, method = 'equi')
  M = k_stat(X, Xk, y)
  selected_index <- analys(M, abs(M), q)
  result = fdp_power(selected_index)
  result
}

MKN = function(X, y, B = 25, gamma = 0.3){
  k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, nlambda=200, family = "poisson")
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