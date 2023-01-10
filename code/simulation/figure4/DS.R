DS <- function(X, y){
  n = dim(X)[1]; p = dim(X)[2]
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  fit1 <- glm(y[sample_index1] ~ X[sample_index1,] - 1, family = negative.binomial(theta = 2))
  fit2 <- glm(y[sample_index2] ~ X[sample_index2,] - 1, family = negative.binomial(theta = 2))
  beta1 <- fit1$coefficients
  beta2 <- fit2$coefficients
  tau1 <- 1/sqrt(diag(solve(t(X[sample_index1, ]) %*% X[sample_index1, ])))
  tau2 <- 1/sqrt(diag(solve(t(X[sample_index2, ]) %*% X[sample_index2, ])))
  M <- beta1*beta2*tau1*tau2
  select_index <- analys(M, abs(M), q)
  num_select <- length(select_index)
  return(list(select_index = select_index, num_select = num_select))
}
