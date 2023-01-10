GM <- function(X, y){
  n = dim(X)[1]; p = dim(X)[2]
  tau <- 1/sqrt(diag(solve(t(X) %*% X/n)))
  M <- rep(0, p)
  for(j in 1:p){
    Z <- rnorm(n)
    P <- X[, -j]%*%solve(t(X[, -j])%*%X[, -j])%*%t(X[, -j])
    c <- sqrt((sum(X[, j]^2) - sum((P%*%X[, j])^2))/(sum(Z^2) - sum((P%*%Z)^2)))
    fit <- glm(y ~ cbind(X[,j] + c*Z, X[,j] - c*Z, X[, -j]) - 1, family = negative.binomial(theta = 2))
    beta1 <- fit$coefficients[1]
    beta2 <- fit$coefficients[2]
    M[j] <- beta1*beta2*(tau[j]^2 + c^2)
  }
  select_index <- analys(M, abs(M), q)
  num_select <- length(select_index)
  return(list(select_index = select_index, num_select = num_select))
}
