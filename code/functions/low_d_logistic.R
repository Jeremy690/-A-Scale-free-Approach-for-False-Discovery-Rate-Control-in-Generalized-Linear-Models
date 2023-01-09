DS <- function(X, y, q, signal_index){
  ## X: design matrix
  ## y: response variable
  ## q: designated fdr level
  ## signal_index: true signal index
  n = dim(X)[1]; p = dim(X)[2]
  ## Split the data into two halves and run MLE
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  fit1 <- glm(y[sample_index1] ~ X[sample_index1,] - 1, family = 'binomial')
  fit2 <- glm(y[sample_index2] ~ X[sample_index2,] - 1, family = 'binomial')
  beta1 <- fit1$coefficients
  beta2 <- fit2$coefficients
  ## Calculate the variance of each estimator(up to a scaling factor)
  tau1 <- 1/sqrt(diag(solve(t(X[sample_index1, ]) %*% X[sample_index1, ])))
  tau2 <- 1/sqrt(diag(solve(t(X[sample_index2, ]) %*% X[sample_index2, ])))
  M <- beta1*beta2*tau1*tau2
  ## Get the selected feature
  select_index <- analys(M, abs(M), q)
  num_select <- length(select_index)
  ## Calculate fdp and power.
  result = fdp_power(select_index, signal_index)
  return(list(fdp = result$fdp, power = result$power, select_index = select_index, num_select = num_select))
}

MDS <- function(X, y, num_split, q, signal_index){
  ## X: design matrix
  ## y: response variable
  ## num_split: The number of split
  ## q: designated fdr level
  ## signal_index: true signal index
  n = dim(X)[1]; p = dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  num_select <- rep(0, num_split)
  ## Calculate the inclusion rate for every split
  for(iter in 1:num_split){
    DS_result <- DS(X, y, q, signal_index)
    inclusion_rate[iter, DS_result$select_index] <- 1/DS_result$num_select
  }
  ## Calculate the inclusion rate
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  ## Sort the features according to inclusion rate 
  feature_rank <- order(inclusion_rate)
  ## Drop features with zero inclusion rate
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  null_feature <- numeric()
  ## Find the cutoff
  for(feature_index in 1:length(feature_rank)){
    if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
      break
    }else{
      null_feature <- c(null_feature, feature_rank[feature_index])
    }
  }
  select_index <- setdiff(feature_rank, null_feature)
  num_select <- length(select_index)
  ## Calculate fdp and power
  result = fdp_power(select_index, signal_index)
  return(list(fdp = result$fdp, power = result$power))
}