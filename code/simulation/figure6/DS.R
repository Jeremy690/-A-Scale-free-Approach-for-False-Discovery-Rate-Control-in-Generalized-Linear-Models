get_M = function(data){
  n = dim(data)[1]
  p = dim(data)[2]
  C = matrix(0, nrow = p, ncol = p)
  tau_square = numeric(p)
  for(j in 1:p){
    y = data[,j]
    X <- data[, -j]
    if(j == 1){
      cvfit <- cv.glmnet(X, y, type.measure = "mse", nfolds = 10, intercept = F, standardize = F)
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
  sample_index1 <- sample(x = c(1:n), size = floor(0.5 * n), replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  n1 = length(sample_index1)
  n2 = length(sample_index2)
  X1 = x[sample_index1, ]
  y1 = y[sample_index1]
  X2 = x[sample_index2, ]
  y2 = y[sample_index2]
  cvfit <- cv.glmnet(X1, y1, family = "poisson", type.measure = "mse", nfolds = 10, intercept = F, standardize = F)
  lambda <- cvfit$lambda.min
  beta1 = as.vector(coef(cvfit, lambda))[-1]
  cvfit <- cv.glmnet(X2, y2, family = "poisson", type.measure = "mse", nfolds = 10, intercept = F, standardize = F)
  lambda <- cvfit$lambda.min
  beta2 = as.vector(coef(cvfit, lambda))[-1]
  w1 = as.vector(sqrt(f(X1%*%beta1)))
  w2 = as.vector(sqrt(f(X2%*%beta2)))
  X1_beta = diag(w1)%*%X1
  X2_beta = diag(w2)%*%X2
  M1 = get_M(X1_beta)
  M2 = get_M(X2_beta)
  beta1_d = as.numeric(beta1 + 1/n1*M1%*%t(x[sample_index1, ])%*%(y[sample_index1] - f(x[sample_index1,] %*% beta1)))
  beta2_d = as.numeric(beta2 + 1/n2*M2%*%t(x[sample_index2, ])%*%(y[sample_index2] - f(x[sample_index2,] %*% beta2)))
  sigma1 = sqrt(diag(M1%*%(t(X1)%*%diag(diag((y1-f(X1%*%beta1))%*%t(y1-f(X1%*%beta1))))%*%X1/n1)%*%t(M1)))
  sigma2 = sqrt(diag(M2%*%(t(X2)%*%diag(diag((y2-f(X2%*%beta2))%*%t(y2-f(X2%*%beta2))))%*%X2/n2)%*%t(M2)))
  M = beta1_d*beta2_d/(sigma1*sigma2)
  selected_index = analys(M, abs(M), q)
  selected_index
}

MDS = function(x, y, num_split){
  n = dim(x)[1]; p = dim(x)[2]
  inclusion_rate_multiple <- matrix(0, nrow = num_split, ncol = p)
  fdr_multiple <- rep(0, num_split)
  power_multiple <- rep(0, num_split)
  num_select <- rep(0, num_split)
  for(iter in 1:num_split){
    current_selected_index = DS_VDG(x, y)
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
  list(DS.fdp = single_split_fdr, DS.power = single_split_power,
       MDS.fdp = result$fdp, MDS.power = result$power)
}
