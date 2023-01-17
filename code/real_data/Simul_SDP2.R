rm(list = ls())
library(mvnfast)
library(glmnet)
library(CVXR)

replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(replicate)

is_posdef = function(A, tol=1e-9) {
  p = nrow(matrix(A))
  
  if (p<500) {
    lambda_min = min(eigen(A)$values)
  }
  else {
    oldw <- getOption("warn")
    options(warn = -1)
    lambda_min = RSpectra::eigs(A, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values
    options(warn = oldw)
    if( length(lambda_min)==0 ) {
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(A)$values)
    }
  }
  return (lambda_min>tol*10)
}

create.multi_sdp = function(Sigma, k, gaptol=1e-6, maxit=1000, verbose=FALSE){
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),] 
  Cs = c((k+1)/k*G)
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef(2*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) {
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")
  
  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}




create.sdp = function(X, k, shrink = T){
  # Estimate the mean vectorand covariance matrix
  mu = colMeans(X)
  
  # Estimate the covariance matrix
  if(!shrink) {
    Sigma = cov(X)
    # Verify that the covariance matrix is positive-definite
    if(!is_posdef(Sigma)) {
      shrink=TRUE
    }
  }
  if(shrink) {
    if (!requireNamespace('corpcor', quietly=T))
      stop('corpcor is not installed', call.=F)
    Sigma = tryCatch({suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X,verbose=F)), nrow=ncol(X)))},
                     warning = function(w){}, error = function(e) {
                       stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
                     }, finally = {})
  }
  
  p = dim(X)[2]
  s = create.multi_sdp(Sigma, k)
  diag_s = diag(s)
  D  = diag_s
  if(all(diag_s==0)) {
    warning("The conditional knockoff covariance matrix is not positive definite. Knockoffs will have no power.")
    return(X)
  }
  SigmaInv_s = solve(Sigma,diag_s)
  mu_k = X - sweep(X,2,mu,"-") %*% SigmaInv_s
  mu = mu_k
  if(k>1){
    for(i in 1:(k-1)){
      mu = cbind(mu, mu_k)
    }
  }
  C = 2*D-D%*%SigmaInv_s
  tilde_Sigma = matrix(0, k*p, k*p)
  for(i in 1:k){
    for(j in 1:k){
      if(i == j){
        tilde_Sigma[1:p+(i-1)*p, 1:p+(j-1)*p] = C
      }else{
        tilde_Sigma[1:p+(i-1)*p, 1:p+(j-1)*p] = C-D
      }
    }
  }
  X_k = mu+matrix(rnorm(n*k*p), n)%*%chol(tilde_Sigma)
  X_k
}

analys = function(kappa, tau, k, q){
  t_set = c(NULL)
  for(t in tau){
    ps = sum((tau>t)*(kappa == 1))
    ng = sum((tau>t)*(kappa != 1))
    rto = (1/k+1/k*ng)/max(ps, 1)
    if(rto<=q){
      t_set = c(t_set, t)
    }
  }
  if(length(t_set)==0){
    nz_est = c(NULL)
  }else{
    thre = min(t_set)
    nz_est = which(tau>=thre& kappa == 1)
  }
  nz_est
}

fdp_power <- function(selected_index){
  if(length(selected_index) == 0){
    return(list(fdp = 0, power = 0))
  }else{
    fdp <- (length(selected_index) - length(intersect(selected_index, signal_index)))/max(length(selected_index), 1)
    power <- length(intersect(selected_index, signal_index))/max(length(signal_index),1)
    return(list(fdp = fdp, power = power))
  }
}

knockoff_sdp = function(y, X, k, q){
  X_k = create.sdp(X, k)
  X_tilde = cbind(X, X_k)
  cvfit = cv.glmnet(X_tilde, y, nfolds = 10, alpha = 0, family = 'binomial')
  fit = glmnet(X_tilde, y, alpha = 0, lambda =cvfit$lambda.min, family = 'binomial')
  T_stat = abs(matrix(fit$beta, ncol = k+1))
  tau = apply(T_stat, 1, sort)[k+1,]-apply(T_stat, 1, sort)[k,]
  kappa = apply(T_stat, 1, which.max)
  nonzero_index  = analys(kappa, tau, k, q)
  nonzero_index
}

data = t(read.table('~/data.txt'))
y <- c(rep(0, 400), rep(1, 2000))
freq <- colSums(data != 0)/nrow(data)
remove_index <- which(freq < 0.1)
X <- data[, -remove_index]
std <- apply(X, 2, sd)
gene_index <- order(std, decreasing = T)
p <- 500
n <- nrow(X)
X <- X[, gene_index[1:p]]
k = 2
nonzero_index = knockoff_sdp(y, X, k, 0.1)
data_save = nonzero_index
save(data_save, file = paste('~/simul_KN/SDP_', k, '_replicate_', replicate,  '.RData', sep = ''))
