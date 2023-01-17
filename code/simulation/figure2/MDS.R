MDS <- function(X, y, num_split){
  ## X: design matrix
  ## y: response variable
  ## num_split: The number of split
  n = dim(X)[1]; p = dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  num_select <- rep(0, num_split)
  ## Calculate the inclusion rate for every split
  for(iter in 1:num_split){
    DS_result <- DS(X, y)
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
  return(list(select_index = select_index, num_select = num_select))
}
