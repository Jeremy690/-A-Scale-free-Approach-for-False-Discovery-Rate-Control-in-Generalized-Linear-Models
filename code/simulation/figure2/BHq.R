BHq <- function(X, y){
  n = dim(X)[1]; p = dim(X)[2]
  fit <- glm(y ~ X - 1, family = 'binomial')
  pvalues <- summary(fit)$coefficients[,4]
  sorted_pvalues = sort(pvalues, decreasing = F, index.return = T)
  if(sum(sorted_pvalues$x <= (1:p)*q/p) > 0){
    BHq_index = max(which(sorted_pvalues$x <= (1:p)*q/p))
    select_index = sorted_pvalues$ix[1:BHq_index]
    num_select <- length(select_index)
    return(list(select_index = select_index, num_select = num_select))
  }else{
    return(list(select_index = NULL, num_select = 0))
  }
}
