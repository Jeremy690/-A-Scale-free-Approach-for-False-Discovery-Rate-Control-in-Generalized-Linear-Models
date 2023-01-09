## Calculate fdp and power given selected feature index and the true feature index
fdp_power <- function(selected_index, signal_index){
  fdp <- (length(selected_index) - length(intersect(selected_index, signal_index)))/max(length(selected_index), 1)
  power <- length(intersect(selected_index, signal_index))/length(signal_index)
  return(list(fdp = fdp, power = power))
}