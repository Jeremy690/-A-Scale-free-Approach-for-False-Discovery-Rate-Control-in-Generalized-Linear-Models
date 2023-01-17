f = function(x){
  exp(x)/(1+exp(x))
}


fdp_power <- function(selected_index, signal_index){
  fdp <- (length(selected_index) - length(intersect(selected_index, signal_index)))/max(length(selected_index), 1)
  power <- length(intersect(selected_index, signal_index))/length(signal_index)
  return(list(fdp = fdp, power = power))
}

analys = function(mm, ww, q){
  t_set = max(ww)
  for(t in ww){
    ps = length(mm[mm>=t])
    ng = length(na.omit(mm[mm<=-t]))
    rto = (ng+1)/max(ps, 1)
    if(rto<=q){
      t_set = c(t_set, t)
    }
  }
  thre = min(t_set)
  nz_est = which(mm>thre)
  nz_est
}
