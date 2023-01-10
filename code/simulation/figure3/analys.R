analys = function(mm, ww, q){
  t_set = max(ww)
  for(t in ww){
    ps = length(mm[mm >= t])
    ng = length(na.omit(mm[mm <= -t]))
    rto = (ng + 1)/max(ps, 1)
    if(rto <= q){
      t_set = c(t_set, t)
    }
  }
  thre = min(t_set)
  nz_est = which(mm >= thre)
  return(as.vector(nz_est))
}