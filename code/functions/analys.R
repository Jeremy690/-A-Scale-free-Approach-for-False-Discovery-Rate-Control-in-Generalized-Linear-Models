analys = function(mm, ww, q){
  ## mm:Tstat defined in the paper
  ## ww:Absolute value of Tstat
  t_set = max(ww)
  ## Choose the maximum value of ww as threshold
  for(t in ww){
    ## Calculate the estimates of FDP.
    ps = length(mm[mm >= t])
    ng = length(na.omit(mm[mm <= -t]))
    rto = (ng + 1)/max(ps, 1)
    if(rto <= q){
      t_set = c(t_set, t)
    }
  }
  ## Choose the minimum threshold
  thre = min(t_set)
  nz_est = which(mm >= thre)
  return(as.vector(nz_est))
}