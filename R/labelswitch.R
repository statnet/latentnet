ergmm.labelswitch <- function(Z.K,per.to)
{
  Z.K.new <- Z.K
  for(i in 1:length(per.to))
    Z.K.new[Z.K == per.to[i]] <- i
  return(Z.K.new)
}
