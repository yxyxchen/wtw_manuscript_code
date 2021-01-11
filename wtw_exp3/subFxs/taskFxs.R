# this scripts contains some supporting functions for runing the simulation
# and transfer vaWaits and vaQuits
drawSample = function(cond){
  # % generates a sample from next quantile of the designated distribution
  # % timing parameters are specified within this function

  k = pareto[['k']]
  mu = pareto[['mu']]
  sigma = pareto[['sigma']]
  
  if(cond == 'HP'){
    sample = runif(1, min = 0, max = tMaxs[1])
  }else{
    sample = min(mu + sigma * (runif(1) ^ (-k) - 1) / k, tMaxs[2])
  }
  return(sample)
}

