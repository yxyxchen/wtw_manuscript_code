# average simualtion results over multiple repetitions 
averageRes = function(sim_, lower_boundaries, upper_boundaries){
  # constants 
  nRep = length(sim_) # number of repetitions
  stepSec = 1 # 
  nStep = length(sim_[[1]]$Qwaits_[,1]) # number of time steps, since here we assume the interval is [ , ) yet in we need to include the end point at t = 20
  nBreak = length(lower_boundaries) # number of analysis time windows 
  
  # initialize the outputs
  wait_reps = array(dim = c(nStep, nBreak, nRep))
  quit_reps = matrix(NA, nrow = nBreak, ncol = nRep)
  wait_minus_quit_reps = array(dim = c(nStep, nBreak, nRep))
  sub_auc_reps = matrix(NA, nrow = nBreak, ncol = nRep)
  # loop over reps
  for(i in 1 : nRep){
    for(j in 1 : nBreak){
      sim = sim_[[i]]
      wait_reps[, j, i] = sim$Qwaits_[, which.min(abs(sim$sellTime - upper_boundaries[j]))]
      quit_reps[j, i] = sim$V_[which.min(abs(sim$sellTime - upper_boundaries[j]))] 
      wait_minus_quit_reps[ , j, i] =wait_reps[, j, i]  - quit_reps[j, i] 
      sim$Qwaits_ = NULL
      sim$Gs_ = NULL
      sim = as.data.frame(sim)
      sub_auc_reps[j,i] = kmsc(sim[lower_boundaries[j] < sim$sellTime & sim$sellTime <= upper_boundaries[j],], min(delayMaxs), F, kmGrid)$auc
    }
  }
  wait_minus_quit_ = apply(wait_minus_quit_reps, MARGIN = c(1,2), FUN = mean)
  sub_auc_ = apply(sub_auc_reps, MARGIN = 1, FUN = mean)
  wait_ = apply(wait_reps, MARGIN = c(1,2), FUN = mean)
  quit_ = apply(quit_reps, MARGIN = 1, FUN = mean)
  output = list(
     sub_auc_ = sub_auc_,
     wait_minus_quit_ = wait_minus_quit_,
     wait_ = wait_,
     quit_ = quit_
   )
  return(output)
}
  

  
