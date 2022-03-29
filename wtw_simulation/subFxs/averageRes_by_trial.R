# average simualtion results over multiple repetitions 
averageRes_by_trial = function(sim_, record_trials){
  # constants 
  nRep = length(sim_) # number of repetitions
  stepSec = 1 # 
  nStep = length(sim_[[1]]$Qwaits_[,1]) # number of time steps, since here we assume the interval is [ , ) yet in we need to include the end point at t = 20
  nRecord = length(record_trials) # how many times we want to record action values
  
  
  # initialize the outputs
  wait_reps = array(dim = c(nStep, nRecord, nRep))
  quit_reps = matrix(NA, nrow = nRecord, ncol = nRep)
  wait_minus_quit_reps = array(dim = c(nStep, nRecord, nRep))
  duration_reps = matrix(nrow = nRecord, ncol = nRep) # record how long on average to reach each recorded trials
  trialWTW_reps = matrix(nrow = max(record_trials), ncol = nRep) # record trial-wise WTW 
  asym_AUC = vector(length = nRep) # 
 
  
  # loop over reps
  for(i in 1 : nRep){
    for(j in 1 : nRecord){
      record_trial = record_trials[j]
      sim = sim_[[i]]
      wait_reps[, j, i] = sim$Qwaits_[, which.min(abs(sim$sellTime - upper_boundaries[j]))]
      quit_reps[j, i] = sim$V_[which.min(abs(sim$sellTime - upper_boundaries[j]))] 
      wait_minus_quit_reps[ , j, i] =wait_reps[, j, i]  - quit_reps[j, i] 
      
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
  

  
