# average simualtion results over multiple repetitions 
averageRes = function(sim_, recorded_timepoints, paras){
  # load experiment parameters
  load('expParas.RData')
  source(file.path("subFxs", "analysisFxs.R"))
  
  # constants 
  nRep = length(sim_) # number of repetitions
  stepSec = 1 # 
  nStep = length(sim_[[1]]$Qwaits_[,1]) # number of time steps, since here we assume the interval is [ , ) yet in we need to include the end point at t = 20
  nRecord = length(recorded_timepoints) # number of analysis time windows 
  tau = paras[2] # inverse temperature 
  
  # initialize the outputs
  # wait_reps = array(dim = c(nStep, nRecord, nRep))
  # quit_reps = matrix(NA, nrow = nRecord, ncol = nRep)
  wait_minus_quit_reps = array(dim = c(nStep, nRecord, nRep))
  sub_auc_reps = matrix(NA, nrow = nRecord, ncol = nRep)
  
  # loop over reps
  for(i in 1 : nRep){
    for(j in 1 : nRecord){
      sim = sim_[[i]]
      # wait_reps[, j, i] = sim$Qwaits_[, which.min(abs(sim$sellTime - recorded[j]))]
      # quit_reps[j, i] = sim$V_[which.min(abs(sim$sellTime - recorded_timepoints[j]))] 
      wait_minus_quit = sim$Qwaits_[, which.min(abs(sim$sellTime - recorded_timepoints[j]))] - sim$V_[which.min(abs(sim$sellTime - recorded_timepoints[j]))] 
      kmsc_res = model_based_kmsc(wait_minus_quit, tau, min(delayMaxs))
      # record these values 
      wait_minus_quit_reps[ , j, i] = wait_minus_quit
      sub_auc_reps[j,i] = kmsc_res$auc
    }
  }
  wait_minus_quit_ = apply(wait_minus_quit_reps, MARGIN = c(1,2), FUN = mean)
  sub_auc_ = apply(sub_auc_reps, MARGIN = 1, FUN = mean)
  # wait_ = apply(wait_reps, MARGIN = c(1,2), FUN = mean)
  # quit_ = apply(quit_reps, MARGIN = 1, FUN = mean)
  output = list(
     sub_auc_ = sub_auc_,
     wait_minus_quit_ = wait_minus_quit_
   )
  return(output)
}
  

  
