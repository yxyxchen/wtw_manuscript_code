library("stringr")
# return learning parameters for each model 
getParaNames = function(modelName){
  if(modelName == "QL1") paraNames = c("alpha", "tau", "gamma", "eta")
  else if(modelName == "QL2") paraNames = c("alpha", "rho", "tau", "gamma", "eta")
  else if(modelName == "RL1") paraNames = c("alpha", "tau", "eta", "beta")
  else if(modelName == "RL2") paraNames = c("alpha", "rho", "tau", "eta", "beta")
  else if(modelName == "naive") paraNames = c("theta")
  else if(modelName == "omni") paraNames = c("tau")
  return(paraNames)
}

# check MCMC fitting results 
checkFit = function(paraNames, expPara){
  ids = expPara$id
  # detect participants with high Rhats 
  RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(paraNames)] # columns recording Rhats
  if(length(RhatCols) > 1){
    high_Rhat_ids = ids[apply(expPara[,RhatCols] >= 1.01, MARGIN = 1, sum) > 0]
  }else{
    high_Rhat_ids = ids[expPara[,RhatCols] >= 1.01 ]
  }
  
  # detect participants with low ESSs
  ESSCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(paraNames)]# columns recording ESSs
  if(length(ESSCols) > 1){
    low_ESS_ids = ids[apply(expPara[,ESSCols] < (4 * 100), MARGIN = 1, sum) > 0]
  }else{
    low_ESS_ids = ids[expPara[,ESSCols] < (4 * 100)]
  }

  # detect divergent transitions
  dt_ids = ids[expPara$nDt > 0]
  
  # identify participants satisifying all three criteria:
  passCheck = !ids %in% unique(c(dt_ids, high_Rhat_ids, low_ESS_ids))
  
  return(passCheck)
}

# resample pair-wise sequences
# inputs:
# ys: y in the original sequence
# xs: x in the original sequence
# Xs: x in the new sequence
# outputs: 
# Ys : y in the new sequence 
resample = function(ys, xs, Xs){
  isBreak = F
  # initialize Ys
  Ys = rep(NA, length = length(Xs))
  for(i in 1 : length(Xs)){
    # for each X in Xs
    X = Xs[i]
    # find the index of cloest x value on the right
    # if closest_right_x_idx exists 
    if(X <= tail(xs,1)) {
      # Y takes the corresonding y value
      closest_right_x_idx = min(which(xs >= X))
      Ys[i] = ys[closest_right_x_idx]
    }else{
      isBreak = T
      lastY = i - 1
      break
    }
  }
  # fill the remaining elements in Ys by the exisiting last element
  if(isBreak){
    Ys[(lastY + 1) : length(Xs)] = Ys[lastY]
  }
  return(Ys)
}

