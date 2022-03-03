# average simualtion results over multiple repetitions 
averageRes = function(HPSim_, LPSim_){
  # check whether the number of repetitions is equal in two environments
  stopifnot(length(HPSim_) == length(LPSim_))
  
  # constants 
  stepSec = 1
  nSim = length(HPSim_)
  chunkBreaks = c(0, 4, 8, 12, 16) * 60
  nBreak = length(chunkBreaks)
  nHPStep = delayMaxs[1] / stepSec + 1 # since here we assume the interval is [ , ) yet in we need to include the end point at t = 20
  nLPStep = delayMaxs[2] / stepSec + 1 
  
  # record and average the relative value of waiting
  if(!is_null(HPSim_[[1]]$Qwaits_)){
    ## initialize value function outputs
    HPRvQwaits_ = array(dim = c(nHPStep, nBreak, nSim))
    LPRvQwaits_ = array(dim = c(nLPStep, nBreak, nSim))
    asymHPRvQwaits_ = matrix(NA, nHPStep, nSim)
    asymLPRvQwaits_ = matrix(NA, nLPStep, nSim)
    ## record values 
    for(i in 1 : nSim){
      HPSim = HPSim_[[i]]
      LPSim = LPSim_[[i]]
      for(j in 1 : nBreak){
        HPRvQwaits_[ , j, i] = HPSim$Qwaits_[1 : nHPStep, which.min(abs(HPSim$sellTime - chunkBreaks[j]))] - HPSim$V_[which.min(abs(HPSim$sellTime - chunkBreaks[j]))]
        LPRvQwaits_[ , j, i] = LPSim$Qwaits_[1 : nLPStep, which.min(abs(LPSim$sellTime - chunkBreaks[j]))] - LPSim$V_[which.min(abs(LPSim$sellTime - chunkBreaks[j]))]
      }
      asymHPRvQwaits_[,i] = HPSim$Qwaits_[ , ncol(HPSim$Qwaits_)] - HPSim$V_[ncol(HPSim$Qwaits_)]
      asymLPRvQwaits_[,i] = LPSim$Qwaits_[ , ncol(LPSim$Qwaits_)] - LPSim$V_[ncol(LPSim$Qwaits_)]
    }
    ## average across simulations
    HPRvQwaits = list(
      mu = apply(HPRvQwaits_, FUN = mean, MARGIN = c(1, 2)),
      se = apply(HPRvQwaits_, FUN = calc_se, MARGIN = c(1, 2))
    )
    LPRvQwaits = list(
      mu = apply(LPRvQwaits_, FUN = mean, MARGIN = c(1, 2)),
      se = apply(LPRvQwaits_, FUN = calc_se, MARGIN = c(1, 2))
    )
    asymHPRvQwaits = list(mu = apply(asymHPRvQwaits_, FUN = mean, MARGIN = 1), se = apply(asymHPRvQwaits_, FUN = calc_se, MARGIN = 1))
    asymLPRvQwaits = list(mu = apply(asymLPRvQwaits_, FUN = mean, MARGIN = 1), se = apply(asymLPRvQwaits_, FUN = calc_se, MARGIN = 1))
  }

  
  # record and average AUC
  # initialize auc outputs
  HPaucs_ = matrix(NA, nrow = nBreak - 1, ncol = nSim)
  LPaucs_ = matrix(NA, nrow = nBreak - 1, ncol = nSim)
  asymHPauc_ = vector(length = nSim)
  asymLPauc_ = vector(length = nSim)
  for(i in 1 : nSim){
    HPSim = HPSim_[[i]]
    HPSim$Qwaits_ = NULL
    HPSim$Gs_ = NULL
    HPSim = as.data.frame(HPSim)
    HPaucs_[,i] =  sapply(1 : (nBreak - 1), function(x) kmsc(HPSim[(x-1) * 4 * 60 <= HPSim$sellTime & HPSim$sellTime < x * 4 * 60,], min(delayMaxs), F, kmGrid)$auc)
    asymHPauc_[i] = kmsc(HPSim[116 * 60 <= HPSim$sellTime,], min(delayMaxs), F, kmGrid)$auc
      
    LPSim = LPSim_[[i]]
    LPSim$Qwaits_ = NULL
    LPSim$Gs_ = NULL
    LPSim = as.data.frame(LPSim)
    LPaucs_[,i] = sapply(1 : (nBreak - 1), function(x) kmsc(LPSim[(x-1) * 4 * 60 <= LPSim$sellTime & LPSim$sellTime < x * 4 * 60,], min(delayMaxs), F, kmGrid)$auc)
    asymLPauc_[i] = kmsc(LPSim[116 * 60 <= LPSim$sellTime,], min(delayMaxs), F, kmGrid)$auc
  }
  HPaucs = list(
    mu = apply(HPaucs_, FUN = mean, MARGIN = 1), se = apply(HPaucs_, FUN = calc_se, MARGIN = 1)
  )
  LPaucs = list(
    mu = apply(LPaucs_, FUN = mean, MARGIN = 1), se = apply(LPaucs_, FUN = calc_se, MARGIN = 1)
  )
  
  asymHPauc = list(mu = mean(asymHPauc_), se = calc_se(asymHPauc_))
  asymLPauc = list(mu = mean(asymLPauc_), se = calc_se(asymLPauc_))
  
  # return the output
 if(!is_null(HPSim_[[1]]$Qwaits_)){
   output = list(
     HPRvQwaits = HPRvQwaits,
     LPRvQwaits = LPRvQwaits,
     asymHPRvQwaits = asymHPRvQwaits,
     asymLPRvQwaits = asymLPRvQwaits,
     HPaucs = HPaucs,
     LPaucs = LPaucs,
     asymHPauc = asymHPauc,
     asymLPauc = asymLPauc
   )
 }else{
   output = list(
     HPaucs = HPaucs,
     LPaucs = LPaucs,
     asymHPauc = asymHPauc,
     asymLPauc = asymLPauc
   )
 }
  return(output)
}
  

  
