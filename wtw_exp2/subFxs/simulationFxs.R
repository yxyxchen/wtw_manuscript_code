getSimFun = function(modelName){
  if(modelName == "QL1") simFun = QL1
  else if(modelName == "QL2") simFun = QL2
  else if(modelName == "RL1") simFun = RL1
  else if(modelName == "RL2") simFun = RL2
  else if(modelName == "BL") simFun = BL
  else{
    return("wrong model name!")
  }
  return(simFun)
}

# change the parameter 
# change the decision rule
simulateUnitSingle = function(paras, nSim, modelName, cb, blockDuration){
  tGrid = seq(0, blockDuration * 60, by = 5)
  
  # get simFun
  simFun = getSimFun(modelName)
  
  # initialize outputs
  auc_ = vector(length = nSim)
  auc1_ = vector(length = nSim)
  auc2_ = vector(length = nSim)
  auc3_ = vector(length = nSim)
  wtw_ = matrix(NA, nrow = length(tGrid) * length(cb), ncol = nSim)
  reRate_ = matrix(NA, nrow = blockDuration * 60 / iti, ncol = nSim)
  # usually can not use foreach to fill a matrix or a vector
  for(i in 1 : nSim){
    set.seed(i)
    thisTrialData = simFun(paras, cb, blockDuration)
    thisTrialData$Qwaits = NULL
    thisTrialData = as.data.frame(thisTrialData)
    kmscResults = kmsc(thisTrialData, min(tMaxs), "", F, kmGrid)
    auc_[i] = kmscResults$auc
    wtwtsResults = wtwTS(thisTrialData, tGrid, min(tMaxs), "", F )
    wtw_[,i] = wtwtsResults$timeWTW
    junk = nrow(thisTrialData)
    reRate_[1 : length(thisTrialData$trialNum),i] = thisTrialData$reRates
    
    tempt = lapply(1:3, function(i) kmsc(thisTrialData[thisTrialData$sellTime < blockDuration /3 * i * 60
                                                       &thisTrialData$sellTime >= blockDuration /3 * (i-1)*60,],
                                         min(tMaxs), "", F, kmGrid))
    auc1_[i] = tempt[[1]]$auc
    auc2_[i] = tempt[[2]]$auc
    auc3_[i] = tempt[[3]]$auc
  }
  # summarise 
  outputs = list(auc = mean(auc_),
                 aucSD = sd(auc_),
                 wtw = apply(wtw_, MARGIN = 1, mean),
                 wtwSD = apply(wtw_, MARGIN = 1, sd),
                 auc1 = mean(auc1_),
                 auc2 = mean(auc2_),
                 auc3 = mean(auc3_),
                 aucSD1 = sd(auc1_),
                 aucSD2 = sd(auc2_),
                 aucSD3 = sd(auc3_),
                 reRate = apply(reRate_, MARGIN = 1, FUN = function(x) mean(x, na.rm = T)))
  return(outputs)
}

simulateUnit = function(paras, nSim, modelName, cb){
  # get simFun
  simFun = getSimFun(modelName)
  
  # initialize outputs
  aucHP_ = vector(length = nSim)
  aucLP_ = vector(length = nSim)
  wtwHP_ = matrix(NA, nrow = length(tGrid), ncol = nSim)
  wtwLP_ = matrix(NA, nrow = length(tGrid), ncol = nSim)
  reRate_ = vector(length = nSim)
  # usually can not use foreach to fill a matrix or a vector
  for(i in 1 : nSim){
    set.seed(i)
    thisTrialData = simFun(paras, cb)
    thisTrialData$Qwaits = NULL
    thisTrialData = as.data.frame(thisTrialData)
    # HP
    kmscResults = kmsc(thisTrialData[thisTrialData$condition == "HP",], min(tMaxs), "", F, kmGrid)
    aucHP_[i] = kmscResults$auc
    wtwtsResults = wtwTS(thisTrialData[thisTrialData$condition == "HP",], tGrid, min(tMaxs), "", F )
    wtwHP_[,i] = wtwtsResults$timeWTW
    # LP
    kmscResults = kmsc(thisTrialData[thisTrialData$condition == "LP",], min(tMaxs), "", F, kmGrid)
    aucLP_[i] = kmscResults$auc
    wtwtsResults = wtwTS(thisTrialData[thisTrialData$condition == "LP",], tGrid, min(tMaxs), "", F )
    wtwLP_[,i] = wtwtsResults$timeWTW
    
    junk = nrow(thisTrialData)
    
    reRate_[i] = mean(thisTrialData$reRates[(junk - 10) : junk])
  }
  # summarise 
  outputs = list(aucHP = mean(aucHP_),
                 aucLP = mean(aucLP_),
                 aucHPSD = sd(aucHP_),
                 aucLPSD = sd(aucLP_),
                 wtwHP = apply(wtwHP_, MARGIN = 1, mean),
                 wtwLP = apply(wtwLP_, MARGIN = 1, mean),
                 wtwHPSD = apply(wtwHP_, MARGIN = 1, sd),
                 wtwLPSD = apply(wtwHP_, MARGIN = 1, sd),
                 reRate = mean(reRate_))
  return(outputs)
}
RL2 = function(paras, cb, blockDuration){
  # parse para
  phi = paras[1]; nega = paras[2]; tau = paras[3]; prior = paras[4]
  beta = paras[5]; 
  
  # prepare inputs
  nTrialMax = blockDuration * 60 / iti *2
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  Viti = 0; reRate = wIni 
  Qwait = prior*0.1 - 0.1*(0 : (nTimeStep - 1)) + Viti
  
  # recording variables
  reRates = vector(length = nTrialMax); reRates[1] = reRate
  Vitis = vector(length = nTrialMax); Vitis[1] = Viti
  Qwaits = matrix(NA, nrow = nTimeStep, ncol = nTrialMax); Qwaits[,1] = Qwait 
  
  # initialize outputs 
  trialEarnings = rep(0, nTrialMax); timeWaited = rep(0, nTrialMax); sellTime = rep(0, nTrialMax);
  condition = rep(0, nTrialMax); scheduledWait = rep(0, nTrialMax);
  
  
  # loop over blocks
  tIdx = 1
  for(bkIdx in 1 : length(cb)){
    # determine distrib
    seq = c()
    distrib = ifelse(cb[bkIdx] == "HP", "unif16", "exp32")
    elapsedTime = 0
    # loop over trials
    while(elapsedTime <= blockDuration * 60) {
      # determine scheduledWait
      junk = drawSample(distrib, seq)
      thisScheduledWait = junk[['delay']]
      seq = junk[['seq']]
      # loop over timesteps
      t = 1
      while(t <= nTimeStep){
        # determine At
        pWait =  1 / sum(1  + exp((Viti - Qwait[t])* tau))
        action = ifelse(runif(1) < pWait, 'wait', 'quit')
        # observe St+1 and Rt+1
        rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
        getReward = (action == 'wait' && rewardOccur);
        nextReward = ifelse(getReward, tokenValue, 0) 
        # dertime whether St+1 is the terminal state
        nextStateTerminal = (getReward || action == "quit")
        if(nextStateTerminal){
          elapsedTime = elapsedTime + ifelse(getReward, thisScheduledWait, t * stepDuration) + iti
          # only record values before the end of the block
          if(elapsedTime<= blockDuration * 60){
            T = t+1
            trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
            timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
            sellTime[tIdx] = elapsedTime - iti
            condition[tIdx] = cb[bkIdx]
            scheduledWait[tIdx] = thisScheduledWait
          }
          break
        }else{
          t = t + 1
        }
      }# end of the loop over timesteps
      # update action values before the end of the block
      if(elapsedTime <= blockDuration * 60){
        returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
        if(getReward){
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }else{
          if(T > 2){
            Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi * nega * (returns[1 : (T-2)] - Qwait[1 : (T-2)])
          }
        }
        # update Viti
        delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
        Viti = ifelse(nextReward > 0, Viti + phi * delta, Viti + phi * nega * delta)
        
        # update reRate 
        reRate = ifelse(nextReward > 0, reRate + beta * delta, reRate + beta * nega * delta)   
        
        # record variables
        reRates[tIdx + 1] = reRate
        Vitis[tIdx + 1] = Viti
        Qwaits[,tIdx + 1] = Qwait
        # update tIdx and go to the next trial
        tIdx = tIdx + 1
      }
    } # end of the loop over trials
  }# end of the loop over blocks
  
  # cut off the last trial and return outputs
  trialNum = as.vector(sapply(1 : length(cb), function(i) 1 : sum(condition == cb[i])))
  outputs = list( 
    "trialNum" = trialNum, "trialEarnings" = trialEarnings[1 : (tIdx - 1)],
    "timeWaited" = timeWaited[1 : (tIdx - 1)], "sellTime" = sellTime[1 : (tIdx - 1)],
    "scheduledWait" = scheduledWait[1 : (tIdx - 1)], "condition" = condition[1 : (tIdx - 1)],
    "reRates" = reRates[1: (tIdx - 1)], "Vitis" = Vitis[1 : (tIdx - 1)],
    "Qwaits" = Qwaits[,1 : (tIdx - 1)]
  )
  return(outputs)
}
