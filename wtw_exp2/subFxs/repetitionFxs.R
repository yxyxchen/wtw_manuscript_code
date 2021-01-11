getRepFun = function(modelName){
  if(modelName == "QL1") repFun = QL1
  else if(modelName == "QL2") repFun = QL2
  else if(modelName == "RL1") repFun = RL1
  else if(modelName == "RL2") repFun = RL2
  else if(modelName == "BL") repFun = BL
  else{
    return("wrong model name!")
  }
  return(repFun)
}

# this function replicates the behavioral data using inividual fitted parameters
modelRepitation = function(modelName, summaryData, trialData,  nComb){
  # determine repFun
  repFun = getRepFun(modelName)
  
  # load inividual fitted parameters
  paraNames = getParaNames(modelName)
  parentDir ="../../genData/wtw_exp2/expModelFitting"; dirName = sprintf("%s/%sdb",parentDir, modelName)
  expPara = loadExpPara(paraNames, dirName)
  ids = expPara$id
  nSub = length(ids)
  
  # initialize outputs
  repTrialData = vector(length = nSub * nComb, mode ='list')
  repNo = matrix(1 : (nSub * nComb), nrow = nComb, ncol = nSub)
  
  # draw nComb parameter combination samples, and simualte nComb times for each participant
  set.seed(231)
  for(sIdx in 1 : nSub){
    # prepare inputs
    id = ids[[sIdx]]
    paras_ = read.table(sprintf("%s/%sdb/s%s.txt", parentDir, modelName, id),sep = ",", row.names = NULL)
    thisTrialData = trialData[[id]] # here we id instead of sIdx
    # excluded some trials
    excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                              thisTrialData$condition == conditions[1])
    excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                              thisTrialData$condition == conditions[2])
    excluedTrials = c(excluedTrials1, excluedTrials2)
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials &
                                    thisTrialData$blockNum <= 2,]
    condtion = thisTrialData$condition
    scheduledWait = thisTrialData$scheduledWait
    # simulate nComb times
    for(cbIdx in 1 : nComb){
      paras = as.double(paras_[sample(1 : nrow(paras_), 1), 1 : length(paraNames)])
      tempt = repFun(paras, condtion, scheduledWait)
      repTrialData[[repNo[cbIdx, sIdx]]] = tempt
    }
  }
  outputs = list(expPara = expPara, repTrialData = repTrialData, repNo = repNo)
  return(outputs)
}


QL1 = function(paras, condtion, scheduledWait){
  # parse paras
  phi = paras[1]; tau = paras[2]; gamma = paras[3]; prior = paras[4]
  
  # prepare inputs
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9) * subOptimalRatio
  Viti = wIni 
  Qwait = prior*0.1 - 0.1*(0 : (nTimeStep - 1)) + Viti
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  
  # initialize variables for recording targets and deltas in updating Qwait
  deltas = matrix(NA, nTimeStep, nTrial)
  targets = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # select actions
    t = 1
    thisScheduledWait = scheduledWait[tIdx]
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
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update action values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        targets[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          targets[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = gamma^(iti / stepDuration) * returns[1] - Viti
      Viti = Viti + phi * delta
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "targets" = targets, "deltas" = deltas,
    "Vitis" = Vitis, "condtion" = condtion
  )
  return(outputs)
}

QL2 = function(paras, condtion, scheduledWait){
  # parse paras
  phi = paras[1]; nega = paras[2]; tau = paras[3]; gamma = paras[4]; prior = paras[5]
  
  # prepare inputs
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9) * subOptimalRatio
  Viti = wIni 
  Qwait = prior*0.1 - 0.1*(0 : (nTimeStep - 1)) + Viti
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  
  # initialize variables for recording targets and deltas in updating Qwait
  deltas = matrix(NA, nTimeStep, nTrial)
  targets = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # select actions
    t = 1
    thisScheduledWait = scheduledWait[tIdx]
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
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update action values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        targets[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          targets[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi * nega * (returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = gamma^(iti / stepDuration) * returns[1] - Viti
      Viti = ifelse(nextReward > 0, Viti + phi * delta, Viti + phi * nega * delta)
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "targets" = targets, "deltas" = deltas,
    "Vitis" = Vitis, "condtion" = condtion
  )
  return(outputs)
}

RL1 = function(paras, condtion, scheduledWait){
  # parse para
  phi = paras[1]; tau = paras[2]; prior = paras[3]; beta = paras[4]
  
  # prepare inputs
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  Viti = 0; reRate = wIni 
  Qwait = prior*0.1 - 0.1*(0 : (nTimeStep - 1)) + Viti
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  reRates = vector(length = nTrial); reRates[1] = reRate
  
  # initialize variables for recording targets and deltas in updating Qwait
  deltas = matrix(NA, nTimeStep, nTrial)
  targets = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # select actions
    t = 1
    thisScheduledWait = scheduledWait[tIdx]
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
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update action values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
      if(getReward){
        targets[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          targets[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
      Viti = Viti + phi * delta
      # update reRate 
      reRate = reRate + beta * delta  
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Vitis[tIdx + 1] = Viti
      reRates[tIdx + 1] = reRate
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "targets" = targets, "deltas" = deltas,
    "Vitis" = Vitis, "reRates" = reRates, "condtion" = condtion
  )
  return(outputs)
}

RL2 = function(paras, condtion, scheduledWait){
  # parse para
  phi = paras[1]; nega = paras[2]; tau = paras[3]; prior = paras[4]
  beta = paras[5]; 
  
  # prepare inputs
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  Viti = 0; reRate = wIni 
  Qwait = prior*0.1 - 0.1*(0 : (nTimeStep - 1)) + Viti
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  reRates = vector(length = nTrial); reRates[1] = reRate
  
  # initialize variables for recording targets and deltas in updating Qwait
  deltas = matrix(NA, nTimeStep, nTrial)
  targets = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # select actions
    t = 1
    thisScheduledWait = scheduledWait[tIdx]
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
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update action values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
      if(getReward){
        targets[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          targets[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi * nega *(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
      Viti = ifelse(nextReward > 0, Viti + phi * delta, Viti + phi * nega * delta)
      # update reRate 
      reRate = ifelse(nextReward > 0, reRate + beta * delta, reRate + beta * nega * delta)   
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Vitis[tIdx + 1] = Viti
      reRates[tIdx + 1] = reRate
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "targets" = targets, "deltas" = deltas,
    "Vitis" = Vitis, "reRates" = reRates, "condtion" = condtion
  )
  return(outputs)
}

BL = function(paras, condtion, scheduledWait){
  # parse 
  pWait = paras[1]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # action selections
    t = 1
    thisScheduledWait = scheduledWait[tIdx]
    while(t <= nTimeStep){
      # determine At
      action = ifelse(runif(1) < pWait, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
  }
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait, "condtion" = condtion
  )
  return(outputs)
}
