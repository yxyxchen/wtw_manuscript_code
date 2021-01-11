modelRep = function(trialData, ids, nRep, isTrct, modelName){
  # 
  nSub = length(ids)
  nBlock = 2
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # normative analysis 
  iti = 2
  normResults = expSchematics(0, iti, F)
  slowPDF = normResults$slowPDF
  fastPDF = normResults$fastPDF
  time = normResults$time
  pFasts = fastPDF / (slowPDF + fastPDF) # prob of being drawn from the fast dist
  
  # initialize outputs
  repTrialData = vector(length = nSub * nRep, mode ='list')
  repNo = matrix(1 : (nSub * nRep), nrow = nRep, ncol = nSub)
  
  # loop over participants
  for(sIdx in 1 : nSub){
    # prepare empirical data 
    id = ids[sIdx]
    thisTrialData = trialData[[id]] 
    # excluded trials at the end of blocks 
    if(isTrct){
      excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(delayMaxs)))
      thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
    }
    thisTrialData = thisTrialData[thisTrialData$blockNum <= 2,]
    # scheduled rewards on non-rewards trials are not recorded and not used in parameter
    # estimation.However, given the scheduled delay, we can infer whether the scheduled 
    # reward is more likely to be 8 or -1. Since scheduled rewards on non-rewards were
    # unknown for participants and didn't influence their choices, we can use the 
    # most likely values as if they were the real scheduled rewards. 
    scheduledWait = thisTrialData$scheduledWait
    condition = thisTrialData$condition
    trialEarnings = thisTrialData$trialEarnings
    
    pFast_ = sapply(1 : length(scheduledWait),
                    function(i) pFasts[which.min(abs(time$HP - scheduledWait[i]))])
    inferredReward= sapply(1 : length(scheduledWait),
                           function(i) ifelse(runif(1) <= pFast_[i],
                                              ifelse(condition[i] == "HP", min(tokenValue), max(tokenValue)),
                                              ifelse(condition[i] == "HP", max(tokenValue), min(tokenValue))))
    scheduledReward = vector(length = length(trialEarnings))
    scheduledReward[trialEarnings == 0] = inferredReward[trialEarnings == 0]
    scheduledReward[trialEarnings != 0] = trialEarnings[trialEarnings != 0]
    # load individually fitted paramters 
    fitSummary = read.table(sprintf("../../genData/wtw_exp3/expModelFit/%s/s%s_summary.txt",  modelName, id),sep = ",", row.names = NULL)
    paras =  fitSummary[1 : nPara,1]
    # simulate nRep times
    for(rIdx in 1 : nRep){
      tempt = gnrModel(paras, condition, scheduledWait, scheduledReward, normResults)
      tempt$blockNum = thisTrialData$blockNum
      repTrialData[[repNo[rIdx, sIdx]]] = tempt
    }
  }
  
  # analysis
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub * nBlock)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub * nBlock)
  # timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub * nBlock) # block duration is not constant
  noIdx = 1
  for(sIdx in 1 : nSub){
    id = ids[sIdx]
    for(bkIdx in 1: 2){
      thisTimeWTW_ = matrix(NA, nrow = length(tGrid), ncol = nRep)
      for(rIdx in 1 : nRep){
        thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
        thisRepTrialData = do.call(cbind.data.frame, thisRepTrialData[c(names(thisRepTrialData)[1:6], "blockNum")])
        thisRepTrialData = thisRepTrialData[thisRepTrialData$blockNum == bkIdx,]
        kmscResults = kmsc(thisRepTrialData, min(delayMaxs), F, kmGrid)
        muWTWRep_[rIdx, noIdx] = kmscResults$auc
        stdWTWRep_[rIdx, noIdx] = kmscResults$stdWTW
        # wtwResults = wtwTS(thisRepTrialData, tGrid, min(delayMaxs), F)
        # thisTimeWTW_[, rIdx] = wtwResults$timeWTW
      }
      # timeWTW_[,noIdx] = apply(thisTimeWTW_, 1, function(x) mean(x, na.rm = F))
      noIdx = noIdx + 1
    }
  }
  # wtw timecourse
  ## summarise WTW across simulations for replicated data 
  muWTWRep_mu = apply(muWTWRep_, MARGIN = 2, mean) # mean of average willingness to wait
  muWTWRep_std = apply(muWTWRep_, MARGIN = 2, sd) # std of average willingess to wait
  stdWTWRep_mu = apply(stdWTWRep_, MARGIN = 2, mean) # mean of std willingness to wait
  stdWTWRep_std = apply(stdWTWRep_, MARGIN = 2, sd) # std of std willingess to wait
  
  outputs = list(
    muWTWRep_mu = muWTWRep_mu,
    muWTWRep_std = muWTWRep_std,
    stdWTWRep_mu = stdWTWRep_mu,
    stdWTWRep_std = stdWTWRep_std,
    # "timeWTW_" = timeWTW_,
    repTrialData = repTrialData,
    'repNo' = repNo
  )
  return(outputs)
}
