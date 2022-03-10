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
  aucs = vector(length = nSub * 2) # auc in both conditions: 
  stdWTWs = vector(length = nSub * 2) # sigma_wtw in both conditions
  sub_auc_ = matrix(NA, nrow  = nSub, ncol = 4) # auc in each half block for each participant 
  timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub * 2) # wtw timecourse for each subj and each condition 
  trialWTW_ = list() # wtw timecourse for each subj and each condition 

  for(sIdx in 1 : nSub){
    ID = ids[sIdx]
    thisTrialData = trialData[[ID]]
    # excluded trials at the end of blocks 
    if(isTrct){
      excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(delayMaxs)))
      if(length(excluedTrials) > 0){
        thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
      }
    }
    # calc for each block 
    for(bkIdx in 1 : 2){
      # select data for this block
      blockData = thisTrialData[thisTrialData$blockNum == bkIdx,]
      condition = blockData$condition[1]
      # initialize data outputs 
      auc_reps = vector(length = nRep)
      stdWTW_reps = vector(length = nRep)
      first_half_auc_reps = vector(length = nRep)
      second_half_auc_reps = vector(length = nRep)
      nTrial = nrow(blockData)
      trialWTW_reps = matrix(NA, nrow = nTrial, ncol = nRep)
      # loop over reps 
      for(rIdx in 1 : nRep){
        # extract each replicated blockdata
        thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
        thisRepTrialData = do.call(cbind.data.frame, thisRepTrialData[1:6])
        repBlockData = thisRepTrialData[thisRepTrialData$condition == condition,]
        # 
        kmscResults = kmsc(repBlockData, min(delayMaxs), F, kmGrid)
        auc_reps[rIdx] = kmscResults$auc
        stdWTW_reps[rIdx] = kmscResults$stdWTW
        # calc auc for the first and second halves 
        kmscResults = kmsc(repBlockData[blockData$sellTime < 300 & blockData$sellTime >= 0,], min(delayMaxs), F, kmGrid)
        first_half_auc_reps[rIdx] = kmscResults$auc
        kmscResults = kmsc(repBlockData[blockData$sellTime < 600 & blockData$sellTime >= 300,], min(delayMaxs), F, kmGrid)
        second_half_auc_reps[rIdx] = kmscResults$auc
        # calc WTW
        wtwResults = wtwTS(repBlockData, tGrid, min(delayMaxs), F)
        trialWTW_reps[, rIdx] = wtwResults$trialWTW
      }
      # average across reps 
      aucs[(sIdx-1) * 2 + bkIdx] = mean(auc_reps)
      stdWTWs[(sIdx-1) * 2 + bkIdx] = mean(stdWTW_reps)
      sub_auc_[sIdx,  bkIdx * 2 - 1] =  mean(first_half_auc_reps)
      sub_auc_[sIdx,  bkIdx * 2] =  mean(second_half_auc_reps)
      trialWTW_[[(sIdx-1) * 2 + bkIdx]] = apply(trialWTW_reps, 1, mean)
      timeWTW_[,(sIdx-1) * 2 + bkIdx] =  resample(trialWTW_[[(sIdx-1) * 2 + bkIdx]], blockData$sellTime, tGrid)
    }
    # deltas[sIdx] = aucs[sIdx * 2]  - aucs[(sIdx-1) * 2 + 1] 
  }  
  
  # aucRep_ = matrix(NA, nrow = nRep , ncol = nSub * nBlock)
  # stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub * nBlock)
  # # timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub * nBlock) # block duration is not constant
  # noIdx = 1
  # for(sIdx in 1 : nSub){
  #   id = ids[sIdx]
  #   for(bkIdx in 1: 2){
  #     thisTimeWTW_ = matrix(NA, nrow = length(tGrid), ncol = nRep)
  #     for(rIdx in 1 : nRep){
  #       thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
  #       thisRepTrialData = do.call(cbind.data.frame, thisRepTrialData[c(names(thisRepTrialData)[1:6], "blockNum")])
  #       thisRepTrialData = thisRepTrialData[thisRepTrialData$blockNum == bkIdx,]
  #       kmscResults = kmsc(thisRepTrialData, min(delayMaxs), F, kmGrid)
  #       aucRep_[rIdx, noIdx] = kmscResults$auc
  #       stdWTWRep_[rIdx, noIdx] = kmscResults$stdWTW
  #       # wtwResults = wtwTS(thisRepTrialData, tGrid, min(delayMaxs), F)
  #       # thisTimeWTW_[, rIdx] = wtwResults$timeWTW
  #     }
  #     # timeWTW_[,noIdx] = apply(thisTimeWTW_, 1, function(x) mean(x, na.rm = F))
  #     noIdx = noIdx + 1
  #   }
  # }
  # # wtw timecourse
  # ## summarise WTW across simulations for replicated data 
  # aucRep_mu = apply(aucRep_, MARGIN = 2, mean) # mean of average willingness to wait
  # aucRep_std = apply(aucRep_, MARGIN = 2, sd) # std of average willingess to wait
  # stdWTWRep_mu = apply(stdWTWRep_, MARGIN = 2, mean) # mean of std willingness to wait
  # stdWTWRep_std = apply(stdWTWRep_, MARGIN = 2, sd) # std of std willingess to wait
  
  # return 
  outputs = list(
    auc = aucs,
    stdWTW = stdWTWs,
    sub_auc_ = sub_auc_,
    timeWTW_ =  timeWTW_,
    trialWTW_ = trialWTW_,
    repTrialData = repTrialData,
    'repNo' = repNo
  )
  return(outputs)
}
