modelRep = function(trialData, ids, nRep, isTrct, modelName){
  # 
  nSub = length(ids)
  iti = 2
  load("expParas.RData")
  source("exp2_expSchematics.R")
  
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # normative analysis 
  normResults = expSchematics(0, iti, F)
  
  # initialize outputs
  repTrialData = vector(length = nSub * nRep, mode ='list')
  repNo = matrix(1 : (nSub * nRep), nrow = nRep, ncol = nSub)
  
  # loop over participants
  for(sIdx in 1 : nSub){
    # prepare empirical data 
    ID = ids[sIdx]
    thisTrialData = trialData[[ID]]
    # excluded trials at the end of blocks 
    if(isTrct){
      excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(delayMaxs)))
      if(length(excluedTrials) > 0){
        thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
      }
    }
    # load individually fitted paramters 
    fitSummary = read.table(sprintf("../../genData/wtw_exp2/expModelFit/%s/s%s_summary.txt",  modelName, ID),sep = ",", row.names = NULL)
    paras =  fitSummary[1 : nPara,1]
    # simulate nRep times
    for(rIdx in 1 : nRep){
      tempt = gnrModel(paras, thisTrialData$condition, thisTrialData$scheduledWait, normResults)
      repTrialData[[repNo[rIdx, sIdx]]] = tempt
    }
  }
  
  # initialize outputs
  aucs = vector(length = nSub * 2) # auc in both conditions: 
  stdWTWs = vector(length = nSub * 2) # sigma_wtw in both conditions
  deltas = vector(length = nSub) # differences in auc_HP and auc_LP
  sub_auc_ = matrix(NA, nrow  = nSub, ncol = 4) # auc in each half block for each participant 
  timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub * 2) # wtw timecourse for each subj and each condition 
  trialWTW_ = list() # wtw timecourse for each subj and each condition 
  
  # timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub * 2)
  # noIdx = 1
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
    deltas[sIdx] = aucs[sIdx * 2]  - aucs[(sIdx-1) * 2 + 1] 
  }
  # return 
  outputs = list(
    auc = aucs,
    stdWTW = stdWTWs,
    delta = deltas,
    sub_auc_ = sub_auc_,
    timeWTW_ =  timeWTW_,
    trialWTW_ = trialWTW_,
    repTrialData = repTrialData,
    'repNo' = repNo
  )
  return(outputs)
}
