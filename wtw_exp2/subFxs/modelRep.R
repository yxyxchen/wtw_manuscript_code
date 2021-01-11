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
  
  # initialize 
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub * 2)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub * 2)
  # timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub * 2)
  noIdx = 1
  for(sIdx in 1 : nSub){
    ID = ids[sIdx]
    for(bkIdx in 1 : 2){
      condition = ifelse(bkIdx== 1, "LP", "HP")
      # timeWTW = matrix(NA, nrow = length(tGrid), ncol = nRep)
      for(rIdx in 1 : nRep){
        thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
        thisRepTrialData = do.call(cbind.data.frame, thisRepTrialData[1:6])
        thisRepTrialData = thisRepTrialData[thisRepTrialData$condition == condition,]
        kmscResults = kmsc(thisRepTrialData, min(delayMaxs), F, kmGrid)
        muWTWRep_[rIdx, noIdx] = kmscResults$auc
        stdWTWRep_[rIdx, noIdx] = kmscResults$stdWTW
        # wtwResults = wtwTS(thisRepTrialData, tGrid, min(delayMaxs), F)
        # timeWTW[,rIdx] = wtwResults$timeWTW
      }
      # timeWTW_[,noIdx] = apply(timeWTW, 1, mean, na.rm = T)
      noIdx = noIdx + 1
    }
  }
  
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
