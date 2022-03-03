modelRep = function(modelName, trialData, ids, nRep, isTrct, aveParas = NULL){
  # 
  nSub = length(ids)
  load("expParas.RData")
  
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # normative analysis
  iti = 2
  normResults = expSchematics(0, iti, F)
  
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
    # load individually fitted paramters 
    if(is.null(aveParas)){
      fitSummary = read.table(sprintf("../../genData/wtw_exp1/expModelFit/%s/s%s_summary.txt",  modelName, id),sep = ",", row.names = NULL)
      paras =  fitSummary[1 : nPara,1]
    }else{
      paras = aveParas
    }
    # simulate nRep times
    for(rIdx in 1 : nRep){
      tempt = gnrModel(paras, thisTrialData$condition, thisTrialData$scheduledWait, normResults)
      repTrialData[[repNo[rIdx, sIdx]]] = tempt
    }
  }
  # initialize 
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub)
  timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub)
  trialWTW_ = list()
  sub_auc_ = matrix(NA, nrow  = nSub, ncol = 6)
  for(sIdx in 1 : nSub){
    id = ids[sIdx]
    thisTrialData = trialData[[id]]
    this_sub_auc_ = matrix(NA, nrow  = nRep, ncol = 6)
    if(isTrct){
      excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(delayMaxs)))
      thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
    }
    thisTrialData = block2session(thisTrialData)
    nTrial =  length(thisTrialData$scheduledWait)
    # timeWTW = matrix(NA, nrow = length(tGrid), ncol = nRep)
    trialWTW =  matrix(NA, nrow = nTrial, ncol = nRep)
    for(rIdx in 1 : nRep){
      thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
      thisRepTrialData$Qwaits_ = NULL
      thisRepTrialData = data.frame(thisRepTrialData)
      kmscResults = kmsc(thisRepTrialData, min(delayMaxs), F, kmGrid)
      muWTWRep_[rIdx,sIdx] = kmscResults$auc
      stdWTWRep_[rIdx, sIdx] = kmscResults$stdWTW
      wtwResults = wtwTS(thisRepTrialData, tGrid, min(delayMaxs), F)
      # timeWTW[,rIdx] = wtwResults$timeWTW
      trialWTW[, rIdx] = wtwResults$trialWTW
      for(k in 1 : 3){
        # first half block
        sub_kmsc_res = kmsc(thisRepTrialData[thisTrialData$sellTime < k * 7 * 60 - 210 & thisTrialData$sellTime >= k * 7 * 60 - 420,], min(delayMaxs), F, kmGrid)
        this_sub_auc_[rIdx, k * 2 - 1] = sub_kmsc_res$auc
        sub_kmsc_res = kmsc(thisRepTrialData[thisTrialData$sellTime >= k * 7 * 60 - 210 & thisTrialData$sellTime < k * 7 * 60,], min(delayMaxs), F, kmGrid)
        this_sub_auc_[rIdx, k * 2]  = sub_kmsc_res$auc
      }
    }

    trialWTW_[[id]] =  apply(trialWTW, 1, mean)
    timeWTW_[,sIdx] =  resample(trialWTW_[[id]], thisTrialData$sellTime, tGrid)
    sub_auc_[sIdx, ] = apply(this_sub_auc_, FUN = mean, MARGIN = 2)
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
    timeWTW_ = timeWTW_,
    trialWTW_ = trialWTW_,
    repTrialData = repTrialData,
    repNo = repNo,
    sub_auc_ = sub_auc_
  )
  return(outputs)
}
