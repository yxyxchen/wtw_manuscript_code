# model-free analysis

# inputs:
# isTrct : logical variable determining whether the last portion in each block is truncated 

# outputs (summarised stats for each participant and each condition, 40 * 4):
# sumStats = {
  # id : [160x1 id]
  # condition : [160x1 fac]
  # nExcl : [160x1 int] # total number of excluded trials 
  # aucs : [160x1 num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
  # stdWTWs : [160x1 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
  # totalEarnings_s :  [160x1 num] 
# }
# timeWTW_ : list(160x1) # wtw timecourse, each element is a vector
# trialWTW_ : list(160x1) # trial-wise WTW, each element is a vector
# survCurve_ : list(160x1) # Kaplan-Meier survival curve, each element is a vector

MFAnalysis = function(isTrct){
  # load libraries
  source('subFxs/loadFxs.R') 
  source('subFxs/analysisFxs.R') 
  library('dplyr')
  library("tidyr")
  
  # create the output directory
  dir.create("../../genData/wtw_exp3")
  dir.create("../../genData/wtw_exp3/MFAnalysis")
  
  # load experiment parameters
  load("expParas.RData")
  nBlock = 4
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = unique(hdrData$id) 
  nSub = length(ids)                    # n
  cat('Analyzing data for',nSub,'subjects.\n')
  
################### block stats ############################
  # initialize output variables 
  nExcls = numeric(length = nSub * nBlock)
  aucs = numeric(length = nSub * nBlock) 
  stdWTWs = numeric(length = nSub * nBlock) 
  totalEarnings_s =  numeric(length = nSub * nBlock) 
  conditions = numeric(length = nSub * nBlock) 
  timeWTW_ = vector(mode = "list", length = nSub * nBlock) 
  trialWTW_ = vector(mode = "list", length = nSub * nBlock) 
  survCurve_ = vector(mode = "list", length = nSub * nBlock) 
  sub_auc_ = matrix(NA, nrow  = nSub, ncol = nBlock * 2)
  
  # loop over inidviduals
  for (sIdx in 1 : nSub) {
    # loop over blocks
    for(bkIdx in 1 : nBlock){
      # load trialData 
      id = ids[sIdx]
      thisTrialData = trialData[[id]]
      # index for elements in trialData
      noIdx = (sIdx - 1) * nBlock + bkIdx # 
      # extract (and truncate) trialData for this block
      thisTrialData = thisTrialData %>% filter(thisTrialData$blockNum == bkIdx)
      if(isTrct){
        trctLine = blockSec - max(delayMaxs)
        # truncate trials completed after tractline in each block
        nExcls[noIdx] = sum(thisTrialData$trialStartTime > trctLine)
        thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine )
      }else{
        nExcls[noIdx] = 0
      }
      
      # determine condition
      conditions[noIdx] = unique(thisTrialData$condition)
      
        
      # calcualte totalEarnings
      totalEarnings_s[noIdx] =  sum(thisTrialData$trialEarnings)
      
      # survival analysis
      kmscResults = kmsc(thisTrialData, min(delayMaxs), F, kmGrid)
      
      # 
      aucs[noIdx] = kmscResults[['auc']]
      stdWTWs[[noIdx]] = kmscResults$stdWTW
      survCurve_[[noIdx]] = kmscResults$kmOnGrid
      
      # WTW timecourse
      wtwtsResults = wtwTS(thisTrialData, tGrid, min(delayMaxs), F)
      timeWTW_[[noIdx]] = wtwtsResults$timeWTW
      trialWTW_[[noIdx]] = wtwtsResults$trialWTW
      
      # calc sub_auc
      # first half block
      sub_kmsc_res = kmsc(thisTrialData[thisTrialData$sellTime < 300,], min(delayMaxs), F, kmGrid)
      sub_auc_[sIdx, bkIdx * 2 - 1] = sub_kmsc_res$auc
      sub_kmsc_res = kmsc(thisTrialData[thisTrialData$sellTime >= 300,], min(delayMaxs), F, kmGrid)
      sub_auc_[sIdx, bkIdx * 2]  = sub_kmsc_res$auc
    }
      
  }
  # return outputs
  blockStats = data.frame(
    id = rep(ids, each = nBlock),
    blockNum = rep(1 : nBlock, nSub),
    condition = conditions,
    cbal = rep(hdrData$cbal, each = nBlock),
    nExcl = nExcls,
    totalEarnings = totalEarnings_s,
    auc = aucs,
    stdWTW = stdWTWs
  )
  outputs = list(
    blockStats = blockStats,
    survCurve_ = survCurve_,
    trialWTW_ = trialWTW_,
    timeWTW_ = timeWTW_ ,
    sub_auc_ = sub_auc_
  )
  return(outputs)
}
