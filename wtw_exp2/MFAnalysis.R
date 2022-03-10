# model-free analysis

# inputs:
# isTrct : logical variable determining whether the last portion in each block is truncated 

# outputs (summarised stats for each participant and each condition, 42 * 2):
# sumStats = {
  # ID : [84x1 ID]
  # condition : [84x1 fac]
  # nExcl : [84x1 int] # total number of excluded trials 
  # aucs : [84x1 num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
  # stdWTWs : [84x1 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
  # totalEarnings_s :  [84x1 num] 
# }
# timeWTW_ : list(84x1) # wtw timecourse, each element is a vector
# trialWTW_ : list(84x1) # trial-wise WTW, each element is a vector
# survCurve_ : list(84x1) # Kaplan-Meier survival curve, each element is a vector

MFAnalysis = function(isTrct){
  # load libraries
  source('subFxs/loadFxs.R') 
  source('subFxs/analysisFxs.R') 
  library('tidyverse')
  
  # create the output directory
  dir.create("../../genData/wtw_exp2")
  dir.create("../../genData/wtw_exp2/MFAnalysis")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id
  nSub = length(ids)                    # n
  cat('Analyzing data for',nSub,'subjects.\n')
  
  # initialize output variables 
  nExcls = numeric(length = nSub * nBlock)
  aucs = numeric(length = nSub * nBlock) 
  stdWTWs = numeric(length = nSub * nBlock) 
  totalEarnings_s =  numeric(length = nSub * nBlock) 
  sub_auc_ = matrix(NA, nrow  = nSub, ncol = nBlock * 2)
  timeWTW_ = vector(mode = "list", length = nSub * nBlock) 
  trialWTW_ = vector(mode = "list", length = nSub * nBlock) 
  survCurve_ = vector(mode = "list", length = nSub * nBlock) 
  print("check")
  # loop over individuals
  for (sIdx in 1 : nSub) {
    # loop over blocks
    for(bkIdx in 1 : nBlock){
      # load trialData 
      ID = ids[sIdx]
      thisTrialData = trialData[[ID]]
      # index for elements in trialData
      noidx = (sIdx - 1) * nBlock + bkIdx # 
      # extract (and truncate) trialData for this block
      thisTrialData = thisTrialData %>% filter(thisTrialData$blockNum == bkIdx)
      if(isTrct){
        trctLine = blockSec - max(delayMaxs)
        # truncate trials completed after tractline in each block
        nExcls[noidx] = sum(thisTrialData$trialStartTime > trctLine)
        thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine )
      }else{
        nExcls[noidx] = 0
      }
      
      # first half block
      sub_kmsc_res = kmsc(thisTrialData[thisTrialData$sellTime < 5 * 60,], min(delayMaxs), F, kmGrid)
      sub_auc_[sIdx, bkIdx * 2 - 1] = sub_kmsc_res$auc
      sub_kmsc_res = kmsc(thisTrialData[thisTrialData$sellTime >= 5 * 60,], min(delayMaxs), F, kmGrid)
      sub_auc_[sIdx, bkIdx * 2]  = sub_kmsc_res$auc

      
      # calcualte totalEarnings
      totalEarnings_s[noidx] =  sum(thisTrialData$trialEarnings)
      
      # survival analysis
      kmscResults = kmsc(thisTrialData, min(delayMaxs), F, kmGrid)
      
      # 
      aucs[noidx] = kmscResults[['auc']]
      survCurve_[[noidx]] = kmscResults$kmOnGrid
      stdWTWs[[noidx]] = kmscResults$stdWTW
      
      # WTW timecourse
      wtwtsResults = wtwTS(thisTrialData, tGrid, min(delayMaxs), F)
      timeWTW_[[noidx]] = wtwtsResults$timeWTW
      trialWTW_[[noidx]] = wtwtsResults$trialWTW
    }
  }
  
  # return outputs
  sumStats = data.frame(
    ID = rep(ids, each = 2),
    condition = rep(c('LP', 'HP'), nSub),
    nExcl = nExcls,
    totalEarnings = totalEarnings_s,
    auc = aucs,
    stdWTW = stdWTWs
  )
  outputs = list(
    sumStats = sumStats,
    survCurve_ = survCurve_,
    trialWTW_ = trialWTW_,
    timeWTW_ = timeWTW_,
    sub_auc_ = sub_auc_
  )
  return(outputs)
}
