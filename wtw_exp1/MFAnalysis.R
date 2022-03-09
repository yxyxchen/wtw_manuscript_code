# model-free analysis

# inputs:
# isTrct : logical variable determining whether the last portion in each block is truncated 

# outputs:
# sumStats = {
  # id : [nSubx1 id]
  # condition : [nSubx1 fac]
  # nExcl : [nSubx1 int] # total number of excluded trials 
  # auc : [nSubx1 num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
  # stdWTW : [nSubx1 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
  # totalEarnings :  [nSubx1 num] 
# }

# blockStats = {
  # id : [nSubx3 id]
  # condition : [nSubx3 fac]
  # nExcl : [nSubx3  int] # total number of excluded trials 
  # auc : [nSubx3  num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
  # stdWTW : [nSubx3 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
  # totalEarnings :  [nSubx3 num] 
# }

# timeWTW_ : list(nSubx1) # wtw timecourse, each element is a vector
# trialWTW_ : list(nSubx1) # trial-wise WTW, each element is a vector
# survCurve_ : list(nSubx1) # Kaplan-Meier survival curve, each element is a vector

MFAnalysis = function(isTrct){
  # load libraries
  source('subFxs/loadFxs.R') 
  source('subFxs/analysisFxs.R') 
  library('dplyr')
  library("tidyr")
  
  # create the output directory
  dir.create("../../genData/wtw_exp1")
  dir.create("../../genData/wtw_exp1/MFAnalysis")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]  
  nSub = length(ids)                    # n
  cat('Analyzing data for',nSub,'subjects.\n')
  
  # calculate demographic stats 
  nFemale = sum(hdrData$gender[hdrData$stress == "no_stress"] == "Female")

  # initialize output variables 
  nExcl = numeric(length = nSub)
  auc = numeric(length = nSub) 
  sub_auc_ = matrix(NA, nrow  = nSub, ncol = nBlock * 2)
  stdWTW = numeric(length = nSub) 
  totalEarnings =  numeric(length = nSub) 
  timeWTW_ = vector(mode = "list", length = nSub) 
  trialWTW_ = vector(mode = "list", length = nSub) 
  survCurve_ = vector(mode = "list", length = nSub) 
  
  # loop over inidviduals
  for (sIdx in 1 : nSub) {
    # load (and truncate) trialData for this individual
    id = ids[sIdx]
    thisTrialData = trialData[[id]]
    # trialPlots(thisTrialData)
    if(isTrct){
      trctLine = blockSec - max(delayMaxs)
      # truncate trials completed after tractline in each block
      nExcl[sIdx] = sum(thisTrialData$trialStartTime > trctLine)
      thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine )
    }else{
      nExcl[sIdx] = 0
    }
    # calcualte totalEarnings
    totalEarnings[sIdx] =  sum(thisTrialData$trialEarnings)

    # divide data into sub blocks and calculate AUCs 
    for(k in 1 : nBlock){
      blockdata = thisTrialData[thisTrialData$blockNum == k,]
      # first half block
      sub_kmsc_res = kmsc(blockdata[blockdata$sellTime < 210,], min(delayMaxs), F, kmGrid)
      sub_auc_[sIdx, k * 2 - 1] = sub_kmsc_res$auc
      sub_kmsc_res = kmsc(blockdata[blockdata$sellTime >= 210,], min(delayMaxs), F, kmGrid)
      sub_auc_[sIdx, k * 2]  = sub_kmsc_res$auc
    }
    
    # intergrate data of 3 blocks
    thisTrialData = block2session(thisTrialData) 

    # survival analysis
    kmscResults = kmsc(thisTrialData, min(delayMaxs), F, kmGrid)
    auc[sIdx] = kmscResults[['auc']]
    survCurve_[[sIdx]] = kmscResults$kmOnGrid
    stdWTW[[sIdx]] = kmscResults$stdWTW
    
    # WTW timecourse
    wtwtsResults = wtwTS(thisTrialData, tGrid, min(delayMaxs), F)
    timeWTW_[[sIdx]] = wtwtsResults$timeWTW
    trialWTW_[[sIdx]] = wtwtsResults$trialWTW
  }
  
  sumStats = data.frame(
    id = ids,
    condition = factor(hdrData$condition[hdrData$stress == "no_stress"], levels = c("HP", "LP")),
    nExcl = nExcl,
    totalEarnings = totalEarnings,
    auc = auc,
    stdWTW = stdWTW
  )
  
  # initialize outputs on the block level
  nExcl = numeric(length = nSub * 3)
  auc = numeric(length = nSub * 3) 
  stdWTW = numeric(length = nSub * 3) 
  totalEarnings =  numeric(length = nSub * 3) 
  
  # loop over blocks 
  noIdx = 1
  for(sIdx in 1 : nSub){
    # load subject ID
    id = ids[sIdx]

    for(bkIdx in 1 : 3){
      # load trialdata
      thisTrialData = trialData[[id]]
      thisTrialData = thisTrialData[thisTrialData$blockNum == bkIdx, ]
      
      if(isTrct){
        trctLine = blockSec - max(delayMaxs)
        # truncate trials completed after tractline in each block
        nExcl[noIdx] = sum(thisTrialData$trialStartTime > trctLine)
        thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine )
      }else{
        nExcl[noIdx] = 0
      }
      
      # calcualte totalEarnings
      totalEarnings[noIdx] =  sum(thisTrialData$trialEarnings)
      
      # survival analysis
      kmscResults = kmsc(thisTrialData, min(delayMaxs), F, kmGrid)
      auc[noIdx] = kmscResults[['auc']]
      stdWTW[[noIdx]] = kmscResults$stdWTW
      
      # update noIdx
      noIdx = noIdx + 1
    }
  }
  blockStats = data.frame(
    id = rep(ids, each = 3),
    condition = rep(factor(hdrData$condition[hdrData$stress == "no_stress"], levels = c("HP", "LP")), each = 3),
    blockNum = rep(1:3, nSub),
    nExcl = nExcl,
    totalEarnings = totalEarnings,
    auc = auc,
    stdWTW = stdWTW
  )
  
  # return outputs
  outputs = list(
    sumStats = sumStats,
    blockStats = blockStats,
    survCurve_ = survCurve_,
    trialWTW_ = trialWTW_,
    timeWTW_ = timeWTW_,
    sub_auc_ = sub_auc_
  )
  return(outputs)
}
