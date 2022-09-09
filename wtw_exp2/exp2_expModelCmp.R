expModelCmp = function(){
  # libraries and scripts
  library("ggplot2")
  library("dplyr")
  library("tidyr")
  source("subFxs/helpFxs.R")
  source("subFxs/loadFxs.R")
  source("subFxs/plotThemes.R")
  load("expParas.RData")
  
  # load data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id             
  nSub = length(ids) 
  
  # check fit
  modelNames = c("QL1", "QL2", "RL1", "RL2", "naive")
  nModel = length(modelNames)
  passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
  for(i in 1 : nModel){
    modelName = modelNames[i]
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp2/expModelFit/%s", modelName))
    passCheck_[,i] = checkFit(paraNames, expPara)
  }
  
  # extract leave-one-out results
  waic_ =  matrix(NA, nSub, nModel)
  for(m in 1 : nModel){
    modelName = modelNames[m]
    for(sIdx in 1 : nSub ){
      id = ids[sIdx]
      fileName = sprintf("../../genData/wtw_exp2/expModelFit/%s/s%s_waic.RData", modelName, id)
      load(fileName)
      waic_[sIdx, m] = WAIC$waic
    }
  }
  outputTable = cbind(waic_, passCheck_,
                      ifelse(hdrData$condition[hdrData$stress == "no_stress"] == "HP", 1, 2))
  write.table(outputTable, "../../genData/wtw_exp2/waic.csv", sep=",",  col.names=FALSE, row.names=FALSE)
  
  # calculate delta waic (using QL2 as the baseline)
  full_waic_ = waic_[apply(passCheck_, MARGIN = 1, all),] 
  full_waic_aves = round(full_waic_ %>% apply(2, mean), 2)
  full_waic_ses = round(full_waic_ %>% apply(2, function(x) sd(x) / sqrt(length(x))), 2)
  
  ################ compare all models ###########
  # calcualted num best fit
  allPass = apply(passCheck_, 1, sum) == nModel
  bestFitNums = sapply(1 : nModel, function(i) sum(apply(waic_[allPass,], 1, which.min) == i))
  bestFitDf = data.frame(
    modelName = factor(modelNames, levels = modelNames),
    bestFitNum = as.vector(bestFitNums)
  )
  bestFitDF = bestFitDf ; bestFitDF$modelName = as.character(bestFitDF$modelName); bestFitDF[nrow(bestFitDF) + 1, ] = c("all", sum(allPass))

  ############# return outputs #############
  outputs = list(
    "bestFit" = bestFitDF,
    "full_waic_aves" = full_waic_aves,
    "full_waic_ses" = full_waic_ses
  )
  
  return(outputs)
}


