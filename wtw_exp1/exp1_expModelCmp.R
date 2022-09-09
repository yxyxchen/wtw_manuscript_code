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
  ids = hdrData$id[hdrData$stress == "no_stress"]                 
  nSub = length(ids) 
  
  # check fit
  modelNames = c("QL1", "QL2", "RL1", "RL2", "naive")
  nModel = length(modelNames)
  passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
  for(i in 1 : nModel){
    modelName = modelNames[i]
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", modelName))
    passCheck_[,i] = checkFit(paraNames, expPara)
  }
  
  # extract leave-one-out results
  waic_ =  matrix(NA, nSub, nModel)
  for(m in 1 : nModel){
    modelName = modelNames[m]
    for(sIdx in 1 : nSub ){
      id = ids[sIdx]
      fileName = sprintf("../../genData/wtw_exp1/expModelFit/%s/s%s_waic.RData", modelName, id)
      load(fileName)
      waic_[sIdx, m] = WAIC$waic
    }
  }
  outputTable = cbind(waic_, passCheck_,
                      ifelse(hdrData$condition[hdrData$stress == "no_stress"] == "HP", 1, 2))
  write.table(outputTable, "../../genData/wtw_exp1/waic.csv", sep=",",  col.names=FALSE, row.names=FALSE)

  # delta_waic_ave_6 = round(waic_ave_6 - waic_ave_6[2], 2)
  full_waic_ = waic_[apply(passCheck_, MARGIN = 1, all),] 
  full_waic_aves = round(full_waic_ %>% apply(2, mean), 2)
  full_waic_ses = round(full_waic_ %>% apply(2, function(x) sd(x) / sqrt(length(x))), 2)
  
  allPass = apply(passCheck_, 1, sum) == nModel
  bestFitNums = sapply(1 : nModel, function(i) sum(apply(waic_[allPass,], 1, which.min) == i))
  bestFitDf = data.frame(
    modelName = factor(modelNames, levels = modelNames),
    bestFitNum = as.vector(bestFitNums)
  )
  bestFitDf$modelName = as.character(bestFitDf$modelName); bestFitDf[nrow(bestFitDf) + 1, ] = c("all", sum(allPass))
  
  # reduced_waic_ = waic_[apply(passCheck_[,1:4], MARGIN = 1, all), 1 : 4] 
  # reduced_waic_aves = reduced_waic_ %>% apply(2, mean)
  # reduced_waic_ses = reduced_waic_ %>% apply(2, function(x) sd(x) / sqrt(length(x)))
  
  #################################################################
  ##                 figures with only RL models                 ##
  #################################################################
  # modelColors = c(
  #   "#a6cee3",
  #   "#1f78b4",
  #   "#b2df8a",
  #   "#33a02c"
  # )
  # num of best explained participants 
  # allPass = apply(passCheck_[,1:4], 1, sum) == 4
  # bestFitNums = sapply(1 : 4, function(i) sum(apply(waic_[allPass, 1 : 4], 1, which.min) == i))
  # bestFitDf = data.frame(
  #   modelName = factor(c("QL1", "QL2", "RL1", "RL2"), levels = c("QL1", "QL2", "RL1", "RL2")),
  #   bestFitNum = as.vector(bestFitNums)
  # )
  # # save for the output
  # bestFitDf4 = bestFitDf ; bestFitDf4$modelName = as.character(bestFitDf4$modelName); bestFitDf4[nrow(bestFitDf4) + 1, ] = c("all", sum(allPass))

  ##################################################################
  ##                 figures with all six models                  ##
  ##################################################################
  # calcualted num best fit
  # modelColors = c(
  #   "#a6cee3",
  #   "#1f78b4",
  #   "#b2df8a",
  #   "#33a02c",
  #   "#fb9a99",
  #   "#e31a1c"
  # )

  
  
  ############# return outputs #############
  outputs = list(
    "bestFit" = bestFitDf,
    "full_waic_aves" = full_waic_aves,
    "full_waic_ses" = full_waic_ses
  )
  return(outputs)
}


