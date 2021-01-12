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
  modelNames = c("QL1", "QL2", "RL1", "RL2", "naive", "omni")
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

  #################################################################
  ##                 figures with only RL models                 ##
  #################################################################
  modelColors = c(
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c"
  )
  # num of best explained participants 
  allPass = apply(passCheck_[,1:4], 1, sum) == 4
  bestFitNums = sapply(1 : 4, function(i) sum(apply(waic_[allPass, 1 : 4], 1, which.min) == i))
  bestFitDf = data.frame(
    modelName = factor(c("QL1", "QL2", "RL1", "RL2"), levels = c("QL1", "QL2", "RL1", "RL2")),
    bestFitNum = as.vector(bestFitNums)
  )
  # save for the output
  bestFitDf4 = bestFitDf ; bestFitDf4$modelName = as.character(bestFitDf4$modelName); bestFitDf4[nrow(bestFitDf4) + 1, ] = c("all", sum(allPass))
  fig4pie = data.frame(
    modelName = factor(c("QL1", "QL2", "RL1", "RL2"), levels = c("QL1", "QL2", "RL1", "RL2")),
    bestFitNum = as.vector(bestFitNums)
  ) %>% arrange(desc(modelName)) %>%
    mutate(prop = bestFitNum / sum(bestFitNum) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop) %>%
    ggplot(aes(x = "", y = prop, fill = modelName)) +
    geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
    scale_fill_manual(values = modelColors) + xlab("") + ylab("") +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=0, face = "bold"),
          axis.title=element_text(size= 20, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(sprintf("Exp.1, N = %d", sum(allPass))) 

  # waic
  waic_[!allPass,] = NA
  fig4WAIC = data.frame(
    modelName = c("QL1", "QL2", "RL1", "RL2", "naive", "omni"),
    mu = apply(waic_, 2, mean, na.rm = T),
    se = apply(waic_, 2, sd, na.rm = T) / sqrt(apply(waic_, 2, function(x) sum(!is.na(x))))
  ) %>% mutate(modelName = factor(modelName, levels = c("QL1", "QL2", "RL1", "RL2", "naive", "omni"))) %>%
    filter(modelName %in% c("QL1", "QL2", "RL1", "RL2")) %>% ggplot(aes(modelName, mu, fill = modelName)) + geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mu - se, ymax = mu + se), width = 0.4) +
    myTheme + xlab("") + ylab("WAIC") +
    scale_fill_manual(values = modelColors) + theme(legend.position = "none") +
    ylim(c(0, 400))
  
  ##################################################################
  ##                 figures with all six models                  ##
  ##################################################################
  # calcualted num best fit
  modelColors = c(
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c"
  )
  allPass = apply(passCheck_, 1, sum) == nModel
  bestFitNums = sapply(1 : nModel, function(i) sum(apply(waic_[allPass,], 1, which.min) == i))
  bestFitDf = data.frame(
    modelName = factor(modelNames, levels = modelNames),
    bestFitNum = as.vector(bestFitNums)
  )
  bestFitDf6 = bestFitDf ; bestFitDf6$modelName = as.character(bestFitDf6$modelName); bestFitDf6[nrow(bestFitDf6) + 1, ] = c("all", sum(allPass))
  fig6pie = data.frame(
    modelName = factor(modelNames, levels = modelNames),
    modelLabel = factor(c("QL1", "QL2", "RL1", "RL2", "naive", "omni"), levels = c("QL1", "QL2", "RL1", "RL2", "naive", "omni")),
    bestFitNum = as.vector(bestFitNums)
  ) %>% arrange(desc(modelName)) %>%
    mutate(prop = bestFitNum / sum(bestFitNum) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop) %>%
    ggplot(aes(x = "", y = prop, fill = modelName)) +
    geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
    scale_fill_manual(values = modelColors) + xlab("") + ylab("") +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=0, face = "bold"),
          axis.title=element_text(size= 20, face = "bold"))  +
    ggtitle(sprintf("Exp.1, N = %d", sum(allPass))) 
  
  # plot WAIC
  waic_[!allPass,] = NA
  fig6WAIC = data.frame(
    modelName = c("QL1", "QL2", "RL1", "RL2", "naive", "omni"),
    mu = apply(waic_, 2, mean, na.rm = T),
    se = apply(waic_, 2, sd, na.rm = T) / sqrt(apply(waic_, 2, function(x) sum(!is.na(x))))
  ) %>% mutate(modelName = factor(modelName, levels = c("QL1", "QL2", "RL1", "RL2", "naive", "omni"))) %>%
    ggplot(aes(modelName, mu, fill = modelName)) + geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mu - se, ymax = mu + se), width = 0.4) +
    myTheme + xlab("") + ylab("WAIC") +
    scale_fill_manual(values = modelColors) + theme(legend.position = "none") +
    ylim(c(0, 650))
  
  ############# return outputs #############
  outputs = list(
    "4pie" = fig4pie,
    "4waic" = fig4WAIC,
    "6pie" = fig6pie,
    "6waic" = fig6WAIC,
    "df4" = bestFitDf4,
    "df6" = bestFitDf6
  )
  return(outputs)
}


