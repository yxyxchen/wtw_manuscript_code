# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName, allData = NULL, MFResults = NULL, repOutputs = NULL){
  set.seed(123)
  # create output directories
  # dir.create("../../figures/wtw_exp3")
  # dir.create("../../figures/wtw_exp3/expModelRep/")
  dir.create("../../genData/wtw_exp3/expModelRep")
  # dir.create(sprintf("../../figures/wtw_exp3/expModelRep/%s",modelName))
  
  # load experiment parameters
  load('expParas.RData')
  nBlock = 2
  
  # load packages and sub functions 
  library("tidyverse")
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  library("patchwork")
  source("exp3_expSchematics.R")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  source("subFxs/modelRep.R")
  
  # num of repetitions 
  nRep = 10
  
  # load empirical data 
  if(is.null(allData)){
    allData = loadAllData()
  }
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id
  nSub = length(ids)

  # check fit
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp3/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  # empirical data 
  if(is.null(MFResults)){
    source("MFAnalysis.R")
    MFResults = MFAnalysis(isTrct = T)
  }
  blockStats = MFResults[['blockStats']]
  blockStats = blockStats[blockStats$blockNum <= 2, ]
  aucEmp = blockStats$auc
  stdWTWEmp = blockStats$stdWTW
  
  # replicated data 
  if(is.null(repOutputs)){
    repOutputs =  modelRep(trialData, ids, nRep, T, modelName)
    save(repOutputs, file = sprintf("../../genData/wtw_exp3/expModelRep/%s_trct.RData", modelName))
  }

  plotData = data.frame(id = blockStats$id, mu =  repOutputs$aucRep_mu, std = repOutputs$stdWTWRep_mu,
                        empMu = aucEmp, empStd = stdWTWEmp,
                        passCheck = rep(passCheck, each = 2), 
                        condition = blockStats$condition) %>% filter(passCheck)
  
  ##########################
  ##      calc RSME      ##
  ##########################
  mu_sqerr = (plotData$empMu - plotData$mu)^2
  std_sqerr = (plotData$empStd - plotData$std)^2
  sqerr_df = data.frame(
    id = plotData$id, 
    condition = plotData$condition,
    mu_sqerr = mu_sqerr,
    std_sqerr = std_sqerr)
  
  #################################################################
  ##      compare observed and replicated AUC and sigma_wtw      ##
  #################################################################
  ## plot to compare average willingess to wait
  aucFig = plotData %>%
    ggplot(aes(empMu, mu)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0)  + 
    ylab("Model-generated AUC (s)") + xlab("Observed AUC (s)") +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 30), limits = c(-1, 31)) + 
    scale_y_continuous(breaks = c(0, 30), limits = c(-1, 31)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + ggtitle(sprintf("%s, N = %d", modelName, sum(passCheck)))
  
  ## plot to compare std willingess to wait
  cipFig = plotData %>%
    ggplot(aes(empStd, std, shape = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated ", sigma["WTW"], " (s"^2,")")))) +
    xlab(expression(bold(paste("Observed ", sigma["WTW"], " (s"^2,")")))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    scale_x_continuous(breaks = c(0, 15), limits = c(0, 15)) + 
    scale_y_continuous(breaks = c(0, 15), limits = c(0, 15)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none")
  
  # rep = grid.arrange(aucFig, cipFig, nrow = 2, top = textGrob(sprintf("%s, n = %d", modelName, sum(passCheck)), gp = gpar(fontsize = 20,font = 1)))
  rep = aucFig / cipFig +  plot_annotation(title = sprintf("%s, n = %d", modelName, sum(passCheck)),
                                     theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  
  #################################################################
  ##                     example participant                     ##
  #################################################################
  # repTrialData = repOutputs$repTrialData
  # repNo = repOutputs$repNo
  # sIdx = 3
  # thisTrialData = trialData[[ids[sIdx]]]
  # thisTrialData  = thisTrialData %>% filter(trialStartTime <=  blockSec - max(delayMaxs) & thisTrialData$blockNum <= 2)
  # thisTrialData = block2session(thisTrialData)
  # thisRepTrialData = repTrialData[[repNo[1, sIdx]]]
  # thisRepTrialData = data.frame(thisRepTrialData[1:6])
  # blockBoundary = min(which(thisTrialData$blockNum == 2)) # the trial boundary between two conditions
  # firstBlockMid = blockBoundary / 2
  # secondBlockMid = nrow(thisTrialData) - (nrow(thisTrialData) - blockBoundary) / 2
  # firstCondition = unique(thisTrialData$condition[thisTrialData$blockNum == 1])
  # secondCondition = unique(thisTrialData$condition[thisTrialData$blockNum == 2])
  # 
  # # model generated 
  # modelFig = trialPlots(thisRepTrialData) +  
  #   ggtitle("Model-generated") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none") + 
  #   geom_vline(xintercept = min(which(thisTrialData$blockNum == 2)) - 0.5, linetype = "dashed", color = "#999999") +
  #   myTheme +
  #   annotate("text", x = firstBlockMid, y = 22, label = firstCondition, size = 6) + 
  #   annotate("text", x = secondBlockMid , y = 22, label = secondCondition, size = 6)
  # 
  # # observed
  # empFig = trialPlots(thisTrialData) + ggtitle("Observed") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none") + 
  #   geom_vline(xintercept = blockBoundary - 0.5, linetype = "dashed", color = "#999999") +
  #   myTheme +
  #   annotate("text", x = firstBlockMid, y = 22, label = firstCondition, size = 6) + 
  #   annotate("text", x = secondBlockMid , y = 22, label = secondCondition, size = 6)
  #   
  # # example = grid.arrange(modelFig, empFig, nrow = 2, top = textGrob(sprintf("%s", modelName), gp = gpar(fontsize = 20, font = 2)))
  # example = modelFig / empFig
  # ggsave(sprintf("../../figures/wtw_exp3/expModelRep/%s/example.eps", modelName), example, width = 6, height = 8) 
  
  ################# return figure outputs ###############
  # outputs = list("rep" = rep, "example" = example)
  outputs = list("rep" = rep, "sqerr_df" = sqerr_df)
  return(outputs)
}

