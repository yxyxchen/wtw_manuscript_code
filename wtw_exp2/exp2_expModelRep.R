expModelRep = function(modelName, allData = NULL, MFResults = NULL, repOutputs = NULL){
  set.seed(123)
  
  # create output directories
  # dir.create("../../figures/wtw_exp2")
  # dir.create("../../figures/wtw_exp2/expModelRep/")
  # dir.create(sprintf("../../figures/wtw_exp2/expModelRep/%s",modelName))
  dir.create("../../genData/wtw_exp2/expModelRep")
  source("exp2_expSchematics.R")
  
  # load experiment parameters
  load('expParas.RData')
  
  # load packages and sub functions 
  library("tidyverse")
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  source("exp2_expSchematics.R")
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
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp2/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  ################ compare AUCs and CIPs from empirical and replicated data ################
  ## AUCs and CIPs from empirical data 
  source("MFAnalysis.R")
  if(is.null(MFResults)){
    MFResults = MFAnalysis(isTrct = T)
  }
  sumStats = MFResults[['sumStats']]
  muWTWEmp = sumStats$muWTW
  stdWTWEmp = sumStats$stdWTW
  ## WTW from empirical data 
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub)
  timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub)
  
  ## AUCs and CIPs from simulated data 
  if(is.null(repOutputs)){
    repOutputs =  modelRep(trialData, ids, nRep, T, modelName)
    save(repOutputs, file = sprintf("../../genData/wtw_exp2/expModelRep/%s_trct.RData", modelName))
  }
  plotData = data.frame(mu =  repOutputs$muWTWRep_mu, std = repOutputs$stdWTWRep_mu,
                        empMu = muWTWEmp, empStd = stdWTWEmp,
                        id = sumStats$ID,
                        passCheck = rep(passCheck, each = 2),
                        condition = sumStats$condition) %>% filter(passCheck)
   ################### observed stats vs model generated stats ####################
   ## plot to compare average willingess to wait
  aucFig = plotData %>%
    ggplot(aes(empMu, mu)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0)  + 
    ylab("Model-generated AUC (s)") + xlab("Observed AUC (s)") +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 16), limits = c(-1, 17)) + 
    scale_y_continuous(breaks = c(0, 16), limits = c(-1, 17)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + 
    ggtitle(sprintf("%s, N = %d", modelName, sum(passCheck)))
  
  ## plot to compare std willingess to wait
  cipFig = plotData %>%
    ggplot(aes(empStd, std, shape = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated ", sigma["WTW"], " (s"^2,")")))) +
    xlab(expression(bold(paste("Observed ", sigma["WTW"], " (s"^2,")")))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    scale_x_continuous(breaks = c(0, 7), limits = c(0, 7)) + 
    scale_y_continuous(breaks = c(0, 7), limits = c(0, 7)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") 
  rep = aucFig / cipFig 
  
  ##################### Example Participant ##############
  # prepare the data
  # repTrialData = repOutputs$repTrialData
  # repNo = repOutputs$repNo
  # sIdx = 31
  # thisTrialData = trialData[[ids[sIdx]]]
  # thisTrialData  = thisTrialData %>% filter(trialStartTime <=  blockSec - max(delayMaxs) )
  # thisTrialData = block2session(thisTrialData)
  # thisRepTrialData = repTrialData[[repNo[1, sIdx]]]
  # thisRepTrialData = data.frame(thisRepTrialData[1:6])
  # blockBoundary = min(which(thisTrialData$blockNum == 2)) # the trial boundary between two conditions
  # firstBlockMid = blockBoundary / 2
  # secondBlockMid = nrow(thisTrialData) - (nrow(thisTrialData) - blockBoundary) / 2
  # firstCondition = unique(thisTrialData$condition[thisTrialData$blockNum == 1])
  # secondCondition = unique(thisTrialData$condition[thisTrialData$blockNum == 2])
  #                         
  #                         
  #                         
  # # model generated 
  # modelFig = trialPlots(thisRepTrialData) +  
  #   ggtitle("Model-generated") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none") + 
  #   geom_vline(xintercept = min(which(thisTrialData$blockNum == 2)) - 0.5, linetype = "dashed", color = "#999999") +
  #   myTheme +
  #   annotate("text", x = firstBlockMid, y = 35, label = firstCondition, size = 6) + 
  #   annotate("text", x = secondBlockMid , y = 35, label = secondCondition, size = 6)
  # 
  # # observed
  # empFig = trialPlots(thisTrialData) + ggtitle("Observed") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none") + 
  #   geom_vline(xintercept = blockBoundary - 0.5, linetype = "dashed", color = "#999999") +
  #   myTheme +
  #   annotate("text", x = firstBlockMid, y = 35, label = firstCondition, size = 6) + 
  #   annotate("text", x = secondBlockMid , y = 35, label = secondCondition, size = 6)
  # 
  # example = modelFig / empFig +  plot_annotation(title = sprintf("%s", modelName),
  #                                      theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  # ggsave(sprintf("../../figures/wtw_exp2/expModelRep/%s/example.eps", modelName), example, width = 6, height = 8) 
  
  ################# return figure outputs ###############
  # outputs = list('rep' = rep, 'example' = example)
  outputs = list('rep' = rep)
  return(outputs)

}

