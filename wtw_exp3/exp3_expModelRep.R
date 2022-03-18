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
  timeWTW_ = MFResults$timeWTW_[blockStats$blockNum <= 2]
  blockStats = blockStats[blockStats$blockNum <= 2, ]
  
  
  # replicated data 
  if(is.null(repOutputs)){
    #repOutputs =  modelRep(trialData, ids, nRep, T, modelName)
    #save(repOutputs, file = sprintf("../../genData/wtw_exp3/expModelRep/%s_trct.RData", modelName))
    load(file = sprintf("../../genData/wtw_exp3/expModelRep/%s_trct.RData", modelName))
  }

  ################ observed WTW vs model generated WTW ############
  repTimeWTW_ = repOutputs$timeWTW_
  plotdf = data.frame(
    rep = as.vector(repTimeWTW_),
    emp = unlist(timeWTW_),
    time = rep(seq(0, blockSec * 2 -1, by = 2) / 60, nSub),
    condition = rep(blockStats$condition, each = length(tGrid)),
    passCheck = rep(passCheck, each = length(tGrid) * 2),
    cbal = rep(ifelse(blockStats$cbal == 1, "HP First", "LP First"), each = length(tGrid))
  ) %>% filter(passCheck) %>%
    gather(key = "type", value = "wtw", -c(condition, passCheck, time, cbal)) %>%
    group_by(condition, time, type, cbal) %>% 
    summarise(mu = mean(wtw),
              se = sd(wtw) / sqrt(length(wtw))) %>%
    mutate(ymin = mu - se,
           ymax = mu + se) 
  
  rectData = data.frame(
    xmin = rep(((1:2) - 1) * blockSec, 2) / 60,
    xmax = rep(1:2 * blockSec, 2) / 60,
    condition = c("HP", "LP", "LP", "HP"),
    cbal = c("HP First", "HP First", "LP First", "LP First")
  )
  figWTW = ggplot(plotdf, aes(time, mu)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, fill = condition),
              data = rectData, ymin = 0, ymax = 20, alpha = 0.75, inherit.aes = F) +
    facet_grid(~cbal) +
    geom_line(aes(time, mu, color = type))  +
    myTheme +
    scale_color_manual(values = c("black", "#b2182b"))+
    scale_fill_manual(values = c(conditionColorBacks))  +
    scale_x_continuous(breaks = c(0, 10, 20)) +
    xlab("Task time (s)") + ylab("WTW (s)") + ylim(c(0, 20)) +
    theme(legend.position = "None") + ylim(c(0, 20))

  # figWTW = ggplot(plotdf, aes(time, mu)) + 
  #   facet_grid(~cbal) + 
  #   geom_line(aes(time, mu, color = condition, linetype = type))  +
  #   myTheme +
  #   scale_color_manual(values = conditionColors)+
  #   scale_x_continuous(breaks = c(0, 10, 20)) + 
  #   xlab("Task time (s)") + ylab("WTW (s)") + ylim(c(0, 20)) +
  #   theme(legend.position = "None")  
  #################################################################
  ##      compare observed and replicated AUC and sigma_wtw      ##
  #################################################################
  ## plot to compare average willingess to wait
  plotData = data.frame(
    emp_auc = blockStats$auc,
    auc = repOutputs$auc,
    emp_std = blockStats$stdWTW,
    std = repOutputs$stdWTW,
    condition = blockStats$condition,
    cbal = blockStats$cbal,
    passCheck = rep(passCheck, each = 2)
  ) %>% filter(passCheck)
  
  aucFig = plotData %>% 
    ggplot(aes(emp_auc, auc)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0)  + 
    ylab("Model-generated (s)") + xlab("Observed (s)") +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 30), limits = c(-1, 31)) + 
    scale_y_continuous(breaks = c(0, 30), limits = c(-1, 31)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + ggtitle("AUC") 
  
  ## plot to compare std willingess to wait
  stdFig = plotData %>%
    ggplot(aes(emp_std, std, shape = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated"," (s"^2,")")))) +
    xlab(expression(bold(paste("Model-generated"," (s"^2,")"))))+
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    scale_x_continuous(breaks = c(0, 15), limits = c(0, 15)) + 
    scale_y_continuous(breaks = c(0, 15), limits = c(0, 15)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + 
    ggtitle(expression(bold(sigma[WTW])))
    
  delta_df = data.frame(
    delta = as.vector(repOutputs$sub_auc_[, c(4,2)] -repOutputs$sub_auc_[, c(3,1)]),
    emp_delta = as.vector(MFResults$sub_auc_[, c(4,2)] - MFResults$sub_auc_[, c(3,1)]),
    condition = blockStats$condition,
    cbal = blockStats$cbal,
    passCheck = rep(passCheck, each = 2)
  ) %>% filter(passCheck)
  
 deltaFig = delta_df %>% ggplot(aes(emp_delta, delta)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) + 
    xlim(c(-15, 15)) + ylim(c(-15, 15)) + 
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + 
    ggtitle(expression(bold(AUC[end] - AUC[start]))) +
    ylab("Model-generated (s)") + xlab("Observed (s)") + myTheme
   
  #################### calc % explained ######## 
  summary(lm(emp_auc~auc, plotData[plotData$condition == "HP",]))$r.squared
  summary(lm(emp_auc~auc, plotData[plotData$condition == "LP",]))$r.squared
  summary(lm(emp_std~std, plotData[plotData$condition == "HP",]))$r.squared
  summary(lm(emp_std~std, plotData[plotData$condition == "LP",]))$r.squared
  summary(lm(emp_delta~delta, delta_df[delta_df$condition == "HP",]))$r.squared
  summary(lm(emp_delta~delta, delta_df[delta_df$condition == "LP",]))$r.squared
  
  # combind
  figStats = aucFig / deltaFig / stdFig +
    plot_annotation(title = sprintf("%s, n = %d", modelName, sum(passCheck)), theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  
  
  ################# return figure outputs ###############
  # outputs = list("rep" = rep, "example" = example)
  outputs = list('figWTW' = figWTW, 'figStats' = figStats)
  return(outputs)
}

