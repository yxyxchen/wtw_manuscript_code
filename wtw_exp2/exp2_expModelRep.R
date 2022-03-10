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
  library(patchwork)
  
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
  
  ## AUCs and CIPs from simulated data 
  if(is.null(repOutputs)){
    #repOutputs =  modelRep(trialData, ids, nRep, T, modelName)
    #save(repOutputs, file = sprintf("../../genData/wtw_exp2/expModelRep/%s_trct.RData", modelName))
    load(file = sprintf("../../genData/wtw_exp2/expModelRep/%s_trct.RData", modelName))
  }

  ################ observed WTW vs model generated WTW ############
  repTimeWTW_ = repOutputs$timeWTW_
  timeWTW_ = MFResults$timeWTW_
  plotdf = data.frame(
    rep = as.vector(repTimeWTW_),
    emp = unlist(timeWTW_),
    time = rep(seq(0, blockSec * nBlock -1, by = 1), nSub),
    condition = rep(rep(c("LP", "HP"), each = length(tGrid))),
    passCheck = rep(passCheck, each = length(tGrid) * 2)
  ) %>% filter(passCheck) %>%
    gather(key = "type", value = "wtw", -c(condition, passCheck, time)) %>%
    group_by(condition, time, type) %>% 
    summarise(mu = mean(wtw),
              se = sd(wtw) / sqrt(length(wtw))) %>%
    mutate(ymin = mu - se,
          ymax = mu + se) 
  
    # plotdf$mu[plotdf$time %in% c((blockSec - max(delayMaxs)) : blockSec, (blockSec*2 - max(delayMaxs)) : (blockSec*2))] = NA
    # plotdf$ymin[plotdf$time %in% c((blockSec - max(delayMaxs)) : blockSec, (blockSec*2 - max(delayMaxs)) : (blockSec*2))]  = NA
    # plotdf$ymax[plotdf$time %in% c((blockSec - max(delayMaxs)) : blockSec, (blockSec*2 - max(delayMaxs)) : (blockSec*2))]  = NA
  figWTW = ggplot(plotdf, aes(time, mu)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax, fill = type), color = NA, alpha = 0.5)  + 
    geom_line(aes(time, mu, color = type)) +
    myTheme +
    scale_color_manual(values = c("black", "#b2182b"))+
    scale_fill_manual(values = c("#969696", "#fa9fb5")) + 
    theme(legend.position = "None") +
      scale_x_continuous(breaks = c(0, 300, 600, 900, 1200)) + 
      xlab("Task time (s)") + ylab("WTW (s)")
    
  ################### observed stats vs model generated stats ####################
  plotData = data.frame(id = sumStats$ID, auc =  repOutputs$auc, std = repOutputs$stdWTW,
                        emp_auc = sumStats$auc, emp_std = sumStats$stdWTW,
                        passCheck = rep(passCheck, each = 2),
                        condition = sumStats$condition) %>% filter(passCheck)
  
  ## plot to compare average willingess to wait
  aucFig = plotData %>%
    ggplot(aes(emp_auc, auc)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0)  + 
    ylab("Model-generated (s)") + xlab("Observed (s)")  + 
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 16), limits = c(-1, 17)) + 
    scale_y_continuous(breaks = c(0, 16), limits = c(-1, 17)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none")  + 
    ggtitle(sprintf("%s, N = %d", modelName, sum(passCheck)))
  
  ## plot to compare std willingess to wait
  stdFig = plotData %>%
    ggplot(aes(emp_std, std, shape = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated ","(s"^2,")")))) +
    xlab(expression(bold(paste("Observed ", "(s"^2,")")))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    scale_x_continuous(breaks = c(0, 7.5), limits = c(0, 7)) + 
    scale_y_continuous(breaks = c(0, 7.5), limits = c(0, 7)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") 
  
  ## plot to compare delta auc
  delta_df = data.frame(
    id = rep(sumStats$ID[sumStats$condition == "HP"], 2), 
    delta =  c(repOutputs$sub_auc_[,2] - repOutputs$sub_auc_[,1], repOutputs$sub_auc_[,4] - repOutputs$sub_auc_[,3]),
    emp_delta = c(MFResults$sub_auc_[,2] - MFResults$sub_auc_[,1], MFResults$sub_auc_[,4] - MFResults$sub_auc_[,3]),
    condition = rep(c("LP", "HP"), each = sum(sumStats$condition == "HP")),
    passCheck = rep(passCheck, 2)
  ) %>% filter(passCheck)
  
  deltaFig = delta_df %>%
    ggplot(aes(emp_delta, delta, color = condition)) + 
    geom_point(size = 4, stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) + 
    ylab("Model-generated (s)") +
    xlab("Observed (s)")  +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    coord_fixed() + xlim(c(-8, 8)) + 
    ylim(c(-8, 8)) + scale_color_manual(values = conditionColors) +
    theme(legend.position = "none")

 
  # combine figures 
  figStats = (aucFig / deltaFig / stdFig)
  
  ################### calc variance explained ####################
  # summary(lm(plotData$auc[plotData$condition == "HP" & plotData$passCheck] ~ plotData$emp_auc[plotData$condition == "HP" & plotData$passCheck]))$r.squared
  # summary(lm(plotData$std[plotData$condition == "HP" & plotData$passCheck] ~ plotData$emp_std[plotData$condition == "HP" & plotData$passCheck]))$r.squared
  # summary(lm(delta_df$delta[passCheck] ~ delta_df$emp_delta[passCheck]))$r.squared
  # summary(lm(delta_df$delta[passCheck] ~ delta_df$emp_delta[passCheck]))$r.squared
  
  outputs = list('figWTW' = figWTW, 'figStats' = figStats)
  return(outputs)

}

