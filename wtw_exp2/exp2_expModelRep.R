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
    time = rep(seq(0, blockSec * nBlock -1, by = 2), nSub),
    condition = rep(rep(c("LP", "HP"), each = length(tGrid))),
    passCheck = rep(passCheck, each = length(tGrid) * 2)
  ) %>% filter(passCheck) %>%
    gather(key = "type", value = "wtw", -c(condition, passCheck, time)) %>%
    group_by(condition, time, type) %>% 
    summarise(mu = mean(wtw),
              se = sd(wtw) / sqrt(length(wtw))) %>%
    mutate(ymin = mu - se,
          ymax = mu + se) 
  
  rectData = data.frame(
    xmin = ((1:2) - 1) * blockSec,
    xmax = (1:2) * blockSec,
    condition = c("LP", "HP")
  )
  
  figWTW = ggplot(plotdf, aes(time, mu))  +
    geom_rect(aes(xmin = xmin, xmax = xmax, fill = condition),
              data = rectData, ymin = 0, ymax = 16, alpha = 0.75, inherit.aes = F) +
    geom_line(aes(time, mu, color = type)) +
    myTheme +
    scale_color_manual(values = c("black", "#b2182b"))+
    scale_fill_manual(values =  conditionColorBacks) +
    theme(legend.position = "None") +
      scale_x_continuous(breaks = c(0, 300, 600, 900, 1200), labels = c(0, 300, 600, 900, 1200) / 60) +
      xlab("Task time (min)") + ylab("WTW (s)") + 
    scale_y_continuous(limits = c(0, 16.5), breaks = c(0, 8, 16), labels = c(0, 8, 16))
  
    
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
    ggtitle('AUC')
  
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
    theme(legend.position = "none") + ggtitle(expression(bold(sigma["WTW"])))
  
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
    theme(legend.position = "none") + ggtitle(expression(bold(AUC[end] - AUC[start]))) 

 
  # combine figures 
  figStats = (aucFig / deltaFig / stdFig)  +
    plot_annotation(title = sprintf("%s, n = %d", modelName, sum(passCheck)), theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  
  ################### calc variance explained ####################
  
  r2_df = data.frame(
    variable = rep(c("auc", "delta", "stdWTW"), each = 2),
    condition = rep(c("HP", "LP"), 3),
    r2 = round(c(
      summary(lm(plotData$auc[plotData$condition == "HP"] ~ plotData$emp_auc[plotData$condition == "HP"]))$r.squared,
      summary(lm(plotData$auc[plotData$condition == "LP"] ~ plotData$emp_auc[plotData$condition == "LP"]))$r.squared,
      summary(lm(delta_df$delta[delta_df$condition == "HP"] ~ delta_df$emp_delta[delta_df$condition == "HP"]))$r.squared,
      summary(lm(delta_df$delta[delta_df$condition == "LP"] ~ delta_df$emp_delta[delta_df$condition == "LP"]))$r.squared,
      summary(lm(plotData$std[plotData$condition == "HP"] ~ plotData$emp_std[plotData$condition == "HP"]))$r.squared,
      summary(lm(plotData$std[plotData$condition == "LP"] ~ plotData$emp_std[plotData$condition == "LP"]))$r.squared
    ) * 100, 2)
  )

  
  outputs = list('figWTW' = figWTW, 'figStats' = figStats, "r2_df" = r2_df)
  return(outputs)

}

