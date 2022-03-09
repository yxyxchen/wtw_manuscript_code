# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName, allData = NULL, MFResults = NULL, repOutputs = NULL){
  set.seed(123)
  # create output directories
  # dir.create("../../figures/wtw_exp1")
  # dir.create("../../figures/wtw_exp1/expModelRep/")
  # dir.create(sprintf("../../figures/wtw_exp1/expModelRep/%s",modelName))
  
  # load experiment parameters
  load('expParas.RData')
  
  # load packages and sub functions 
  library("tidyverse")
  library("latex2exp")
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library("patchwork")
  library(lattice)
  source("exp1_expSchematics.R")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  source("subFxs/modelRep.R")
  
  
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  
  # load empirical data 
  if(is.null(allData)){
    allData = loadAllData()
  }
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]
  nSub = length(ids)

  # check fit
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  # compare willingness to wait (WTW) from empirical and replicated data
  ## WTW from empirical data 
  if(is.null(MFResults)){
    source("MFAnalysis.R")
    MFResults = MFAnalysis(isTrct = T)
  }
  sumStats = MFResults[['sumStats']]
  timeWTW_ = matrix(NA, nrow = nSub, ncol = length(tGrid))
  for(i in 1 : 60){
    timeWTW_[i,] = MFResults[['timeWTW_']][[i]]
  }
  
  ## replicate data
  if(is.null(repOutputs)){
    
    #repOutputs =  modelRep(modelName, trialData, ids, nRep, T)
    #save(repOutputs, file = sprintf("../../genData/wtw_exp1/expModelRep/%s_trct.RData", modelName))
    load(file = sprintf("../../genData/wtw_exp1/expModelRep/%s_trct.RData", modelName))
  }

  #################################################################
  ##      compare observed and replicated WTW time courses    ##
  #################################################################
  repTimeWTW_ = repOutputs$timeWTW_
  timeWTW_ = MFResults$timeWTW_
  plotdf = data.frame(
    rep = as.vector(repTimeWTW_),
    emp = unlist(timeWTW_),
    time = rep(tGrid, nSub),
    condition = rep(sumStats$condition, each = length(tGrid)),
    passCheck = rep(passCheck, each = length(tGrid))
  ) %>% filter(passCheck) %>%
    gather(key = "type", value = "wtw", -c(condition, passCheck, time)) %>%
    group_by(condition, time, type) %>% 
    summarise(mu = median(wtw),
              se = sd(wtw) / sqrt(length(wtw))) %>%
    mutate(ymin = mu - se,
           ymax = mu + se) %>% ungroup()
  
    # 
    figWTW = ggplot(plotdf, aes(time, mu, color = type)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax, fill = type), color = NA, alpha = 0.5)  + 
    geom_line(aes(time, mu, color = type)) +
      facet_grid(~condition) + myTheme +
    scale_color_manual(values = c("black", "#b2182b"))+
    scale_fill_manual(values = c("#969696", "#fa9fb5")) + 
    theme(legend.position = "None") +
    scale_x_continuous(breaks = 0:3 * 60 * 7, labels = 0:3 * 7) + 
    xlab("Task time (min)") + ylab("WTW (s)")
  
  #################################################################
  ##      compare observed and replicated AUC and sigma_wtw      ##
  #################################################################
  # plot to compare delta wtw 
  delta_df = data.frame(
    delta =  repOutputs$sub_auc_[,6] - repOutputs$sub_auc_[,1],
    emp_delta = MFResults$sub_auc_[,6] - MFResults$sub_auc_[,1],
    condition = sumStats$condition,
    passCheck = passCheck
  ) %>% filter(passCheck)
  
  deltaFig = delta_df %>%
    ggplot(aes(emp_delta, delta, color = condition)) + 
    geom_point(size = 4, stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) + 
    ylab("Model-generated (s)") +
    xlab("Observed (s)")  +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    coord_fixed() + xlim(c(-20, 8)) + 
    ylim(c(-20, 8)) + scale_color_manual(values = conditionColors) +
    theme(legend.position = "none")
    
  # plot to compare average willingess to wait
  plotData = data.frame(id = ids, mu =  repOutputs$auc, std = repOutputs$stdWTW,
                        empMu = sumStats$auc, empStd = sumStats$stdWTW,
                        passCheck = passCheck, 
                        condition = sumStats$condition) %>% filter(passCheck)
  
  aucFig = plotData %>%
  ggplot(aes(empMu, mu)) + 
  geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
  geom_abline(slope = 1, intercept = 0)  + 
  ylab("Model-generated (s)") + xlab("Observed (s)") +
  myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(breaks = c(0, 20), limits = c(-1, 21)) + 
  scale_y_continuous(breaks = c(0, 20), limits = c(-1, 21)) +
  scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + 
    ggtitle("AUC")
  summary(lm(empMu~mu, plotData[plotData$passCheck & plotData$condition == "HP",]))$r.squared
  summary(lm(empMu~mu, plotData[plotData$passCheck & plotData$condition == "LP",]))$r.squared

  # plot to compare std willingess to wait
  stdFig = plotData %>%
  ggplot(aes(empStd, std, shape = condition)) + 
  geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
  geom_abline(slope = 1, intercept = 0) +
  ylab(expression(bold(paste("Model-generated", " (s"^2,")")))) +
    xlab(expression(bold(paste("Observed", " (s"^2,")")))) +
  myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_x_continuous(breaks = c(0, 10), limits = c(0, 10)) + 
  scale_y_continuous(breaks = c(0, 10), limits = c(0, 10)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") +
    ggtitle(expression(bold(sigma["WTW"])))

  # combind
  figStats = aucFig / deltaFig / stdFig +
    plot_annotation(title = sprintf("%s, n = %d", modelName, sum(passCheck)), theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  
  ################# return figure outputs ###############
  # outputs = list("rep" = rep, "example" = example)
  outputs = list('figWTW' = figWTW, 'figStats' = figStats, "rep_wtw_df" = plotdf)
  return(outputs)
  
}

