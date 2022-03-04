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
  # muWTWEmp = sumStats$auc
  # stdWTWEmp = sumStats$stdWTW
  ## WTW from empirical data 
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub)
  timeWTW_ =  matrix(NA, nrow = length(tGrid), ncol = nSub)
  
  ## AUCs and CIPs from simulated data 
  if(is.null(repOutputs)){
    repOutputs =  modelRep(trialData, ids, nRep, T, modelName)
    save(repOutputs, file = sprintf("../../genData/wtw_exp2/expModelRep/%s_trct.RData", modelName))
  }
  plotData = data.frame(id = sumStats$ID, auc =  repOutputs$auc, std = repOutputs$stdWTW,
                        emp_auc = sumStats$auc, emp_std = sumStats$stdWTW,
                        passCheck = rep(passCheck, each = 2),
                        condition = sumStats$condition) %>% filter(passCheck)
  
  
  ################### observed stats vs model generated stats ####################
  ## plot to compare average willingess to wait
  aucFig = plotData %>%
    ggplot(aes(emp_auc, auc)) + 
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
  stdFig = plotData %>%
    ggplot(aes(emp_std, std, shape = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated ", sigma["WTW"], " (s"^2,")")))) +
    xlab(expression(bold(paste("Observed ", sigma["WTW"], " (s"^2,")")))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    scale_x_continuous(breaks = c(0, 7.5), limits = c(0, 7)) + 
    scale_y_continuous(breaks = c(0, 7.5), limits = c(0, 7)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") 
  
  ## plot to compare delta auc 
 delta_df = data.frame(
    id = sumStats$ID, delta =  repOutputs$delta, 
    emp_delta = sumStats[sumStats$condition == "HP", "auc"] - sumStats[sumStats$condition == "LP", "auc"],
    passCheck = passCheck
  ) 
 delta_df %>%
    ggplot(aes(emp_delta, delta)) + 
    geom_point(size = 4, stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) + 
    ylab(expression(bold(paste("Model-generated ", AUC[HP] - AUC_[LP])))) +
    xlab(expression(bold(paste("Observed", AUC[HP] - AUC_[LP]))))  +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) 
  
 ## check another 
 delta_df = data.frame(
   id = sumStats$ID, delta =  repOutputs$sub_auc_[,4] - repOutputs$sub_auc_[,3], 
   emp_delta = sub_auc_[,4] - sub_auc_[,3],
   passCheck = passCheck
 ) 
 delta_df %>%
   ggplot(aes(emp_delta, delta)) + 
   geom_point(size = 4, stroke = 1, shape = 21)  + 
   geom_abline(slope = 1, intercept = 0) + 
   ylab(expression(bold(paste("Model-generated ", AUC[HP] - AUC[LP])))) +
   xlab(expression(bold(paste("Observed", AUC[HP] - AUC[LP]))))  +
   myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) 
 
  # combine figures 
  rep = aucFig / deltaFig / stdFig 
  
  ################### calc variance explained ####################
  summary(lm(plotData$auc[plotData$condition == "HP" & plotData$passCheck] ~ plotData$emp_auc[plotData$condition == "HP" & plotData$passCheck]))$r.squared
  summary(lm(plotData$std[plotData$condition == "HP" & plotData$passCheck] ~ plotData$emp_std[plotData$condition == "HP" & plotData$passCheck]))$r.squared
  summary(lm(delta_df$delta[passCheck] ~ delta_df$emp_delta[passCheck]))$r.squared
  summary(lm(delta_df$delta[passCheck] ~ delta_df$emp_delta[passCheck]))$r.squared
  
  outputs = list('rep' = rep)
  return(outputs)

}

