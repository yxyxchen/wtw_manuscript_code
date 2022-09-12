expModelRepInd = function(){
  # load libraries
  source('subFxs/loadFxs.R') 
  source('subFxs/analysisFxs.R') 
  source("subFxs/plotThemes.R")
  library('dplyr')
  library("tidyr")
  library("ggplot2")
  source("subFxs/helpFxs.R")
  source("subFxs/modelRep.R")
  source("exp1_expSchematics.R")
  
  # create the output directory
  dir.create("../../genData/wtw_exp1")
  dir.create("../../genData/wtw_exp1/MFAnalysis")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]  
  nSub = length(ids)                    # n
  cat('Analyzing data for',nSub,'subjects.\n')
  
  # conduct wtw analysis for several example participants
  cut =1260
  selectedIds = c("11","20", "79")
  outputs_ = list()
  sellTimes_ = list()
  i = 1
  for (id in selectedIds) {
    thisTrialData = trialData[[id]] 
    if(unique(thisTrialData$condition) == "LP"){
      trctLine = blockSec - max(delayMaxs)
      # truncate trials completed after tractline in each block
      thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine) 
      thisTrialData = block2session(thisTrialData)
      # thisTrialData = thisTrialData[thisTrialData$trialStartTime < cut,]
      # p = trialPlots(thisTrialData) + ylim(c(0, 40))
      outputs = wtwTS(thisTrialData, tGrid[1:cut], 20, plotWTW = F) 
      # outputs_[[i]] = outputs$timeWTW
      outputs_[[i]] = outputs$trialWTW
      sellTimes_[[i]] = thisTrialData$sellTime
      i = i + 1
    }
  }
  
  
  source("MFAnalysis.R")
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults$sumStats
  sumStats$delta = MFResults$sub_auc_[,6] - MFResults$sub_auc_[,1]
  ############### plot empirical results of three example participants #######
  figEmp = data.frame(
    wtw = unlist(outputs_),
    time = unlist(sellTimes_),
    sub = rep(selectedIds, sapply(sellTimes_, length))
  ) %>% ggplot(aes(time, wtw,  color = sub)) + geom_line(size = 1) + geom_point() +
    scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
    myTheme + scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
    xlab("Task time (min)") + ylab("WTW (s)") +
    theme(legend.position = "None") +
    geom_hline(aes(yintercept = 2.2), color = conditionColors[2], size = 1, linetype = "dashed") +
    ggtitle("Observed") + theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(c(0, 20))
  
  ############### replicate results using individual parameters#######
  ################# plot on trial level #######
  set.seed(123)
  source("subFxs/modelRep.R")
  source('subFxs/helpFxs.R')
  source('exp1_expSchematics.R')
  load(file = sprintf("../../genData/wtw_exp1/expModelRep/%s_trct.RData", "QL2"))
  timeWTW_ = repOutputs$timeWTW_[c("11","20", "79")]
  trialWTW_ = repOutputs$trialWTW_[c("11","20", "79")]
  figModelInd = data.frame(
    time =  unlist(sellTimes_),
    wtw = unlist(trialWTW_),
    id = rep(selectedIds, sapply(sellTimes_, length))
  ) %>% ggplot(aes(time, wtw, color = id)) + geom_line(size = 1) + geom_point() +
    myTheme + scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
    xlab("Task time (min)") + 
    scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
    ylab("WTW (s)") +
    theme(legend.position = "None") +
    geom_hline(aes(yintercept = 2.2), color = conditionColors[2], size = 1, linetype = "dashed") +
    ggtitle("Model-generated\nwith individually-fitted parameters") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(c(0, 20))

  
  ############### print parameters ##############
  paraNames = getParaNames("QL2")
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", "QL2"))
  # get the individual parameters for each 
  expPara[expPara$id %in% selectedIds,c(paraNames, "id")]
  
  ############### resimulate with group average ############
  hdrData = hdrData[hdrData$stress == "no_stress",]
  indParas = expPara[expPara$id %in% selectedIds,]
  indStats = sumStats[sumStats$id %in% selectedIds,]
  aveParas = as.numeric(apply(expPara[,1:5], mean, MARGIN = 2)) 
  aveRepOutputs = modelRep("QL2", trialData[selectedIds], selectedIds, nRep = 10, T, aveParas) 
  figModelGroup = data.frame(
    time = unlist(sellTimes_),
    wtw = unlist(aveRepOutputs$trialWTW_),
    id = rep(selectedIds, sapply(sellTimes_, length))
  ) %>% ggplot(aes(time, wtw, color = id)) + geom_line(size = 1) + geom_point() +
    myTheme + scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
    xlab("Task time (min)") + 
    scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
    ylab("WTW (s)") +
    theme(legend.position = "None") +
    geom_hline(aes(yintercept = 2.2), color = conditionColors[2], size = 1, linetype = "dashed") +
    ggtitle("Model-generated\nwith group-averaged parameters") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(c(0, 20))
  
  
  ############# return outputs #############
  outputs = list(
    "emp" = figEmp,
    "modelInd" = figModelInd,
    "modelGroup" = figModelGroup,
    "indParas" = indParas,
    "indStats" = indStats
  )
  return(outputs)
}


######### replicate
# iti = 2
# source("exp1_expSchematics.R")
# normResults = expSchematics(0, iti, F)
# optimWaitThresholds = normResults$optimWaitThresholds
# data.frame(
#   wtw = as.vector(timeWTW_),
#   id = rep(ids, each = length(tGrid)),
#   time = rep(tGrid, length(ids)),
#   condition = rep(hdrData$condition, each = length(tGrid))) %>%
#   group_by(condition, time) %>%
#   summarise(mu = mean(wtw), se = sd(wtw) / sqrt(30)) %>% ungroup() %>%
#   ggplot(aes(time, mu, color = condition)) + geom_line() + myTheme +
#   scale_x_continuous(breaks = 0:3 * 420, labels = 0:3 * 7) + 
#   geom_ribbon(aes(ymin= mu - se, ymax= mu + se, fill = condition, color = NA), alpha = 0.5) +
#   theme(legend.position = "none") +
#   scale_fill_manual(values = conditionColors) +
#   scale_color_manual(values = conditionColors) +
#   geom_hline(aes(yintercept = optimWaitThresholds$HP), color = "red", size = 2,  linetype = "dashed") +
#   geom_hline(aes(yintercept = optimWaitThresholds$LP), color = "red", size = 2, linetype = "dashed") +
#   ylab("Willingness to wait (s)") + xlab("Simulation time (min)") + ylim(c(0, 20))
# ggsave("tmp/fig_WTW_rep.pdf", width = 4, height = 4)
