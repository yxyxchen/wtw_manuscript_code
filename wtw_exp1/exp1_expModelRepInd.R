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
  
  # a function to plot trial-wise behavior for each participant one by one
  # for (sIdx in 1 : nSub) {
  #   id = ids[sIdx]
  #   thisTrialData = trialData[[id]] 
  #   if(unique(thisTrialData$condition) == "LP"){
  #     trctLine = blockSec - max(delayMaxs)
  #     # truncate trials completed after tractline in each block
  #     thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine)
  #     thisTrialData = block2session(thisTrialData)
  #     # p = trialPlots(thisTrialData) + ylim(c(0, 40))
  #     wtwTS(thisTrialData, tGrid, 20, plotWTW = T) 
  #     print(id)
  #     # ggsave("p26.eps", width = 8, height = 4)
  #     readline(prompt="Press [enter] to continue")
  #   }
  # }
  
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
      thisTrialData = thisTrialData[thisTrialData$sellTime < cut,]
      # p = trialPlots(thisTrialData) + ylim(c(0, 40))
      outputs = wtwTS(thisTrialData, tGrid[1:cut], 20, plotWTW = F) 
      # outputs_[[i]] = outputs$timeWTW
      outputs_[[i]] = outputs$trialWTW
      sellTimes_[[i]] = thisTrialData$sellTime
      i = i + 1
    }
  }
  
  ############### plot empirical results of three example participants #######
  figEmp = data.frame(
    wtw = unlist(outputs_),
    time = unlist(sellTimes_),
    sub = rep(selectedIds, sapply(sellTimes_, length))
  ) %>% ggplot(aes(time, wtw,  color = sub)) + geom_line(size = 1) + geom_point() +
    scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
    myTheme + scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
    xlab("Task time (min)") + ylab("Willingness to wait (s)") +
    theme(legend.position = "None") +
    geom_hline(aes(yintercept = 2.2), color = conditionColors[2], size = 1, linetype = "dashed") +
    ggtitle("Observed") + theme(plot.title = element_text(hjust = 0.5))
  
  # figEmp = data.frame(
  #   wtw = unlist(outputs_),
  #   time = rep(tGrid[1:cut], 3),
  #   sub = rep(selectedIds, each = cut)
  # ) %>% ggplot(aes(time, wtw,  color = sub)) + geom_line(size = 2) +
  #   scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) + 
  #   myTheme + scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
  #   xlab("Task time (min)") + ylab("Willingness to wait (s)") +
  #   theme(legend.position = "None") + 
  #   geom_hline(aes(yintercept = 2.2), color = conditionColors[2], size = 1, linetype = "dashed") +
  #   ggtitle("Observed") + theme(plot.title = element_text(hjust = 0.5))
  
  
  ############### replicate results using individual parameters#######
  ################# plot on trial level #######
  set.seed(123)
  source("subFxs/modelRep.R")
  source('subFxs/helpFxs.R')
  source('exp1_expSchematics.R')
  repOutputs = modelRep("QL2", trialData[c("11","20", "79")], c("11","20", "79"), nRep = 10, T)
  timeWTW_ = repOutputs$timeWTW_
  trialWTW_ = repOutputs$trialWTW_
  figModelInd = data.frame(
    time =  unlist(sellTimes_),
    wtw = unlist(trialWTW_),
    id = rep(selectedIds, sapply(sellTimes_, length))
  ) %>% ggplot(aes(time, wtw, color = id)) + geom_line(size = 1) + geom_point() +
    myTheme + scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
    xlab("Task time (min)") + 
    scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
    ylab("Willingness to wait (s)") +
    theme(legend.position = "None") +
    geom_hline(aes(yintercept = 2.2), color = conditionColors[2], size = 1, linetype = "dashed") +
    ggtitle("Model-generated\nwith individually-fitted parameters") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path("../figures/cmb","exp1_ind_rep.eps"), width = 4, height = 4)  
  
  # obs =  data.frame(
  #   wtw = unlist(outputs_),
  #   time = unlist(sellTimes_),
  #   id = rep(selectedIds, sapply(sellTimes_, length))
  # )
  # obs[['type']] = "obs"
  # rep = data.frame(
  #   time =  unlist(sellTimes_),
  #   wtw = unlist(trialWTW_),
  #   id = rep(selectedIds, sapply(sellTimes_, length)))
  # rep[['type']] = "rep"
  # plotdf = rbind(obs, rep)
  # plotdf %>% 
  #   ggplot(aes(time, wtw, color = type)) +
  #   geom_line() + myTheme +
  #   geom_point() + facet_grid(id~.)
      
    
  
  
  ######## plot with time #######
  # data.frame(
  #   time =  rep(tGrid, length(selectedIds)),
  #   wtw = as.vector(timeWTW_),
  #   id = rep(selectedIds, each = length(tGrid))
  # ) %>% ggplot(aes(time, wtw, color = id)) + geom_line(size = 1) +
  #   myTheme + scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
  #   xlab("Task time (min)") + 
  #   scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
  #   ylab("Willingness to wait (s)") +
  #   theme(legend.position = "None") +
  #   geom_hline(aes(yintercept = 2.2), color = conditionColors[2], size = 1, linetype = "dashed") +
  #   ggtitle("Model-generated\nwith individually-fitted parameters") +
  #   theme(plot.title = element_text(hjust = 0.5))
  
  ############### print parameters ##############
  paraNames = getParaNames("QL2")
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", "QL2"))
  # get the individual parameters for each 
  expPara[expPara$id %in% selectedIds,c(paraNames, "id")]
  
  ############### resimulate with group average ############
  hdrData = hdrData[hdrData$stress == "no_stress",]
  aveLPParas = as.numeric(apply(expPara[hdrData$condition == "LP",  1:5], mean, MARGIN = 2)) 
  aveRepOutputs = modelRep("QL2", trialData[selectedIds], selectedIds, nRep = 10, T, aveLPParas) 
  
  figModelGroup = data.frame(
    time = rep(tGrid, length(selectedIds)),
    wtw = as.vector(aveRepOutputs$timeWTW_),
    id = rep(selectedIds, each = length(tGrid))
  ) %>% ggplot(aes(time, wtw, color = id)) + geom_line(size = 2) +
    myTheme + scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
    xlab("Simulation time (min)") + 
    scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
    ylab("Willingness to wait (s)") +
    theme(legend.position = "None") +
    geom_hline(aes(yintercept = 2.2), color = conditionColors[2], linetype = "dashed", size = 1) +
    ggtitle("Model-generated\nwith group-averaged parameters") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ############# return outputs #############
  outputs = list(
    "emp" = figEmp,
    "modelInd" = figModelInd,
    "modelGroup" = figModelGroup
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
