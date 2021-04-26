# load libraries
source('subFxs/loadFxs.R') 
source('subFxs/analysisFxs.R') 
source("subFxs/plotThemes.R")
library('dplyr')
library("tidyr")
library("ggplot2")

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
for (sIdx in 1 : nSub) {
  id = ids[sIdx]
  thisTrialData = trialData[[id]] 
  if(unique(thisTrialData$condition) == "LP"){
    trctLine = blockSec - max(delayMaxs)
    # truncate trials completed after tractline in each block
    thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine)
    thisTrialData = block2session(thisTrialData)
    # p = trialPlots(thisTrialData) + ylim(c(0, 40))
    wtwTS(thisTrialData, tGrid, 20, plotWTW = T) 
    print(id)
    # ggsave("p26.eps", width = 8, height = 4)
    readline(prompt="Press [enter] to continue")
  }
}

# conduct wtw analysis for several example participants
cut =1260
selectedIds = c("11","20", "79")
outputs_ = list()
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
    outputs = wtwTS(thisTrialData, tGrid[1:cut], 20, plotWTW = T) 
    outputs_[[i]] = outputs$timeWTW
    i = i + 1
    print(id)
    ggsave(sprintf("tmp/wtw_%s.eps", id), width = 4, height = 4)
    # ggsave("p26.eps", width = 8, height = 4)
    # readline(prompt="Press [enter] to continue")
  }
}

############### plot empirical results of three example participants #######
# maybe I want to switch back to the non-trct version 
data.frame(
  wtw = unlist(outputs_),
  time = rep(tGrid[1:cut], 3),
  sub = rep(selectedIds, each = cut)
) %>% ggplot(aes(time, wtw,  color = sub)) + geom_line(size = 2) +
  scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) + 
  myTheme + scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
  xlab("Task time (min)") + ylab("Willingness to wait (s)") +
  theme(legend.position = "None") + 
  geom_hline(aes(yintercept = 2.2), color = "red")
ggsave("tmp/fig_ind.eps", width = 4, height = 4) 
  

############### replicate results using individual parameters#######
set.seed(123)
source("subFxs/modelRep.R")
source('subFxs/helpFxs.R')
source('exp1_expSchematics.R')
repOutputs = modelRep("QL2", trialData, ids, nRep = 10, T)
timeWTW_ = repOutputs$timeWTW_
data.frame(
  time = rep(tGrid, length(selectedIds)),
  wtw = as.vector(timeWTW_[,which(ids %in% selectedIds)]),
  id = rep(selectedIds, each = length(tGrid))
) %>% ggplot(aes(time, wtw, color = id)) + geom_line(size = 2) +
  myTheme + scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
  xlab("Simulation time (min)") + 
  scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
  ylab("Willingness to wait (s)") +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept = 2.2), color = "red")
ggsave("tmp/fig_rep_ind.pdf", width = 4, height = 4)

############### print parameters ##############
paraNames = getParaNames("QL2")
expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", "QL2"))
# get the individual parameters for each 
expPara[expPara$id %in% selectedIds,c(paraNames, "id")]

############### resimulate with group average ############
hdrData = hdrData[hdrData$stress == "no_stress",]
aveLPParas = as.numeric(apply(expPara[hdrData$condition == "LP",  1:5], mean, MARGIN = 2)) 
aveRepOutputs = modelRep("QL2", trialData[selectedIds], selectedIds, nRep = 10, T, aveLPParas) 


data.frame(
  time = rep(tGrid, length(selectedIds)),
  wtw = as.vector(aveRepOutputs$timeWTW_),
  id = rep(selectedIds, each = length(tGrid))
) %>% ggplot(aes(time, wtw, color = id)) + geom_line(size = 2) +
  myTheme + scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00")) +
  xlab("Simulation time (min)") + 
  scale_x_continuous(breaks = c(0, 420, 840, 1260), labels = c(0, 420, 840, 1260)/60) +
  ylab("Willingness to wait (s)") +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept = 2.2), color = "red")
ggsave("tmp/fig_rep_ind.pdf", width = 4, height = 4)


######### rep on the groupLevel
iti = 2
source("exp1_expSchematics.R")
normResults = expSchematics(0, iti, F)
optimWaitThresholds = normResults$optimWaitThresholds
data.frame(
  wtw = as.vector(timeWTW_),
  id = rep(ids, each = length(tGrid)),
  time = rep(tGrid, length(ids)),
  condition = rep(hdrData$condition, each = length(tGrid))) %>%
  group_by(condition, time) %>%
  summarise(mu = mean(wtw), se = sd(wtw) / sqrt(30)) %>% ungroup() %>%
  ggplot(aes(time, mu, color = condition)) + geom_line() + myTheme +
  scale_x_continuous(breaks = 0:3 * 420, labels = 0:3 * 7) + 
  geom_ribbon(aes(ymin= mu - se, ymax= mu + se, fill = condition, color = NA), alpha = 0.5) +
  theme(legend.position = "none") +
  scale_fill_manual(values = conditionColors) +
  scale_color_manual(values = conditionColors) +
  geom_hline(aes(yintercept = optimWaitThresholds$HP), color = "red", size = 2,  linetype = "dashed") +
  geom_hline(aes(yintercept = optimWaitThresholds$LP), color = "red", size = 2, linetype = "dashed") +
  ylab("Willingness to wait (s)") + xlab("Simulation time (min)") + ylim(c(0, 20))
ggsave("tmp/fig_WTW_rep.pdf", width = 4, height = 4)
