library('dplyr')
library("tidyr")
library("ggplot2")
library("ggpubr")
library("lme4")
source("subFxs/plotThemes.R")
source("MFAnalysis.R")
library("lmerTest")
library("latex2exp")
source('subFxs/loadFxs.R') 
source('subFxs/analysisFxs.R') 


# output dir
dir.create('../../figures/wtw_exp1/MFplot')

# load experiment parameters
load("expParas.RData")

# load exp data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$id[hdrData$stress == "no_stress"]  


b = kmsc(trialData[['31']], min(delayMaxs), F, grid = seq(0, 20, length.out = 20))$kmOnGrid
c = kmsc(trialData[["17"]], min(delayMaxs), F, grid = seq(0, 20, length.out = 20))$kmOnGrid

figAUC = data.frame(
  ts = rep(seq(0, 12, length.out = 20), 2),
  st = c(b, c),
  id = rep(c("31", "17"), each = 20)
) %>% ggplot(aes(ts, st, color = id)) + geom_line(size = 3) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
  theme(legend.position = 'None', text=element_text(face = "bold", size = 22)) + 
  scale_color_manual(values = c("#5ab4ac", "#d8b365")) +
  scale_x_continuous(limits = c(-0.5, 12.5), breaks = c(0, 3, 6, 9, 12))
ggsave("~/Downloads/figAUC.eps", figAUC, width = 4, height = 4)

###################################
b = c(1, 1, 1, 1, 0.98809524, 
      0.96222349, 0.93238164, 0.86989942, 0.81793093, 0.74389520, 
      0.56346534, 0.45293945, 0.34756812, 0.15798551, 0.07899275, 
      0.02899275, 0.01899275, 0.012, 0.011, 0.011)
c = c(1, 0.985, 0.953, 0.850, 0.820,
      0.810, 0.760, 0.742, 0.641, 0.531,
      0.521, 0.450, 0.430, 0.390, 0.330, 
      0.252, 0.202, 0.153, 0.053, 0.042)
figSigma = data.frame(
  ts = rep(seq(0, 12, length.out = 20), 2),
  st = c(b, c),
  id = rep(c("31", "17"), each = 20)
) %>% ggplot(aes(ts, st, color = id)) + geom_line(size = 3) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
  theme(legend.position = 'None', text=element_text(face = "bold", size = 22)) + 
  scale_color_manual(values = c("#d8b365", "#5ab4ac")) +
  scale_x_continuous(limits = c(-0.5, 12.5), breaks = c(0, 3, 6, 9, 12))
ggsave("~/Downloads/figSigma.eps", figSigma, width = 4, height = 4)

###################################
thisTrialData = trialData[['1']]
b_first = kmsc(thisTrialData[thisTrialData$sellTime < 210 & thisTrialData$blockNum == 1,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid
b_second = kmsc(thisTrialData[thisTrialData$sellTime >= 210 & thisTrialData$blockNum == 3,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid


thisTrialData = trialData[['2']]
c_first = kmsc(thisTrialData[thisTrialData$sellTime < 210 & thisTrialData$blockNum == 1,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid
c_first[4:20]  = c_first[4:20]  + 0.1
c_second = kmsc(thisTrialData[thisTrialData$sellTime >= 210 & thisTrialData$blockNum == 3,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid
c_second = c(1, 0.923, 0.850, 0.723, 0.711,
             0.702, 0.620, 0.612, 0.610, 0.609,
             0.609, 0.421, 0.421,0.421, 0.311,
             0.301, 0.273, 0.273, 0.273, 0.273)
figDelta = data.frame(
  ts = rep(seq(0, 12, length.out = 20), 4),
  st = c(b_first, b_second, c_first, c_second),
  chunck = factor(rep(rep(c("HP", "LP"), each = 20), 2), levels = c("HP", "LP")),
  id = rep(c("1", "2"), each = 40),
  color = rep(c("1", "2", "3", "4"), each = 20)
) %>% ggplot(aes(ts, st, color = id, linetype = chunck), alpha = 0.8) + geom_line(size = 3) + myTheme +
  scale_color_manual(values = c("#d8b365", "#5ab4ac")) +
  theme(legend.title = element_blank(), text=element_text(face = "bold", size = 22)) +
  xlab("Elapsed time (s)") + ylab("Survival rate") +
  scale_x_continuous(limits = c(-0.5, 12.5), breaks = c(0, 3, 6, 9, 12)) +
  scale_linetype_manual( values = c("solid", "dashed")) +
  theme(legend.position = "None")
# c("#d8b365", "#8c510a", "#5ab4ac", "#01665e")
ggsave("~/Downloads/figDelta.eps", figDelta, width = 4, height = 4)
  

