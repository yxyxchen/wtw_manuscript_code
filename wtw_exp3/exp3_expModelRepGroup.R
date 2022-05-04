expModelRepGroup = function(){
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
source("exp3_expSchematics.R")
source("subFxs/plotThemes.R")
source("subFxs/helpFxs.R") 
source("subFxs/loadFxs.R") # 
source("subFxs/analysisFxs.R") 
source("subFxs/modelRep.R")
source("MFAnalysis.R")
# load data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$id             
nSub = length(ids) 
# all passcheck
modelNames = c("QL1", "QL2", "RL1", "RL2", "naive", "omni")
nModel = length(modelNames)
passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
for(k in 1 : nModel){
  modelName = modelNames[k]
  paraNames = getParaNames(modelName)
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp3/expModelFit/%s", modelName))
  passCheck_[,k] = checkFit(paraNames, expPara)
}
apply(passCheck_, MARGIN = 2, sum)
passCheck = apply(passCheck_, MARGIN = 1, all)

# empirical data
MFResults = MFAnalysis(isTrct = T)
blockStats = MFResults$blockStats
plotdf = data.frame(
  wtw = unlist(MFResults$timeWTW_[blockStats$blockNum <= 2]),
  time = rep(seq(0, blockSec * 2 -1, by = 2), nSub),
  condition = rep(blockStats$condition[blockStats$blockNum <= 2], each = length(tGrid)),
  cbal = rep(blockStats$cbal[blockStats$blockNum <= 2], each = length(tGrid)),
  passCheck = rep(passCheck, each = length(tGrid) * 2),
  type = "emp"
) %>% filter(passCheck)  %>%
  group_by(cbal, time, type) %>% 
  summarise(mu = mean(wtw),
            se = sd(wtw) / sqrt(length(wtw)),
            condition = condition) %>%
  mutate(ymin = mu - se,
         ymax = mu + se) %>% ungroup()
# plotdf %>% ggplot(aes(time, mu)) + facet_grid(~condition) + geom_line()

for(j in 1 : nModel){
  model = modelNames[j]
  load(file = sprintf("../../genData/wtw_exp3/expModelRep/%s_trct.RData",  model))
  added_data = data.frame(
    wtw = as.vector(repOutputs$timeWTW_),
    time = rep(seq(0, blockSec * 2 -1, by = 2), nSub),
    condition = rep(blockStats$condition[blockStats$blockNum <= 2], each = length(tGrid)),
    cbal = rep(blockStats$cbal[blockStats$blockNum <= 2], each = length(tGrid)),
    passCheck = rep(passCheck, each = length(tGrid) * 2),
    type = model
  ) %>% filter(passCheck)  %>%
    group_by(cbal, time, type) %>% 
    summarise(mu = mean(wtw),
              se = sd(wtw) / sqrt(length(wtw)),
              condition = condition) %>%
    mutate(ymin = mu - se,
           ymax = mu + se) %>% ungroup()
  
  plotdf = rbind(plotdf, added_data)
}

rectData = data.frame(
  xmin = rep(((1:2) - 1) * blockSec, 2) ,
  xmax = rep(1:2 * blockSec, 2) ,
  condition = c("HP", "LP", "LP", "HP"),
  cbal = c("HP first", "HP first", "LP first", "LP first")
)

figWTW = plotdf %>%
  filter(type %in% c("QL1", "QL2", "RL1", "RL2", "emp")) %>%
  mutate(type = factor(type, levels = c("QL1", "QL2", "RL1", "RL2", "emp")),
         cbal = ifelse(cbal == 1, "HP first", "LP first")) %>%
  ggplot(aes(time, mu, color = type)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, fill = condition),
            data = rectData, ymin = 0, ymax = 20, alpha = 0.75, inherit.aes = F)  + 
  geom_segment(x = blockSec, xend = blockSec, y = 0, yend = 20, color = "#878787", linetype = "dashed") +
  geom_line(alpha = 0.75, size = 1) + myTheme +
  scale_color_manual(values = c("#d6604d", "#b2182b", "#4393c3", "#2166ac", "black"))  + xlab("Task time (min)") + 
  scale_x_continuous(breaks = 0 : 2 * 10 * 60, labels = 0 : 2 * 10) + 
  ylab("WTW (s)") + facet_grid(.~cbal) + theme(legend.position =  "None") +
  ylim(c(0, 20)) + scale_fill_manual(values = conditionColorBacks)

return(figWTW)
}