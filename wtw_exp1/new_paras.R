load("expParas.RData")
library(latex2exp)
library("ggplot2");
library("dplyr"); library("tidyr")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getparaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
source('MFAnalysis.R')
# library("PerformanceAnalytics")

# model Name
modelName = "RL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)

# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
parentDir = "../../genData/wtw_exp1/expModelFit"
dirName = sprintf("%s/%s",parentDir, modelName)
expPara = loadExpPara(paraNames, dirName)
passCheck = checkFit(paraNames, expPara)
hist(expPara$tau[passCheck])
max(expPara$tau)

hist(expPara$eta[passCheck])
max(expPara$eta[passCheck])

# plot hist
plotData = expPara %>% filter(passCheck) %>% select(c(paraNames)) 
tmp = plotData %>% gather(key = 'paraname', value = 'value')
# scales_x <- list(
#   'alpha' = scale_x_continuous(limits =  c(-0.05, 0.35), breaks = c(0,  0.3), labels = c(0,  0.3)),
#   "nu" =  scale_x_continuous(limits =  c(-0.5, 5.5), breaks = c(0, 5), labels = c(0, 5)),
#   "tau" =  scale_x_continuous(limits =  c(-0.5, 22.5), breaks = c(0.1, 22), labels = c(0.1, 22)),
#   "eta" = scale_x_continuous(limits =  c(-0.5, 7), breaks = c(0, 6.5), labels = c(0, 6.5)),
#   "beta" = scale_x_continuous(limits =  c(-0.05, 0.35), breaks = c(0,  0.3), labels = c(0,  0.3))
# )
scales_x <- list(
  'alpha' = scale_x_continuous(limits =  c(-0.05, 0.35), breaks = c(0,  0.3), labels = c(0,  0.3)),
  "tau" =  scale_x_continuous(limits =  c(-0.5, 22.5), breaks = c(0.1, 22), labels = c(0.1, 22)),
  "eta" = scale_x_continuous(limits =  c(-0.5, 15), breaks = c(0, 15), labels = c(0, 15)),
  "beta" = scale_x_continuous(limits =  c(-0.05, 0.35), breaks = c(0,  0.3), labels = c(0,  0.3))
)
tmp$paraname = factor(tmp$paraname, levels = paraNames)
priors = data.frame(
  paraname = paraNames,
  lower = c(0, 0.1, 0, 0),
  upper = c(0.3, 22, 12, 0.3)
)

# priors = data.frame(
#   paraname = paraNames,
#   lower = c(0, 0, 0.1, 0, 0),
#   upper = c(0.3, 5, 22, 6.5, 0.3)
# )

library(scales)
library(facetscales)
ggplot(tmp) + geom_histogram(aes(x=value),  bins = 10) +
  facet_grid_sc(cols = vars(paraname), scales = list(x = scales_x), labeller = label_parsed) + theme(legend.position = "none") + myTheme +
  geom_segment(data = priors, aes(x = lower, xend = upper, y = -0.5, yend = -0.5, color = "red")) + 
  xlab("") + ylab("Count")

