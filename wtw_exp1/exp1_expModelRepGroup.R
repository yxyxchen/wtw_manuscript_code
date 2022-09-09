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
  source("exp1_expSchematics.R")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  source("subFxs/modelRep.R")
  
  allData = loadAllData()
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]
  nSub = length(ids)
  # all passcheck
  modelNames = c("QL1", "QL2", "RL1", "RL2", "naive", "omni")
  nModel = length(modelNames)
  passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
  for(k in 1 : nModel){
    modelName = modelNames[k]
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", modelName))
    passCheck_[,k] = checkFit(paraNames, expPara)
  }
  apply(passCheck_, MARGIN =  2, sum)
  passCheck = apply(passCheck_, MARGIN = 1, all)
  
  # emp MF 
  source("MFAnalysis.R")
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults$sumStats
  plotdf = data.frame(
    wtw = unlist(MFResults$timeWTW_),
    time = rep(tGrid, nSub),
    condition = rep(sumStats$condition, each = length(tGrid)),
    passCheck = rep(passCheck, each = length(tGrid)),
    type = "emp"
  ) %>% filter(passCheck)  %>%
    group_by(condition, time, type) %>% 
    summarise(mu = mean(wtw),
              se = sd(wtw) / sqrt(length(wtw))) %>%
    mutate(ymin = mu - se,
           ymax = mu + se) %>% ungroup()
  # plotdf %>% ggplot(aes(time, mu)) + facet_grid(~condition) + geom_line()
  
  models = c("QL1", "QL2", "RL1", "RL2", "naive", "omni")
  for(j in 1 : nModel){
    # aasdads
    model = models[j]
    load(file = sprintf("../../genData/wtw_exp1/expModelRep/%s_trct.RData", model))
    added_data = data.frame(
      wtw = as.vector(repOutputs$timeWTW_),
      time = rep(tGrid, nSub),
      condition = rep(sumStats$condition, each = length(tGrid)),
      passCheck = rep(passCheck, each = length(tGrid)),
      type = model
    ) %>% filter(passCheck)  %>%
      group_by(condition, time, type) %>% 
      summarise(mu = mean(wtw),
                se = sd(wtw) / sqrt(length(wtw))) %>%
      mutate(ymin = mu - se,
             ymax = mu + se) %>% ungroup()
    
    plotdf = rbind(plotdf, added_data)
  }
  
  figWTW = plotdf %>%
    filter(type %in% c("QL1", "QL2", "RL1", "RL2", "emp")) %>%
    mutate(type = factor(type, levels = c("QL1", "QL2", "RL1", "RL2", "emp"))) %>%
    ggplot(aes(time, mu, color = type)) +
    geom_rect(aes(fill = condition), xmin = 0, xmax = blockSec * 3 , ymin = 0, ymax = 20, color = NA)  +
    geom_segment(x = blockSec, xend = blockSec, y = 0, yend = 20, color = "#878787", linetype = "dashed") + 
    geom_segment(x = blockSec * 2, xend = blockSec * 2, y = 0, yend = 20, color = "#878787", linetype = "dashed") +
    geom_line(alpha = 0.65, size = 1) + myTheme +   facet_grid(~condition) +
    scale_color_manual(values = c("#f4a582", "#b2182b", "#4393c3", "#2166ac", "black"))  +
    xlab("Task time (min)") + 
    scale_x_continuous(breaks = 0 : 3 * 7 * 60, labels = 0 : 3 * 7) + 
    ylab("WTW (s)") + theme(legend.position =  "None") +
    scale_fill_manual(values = conditionColorBacks) + ylim(c(0, 20))
  
  return(figWTW)
}

  
