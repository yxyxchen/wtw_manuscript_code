expParaAnalysis = function(){
  load("expParas.RData")
  library("ggplot2"); 
  library("dplyr"); library("tidyr")
  source("subFxs/plotThemes.R")
  source("subFxs/loadFxs.R") # load blockData and expPara
  source("subFxs/helpFxs.R") # getparaNames
  source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
  source('MFAnalysis.R')
  library(latex2exp)
  
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData        
  trialData = allData$trialData       
  condition = hdrData$condition[hdrData$stress == "no_stress"]
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  
  # model Name
  modelName = "QL2"
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # output directories
  figDir = file.path("..", "..", "figures", "cmb")
  dir.create(figDir)
  
  # load expPara
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  parentDir = "../../genData/wtw_exp1/expModelFit"
  dirName = sprintf("%s/%s",parentDir, modelName)
  expPara = loadExpPara(paraNames, dirName)
  passCheck = checkFit(paraNames, expPara)
  # expPara = merge(x=tempt[,c(paraNames, "id")],y=sumStats, by="id",all.x=TRUE)
  
  #################################################################
  ##                 plot parameter distribution                 ##
  #################################################################
  expPara$condition = condition
  plotData = expPara %>% filter(passCheck) %>% select(c(paraNames, "condition")) 
  
  tmp = gather(plotData, key = 'paraname', value = 'value', -"condition")
  tmp$paraname = factor(tmp$paraname, levels = paraNames)
  tmp %>% group_by(paraname) %>%
    summarise(q1 = quantile(value, 0.1),
              q2 = quantile(value, 0.3),
              q3 = quantile(value, 0.5),
              q4 = quantile(value, 0.7),
              q5 = quantile(value, 0.9))
  scales_x <- list(
    'alpha' = scale_x_continuous(limits =  c(-0.05, 0.35), breaks = c(0,  0.3), labels = c(0,  0.3)),
    "nu" =  scale_x_continuous(limits =  c(-0.5, 5.5), breaks = c(0, 5), labels = c(0, 5)),
    "tau" =  scale_x_continuous(limits =  c(-0.5, 22.5), breaks = c(0.1, 22), labels = c(0.1, 22)),
    "gamma" = scale_x_continuous(limits =  c(0.65, 1.05), breaks = c(0.7, 1), labels = c(0.7, 1)),
    "eta" = scale_x_continuous(limits =  c(-0.5, 16), breaks = c(0, 15), labels = c(0, 15))
  )
  priors = data.frame(
      paraname = paraNames,
      lower = c(0, 0, 0.1, 0.7, 0),
      upper = c(0.3, 5, 22, 1, 15)
    )


  library(scales)
  library(facetscales)

  hist = ggplot(tmp) + geom_histogram(aes(x=value),  bins = 20) +
    facet_grid_sc(cols = vars(paraname), scales = list(x = scales_x), labeller = label_parsed) +
    scale_fill_manual(values = conditionColors) + theme(legend.position = "none") + myTheme +
    geom_segment(data = priors, aes(x = lower, xend = upper, y = -0.5, yend = -0.5, color = "red")) + 
    xlab("") + ylab("Count") + 
    geom_segment(x = 1, xend = 1, y = 0, yend  = 30, color = "black", linetype = "dashed")
  
  para_summary = tmp %>% group_by(paraname) %>% 
    summarise(median = round(median(value), 3),
              q1 = round(quantile(value, 0.25), 3),
              q3 = round(quantile(value, 0.75), 3),
              IQR = round(IQR(value), 3)) %>% ungroup()
  #################################################################
  ##                correlations among parameters                ##
  #################################################################
  pdf("../../figures/cmb/exp1_corr.pdf", width = 7.5, height = 7.5) 
  pairs(plotData[passCheck,1:5], condition = plotData$condition[passCheck], gap=0, lower.panel = my.reg.HP, upper.panel = my.reg.LP, main= "Exp.1") 
  dev.off()
  
  #################################################################
  ##                        optimism bias                        ##
  #################################################################
  # nuHPTest = wilcox.test(expPara$nu[passCheck & condition == "HP"] - 1)
  # nuLPTest = wilcox.test(expPara$nu[passCheck & condition == "LP"] - 1)
  # nuTest = wilcox.test(expPara$nu[passCheck] - 1)
  # numOptim = sum(expPara$nu[passCheck] < 1)
  optim_summary = expPara[passCheck,] %>% group_by(condition) %>%
    summarise(nOptim = sum(nu < 1), nTotal = length(nu))
  
  ##################################################################
  ##             correlations with personality traits             ##
  ##################################################################
  # load and merge trait data
  personality = read.csv("data/hdrData.csv")
  personality = personality[personality$id %in% expPara$id,]
  personality$id = expPara$id
  traits = c("BDI","IUS","DoG","BIS.11", "STAI_T", "PSS") #BIS_11 trait anxiety 
  nTrait = length(traits)
  expPara = merge(x= expPara,y = personality[,c(traits, "id")], by="id",all.x=TRUE)
  
  # plot
  # I am not a big fan of this figure
  pdf("../../figures/cmb/exp1_trait_corr.pdf", width = 13.5, height = 13.5) 
  pairs(expPara[passCheck, c(paraNames, traits)], condition = expPara$condition[passCheck], gap=0, lower.panel = my.reg.HP, upper.panel = my.reg.LP, main= "Exp.1", nCmp = length(traits)) 
  dev.off()
  

  
  ############### return outputs #############
  outputs = list(
    "hist" = hist,
    "para_summar" = para_summary,
    "optim_summary" = optim_summary
  )
  return(outputs)
}


