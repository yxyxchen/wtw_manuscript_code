expParaAnalysis = function(){
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
  modelName = "QL2"
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  
  # load expPara
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  parentDir = "../../genData/wtw_exp2/expModelFit"
  dirName = sprintf("%s/%s",parentDir, modelName)
  expPara = loadExpPara(paraNames, dirName)
  passCheck = checkFit(paraNames, expPara)
  
  
  # plot hist
  plotData = expPara %>% filter(passCheck) %>% select(c(paraNames)) 
  tmp = plotData %>% gather(key = 'paraname', value = 'value')
  scales_x <- list(
    'alpha' = scale_x_continuous(limits =  c(-0.05, 0.35), breaks = c(0,  0.3), labels = c(0,  0.3)),
    "nu" =  scale_x_continuous(limits =  c(-0.5, 5.5), breaks = c(0, 5), labels = c(0, 5)),
    "tau" =  scale_x_continuous(limits =  c(-0.5, 22.5), breaks = c(0.1, 22), labels = c(0.1, 22)),
    "gamma" = scale_x_continuous(limits =  c(0.65, 1.05), breaks = c(0.7, 1), labels = c(0.7, 1)),
    "eta" = scale_x_continuous(limits =  c(-0.5, 7), breaks = c(0, 6.5), labels = c(0, 6.5))
  )
  tmp$paraname = factor(tmp$paraname, levels = paraNames)
  priors = data.frame(
    paraname = paraNames,
    lower = c(0, 0, 0.1, 0.7, 0),
    upper = c(0.3, 5, 22, 1, 6.5)
  )
  
  library(scales)
  library(facetscales)
  hist = ggplot(tmp) + geom_histogram(aes(x=value),  bins = 20) +
    facet_grid_sc(cols = vars(paraname), scales = list(x = scales_x), labeller = label_parsed) + theme(legend.position = "none") + myTheme +
    geom_segment(data = priors, aes(x = lower, xend = upper, y = -0.5, yend = -0.5, color = "red")) + 
    xlab("") + ylab("Count")
  
  ################## plot correlations ####################
  plotData = expPara %>% filter(passCheck ) %>% select(c(paraNames)) %>%
    gather(key = "para", value = "value") %>%
    mutate(para = factor(para, levels = paraNames, labels = paraNames))
  plotData = expPara %>% filter(passCheck ) %>% select(c(paraNames)) 
  
  # chart.Correlation(expPara[passCheck,1:5], histogram=TRUE, pch=19, method = "kendall")
  
  pdf("../../figures/cmb/exp2_corr.pdf", width = 7.5, height = 7.5) 
  pairs(plotData[passCheck,1:5], gap=0, lower.panel = my.reg, upper.panel = my.panel.cor, main= "Exp.2") 
  dev.off()
  
  ########## ptimism bias ########
  nuTest = wilcox.test(expPara$nu - 1)
  numOptim = sum(expPara$nu[passCheck] < 1)
  ###### return outputs ########
  outputs = list(
    "hist" = hist,
    "nuTest" = nuTest,
    "numOptim" = numOptim
  )
  
  return(outputs)
}



