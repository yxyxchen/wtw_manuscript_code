expParaAnalysis = function(){
  load("expParas.RData")
  library("ggplot2"); library("Hmisc"); library("coin")
  library("dplyr"); library("tidyr")
  source("subFxs/plotThemes.R")
  source("subFxs/loadFxs.R") # load blockData and expPara
  source("subFxs/helpFxs.R") # getparaNames
  source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
  source('MFAnalysis.R')
  library(latex2exp)
  library("PerformanceAnalytics")
  
  # model Name
  modelName = "QL2"
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)

  
  # load expPara
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  parentDir = "../../genData/wtw_exp3/expModelFit"
  dirName = sprintf("%s/%s",parentDir, modelName)
  expPara = loadExpPara(paraNames, dirName)
  passCheck = checkFit(paraNames, expPara)
  
  # plot hist 
  plotData = expPara %>% filter(passCheck) %>% select(c(paraNames))
  lims = list(
    c(0, 0.15),
    c(0, 5.1),
    c(0, 15),
    c(0.65, 1.05),
    c(0, 8)
  )
  # bin_nums = c(10, 15, )
  outPs = list()
  for(i in 1 : nPara){
    paraName = paraNames[i]
    thisTitle = parse(text = paraName)
    thisPlotData = data.frame(
      value =  plotData[,paraName]
    )
    outPs[[i]] = thisPlotData %>% ggplot(aes(value)) +
      geom_histogram(bins = 15, color = "black", fill = "grey") + myTheme +
      ggtitle(thisTitle) +
      theme(legend.position = "None",
            plot.title = element_text(hjust = 0.5)) +
      xlab("") + scale_y_continuous(breaks = c(0, 18), limits = c(0, 20)) +
      xlim(lims[[i]]) + ylab("")
  }
  
  library(patchwork)
  figHist = (outPs[[1]] | outPs[[2]] | outPs[[3]] | outPs[[4]] | outPs[[5]]) 
  
  ################## plot correlations ####################
  plotData = expPara %>% filter(passCheck ) %>% select(c(paraNames)) %>%
    gather(key = "para", value = "value") %>%
    mutate(para = factor(para, levels = paraNames, labels = paraNames))
  plotData = expPara %>% filter(passCheck ) %>% select(c(paraNames)) 
  
  # chart.Correlation(expPara[passCheck,1:5], histogram=TRUE, pch=19, method = "spearman")
  pdf("../../figures/cmb/exp3_corr.pdf", width = 7.5, height = 7.5) 
  pairs(plotData[passCheck,1:5], gap=0, lower.panel = my.reg, upper.panel = my.panel.cor, main= "Exp.3") 
  dev.off()
  
  
  ########## ptimism bias ########
  rhoTest = wilcox.test(expPara$rho[passCheck] - 1)
  numOptim = sum(expPara$rho[passCheck] < 1)
  ######## return outputs ##########
  outputs = list(
    "rhoTest" = rhoTest,
    "hist" = figHist,
    "numOptim" = numOptim
  )
  
  return(outputs)
}

