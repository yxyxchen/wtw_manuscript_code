expParaAnalysis = function(){
  load("expParas.RData")
  library(latex2exp)
  library("ggplot2"); library("Hmisc"); library("coin")
  library("dplyr"); library("tidyr")
  source("subFxs/plotThemes.R")
  source("subFxs/loadFxs.R") # load blockData and expPara
  source("subFxs/helpFxs.R") # getparaNames
  source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
  source('MFAnalysis.R')
  library("PerformanceAnalytics")
  
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
    outPs[[i]] = thisPlotData %>% ggplot(aes(value, fill = condition)) +
      geom_histogram(bins = 15, color = "black", fill = "grey") + myTheme +
      ggtitle(thisTitle) +
      theme(legend.position = "None",
            plot.title = element_text(hjust = 0.5)) +
      xlab("") +
      scale_y_continuous(breaks = c(0, 12), limits = c(0, 13)) +
      xlim(lims[[i]]) + ylab("")
  }
  
  library(patchwork)
  hist = (outPs[[1]] | outPs[[2]] | outPs[[3]] | outPs[[4]] | outPs[[5]]) 
  
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
  rhoTest = wilcox.test(expPara$rho - 1)
  numOptim = sum(expPara$rho[passCheck] < 1)
  ###### return outputs ########
  outputs = list(
    "hist" = hist,
    "rhoTest" = rhoTest,
    "numOptim" = numOptim
  )
  
  return(outputs)
}



