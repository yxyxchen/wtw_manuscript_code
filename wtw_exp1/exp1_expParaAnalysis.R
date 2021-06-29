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
  library("corrplot")
  
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData        
  trialData = allData$trialData       
  condition = hdrData$condition[hdrData$stress == "no_stress"]
  MFResults = MFAnalysis(isTrct = F)
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
      value =  plotData[,paraName],
      condition = plotData$condition
    )
    outPs[[i]] = thisPlotData %>% ggplot(aes(value, fill = condition)) +
      geom_histogram(bins = 10, color = "black") +
      facet_grid(condition ~.)  +
      scale_fill_manual(values = conditionColors) + myTheme +
      ggtitle(thisTitle) +
      theme(legend.position = "None",
            plot.title = element_text(hjust = 0.5),
            strip.background = element_blank(),
            strip.text.y = element_blank()) +
      xlab("") + scale_y_continuous(breaks = c(0, 12), limits = c(0, 13)) +
      xlim(lims[[i]]) + ylab("")
  }
  
  library(patchwork)
  hist = (outPs[[1]] | outPs[[2]] | outPs[[3]] | outPs[[4]] | outPs[[5]]) 
  
  #################################################################
  ##                correlations among parameters                ##
  #################################################################
  # chart.Correlation(plotData[plotData$condition == "HP",1:5], histogram=TRUE, pch=19, method = "kendall") 
  # chart.Correlation(plotData[plotData$condition == "LP",1:5], histogram=TRUE, pch=19, method = "kendall") 
  # chart.Correlation(plotData[,1:5], histogram=F, pch=19, method = "kendall")
  pdf("../../figures/cmb/exp1_HP_corr.pdf", width = 7.5, height = 7.5) 
  pairs(plotData[plotData$condition == "HP",1:5], gap=0, lower.panel = my.reg, upper.panel = my.panel.cor, main= "Exp.1 HP") 
  dev.off()

  pdf("../../figures/cmb/exp1_LP_corr.pdf", width = 7.5, height = 7.5) 
  pairs(plotData[plotData$condition == "LP",1:5], gap=0, lower.panel = my.reg, upper.panel = my.panel.cor, main= "Exp.1 LP") 
  dev.off()
  
  #################################################################
  ##                        optimism bias                        ##
  #################################################################
  rhoHPTest = wilcox.test(expPara$rho[passCheck & condition == "HP"] - 1)
  rhoLPTest = wilcox.test(expPara$rho[passCheck & condition == "LP"] - 1)
  rhoTest = wilcox.test(expPara$rho[passCheck] - 1)
  numOptim = sum(expPara$rho[passCheck] < 1)
  
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
  # chart.Correlation(expPara[expPara$condition == "HP" & passCheck, c(paraNames, traits)], hist = F, method = "kendall")
  # chart.Correlation(expPara[expPara$condition == "LP" & passCheck, c(paraNames, traits)], method = "kendall")
  # chart.Correlation(expPara[passCheck, c(paraNames, traits)], method = "kendall")

  
  pdf("../../figures/cmb/HP_trait_corr.pdf", width = 13.5, height = 13.5) 
  pairs(expPara[expPara$condition == "HP" & passCheck, c(paraNames, traits)], gap=0, lower.panel = my.reg, upper.panel = my.panel.cor, main= "Exp.1 HP", nCmp = length(traits)) 
  dev.off()
  
  pdf("../../figures/cmb/LP_trait_corr.pdf", width = 13.5, height = 13.5) 
  pairs(expPara[expPara$condition == "LP", c(paraNames, traits)], gap=0, lower.panel = my.reg, upper.panel = my.panel.cor, main= "Exp.1 LP", nCmp = length(traits)) 
  dev.off()
  
  
  ##################################################################
  ##                     AUC  and self-report                     ##
  ##################################################################
  sumStats$BDI = personality$BDI
  sumStats$IUS = personality$IUS
  sumStats$DoG = personality$DoG
  sumStats$BIS.11 = personality$BIS.11
  sumStats$STAI_T = personality$STAI_T
  sumStats$PSS = personality$PSS
  
  plotData= sumStats[,c("condition", traits , "muWTW", "stdWTW")]
  plotData = gather(plotData, key = "trait", value = "value", -"condition", -"muWTW", -"stdWTW")
  plotData$value = as.numeric(plotData$value)
  
  # compose the correlation equations
  lm_eqn = function(muWTW, value){
    corrTest = cor.test(muWTW, value)
    eq = substitute(italic(r)~"="~corrCoef~","~italic(p)~"="~pvalue, 
                    list(corrCoef = format(as.numeric(corrTest$estimate), digits = 2),
                         pvalue = format(corrTest$p.value, digits = 3)))
    return(as.character(as.expression(eq)))  
  }
  
  # plot the first layer
  p = plotData %>% ggplot(aes(value, muWTW, color = condition)) +
    facet_grid(condition ~ trait, scales = "free") + geom_point(alpha = 0.8) +
    ylab("AUC (s)") + xlab("") + scale_color_manual(values = conditionColors) +
    myTheme +
    theme(legend.position = "None") +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) 
  
  # add the correlation stats
  eqDf = plotData %>% group_by(condition, trait) %>% summarise(eq = lm_eqn(muWTW, value)) %>% ungroup()
  tempt =  ggplot_build(p)$layout$panel_scales_x
  rangeDf = data.frame(
    low = sapply(tempt, function(x) x$range$range[1]),
    up = sapply(tempt, function(x) x$range$range[2])
  ) %>% mutate(pos = low + (up - low) * 0.5)
  eqDf = cbind(eqDf, sapply(rangeDf, rep.int, times = 2))
  figTraitAUC = p + geom_text(data=eqDf, aes(x = pos, y = 21, label = eq), parse = TRUE, inherit.aes=FALSE) +
    ylim(c(0, 22)) +
    ggtitle("AUC vs self-report measurments") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ##################################################################
  ##                     sigma_WTW and self-report                     ##
  ##################################################################
  # plot the first layer
  p = plotData %>% ggplot(aes(value, stdWTW, color = condition)) +
    facet_grid(condition ~ trait, scales = "free") + geom_point(alpha = 0.8) +
    ylab(expression(bold(sigma[WTW])~"("~"s"^2~")")) + xlab("") +
    scale_color_manual(values = conditionColors) +
    myTheme +
    theme(legend.position = "None") +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) 
  
  eqDf = plotData %>% group_by(condition, trait) %>% summarise(eq = lm_eqn(stdWTW, value)) %>% ungroup()
  tempt =  ggplot_build(p)$layout$panel_scales_x
  rangeDf = data.frame(
    low = sapply(tempt, function(x) x$range$range[1]),
    up = sapply(tempt, function(x) x$range$range[2])
  ) %>% mutate(pos = low + (up - low) * 0.5)
  eqDf = cbind(eqDf, sapply(rangeDf, rep.int, times = 2))
  figTraitSigma = p + geom_text(data=eqDf, aes(x = pos, y = 9, label = eq), parse = TRUE, inherit.aes=FALSE) +
    ggtitle(expression(bold(sigma[WTW]~" vs self-report measurements"))) +
    theme(plot.title = element_text(hjust = 0.5))
  
  #################################################################
  ##                  Correlations among traits                  ##
  #################################################################
  pdf(file.path(figDir, "exp1_trait.pdf"), width = 8, height = 8)
  pairs(personality[,traits], gap=0, lower.panel = my.reg, upper.panel = my.panel.cor,
        main= "Self-report measurement") 
  dev.off()
  ####### combine figTraitAUC and figTriatSigma #######
  figTraitPerform = figTraitAUC & figTraitSigma 
  ggsave(file.path(figDir, "exp1_trait_perform.pdf"),  figTraitPerform, width = 24, height = 8)
  
  ############### return outputs #############
  outputs = list(
    "hist" = hist,
    "rhoHPTest" = rhoHPTest,
    "rhoLPTest" = rhoLPTest,
    "rhoTest" = rhoTest,
    "numOptim" = numOptim
  )
  return(outputs)
}


