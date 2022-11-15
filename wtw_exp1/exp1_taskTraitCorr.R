taskTraitCorr = function(){

  load("expParas.RData")
  library("ggplot2"); 
  library("dplyr"); library("tidyr")
  source("subFxs/plotThemes.R")
  source("subFxs/loadFxs.R") # load blockData and expPara
  source("subFxs/helpFxs.R") # getparaNames
  source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
  source('MFAnalysis.R')
  library(latex2exp)
  # library("corrplot")
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData        
  trialData = allData$trialData       
  condition = hdrData$condition[hdrData$stress == "no_stress"]
  ids = hdrData$id[hdrData$stress == "no_stress"]  
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  
  ##################################################################
  ##                     AUC  and self-report                     ##
  ##################################################################
  # load and merge trait data
  personality = read.csv("data/hdrData.csv")
  personality = personality[personality$id %in% ids,]
  traits = c("BDI","IUS","DoG","BIS.11", "STAI_T", "PSS") #BIS_11 trait anxiety 
  nTrait = length(traits)
  sumStats$BDI = personality$BDI
  sumStats$IUS = personality$IUS
  sumStats$DoG = personality$DoG
  sumStats$BIS.11 = personality$BIS.11
  sumStats$STAI_T = personality$STAI_T
  sumStats$PSS = personality$PSS
  
  plotData= sumStats[,c("condition", traits , "auc", "stdWTW")]
  plotData['delta'] = MFResults$sub_auc_[,6] - MFResults$sub_auc_[,1]
  plotData = gather(plotData, key = "trait", value = "value", -"condition", -"auc", -"stdWTW", -"delta")
  plotData$value = as.numeric(plotData$value)
  
  ##################################################################
  ##                     AUC and self-report                     ##
  ##################################################################
  # compose the correlation equations
  lm_eqn = function(task_measure, trait_measure){
    corrTest = cor.test(task_measure, trait_measure, method = "kendall")
    eq = substitute(italic(r)~"="~corrCoef~","~italic(p)~"="~pvalue, 
                    list(corrCoef = format(as.numeric(corrTest$estimate), digits = 2),
                         pvalue = format(corrTest$p.value, digits = 3)))
    return(as.character(as.expression(eq)))  
  }
  
  # plot the first layer
  p = plotData %>% ggplot(aes(value, auc, color = condition)) +
    facet_grid(condition ~ trait, scales = "free") + geom_point(alpha = 0.8) +
    ylab("AUC (s)") + xlab("") + scale_color_manual(values = conditionColors) +
    myTheme +
    theme(legend.position = "None") +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) 
  
  # add the correlation stats
  eqDf = plotData %>% group_by(condition, trait) %>% summarise(eq = lm_eqn(auc, value)) %>% ungroup()
  tempt =  ggplot_build(p)$layout$panel_scales_x
  rangeDf = data.frame(
    low = sapply(tempt, function(x) x$range$range[1]),
    up = sapply(tempt, function(x) x$range$range[2])
  ) %>% mutate(pos = low + (up - low) * 0.5)
  eqDf = cbind(eqDf, sapply(rangeDf, rep.int, times = 2))
  figTraitAUC = p + geom_text(data=eqDf, aes(x = pos, y = 21, label = eq), parse = TRUE, inherit.aes=FALSE) +
    ylim(c(0, 22)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # adjust p values
  pvalue_df = plotData %>% group_by(condition, trait) %>% summarise(p = cor.test(auc, value, method = "kendall")$p.value) %>% ungroup()
  p.adjust(pvalue_df$p[pvalue_df$condition == "HP"], method = "fdr")
  p.adjust(pvalue_df$p[pvalue_df$condition == "LP"], method = "fdr")
  ##################################################################
  ##                     sigma_WTW and self-report                     ##
  ##################################################################
  # plot the first layer
  p = plotData %>% ggplot(aes(value, stdWTW, color = condition)) +
    facet_grid(condition ~ trait, scales = "free") + geom_point(alpha = 0.8) +
    ylab(expression(bold(sigma[WTW]~"(s)"))) + xlab("") +
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
  figTraitSigma = p + geom_text(data=eqDf, aes(x = pos, y = 9, label = eq), parse = TRUE, inherit.aes=FALSE)
  
  pvalue_df = plotData %>% group_by(condition, trait) %>% summarise(p = cor.test(stdWTW, value, method = "kendall")$p.value) %>% ungroup()
  p.adjust(pvalue_df$p[pvalue_df$condition == "HP"], method = "fdr")
  p.adjust(pvalue_df$p[pvalue_df$condition == "LP"], method = "fdr")
  
  ##################################################################
  ##                     delta_WTW and self-report                     ##
  ##################################################################
  # plot the first layer
  p = plotData %>% ggplot(aes(value, delta, color = condition)) +
    facet_grid(condition ~ trait, scales = "free") + geom_point(alpha = 0.8) +
    ylab(expression(bold(paste(AUC[end]-AUC[start], " (s)")))) + xlab("") +
    scale_color_manual(values = conditionColors) +
    myTheme +
    theme(legend.position = "None") +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) 
  
  eqDf = plotData %>% group_by(condition, trait) %>% summarise(eq = lm_eqn(delta, value)) %>% ungroup()
  tempt =  ggplot_build(p)$layout$panel_scales_x
  rangeDf = data.frame(
    low = sapply(tempt, function(x) x$range$range[1]),
    up = sapply(tempt, function(x) x$range$range[2])
  ) %>% mutate(pos = low + (up - low) * 0.5)
  eqDf = cbind(eqDf, sapply(rangeDf, rep.int, times = 2))
  figTraitDelta = p + geom_text(data=eqDf, aes(x = pos, y = 9, label = eq), parse = TRUE, inherit.aes=FALSE)  +
    theme(plot.title = element_text(hjust = 0.5))
  
  pvalue_df = plotData %>% group_by(condition, trait) %>% summarise(p = cor.test(delta, value, method = "kendall")$p.value) %>% ungroup()
  p.adjust(pvalue_df$p[pvalue_df$condition == "HP"], method = "fdr")
  p.adjust(pvalue_df$p[pvalue_df$condition == "LP"], method = "fdr")
  
  
  # return 
  output = list(
    "figTraitAUC" = figTraitAUC,
    "figTraitSigma" = figTraitSigma,
    "figTraitDelta" = figTraitDelta
  )
return(output)
}