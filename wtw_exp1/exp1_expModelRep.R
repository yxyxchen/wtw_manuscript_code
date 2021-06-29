# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName, allData = NULL, MFResults = NULL, repOutputs = NULL){
  set.seed(123)
  # create output directories
  # dir.create("../../figures/wtw_exp1")
  # dir.create("../../figures/wtw_exp1/expModelRep/")
  # dir.create(sprintf("../../figures/wtw_exp1/expModelRep/%s",modelName))
  
  # load experiment parameters
  load('expParas.RData')
  
  # load packages and sub functions 
  library("tidyverse")
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
  
  
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  
  # load empirical data 
  if(is.null(allData)){
    allData = loadAllData()
  }
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]
  nSub = length(ids)

  # check fit
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  # compare willingness to wait (WTW) from empirical and replicated data
  ## WTW from empirical data 
  if(is.null(MFResults)){
    source("MFAnalysis.R")
    MFResults = MFAnalysis(isTrct = T)
  }
  sumStats = MFResults[['sumStats']]
  muWTWEmp = sumStats$muWTW
  stdWTWEmp = sumStats$stdWTW
  
  ## replicate data
  if(is.null(repOutputs)){
    repOutputs =  modelRep(modelName, trialData, ids, nRep, T)
    save(repOutputs, file = sprintf("../../genData/wtw_exp1/expModelRep/%s_trct.RData", modelName))
  }

  plotData = data.frame(mu =  repOutputs$muWTWRep_mu, std = repOutputs$stdWTWRep_mu,
                        empMu = muWTWEmp, empStd = stdWTWEmp,
                        passCheck = passCheck, 
                        condition = sumStats$condition) %>% filter(passCheck)
  
  #################################################################
  ##      compare observed and replicated AUC and sigma_wtw      ##
  #################################################################
  
  # plot to compare average willingess to wait
  aucFig = plotData %>%
  ggplot(aes(empMu, mu)) + 
  geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
  geom_abline(slope = 1, intercept = 0)  + 
  ylab("Model-generated AUC (s)") + xlab("Observed AUC (s)") +
  myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(breaks = c(0, 20), limits = c(-1, 21)) + 
  scale_y_continuous(breaks = c(0, 20), limits = c(-1, 21)) +
  scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + 
    ggtitle(sprintf("%s, N = %d", modelName, sum(passCheck)))
  
  ## plot to compare std willingess to wait
  cipFig = plotData %>%
  ggplot(aes(empStd, std, shape = condition)) + 
  geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
  geom_abline(slope = 1, intercept = 0) +
  ylab(expression(bold(paste("Model-generated ", sigma["WTW"], " (s"^2,")")))) +
    xlab(expression(bold(paste("Observed ", sigma["WTW"], " (s"^2,")")))) +
  myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_x_continuous(breaks = c(0, 10), limits = c(0, 10)) + 
  scale_y_continuous(breaks = c(0, 10), limits = c(0, 10)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none")

  # combind
  rep = aucFig / cipFig +  plot_annotation(title = sprintf("%s, n = %d", modelName, sum(passCheck)),
                                     theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  fileName = sprintf("../../figures/wtw_exp1/expModelRep/%s/rep.eps", modelName) 
  
  ##################################################################
  ##                     example participant                      ##
  ##################################################################
  
  # model generated 
  # repTrialData = repOutputs$repTrialData
  # repNo = repOutputs$repNo
  # sIdx = 1
  # thisRepTrialData = repTrialData[[repNo[1, sIdx]]]
  # thisRepTrialData = data.frame(thisRepTrialData[1:6])
  # modelFig = trialPlots(thisRepTrialData) +  
  #   ggtitle("Model-generated") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none")
  # 
  # # observed
  # thisTrialData = trialData[[ids[sIdx]]]
  # thisTrialData  = thisTrialData %>% filter(trialStartTime <=  blockSec - max(delayMaxs))
  # thisTrialData = block2session(thisTrialData)
  # empFig = trialPlots(thisTrialData) + ggtitle("Observed") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none")
  # example = modelFig / empFig +  plot_annotation(title = sprintf("%s", modelName),
  #                                      theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  # ggsave(sprintf("../../figures/wtw_exp1/expModelRep/%s/example.eps", modelName), example, width = 6, height = 8)
  # 
  
  ################# return figure outputs ###############
  # outputs = list("rep" = rep, "example" = example)
  outputs = list("rep" = rep, "example" = example)
  return(outputs)
  
}

