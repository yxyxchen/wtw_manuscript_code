simPostHoc = function(modelName, paraLabels, paraSamples){
  # default settings 
  smallReward = 0
  iti = 2
  
  # random seed
  set.seed(123)

  
  # load experiment parameters
  load('expParas.RData')
  
  # normative analysis 
  normResults = expSchematics(smallReward, iti)
  optimRewardRates = normResults$optimRewardRates
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # load packages and sub functions 
  library("tidyverse")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  
  # get the generative model 
  source(sprintf("./subFxs/simModels/default/%s.R", modelName))
  simModel = get(modelName)
  paraNames = getParaNames(modelName)
  paraNames = factor(paraNames, levels = paraNames)
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  duration = 120 * 60
  nPeriod = 4
  nCut = 4
  periodBreaks = seq(0, 16 * 60, length.out = nPeriod + 1)
  
  ######################### simulate with multiple parameter combinations #####################
  # parameter combinations 
  paraSampleHPs = paraSamples[['HP']]
  paraSampleLPs = paraSamples[['LP']]

  
  for(condition in conditions){
    if(condition == 'HP'){
      paraSamples = paraSampleHPs
    }else{
      paraSamples = paraSampleLPs
    }
    paraCombs = t(matrix(rep(as.numeric(paraSamples[3, ]), nCut * nPara), nrow = nPara))
    for(pIdx in 1 : nPara){
      paraCombs[(nCut * (pIdx - 1) + 1) : (nCut * pIdx) ,pIdx] = paraSamples[,pIdx] 
    }
    nComb = nrow(paraCombs)
    
    # initialize outputs
    auc_ = matrix(NA, nrow = nComb, ncol = nRep)
    cip_ = matrix(NA, nrow = nComb, ncol = nRep)
    periodAuc_ = matrix(NA, nrow = nComb * nPeriod, ncol = nRep)
    periodCip_ = matrix(NA, nrow = nComb * nPeriod, ncol = nRep)
    asymAuc_ = matrix(NA, nrow = nComb, ncol = nRep)
    asymCip_ = matrix(NA, nrow = nComb, ncol = nRep)
    
    # loop
    for(i in 1 : nComb){
      paras = paraCombs[i, ]
      for(j in 1 : nRep){
        tempt = simModel(as.numeric(paras), condition, duration, normResults)
        # kmscResults = kmsc(tempt, min(delayMaxs), F, kmGrid)
        # auc_[i, j] = kmscResults$auc
        # cip_[i, j] = kmscResults$stdWTW
        # calculate asymmetric auc
        tempt$Qwaits_ = NULL
        tempt$Gs_ = NULL
        tempt = as.data.frame(tempt)
        kmscResults = kmsc(tempt[tempt$sellTime >= (120 - 4) * 60,], min(delayMaxs), F, kmGrid)
        asymAuc_[i, j] = kmscResults$auc
        asymCip_[i, j] = kmscResults$stdWTW       
        junk = lapply(1 : nPeriod, function(x) kmsc(tempt[periodBreaks[x] <= tempt$sellTime & tempt$sellTime < periodBreaks[x+1],], min(delayMaxs), F, kmGrid))
        periodAuc_[((i - 1) * nPeriod + 1) : (i * nPeriod), j] = sapply(1 : nPeriod, function(i) junk[[i]]$auc)
        periodCip_[((i - 1) * nPeriod + 1) : (i * nPeriod), j] = sapply(1 : nPeriod, function(i) junk[[i]]$stdWTW)
        }
      sumDf = data.frame(
        auc = apply(asymAuc_, 1, mean),
        cip = apply(asymCip_, 1, mean),
        paraLabel = rep(paraLabels, each = nCut),
        paraRank = factor(rep(1 : nCut, nPara)),
        condition = condition
      )
      df = data.frame(
        auc = apply(periodAuc_, 1, mean),
        cip = apply(periodCip_, 1, mean),
        period = rep(1 : nPeriod, nComb),
        t = rep((periodBreaks[1 : nPeriod] + periodBreaks[2 : (nPeriod+1)] )/ 120, nComb),
        paraName = rep(paraNames, each = nCut * nPeriod),
        paraLabel = rep(paraLabels, each = nCut * nPeriod),
        paraRank = factor(rep(rep(1 : nCut, each = nPeriod), nPara)),
        condition = condition
      ) 
      if(condition == "HP"){
        HPdf = df
        HPSumDf = sumDf
      }else{
        LPdf = df
        LPSumDf = sumDf
      }
    }
  }
  
  # plot AUC and CIP
  cutValues = c(
    "#cccccc",
    "#969696",
    "#636363",
    "#252525"
  )
  
  sumDf = rbind(HPSumDf, LPSumDf)
  #  
  rbind(HPSumDf, LPSumDf) %>%
    mutate(paraLabel = factor(paraLabel, levels = paraLabels)) %>%
    ggplot(aes(paraRank, cip, color = paraRank)) + geom_point() +
    geom_line(group = 1) + facet_grid(condition~paraLabel, labeller = label_parsed) +
    scale_color_manual(values = cutValues) + myTheme 
  
  
  # plot AUC time courses 
  optimDf = data.frame(
    condition = rep(conditions, each = nPara),
    optimThreshold = rep(unlist(optimWaitThresholds), each = nPara),
    paraLabel= rep(paraLabels, 2)
  )
  fig = rbind(HPdf, LPdf) %>% 
    mutate(paraLabel = factor(paraLabel, levels = paraLabels)) %>%
    ggplot(aes(t, auc, color = paraRank)) +
    geom_line() + geom_point() +
    facet_grid(condition~paraLabel, labeller = label_parsed) +
    geom_hline(data = optimDf, aes(yintercept = optimThreshold), inherit.aes = F, linetype = "dashed", color = rep(conditionColors, each = nPara)) +
    scale_color_manual(values = cutValues) +
    geom_point(aes(color = paraRank, y = auc), x = 22, data = sumDf, inherit.aes = F, shape = 8) + myTheme + xlab("Task time (min)") + ylab("AUC (s)") +
    scale_x_continuous(breaks = c(periodBreaks / 60, 22), limits = c(0, 24), labels = c(periodBreaks / 60, 20)) +
    ylim(c(0, 20))                                                                                
  
  
}





