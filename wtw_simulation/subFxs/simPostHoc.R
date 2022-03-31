simPostHoc = function(modelName, paraLabels, paraSamples_, delays_){
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
  source("subFxs/averageRes.R")
  
  # get the generative model 
  source(sprintf("./subFxs/simModels/default/%s.R", modelName))
  simModel = get(modelName)
  paraNames = getParaNames(modelName)
  paraNames = factor(paraNames, levels = paraNames)
  nPara = length(paraNames)
  
  
  # num of repetitions 
  exist = length(delays_) > 0
  if(exist){
    nRep = length(delays_)
  }else{
    nRep = 10
  }
  duration = 120 * 60
  
  # boundaries of analyses windows
  recorded_timepoints = c(0, 4 * 60, 8 * 60, 12 * 60, 16 * 60, duration)
  nRecord = length(recorded_timepoints)
  cutValues = c(
    "#cccccc",
    "#969696",
    "#636363",
    "#252525"
  )
  
  ######################### simulate with multiple parameter combinations #####################
  # loop over conditions
  for(condition in conditions){
    paraSamples_ = 
      list("HP" = data.frame(
        alpha = c(0.000692, 0.00882, 0.0408, 0.0952),
        nu = c(0.0153, 1.26, 2.23, 2.98),
        tau = c(0.464, 2.17, 3.10, 5.88),
        gamma = c(0.780, 0.841, 0.903, 0.968),
        prior = c(0.924, 2.81, 5.13, 10.4)
      ), "LP" = data.frame(
        alpha = c(0.000692, 0.00882, 0.0408, 0.0952),
        nu = c(0.0153, 1.26, 2.23, 2.98),
        tau = c(0.464, 2.17, 3.10, 5.88),
        gamma = c(0.780, 0.841, 0.903, 0.968),
        prior = c(0.924, 2.81, 5.13, 10.4)
      ))
    default_paras = as.numeric(paraSamples[2, ])
    # default_paras = c(0.01, 1, 5, 0.85, 1)
    # default_paras = c(0.01, 1, 5, 0.85, 4)
    # generate parameter combinations 
    paraSamples = paraSamples_[[condition]]
    paraCombs = t(matrix(rep(default_paras, nCut * nPara), nrow = nPara))
    for(pIdx in 1 : nPara){
      paraCombs[(nCut * (pIdx - 1) + 1) : (nCut * pIdx) ,pIdx] = paraSamples[,pIdx] 
    }
    nComb = nrow(paraCombs)
    # initialize outputs
    sub_auc_ = matrix(NA, nrow = nRecord, ncol = nComb)
    if(condition == "HP"){
      nStep = delayMaxs[1] + 1
    }else{
      nStep = delayMaxs[2] + 1
    }
    shortterm_rv_ = matrix(NA, nrow = nStep, ncol = nComb)
    shortterm_quit_ =  vector( length = nComb)
    shortterm_wait_ = matrix(NA, nrow = nStep, ncol = nComb)
    longterm_rv_ = matrix(NA, nrow = nStep, ncol = nComb)
    longterm_quit_ =  vector( length = nComb)
    longterm_wait_ = matrix(NA, nrow = nStep, ncol = nComb)
    # loop over combinations 
    for(i in 1 : nComb){
      paras = paraCombs[i,]
      sim_ = vector(mode = "list", length = nRep)
      for(j in 1 : nRep){
        if(length(delays_) > 0){
          sim_[[j]] = simModel(paras, condition, duration, normResults, delays_[[j]][[condition]])
        }else{
          sim_[[j]] = simModel(paras, condition, duration, normResults)
        }
      }
      aveRes = averageRes(sim_, lower_boundaries, upper_boundaries)
      sub_auc_[,i] = aveRes$sub_auc_
      shortterm_rv_[,i] = aveRes$wait_minus_quit_[,nRecord - 1]
      shortterm_wait_[,i] = aveRes$wait_[, nRecord - 1]
      shortterm_quit_[i] = aveRes$quit_[nRecord-1]
      longterm_wait_[,i] = aveRes$wait_[, nRecord]
      longterm_rv_[,i] = aveRes$wait_minus_quit_[,nRecord]
      longterm_quit_[i] = aveRes$quit_[nRecord]
    }
    # plot
    auc_df = data.frame(
      auc = as.vector(sub_auc_),
      time = 1 : nRecord,
      parameter = factor(rep(paraLabels, each = nRecord * nCut), levels = paraLabels),
      rank = rep(factor(1 : nCut), each = nRecord)
    )
    auc_df %>% ggplot(aes(time, auc, color = rank)) + facet_grid(~parameter) +
      geom_line() + 
      scale_color_manual(values = cutValues) +
      ylim(c(0, min(delayMaxs)))
    
    longterm_rv_df = data.frame(
      rv = as.vector(longterm_rv_),
      time = 1 : nStep,
      parameter = factor(rep(paraLabels, each = nCut * nStep), levels = paraLabels),
      rank = rep(rep(factor(1 : nCut), each = nStep), nPara)
    )
    
    longterm_rv_df %>% filter(time < nStep) %>% ggplot(aes(time, rv, color = rank)) + facet_grid(~parameter) +
      geom_line()
    
    # combine data 
    shortterm_rv_df = data.frame(
      rv = as.vector(shortterm_rv_),
      time = 1 : nStep,
      parameter = factor(rep(paraLabels, each = nCut * nStep), levels = paraLabels),
      rank = rep(factor(1 : nCut), each = nStep)
    )
    shortterm_rv_df %>% filter(time < nStep) %>% ggplot(aes(time, rv, color = rank)) + facet_grid(~parameter) +
      geom_line()
    

    shortterm_quit_df = data.frame(
      quit = as.vector(shortterm_quit_),
      parameter = factor(rep(paraLabels, each = nCut), levels = paraLabels),
      rank = factor(1 : nCut)
    )
    shortterm_quit_df %>% ggplot(aes(rank, quit)) + facet_grid(~parameter) + 
      geom_point()
  }
    
    


    

    

  
  # plot AUC and CIP

  
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





