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
  lower_boundaries = c(0, 4 * 60, 8 * 60, 12 * 60, duration - 4 * 60)
  upper_boundaries = lower_boundaries + 4 * 60
  
  ######################### simulate with multiple parameter combinations #####################
  # loop over conditions

  for(condition in conditions){
    paraSamples_ = 
      list("HP" = data.frame(
        alpha = seq(0.05, 0.1, length.out = nCut),
        nu = c(0.25, 0.5, 1, 2),
        tau = exp(seq(log(0.5), log(8), length.out = nCut)),
        gamma = seq(0.85, 0.99, length.out = nCut),
        prior = seq(0, 1, length.out = nCut)
      ), "LP" = data.frame(
        alpha = seq(0.05, 0.1, length.out = nCut),
        nu = c(0.25, 0.5, 1, 2),
        tau = exp(seq(log(0.5), log(8), length.out = nCut)),
        gamma = seq(0.85, 0.99, length.out = nCut),
        prior = seq(2, 6, length.out = nCut)
      ))
    # generate parameter combinations 
    paraSamples = paraSamples_[[condition]]
    paraCombs = t(matrix(rep(as.numeric(paraSamples[3, ]), nCut * nPara), nrow = nPara))
    for(pIdx in 1 : nPara){
      paraCombs[(nCut * (pIdx - 1) + 1) : (nCut * pIdx) ,pIdx] = paraSamples[,pIdx] 
    }
    nComb = nrow(paraCombs)
    # initialize outputs
    sub_auc_ = matrix(NA, nrow = nBreak, ncol = nComb)
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
      shortterm_rv_[,i] = aveRes$wait_minus_quit_[,nBreak - 1]
      shortterm_wait_[,i] = aveRes$wait_[, nBreak - 1]
      shortterm_quit_[i] = aveRes$quit_[nBreak-1]
      longterm_wait_[,i] = aveRes$wait_[, nBreak]
      longterm_rv_[,i] = aveRes$wait_minus_quit_[,nBreak]
      longterm_quit_[i] = aveRes$quit_[nBreak]
    }
    # combine data 
    shortterm_rv_df = data.frame(
      rv = as.vector(shortterm_rv_),
      time = 1 : nStep,
      parameter = factor(rep(paraLabels, each = nCut * nStep), levels = paraLabels),
      rank = rep(factor(1 : nCut), each = nStep)
    )
    longterm_rv_df = data.frame(
      rv = as.vector(longterm_rv_),
      time = 1 : nStep,
      parameter = factor(rep(paraLabels, each = nCut * nStep), levels = paraLabels),
      rank = rep(factor(1 : nCut), each = nStep)
    )
    auc_df = data.frame(
      auc = as.vector(sub_auc_),
      time = 1 : nBreak,
      parameter = factor(rep(paraLabels, each = nBreak * nCut), levels = paraLabels),
      rank = rep(factor(1 : nCut), each = nBreak)
    )
    shortterm_quit_df = data.frame(
      quit = as.vector(shortterm_quit_),
      parameter = factor(rep(paraLabels, each = nCut), levels = paraLabels),
      rank = factor(1 : nCut)
    )
    shortterm_quit_df %>% ggplot(aes(rank, quit)) + facet_grid(~parameter) + 
      geom_point()
    
    longterm_quit_df = data.frame(
      quit = as.vector(longterm_quit_),
      parameter = factor(rep(paraLabels, each = nCut), levels = paraLabels),
      rank = factor(1 : nCut)
    )
    
    longterm_quit_df %>% ggplot(aes(rank, quit)) + facet_grid(~parameter) + 
      geom_point()
    
    shortterm_rv_df %>% filter(time < nStep) %>% ggplot(aes(time, rv, color = rank)) + facet_grid(~parameter) +
      geom_line()
    
    longterm_rv_df %>% filter(time < nStep) %>% ggplot(aes(time, rv, color = rank)) + facet_grid(~parameter) +
      geom_line()
    
    auc_df %>% ggplot(aes(time, auc, color = rank)) + facet_grid(~parameter) +
      geom_line() + 
      scale_color_manual(values = cutValues) +
       ylim(c(0, min(delayMaxs)))
    

  
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





