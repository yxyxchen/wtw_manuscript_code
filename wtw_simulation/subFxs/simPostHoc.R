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
    "#737373",
    "#525252",
    "#252525"
  )
  nCut = length(cutValues)
  
  ######################### simulate with multiple parameter combinations #####################
  # loop over conditions
  auc_df_ = list()
  asym_df_ = list()
  for(condition in conditions){
    # default_paras = c(0.01, 1, 5, 0.85, 1)
    # default_paras = c(0.01, 1, 5, 0.85, 4)
    # generate parameter combinations 
    paraSamples = paraSamples_[[condition]]
    default_paras = as.numeric(paraSamples[3, ])
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
    # shortterm_rv_ = matrix(NA, nrow = nStep, ncol = nComb)
    # shortterm_quit_ =  vector( length = nComb)
    # shortterm_wait_ = matrix(NA, nrow = nStep, ncol = nComb)
    # longterm_rv_ = matrix(NA, nrow = nStep, ncol = nComb)
    # longterm_quit_ =  vector( length = nComb)
    # longterm_wait_ = matrix(NA, nrow = nStep, ncol = nComb)
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
      aveRes = averageRes(sim_, recorded_timepoints, paras[3])
      sub_auc_[,i] = aveRes$sub_auc_
    }
    # plot
    auc_df_[[condition]] = data.frame(
      auc = as.vector(sub_auc_[1 : (nRecord - 1),]),
      time = recorded_timepoints[1 : (nRecord - 1)] / 60,
      parameter = factor(rep(paraLabels, each = (nRecord - 1) * nCut), levels = paraLabels),
      rank = rep(factor(1 : nCut), each = (nRecord - 1)),
      condition = condition
    )
    asym_df_[[condition]] = data.frame(
      auc = as.vector(sub_auc_[nRecord,]),
      time = 20,
      parameter = factor(rep(paraLabels, each = nCut), levels = paraLabels),
      rank = factor(1 : nCut),
      condition = condition
    )
    # I understand that the behavior can look wierd ...
    # auc_df = data.frame(
    #   auc = as.vector(sub_auc_),
    #   time = 1 : nRecord,
    #   parameter = factor(rep(paraLabels, each = nRecord * nCut), levels = paraLabels),
    #   rank = rep(factor(1 : nCut), each = nRecord)
    # )
    # auc_df %>% ggplot(aes(time, auc, color = rank)) +
    #   facet_grid(~parameter) + geom_line() + 
    #   geom_point(aes(time, auc, color = rank), data = asym_df, inherit.aes = F, shape=8) +
    #   scale_color_manual(values = cutValues) +
    #   ylim(c(0, min(delayMaxs))) + myTheme 
  }
  auc_df = rbind(auc_df_[["HP"]], auc_df_[["LP"]])
  asym_df = rbind(asym_df_[["HP"]], asym_df_[["LP"]])
  
  auc_df %>% ggplot(aes(time, auc, color = rank)) +
    facet_grid(condition~parameter) + geom_line() + 
    geom_point() +
    scale_color_manual(values = cutValues) +
    ylim(c(0, min(delayMaxs))) + myTheme +
    xlab("Task time (min)") + ylab("AUC (s)") +
    scale_x_continuous(breaks = recorded_timepoints[1 : (nRecord - 1)] / 60) + 
    geom_point(aes(x = time, y = auc, color = rank), data = asym_df, inherit.aes = F, shape = 8)
    
  return(fig)
  
  
  
}





