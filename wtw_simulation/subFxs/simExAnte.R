simExAnte = function(modelName, modelLabel, paras, delays_ = list()){
  set.seed(123)
  load("expParas.RData")
  # default settings 
  smallReward = 0 
  iti = 2
  stepSec = 1
  nHPStep = delayMaxs[1] / stepSec + 1 # since here we assume the interval is [ , ) yet in we need to include the end point at t = 20
  nLPStep = delayMaxs[2] / stepSec + 1 
  
  # random seed
  set.seed(123)
  
  # load experiment parameters
  load('expParas.RData')
  
  # normative analysis 
  source("expSchematics.R")
  normResults = expSchematics(smallReward, iti)
  optimRewardRates = normResults$optimRewardRates
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # load packages and sub functions 
  library("tidyverse")
  source(file.path("subFxs", "plotThemes.R"))
  source(file.path("subFxs", "helpFxs.R"))
  source(file.path("subFxs", "loadFxs.R")) 
  source(file.path("subFxs", "analysisFxs.R"))
  source(file.path("subFxs", "averageRes.R"))
  
  # get the generative model 
  source(file.path("subFxs", "simModels", "default", sprintf("%s.R", modelName)))
  simModel = get(modelName)

  # num of repetitions 
  exist = length(delays_) > 0
  if(exist){
    nRep = length(delays_)
  }else{
    nRep = 10
  }
  duration = 120 * 60
  
  HPSim_ = list()
  LPSim_ = list()
  HPauc_ = matrix(NA, nrow =  5, ncol = nRep)
  LPauc_ = matrix(NA, nrow =  5, ncol = nRep)
  for(i in 1 : nRep){
    if(exist){
      HPSim_[[i]] = simModel(paras$HP, "HP", duration, normResults, delays_[[i]]$HP)
      LPSim_[[i]] = simModel(paras$LP, "LP", duration, normResults, delays_[[i]]$LP)
      print(i)
    }else{
      HPSim_[[i]] = simModel(paras$HP, "HP", duration, normResults)
      LPSim_[[i]] = simModel(paras$LP, "LP", duration, normResults) 
    }
  }
  
  # 
  # average the results across simulations
  recorded_timepoints = c(0, 4 * 60, 8 * 60, 12 * 60, 16 * 60, duration)
  nRecord = length(recorded_timepoints)
  aveResHP = averageRes(HPSim_, recorded_timepoints, paras$HP[2])
  aveResLP = averageRes(LPSim_, recorded_timepoints, paras$LP[2])
  
  condData = data.frame(
    condition = conditions,
    asymAUC = c(aveResHP$sub_auc_[nRecord], aveResLP$sub_auc_[nRecord]),
    optim = unlist(normResults$optimWaitThresholds)
  )
  
  figAUC = data.frame(
    mu = c(aveResHP$sub_auc_[1 : (nRecord-1)], aveResLP$sub_auc_[1 : (nRecord - 1)]),
    condition = c(rep("HP", nRecord - 1), rep("LP", nRecord - 1)),
    t = rep(recorded_timepoints[1 : (nRecord - 1)], 2) / 60
  ) %>% 
    ggplot(aes(t, mu, color = condition)) + geom_point() + geom_line() +
    myTheme + scale_color_manual(values = conditionColors) +
    scale_x_continuous(breaks = c(recorded_timepoints[1 : (nRecord - 1)] / 60, 20),
                       labels = c(recorded_timepoints[1 : (nRecord - 1)] / 60, 120),
                       limits = c(0,22)) +
    xlab("Simulation time (min)") +
    ylab("AUC (s)") + theme(legend.position = "none")  + ylim(c(0, 20)) + 
    geom_hline(aes(yintercept = optim, color = condition), linetype = "dashed", data = condData) +
    geom_point(aes(y = asymAUC, color = condition), x = 20, shape = 8, data = condData) +
    ggtitle(modelLabel) + 
    theme(plot.title = element_text(hjust = 0.5))
  
    
  ##################################################################
  ##                     plot value functions                     ##
  ##################################################################
  

  condData = data.frame(
      rv = c(aveResHP$wait_minus_quit_[1 : (nHPStep - 1), nRecord], aveResLP$wait_minus_quit_[1 : (nLPStep - 1), nRecord]),
      t = c(1 : (nHPStep - 1) - 1, 1 : (nLPStep - 1)) - 1,
      condition = c(rep("HP", nHPStep - 1), rep("LP", nLPStep - 1))
    ) 
  figRV = data.frame(
      rvWait = c(as.vector(aveResHP$wait_minus_quit_[1 : (nHPStep - 1), 1 : (nRecord - 1)]),
                 as.vector(aveResLP$wait_minus_quit_[1 : (nLPStep - 1), 1 : (nRecord - 1)])),
      t = c(rep(1 : (nHPStep - 1) - 1, nRecord - 1), rep(1 : (nLPStep - 1) - 1, nRecord - 1)),
      taskTime = c(rep(1 : (nRecord - 1), each = nHPStep - 1), rep(1 : (nRecord - 1), each = nLPStep - 1)),
      condition = c(rep("HP", (nRecord - 1) * (nHPStep - 1)), rep("LP", (nRecord - 1) * (nLPStep - 1)))
    )  %>%
      mutate(color = factor(c(rep(1 : (nRecord - 1), each = nHPStep - 1), rep(1 : (nRecord - 1) + (nRecord - 1), each = nLPStep - 1)))) %>%
      filter( t< 20) %>% 
      ggplot(aes(t, rvWait, color = color)) + geom_point() + geom_line() +
      facet_grid(. ~ condition) +
      scale_color_manual(values = c("#c7e9c0", "#41ab5d", "#238b45",  "#006d2c", "#00441b", 
                                    "#bcbddc", "#807dba", "#6a51a3", "#3f007d", "#4d004b"))  + 
      geom_point(aes(t, rv), data = condData[condData$t < 20,], color = "black", inherit.aes = F) +
      geom_line(aes(t, rv), data = condData[condData$t < 20,], color = "black", inherit.aes = F) + 
      ylab("Relative value of waiting") + xlab("Elapsed time (s)") + myTheme +
      theme(legend.position = "None",
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(modelLabel) + scale_x_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(0, 20))

 
  
  #################################################################
  ##                       plot a snippet                        ##
  #################################################################
  if(modelName == "QL1"){
    set.seed(123)
    paras = c(0.1, 3.5, 0.70, 2)
    # paras =  c(0.01, 5, 0.85, 1.5)
    exampleSim = simModel(paras, "HP", duration, normResults, delays_[[3]]$HP)
    exampleQwaits = exampleSim$Qwaits_
    exampleGs = exampleSim$Gs_
    exampleV = exampleSim$V_
    exampleRs = exampleSim$trialEarnings
    print(exampleRs[1:4])

    exampleTs = exampleSim$timeWaited
    exampleTs[exampleRs == 10] = exampleSim$scheduledWait[exampleRs== 10]
    print(exampleTs[1:4])
    # initialize the outputs 
    figQwaits_ = vector(mode = "list", length = 4)
    figGs_ = vector(mode = "list", length = 4)
    ## plot value functions
    for(i in 1 : 4){
      if(i > 1){
        figQwaits_[[i]] = data.frame(
          time = -2 : min(delayMaxs),
          type = c(0, rep(1, min(delayMaxs) + 2)),
          value = c(exampleV[i], NA, as.vector(exampleQwaits[,i]))
        ) %>% ggplot(aes(time, value, color = factor(type))) + 
          geom_rect(xmin = -3, xmax = exampleTs[i-1], ymin = 4, ymax = 8, fill ='#fee0d2', color = NA) + 
          geom_point(size = 4) +
          myTheme + ylab('Action value') + xlab("Elapsed time t (s)") +
          ylim(c(4, 8)) + scale_x_continuous(breaks = c(-2, seq(0, min(delayMaxs), by = 5)), limits = c(-3, 21)) +
          geom_vline(xintercept = exampleTs[i-1], color =  '#fc9272', linetype = "dashed", size = 1) +
          scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
          theme(legend.position =  "None",
                panel.grid.major.y = element_line(color = "grey", size = 0.5,linetype = 2))
      }else{
        figQwaits_[[i]] = data.frame(
          time = -2 : min(delayMaxs),
          type = c(0, rep(1, min(delayMaxs) + 2)),
          Qwait = c(exampleV[i], NA, as.vector(exampleQwaits[,i]))
        ) %>% ggplot(aes(time, Qwait, color = factor(type))) + geom_point(size = 4) +
          myTheme + ylab('Action value')  + xlab("Elapsed time t (s)") +
          ylim(c(4, 8)) + scale_x_continuous(breaks = c(-2, seq(0, min(delayMaxs), by = 5)), limits = c(-3, 21)) +
          scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
          theme(legend.position =  "None",  panel.grid.major.y = element_line(color = "grey",
                                                                               size = 0.5,linetype = 2))
      }
    }
    
    ## plot Gs
    for(i in 1 : 4){
      G = exampleSim$Gs_[,i]
      G[G == 0] = NA
      exampleT = round(exampleTs[i], 1)
      figGs_[[i]] =
        data.frame(
        time = -2 : min(delayMaxs),
        G = c(G[1] * paras[3] ^ 2, NA,  G),
        type = c(0, rep(1, min(delayMaxs) + 2))
      ) %>% ggplot(aes(time, G, color = factor(type))) + geom_point(size = 4) +
        myTheme  + ylab('Feedback signal') + xlab("Elpased time t (s)") +
        scale_x_continuous(limits = c(-3, 21), breaks = c(-2, seq(0, min(delayMaxs), by = 5))) +
        ylim(0, 15) +
        geom_vline(xintercept = exampleTs[i], color =  '#fc9272', linetype = "dashed", inherit.axis = F, size = 1) +
        scale_color_manual(values = c("#d7191c", "#2c7bb6")) + theme(legend.position = "None") +
          theme(plot.title = element_text(hjust = 0.5), legend.position = "None")
    }
  }
  
  ## plot priors 
  if(modelName == 'QL1'){
    V = mean(unlist(optimRewardRates)) / (1 - 0.85) # state value for t = 0
    tWaits = seq(0, 20, by = stepSec) # decision points 
    tMax = max(tWaits) #  time point for the last decision point
    Qwaits_prior_equal_1 = -0.1 * (tWaits) + 1 + V 
    Qwaits_prior_equal_2 = -0.1 * (tWaits) + 1.5 + V 
    
    figPrior = data.frame(
      time = c(-2, tWaits, tWaits),
      type = c(0, rep(1, length(tWaits)),rep(2, length(tWaits)) ),
      action_value = c(V, Qwaits_prior_equal_1, Qwaits_prior_equal_2)
    ) %>% ggplot(aes(time, action_value, color = factor(type))) + geom_point(size = 10) +
      myTheme + ylab('Initial action value')  + xlab("Elapsed time t") +
      ylim(c(4.5, 8)) + scale_x_continuous(breaks = c(-2, seq(0, min(delayMaxs), by = 5)), limits = c(-3, 21)) +
      scale_color_manual(values = c("#d7191c", "#2c7bb6", "#9ecae1")) +
      theme(legend.position =  "None")
    

    # add new code to for test-retest studies
    figPrior_testretest_intercept = 
      data.frame(
      time = seq(-1.5, 12, by = 0.5),
      Qquit = rep(3, length = (12+1.5)/0.5 + 1),
      Qwaits_prior_equal_1 = c(NaN, NaN, NaN, -0.1 * seq(0, 12, by = 0.5) + 0.75 + 3),
      Qwaits_prior_equal_2 = c(NaN, NaN, NaN, -0.1 * seq(0, 12, by = 0.5) + 1 + 3)
      ) %>% ggplot(aes(time, Qwaits_prior_equal_1))+ 
      geom_line(aes(time, Qquit), color = "#d7191c", size = 2) +
      geom_point(color = "#9ecae1", size = 3) + 
      geom_point(aes(time, Qwaits_prior_equal_2), size = 3, color = "#2c7bb6", inherit.aes = F)  +
      theme(legend.position =  "None") + myTheme +
      scale_y_continuous(breaks = c(3, 3.75, 4), limits = c(2.5, 4)) +
      scale_x_continuous(breaks = seq(0, 12, by = 3), limits = c(-2, 13)) +
      ylab('Initial action value')  + xlab("Elapsed time t")
      
    ggsave(file.path("~/Downloads/prior_intercept.eps"), figPrior_testretest_intercept, width = 5, height = 4)
    
    figPrior_testretest_slope = 
      data.frame(
        time = seq(-1.5, 12, by = 0.5),
        Qquit = rep(3, length = (12+1.5)/0.5 + 1),
        Qwaits_prior_equal_1 = c(NaN, NaN, NaN, -0.2 * seq(0, 12, by = 0.5) + 1 + 3),
        Qwaits_prior_equal_2 = c(NaN, NaN, NaN, -0.1 * seq(0, 12, by = 0.5) + 1 + 3)
      ) %>% ggplot(aes(time, Qwaits_prior_equal_1))+ 
      geom_line(aes(time, Qquit), color = "#d7191c", size = 2) +
      geom_point(color = "#9ecae1") + 
      geom_point(aes(time, Qwaits_prior_equal_2), color = "#2c7bb6", inherit.aes = F)  +
      theme(legend.position =  "None") + myTheme +
      scale_y_continuous(breaks = c(3, 3.75, 4), limits = c(2.5, 4)) + 
      scale_x_continuous(breaks = seq(0, 12, by = 3), limits = c(-2, 13)) +
      ylab('Initial action value')  + xlab("Elapsed time t")
    
    ggsave(file.path("~/Downloads/prior_slope.eps"), figPrior_testretest_slope, width = 5, height = 4)
  }
  ############ return outputs ###########
  outputs = list(
    "learn" = figAUC
  )
  
  if(modelName %in% c("QL1", "RL1", "omni")){
    outputs[["rv"]] = figRV
  }
  
  if(modelName == "QL1"){
    outputs[['values_']]  = figQwaits_
    outputs[['Gs_']] = figGs_
    outputs[['prior']] = figPrior
  }

  return(outputs)
}





