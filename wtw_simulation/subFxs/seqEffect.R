simExample = function(){
  set.seed(123)
  # default settings 
  modelName = "QL2"
  smallReward = 0 
  iti = 2
  
  # random seed
  set.seed(123)
  
  # create output directories
  dir.create("figures")
  dir.create("figures/simulation/")
  dir.create("figures/simulation/example")
  
  # load experiment parameters
  load('expParas.RData')
  
  # normative analysis 
  normResults = expSchematics(smallReward, iti, F)
  optimRewardRates = normResults$optimRewardRates
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # load packages and sub functions 
  library("tidyverse")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  
  
  # get the generative model 
  source(sprintf("subFxs/simModels/default/%s.R", modelName))
  simModel = get(modelName)
  paraNames = getParaNames(modelName)
  paraNames = factor(paraNames, levels = paraNames)
  nPara = length(paraNames)

  # simulation parameters
  nSeq = 8 # num of sequences
  nRep = 50 # num of repetitions 
  duration = 20 * 60
  nCut = 10 # number of measurements 
  nTrialMax = duration / iti
  
  HPSim_ = list()
  LPSim_ = list()
  HPauc_ = array(rep(NA, nSeq * nCut * nRep), dim = c(nCut, nRep, nSeq))     
  LPauc_ = array(rep(NA, nSeq * nCut * nRep), dim = c(nCut, nRep, nSeq))
  
  
  for(j in 1 : nSeq){
    HPdelays = replicate(nTrialMax, drawSample("HP"))
    LPdelays = replicate(nTrialMax, drawSample("LP"))
    
    for(i in 1 : nRep){
      # HP 
      tempt = simModel(c(0.05, 0.05, 5, 0.85, 6), "HP", duration, normResults, HPdelays)
  
      
      # tempt = simModel(c(0.1, 0.1, 2, 0.85, 3), "HP", duration, normResults)
      tempt$Qwaits_ = NULL
      tempt$Gs_ = NULL
      tempt = as.data.frame(tempt)
      HPauc_[1 : nCut, i, j] = sapply(1 : nCut, function(i) kmsc(tempt[(i-1) * 2 * 60 <= tempt$sellTime & tempt$sellTime < i * 2 * 60,], min(delayMaxs), F, kmGrid)$auc)

      
      # LP
      tempt  = simModel(c(0.05, 0.05, 5, 0.85, 6), "LP", duration, normResults, LPdelays)
      tempt$Qwaits_ = NULL;
      tempt$Gs_ = NULL
      tempt = as.data.frame(tempt)
      LPauc_[1 : nCut, i, j] = sapply(1 : nCut, function(i) kmsc(tempt[(i-1) * 2 * 60 <= tempt$sellTime & tempt$sellTime < i * 2 * 60,], min(delayMaxs), F, kmGrid)$auc)
    }
  }

  
  # plot
  HPdf = data.frame(
    mean = as.vector(apply(HPauc_, FUN = mean, MARGIN = c(1,3))),
    std = as.vector(apply(LPauc_, FUN = sd, MARGIN = c(1,3))),
    seq = rep(1:nSeq, each = nCut),
    time = rep(seq(1, 20, by = 2), nSeq),
    condition = rep("HP", nSeq * nCut)
  ) 
  LPdf = data.frame(
    mean = as.vector(apply(LPauc_, FUN = mean, MARGIN = c(1,3))),
    std = as.vector(apply(LPauc_, FUN = function(x) sd(x), MARGIN = c(1,3))),
    seq = rep(1:nSeq, each = nCut),
    time = rep(seq(1, 20, by = 2), nSeq),
    condition = rep("LP", nSeq * nCut)
  ) 
  
  rbind(HPdf, LPdf) %>% ggplot(aes(time, mean, color = factor(seq))) + geom_line() +
    geom_ribbon(aes(x = time, ymin = mean - std, ymax = mean + std, fill = factor(seq)), alpha = 0.2, color = NA) +
    facet_grid(~condition) + xlab("Task time (min)") + ylab("AUC (s)") + myTheme + 
    theme(
      legend.position = "none"
    )
  ggsave("figures/simulation/seqEffect.png", width = 6, height = 3)
}






