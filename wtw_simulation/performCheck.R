performCheck = function(){
  # libraries 
  source("expSchematics.R")
  source(file.path("subFxs", "simExAnte.R"))
  source(file.path("subFxs", "simPostHoc.R"))
  source(file.path("subFxs", "plotThemes.R"))
  source(file.path("subFxs", "helpFxs.R"))
  source(file.path("subFxs", "taskFxs.R"))
  library(ggpubr)
  library("latex2exp")
  library(patchwork)
  
  # load expParas
  load("expParas.RData")
  
  # model names
  modelName_ = c("QL1", "RL1")
  modelLabel_ = c("QL1", "RL1")
  
  # it is good to use different seqs acoss simulations yet keep them the same across models
  nSim = 10
  set.seed(123)
  delays_ = vector(mode = "list", length = nSim)
  for(i in 1 : nSim){
    delays_[[i]] = list(
      HP = replicate(120 * 60  / 2, drawSample("HP")),
      LP = replicate(120 * 60  / 2, drawSample("LP"))
    )
  }
  
  # simulation parameters 
  paras_ = list(
    list("HP" = c(0.05, 5, 0.85, 4), "LP" = c(0.05, 5, 0.85, 4)), 
    list("HP" = c(0.05, 5, 4, 0.005), "LP" = c(0.05, 5, 4, 0.005))
  )
  
  # simulate and plot 
  figs_ = vector(mode = "list", length = length(modelName_))
  for(i in 1 : length(modelName_)){
    modelName = modelName_[i]
    modelLabel = modelLabel_[i]
    paras = paras_[[i]]
    outs = simExAnte(modelName_[i], modelLabel_[i], paras_[[i]], delays_)
    figs_[[i]] = outs
  }
  
  return(figs_)

}

