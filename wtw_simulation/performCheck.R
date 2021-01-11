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
  modelName_ = c("QL1", "RL1", "omni", "naive", "naive")
  modelLabel_ = c("QL1", "RL1", "omni", "naive75", "naive1")
  
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
    list("HP" = c(0.05, 5, 0.85, 6), "LP" = c(0.05, 5, 0.85, 6)), # "QL2"
    list("HP" = c(0.05, 5, 10, 0.04), "LP" = c(0.05, 5, 10, 0.04)), # "RL2"
    list("HP" = c(5), "LP" = c(5)), # omni
    list("HP" = c(0.75), "LP" = c(0.75)), # naive p = 0.75
    list("HP" = c(1), "LP" = c(1)) # naive p = 0.75
  )
  
  modelName = modelName_[1]
  modelLabel = modelLabel_[1]
  paras = paras_[[1]]
  
  # simulate and plot 
  figs_ = vector(mode = "list", length = length(modelName_))
  for(i in 1 : length(modelName_)){
    outs = simExAnte(modelName_[i], modelLabel_[i], paras_[[i]], delays_)
    figs_[[i]] = outs
  }
  
  return(figs_)
  
}

