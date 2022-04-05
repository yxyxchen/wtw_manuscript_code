paraEffect = function(){
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
  
  # fix the delay sequence 
  nRep = 10
  set.seed(123)
  
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
  
  # simulate with different parameter combinations
  modelName_ = c("QL2")
  paraLabels_ = list(
    c("alpha", "nu", "tau", "gamma", "eta")
  )
  nCut = 3
  paraSamples_options_ = list(
    list("HP" = data.frame(
      alpha = c(0.005, 0.01, 0.02),
      nu = c(0.33, 1, 3),
      tau = c(1, 3, 5),
      gamma = c(0.8, 0.85, 0.9),
      prior = c(1, 2, 4)
    ), "LP" = data.frame(
      alpha = c(0.005, 0.01, 0.02),
      nu = c(0.33, 1, 3),
      tau = c(1, 3, 5),
      gamma = c(0.8, 0.85, 0.9),
      prior = c(1, 2, 4)
    ))
  )
  
  # loop over models
  figs = vector(mode = "list", length = length(modelName_))
  for(i in 1 : length(modelName_)){
    modelName = modelName_[1]
    paraLabels = paraLabels_[[i]]
    paraSamples_ = paraSamples_options_[[i]]
    # figs[[i]] = simPostHoc(modelName, paraLabels, paraSamples_)
    figs_[[i]] =simPostHoc(modelName, paraLabels, paraSamples_, delays_)
  }
  return(figs)
}
