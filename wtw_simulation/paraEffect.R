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
  
  # simulate with different parameter combinations
  modelName_ = c("QL2", "RL2")
  paraLabels_ = list(
    c("alpha", "rho", "tau", "gamma", "eta"),
    c("alpha", "rho", "tau", "eta", "beta")
  )
  nCut = 4
  paraSamples_ = list(
    list("HP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      rho = seq(0.03, 3, length.out = nCut),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      gamma = seq(0.85, 0.95, length.out = nCut),
      prior = seq(0, 1, length.out = nCut)
    ), "LP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      rho = seq(0.03, 3, length.out = nCut),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      gamma = seq(0.85, 0.95, length.out = nCut),
      prior = seq(2, 6, length.out = nCut)
    )),
    list("HP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      rho = seq(0.03, 3, length.out = nCut),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      prior = seq(0, 1, length.out = nCut),
      beta = seq(0.0105, 0.025, length.out = nCut)
    ), "LP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      rho = seq(0.03, 3, length.out = nCut),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      prior = seq(2, 6, length.out = nCut),
      beta = seq(0.0105, 0.025, length.out = nCut)
    ))
  )
  
  figs_ = vector(mode = "list", length = length(modelName_))
  for(i in 1 : length(modelName_)){
    figs_[[i]] =simPostHoc(modelName_[i], paraLabels_[[i]], paraSamples_[[i]])
  }
  
  return(figs_)
}
