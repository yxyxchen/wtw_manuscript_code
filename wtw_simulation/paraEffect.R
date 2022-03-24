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
  delays_ = vector(mode = "list", length = nRep)
  for(i in 1 : nRep){
    delays_[[i]] = list(
      HP = replicate(240 * 60  / 2, drawSample("HP")),
      LP = replicate(240 * 60  / 2, drawSample("LP"))
    )
  }
  
  # simulate with different parameter combinations
  modelName_ = c("QL2", "RL2")
  paraLabels_ = list(
    c("alpha", "nu", "tau", "gamma", "eta"),
    c("alpha", "nu", "tau", "eta", "beta")
  )
  nCut = 4
  # paraSamples_ = list(
  #   list("HP" = data.frame(
  #     alpha = seq(0.01, 0.05, length.out = nCut),
  #     nu = seq(0.03, 3, length.out = nCut),
  #     tau = exp(seq(log(0.5), log(8), length.out = nCut)),
  #     gamma = seq(0.85, 0.95, length.out = nCut),
  #     prior = seq(0, 1, length.out = nCut)
  #   ), "LP" = data.frame(
  #     alpha = seq(0.01, 0.05, length.out = nCut),
  #     nu = seq(0.03, 3, length.out = nCut),
  #     tau = exp(seq(log(0.5), log(8), length.out = nCut)),
  #     gamma = seq(0.85, 0.95, length.out = nCut),
  #     prior = seq(2, 6, length.out = nCut)
  #   )),
  #   list("HP" = data.frame(
  #     alpha = seq(0.01, 0.05, length.out = nCut),
  #     nu = seq(0.03, 3, length.out = nCut),
  #     tau = exp(seq(log(0.5), log(8), length.out = nCut)),
  #     prior = seq(0, 1, length.out = nCut),
  #     beta = seq(0.0105, 0.025, length.out = nCut)
  #   ), "LP" = data.frame(
  #     alpha = seq(0.01, 0.05, length.out = nCut),
  #     nu = seq(0.03, 3, length.out = nCut),
  #     tau = exp(seq(log(0.5), log(8), length.out = nCut)),
  #     prior = seq(2, 6, length.out = nCut),
  #     beta = seq(0.0105, 0.025, length.out = nCut)
  #   ))
  # )
  
  paraSamples_ = list(
    list("HP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      nu = c(0.25, 0.5, 1, 2),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      gamma = seq(0.85, 0.95, length.out = nCut),
      prior = seq(0, 1, length.out = nCut)
    ), "LP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      nu = c(0.25, 0.5, 1, 2),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      gamma = seq(0.85, 0.95, length.out = nCut),
      prior = seq(2, 6, length.out = nCut)
    )),
    list("HP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      nu = seq(0.03, 3, length.out = nCut),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      prior = seq(0, 1, length.out = nCut),
      beta = seq(0.0105, 0.025, length.out = nCut)
    ), "LP" = data.frame(
      alpha = seq(0.01, 0.05, length.out = nCut),
      nu = seq(0.03, 3, length.out = nCut),
      tau = exp(seq(log(0.5), log(8), length.out = nCut)),
      prior = seq(2, 6, length.out = nCut),
      beta = seq(0.0105, 0.025, length.out = nCut)
    ))
  )
  # loop over models
  figs_ = vector(mode = "list", length = length(modelName_))
  for(i in 1 : length(modelName_)){
    figs_[[i]] =simPostHoc(modelName_[i], paraLabels_[[i]], paraSamples_[[i]], delays_)
  }
  
  return(figs_)
}
