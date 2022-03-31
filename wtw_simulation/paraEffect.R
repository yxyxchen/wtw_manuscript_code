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
  # HP_seq = replicate(240 * 60  / 2, drawSample("HP"))
  # LP_seq = replicate(240 * 60  / 2, drawSample("LP"))
  for(i in 1 : nRep){
    delays_[[i]] = list(
      HP = replicate(240 * 60  / 2, drawSample("HP")),
      LP = replicate(240 * 60  / 2, drawSample("LP")),
      HP = HP_seq,
      LP = LP_seq
    )
  }
  
  # simulate with different parameter combinations
  modelName_ = c("QL2")
  paraLabels_ = list(
    c("alpha", "nu", "tau", "gamma", "eta")
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
  
  paraSamples_options_ = list(
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
  )
  
  # loop over models
  figs_ = vector(mode = "list", length = length(modelName_))
  for(i in 1 : length(modelName_)){
    modelName = modelName_[1]
    paraLabels = paraLabels_[[i]]
    paraSamples_ = paraSamples_options_[[i]]
    simPostHoc(modelName, paraLabels, paraSamples_, delays_)
    # figs_[[i]] =simPostHoc(modelName, paraLabels, paraSamples_, delays_)
  }
  
  return(figs_)
}
