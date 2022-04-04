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
      HP = replicate(120 * 60  / 2, drawSample("HP")),
      LP = replicate(120 * 60  / 2, drawSample("LP"))
      # HP = HP_seq,
      # LP = LP_seq
    )
  }
  
  # simulate with different parameter combinations
  modelName_ = c("QL2")
  paraLabels_ = list(
    c("alpha", "nu", "tau", "gamma", "eta")
  )
  nCut = 5
  paraSamples_options_ = list(
    list("HP" = data.frame(
      alpha = c(0.00141, 0.00522, 0.00882, 0.0336, 0.0804),
      nu = c(0.0254, 0.305, 1.26, 2.07, 2.61),
      tau = c(0.670, 1.46, 2.17, 2.85, 4.60),
      gamma = c(0.797, 0.818, 0.841, 0.885, 0.954),
      prior = c(1.11, 1.85, 2.81, 4.59, 8.92)
    ), "LP" = data.frame(
      alpha = c(0.00141, 0.00522, 0.00882, 0.0336, 0.0804),
      nu = c(0.0254, 0.305, 1.26, 2.07, 2.61),
      tau = c(0.670, 1.46, 2.17, 2.85, 4.60),
      gamma = c(0.797, 0.818, 0.841, 0.885, 0.954),
      prior = c(1.11, 1.85, 2.81, 4.59, 8.92)
    ))
  )
  
  # loop over models
  figs = vector(mode = "list", length = length(modelName_))
  for(i in 1 : length(modelName_)){
    modelName = modelName_[1]
    paraLabels = paraLabels_[[i]]
    paraSamples_ = paraSamples_options_[[i]]
    figs[[i]] = simPostHoc(modelName, paraLabels, paraSamples_, delays_)
    # figs_[[i]] =simPostHoc(modelName, paraLabels, paraSamples_, delays_)
  }
  return(figs)
}
