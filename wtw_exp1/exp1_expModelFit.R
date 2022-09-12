expModelFit = function(modelName, isFirstFit, batchIdx = NULL, parallel = FALSE){
  # generate output directories
  dir.create("../../genData/wtw_exp1")
  dir.create("../../genData/wtw_exp1/expModelFit")
  dir.create("stanWarnings")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load sub-functions and packages
  library("dplyr"); library("tidyr")
  source("subFxs/loadFxs.R")
  source("subFxs/helpFxs.R")
  source('subFxs/modelFitGroup.R')
  source("exp1_expSchematics.R")
  
  # prepare inputs
  allData = loadAllData()
  hdrData = allData$hdrData
  trialData = allData$trialData
  trialData = trialData[hdrData$stress == "no_stress"]
  outputDir = sprintf("../../genData/wtw_exp1/expModelFit/%s", modelName)

  if(isFirstFit){
    config = list(
      nChain = 4,
      nIter = 8000,
      adapt_delta = 0.99,
      max_treedepth = 11,
      warningFile = sprintf("stanWarnings/exp_%s.txt", modelName)
    )
    # divide data into small batches if batchIdx exists 
    if(!is.null(batchIdx)){
      if(batchIdx == 1){
        trialData = trialData[1 : 30]
      }else if(batchIdx == 2){
        trialData = trialData[31 : 60]
      }
    }
  }
  # if it is the first time to fit the model, fit all participants
  # otherwise, check model fitting results and refit those that fail any of the following criteria
  ## no divergent transitions 
  ## Rhat < 1.01 
  ## Effective Sample Size (ESS) > nChain * 100
  if(!isFirstFit){
    ids = names(trialData)
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, outputDir)
    passCheck = checkFit(paraNames, expPara)
    trialData = trialData[!passCheck]
    
    # increase the num of Iterations 
    config = list(
      nChain = 4,
      nIter = 20000,
      adapt_delta = 0.99,
      max_treedepth = 11,
      warningFile = sprintf("stanWarnings/exp_refit_%s.txt", modelName)
    )
  }

  # fit the model for all participants
  modelFitGroup(modelName, trialData, config, outputDir, parallel = parallel, isTrct = T)
}

if (sys.nframe() == 0){
  args = commandArgs(trailingOnly = T)
  expModelFit(args[1], as.logical(args[2]), as.numeric(args[3]))
}

