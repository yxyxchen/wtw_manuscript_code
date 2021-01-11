# fit a model for multiple participants in Rstan 

# inputs:
# modelName : the model name
# trialData : a nSubx1 list, each element containing behavioral data for one participant
# config : a list containing the Rstab configuration 
# outputDir: the directory to save parameter estimations

# the config variable contains:
# nChain : number of chains 
# nIter : number of interations on each chain 
# adapt_delta: real number from 0 to 1, and increaseing it forces stan to explore the target distribution in a smaller step
# max_treedepth: maximal depth of the trees that stan evaluates during each iteration
# warningFile : file for saving warnings generated Rstan

modelFitGroup = function(modelName, trialData, config, outputDir, isTrct = T){
  # create the output directory 
  dir.create(outputDir)
 
   # create the file for Rstan warnings and erros
  writeLines("", config[['warningFile']])
  
  #load libraries
  library('plyr'); 
  library('rstan');library("loo");library("coda") 
  library("doMC");library("foreach")
  source('subFxs/modelFitSingle.R') # for fitting each single participant
  load("expParas.RData")
  
  
  # constants 
  iti = 2
  stepSec = 1
  normResults = expSchematics(0, iti, F)
  # compile the Rstan model 
  options(warn= 1) 
  Sys.setenv(USE_CXX14=1) # settings for the local PC
  rstan_options(auto_write = TRUE) 
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  modelNames = c("QL1", "QL2", "RL1", "RL2", "naive", "omni")
  for(i in 1 : length(modelNames)){
    modelName = modelNames[i]
    model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
    outputDir = sprintf("../../genData/wtw_exp3/expModelFit/%s", modelName)
    dir.create(outputDir)
    paraNames = getParaNames(modelName)
    outputFile = sprintf("%s/s%s", outputDir, id)
    modelFitSingle(id, thisTrialData, modelName, paraNames, model, config, outputFile, normResults)
  }

  
  # determine parameters 
  paraNames = getParaNames(modelName)
  
  # load expData
  ids = names(trialData)
  nSub = length(ids)                    
  
  # parallel compuation settings
  nCore = as.numeric(Sys.getenv("NSLOTS")) # settings for SCC
  if(is.na(nCore)) nCore = 1 # settings for SCC
  # nCore = parallel::detectCores() -1 # settings for the local PC
  # registerDoMC(nCore) # settings for the local PC
  
  foreach(i = 1 : nSub) %dopar% {
      id = ids[[i]]
      thisTrialData = trialData[[id]]
      # truncate the last portion in each block 
      if(isTrct){
        excludedTrials = which(thisTrialData$trialStartTime > (blockSec - max(delayMaxs)))
        thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excludedTrials & thisTrialData$blockNum <= 2,]
      }
      outputFile = sprintf("%s/s%s", outputDir, id)
      modelFitSingle(id, thisTrialData, modelName, paraNames, model, config, outputFile, normResults)
  }
}
