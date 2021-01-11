# fit a reinforcement learning model for a single participant in Rstan 

# inputs:
# thisTrialData: behavioral data for this participant
# fileName: the name of the output file
# modelName: the name of   the model 
# paraNames: parameters for the model
# model: the Bayesian model 
# config: a list containing the Rstab configuration 
modelFitSingle = function(id, thisTrialData, modelName, paraNames, model, config, outputFile, normResults){
    # load experiment paras
    load('expParas.RData')
    
    # parse the stan configuration
    nChain = config[['nChain']] # number of MCMC chains
    nIter = config[['nIter']] # number of iterations on each chain
    controlList = list(adapt_delta = config[['adapt_delta']],
                       max_treedepth = config[['max_treedepth']] )
    warningFile = config[['warningFile']] # output file for stan warnings and errors
    
    # analysis constants 
    stepSec = 1  # duration of one time step (namely one temporal state) 
    iti = 2  # duration of iti
    
    # normative analysis
    optimRewardRates = normResults$optimRewardRates
    optimWaitThresholds = normResults$optimWaitThresholds
    subjectValues = normResults$subjectValues
    
    ## ensure timeWaited = scheduledWait on rewarded trials
    thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= 0] = scheduledWait[trialEarnings!= 0]})
    # prepare inputs for fitting the model
    condition = unique(thisTrialData$condition)
    ## number of possible decision points in a trial
    delayMax = max(delayMaxs)
    tWaits = seq(0, delayMax, by = stepSec)
    ## number of wait-or-quit decision time points in a trial
    nWaitOrQuit = length(tWaits) 
    ## number of made actions in each trial 
    nMadeActions = floor(thisTrialData$timeWaited / stepSec) + 1
    ## when a trial ends 
    Ts = thisTrialData$timeWaited 
    Ts[thisTrialData$trialEarnings == 0] = floor(Ts[thisTrialData$trialEarnings == 0] / stepSec) * stepSec # around to the nearest left decision time point
    ## condition idx
    cIdxs = ifelse(thisTrialData$condition == "HP", 1, 2)
    ## the theoretic present value of the awaited reward sampled at 1 hz
    HPtheoreticValues = c(normResults$coarseSubjectValues[['HP']], rep(0, nWaitOrQuit - length(normResults$coarseSubjectValues[['HP']])))
    LPtheoreticValues = normResults$coarseSubjectValues[['LP']]
    
    ## orgianze inputs into a list
    inputs <- list(
      iti = iti,
      stepSec = stepSec,
      nWaitOrQuit = nWaitOrQuit, # number of all possible decision time points
      tWaits = tWaits, # decision time points
      N = length(thisTrialData$trialEarnings), # number of trials
      Rs = thisTrialData$trialEarnings, # rewards on each trial
      Ts = Ts, # time spent on waiting 
      nMadeActions = nMadeActions,
      HPvalues = HPtheoreticValues,
      LPvalues = LPtheoreticValues,
      cIdxs = cIdxs
    )
    if(modelName %in% c("QL1", "QL2")){
      V0_ini = mean(unlist(optimRewardRates)) * stepSec / (1 - 0.85)
      inputs$V0_ini = V0_ini
    }else{
      rewardRate_ini = mean(unlist(optimRewardRates)) * stepSec 
      inputs$rewardRate_ini = rewardRate_ini
    }
   
   # get the path in outputFile
   subName = sub(pattern = sprintf("../../genData/wtw_exp1/(exp|sim)*ModelFit(CV)*/[A-Z0-9]*/*%s/", modelName),
                      replacement = "", outputFile)
    
   # fit the model
    withCallingHandlers({
      fit = sampling(object = model, data = inputs, cores = 1, chains = nChain,
                     iter = nIter, control = controlList) 
      print(sprintf("Finish %s %s !", modelName, subName))
      write(sprintf("Finish %s %s !", modelName, subName), warningFile, append = T, sep = "\n")
    }, warning = function(w){
      write(paste(modelName, subName, w), warningFile, append = T, sep = "\n")
    })
  
  # extract posterior samples
  samples = fit %>% rstan::extract(permuted = F, pars = c(paraNames, "totalLL")) %>%
    adply(2, function(x) x) %>% dplyr::select(-chains) 
  
  # calculate WAIC and Efficient approximate leave-one-out cross-validation (LOO)
  log_lik = extract_log_lik(fit) 
  WAIC = waic(log_lik)
  LOO = loo(log_lik)
  save("WAIC", "LOO", file = sprintf("%s_waic.RData", outputFile))
  
  # summarise posterior parameters and total log likelihood
  fitSummary <- summary(fit, pars = c(paraNames, "totalLL"), use_cache = F)$summary

  # detect participants with low ESSs and high Rhats 
  ESSCols = which(str_detect(colnames(fitSummary), "Effe")) # columns recording ESSs
  if(any(fitSummary[,ESSCols] < nChain * 100)){
    write(paste(modelName, subName, "Low ESS"), warningFile, append = T, sep = "\n")
  }
  RhatCols = which(str_detect(colnames(fitSummary), "Rhat")) # columns recording ESSs
  if(any(fitSummary[,RhatCols] > 1.01)){
    write(paste(modelName, subName, "High Rhat"), warningFile, append = T, sep = "\n")
  } 
  
  # check divergent transitions
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  nDt = sum(divergent)
  fitSummary = cbind(fitSummary, nDt = rep(nDt, nrow(fitSummary)))

  # write outputs  
  write.table(fitSummary, file = sprintf("%s_summary.txt", outputFile), 
              sep = ",", col.names = F, row.names=FALSE)
}

