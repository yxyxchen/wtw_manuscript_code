# fit a reinforcement learning model for a single participant in Rstan 

# inputs:
# thisTrialData: behavioral data for this participant
# fileName: the name of the output file
# modelName: the name of   the model 
# paraNames: parameters for the model
# model: the Bayesian model 
# config: a list containing the Rstab configuration 
modelFitSingle = function(id, thisTrialData, modelName, paraNames, model, config, outputFile){
    # load experiment paras
    load('expParas.RData')
  
    # parse the stan configuration
    nChain = config[['nChain']] # number of MCMC chains
    nIter = config[['nIter']] # number of iterations on each chain
    controlList = list(adapt_delta = config[['adapt_delta']],
                       max_treedepth = config[['max_treedepth']] )
    warningFile = config[['warningFile']] # output file for stan warnings and errors
    
    # duration of one time step (namely one temporal state) 
    stepSec = 1
    
    # prepare inputs for fitting the model
    condition = unique(thisTrialData$condition)
    ## maximal number of steps in a trial
    nStepMax = max(tMaxs) / stepSec
    ## ensure timeWaited = scheduledWait on rewarded trials
    thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= 0] = scheduledWait[trialEarnings!= 0]})
    ## terminal state in each trial
    ## Noticeably, ceiling(timeWaited / stepSec) gives the num of steps
    ## and adding 1 upon it gives the terminal state, namely the state following the final action. 
    Ts = with(thisTrialData, {round(ceiling(timeWaited / stepSec) + 1)}) 
    ## orgianze inputs into a list
    inputs <- list(
      iti = iti,
      stepSec = stepSec,
      nStepMax = nStepMax,
      N = length(thisTrialData$trialEarnings), # number of trials
      Rs = thisTrialData$trialEarnings, # rewards on each trial
      Ts = Ts)
    if(modelName %in% c("QL1", "QL2", "QL1_prime", "QL2_prime")){
      ## in Q-learning, the initial value of the iti state is proportional to 
      ## the discounted total rewards averaged across two conditions
      ## discount factor for one step is 0.85
      VitiIni = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
      inputs$VitiIni = VitiIni
    }else{
      ## in R-learning, the initial reward rate is proportional to
      ## the optimal reward rates averaged across two conditions
      reRateIni = 0.9 * mean(unlist(optimRewardRates)) * stepSec;
      inputs$reRateIni = reRateIni     
    }
   
   # strip the path in outputFile
   outputFile_clean = sub(pattern = sprintf("../../genData/wtw_exp2/(exp|sim)*ModelFit(CV)*/[A-Z0-9]*/*%s/", modelName),
                      replacement = "", outputFile)
  
   # fit the model
    withCallingHandlers({
      fit = sampling(object = model, data = inputs, cores = 1, chains = nChain,
                     iter = nIter, control = controlList) 
      print(sprintf("Finish %s !", outputFile_clean))
      write(sprintf("Finish %s !", outputFile_clean), warningFile, append = T, sep = "\n")
    }, warning = function(w){
      warnText = paste(modelName, outputFile_clean, w)
      write(warnText, warningFile, append = T, sep = "\n")
    })
  
  # extract posterior samples
  samples = fit %>%
    rstan::extract(permuted = F, pars = c(paraNames, "LL_all"))
  # save posterior samples
  samples = samples %>% adply(2, function(x) x) %>% dplyr::select(-chains) 
  write.table(matrix(unlist(samples), ncol = length(paraNames) + 1), file = sprintf("%s.txt", outputFile), sep = ",",
              col.names = F, row.names=FALSE) 
  
  # calculate WAIC and Efficient approximate leave-one-out cross-validation (LOO)
  log_lik = extract_log_lik(fit) 
  WAIC = waic(log_lik)
  LOO = loo(log_lik)
  save("WAIC", "LOO", file = sprintf("%s_waic.RData", outputFile))
  
  # summarise posterior parameters and LL_all
  fitSummary <- summary(fit, pars = c(paraNames, "LL_all"), use_cache = F)$summary
  
  # check ESS and Rhat
  # detect participants with low ESSs and high Rhats 
  ESSCols = which(str_detect(colnames(fitSummary), "Effe")) # columns recording ESSs
  if(any(fitSummary[,ESSCols] < nChain * 100)){
    warnText = paste(modelName, id, "Low ESS")
    write(warnText, warningFile, append = T, sep = "\n")
  }
  RhatCols = which(str_detect(colnames(fitSummary), "Rhat")) # columns recording ESSs
  if(any(fitSummary[,RhatCols] > 1.01)){
    warnText = paste(modelName, id, "High Rhat")
    write(warnText, warningFile, append = T, sep = "\n")
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

