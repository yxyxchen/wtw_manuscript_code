modelFitting = function(thisTrialData, fileName, paraNames, model, modelName){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 50
  
  # determine wIni, the first change
  subOptimalRatio = 0.9 
  if(any(paraNames  == "gamma") || modelName == "BL" ){
    wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9) * subOptimalRatio
  }else{
    wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings !=0 ] = scheduledWait[trialEarnings != 0] # third change
  tMax = max(tMaxs) # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    # real data
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    #
                    iti = iti,
                    stepDuration = stepDuration)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(paraNames, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    dplyr::select(-chains) 
  write.table(matrix(unlist(tempt), ncol = length(paraNames) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE) 
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming
  WAIC = waic(log_lik)
  looStat = loo(log_lik)
  save("WAIC", "looStat", file = sprintf("%s_waic.RData", fileName))
  fitSummary <- summary(fit,pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
}

modelFittingCV = function(thisTrialData, fileName, paraNames, model, modelName){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni, the first change
  subOptimalRatio = 0.9 
  if(any(paraNames  == "gamma") || modelName == "BL" ){
    wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9) * subOptimalRatio
  }else{
    wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  }
  
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  tMax = max(tMaxs) # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    # real data
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    #
                    iti = iti,
                    stepDuration = stepDuration)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  fitSummary <- summary(fit,pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
}

modelFittingdb = function(thisTrialData, fileName, paraNames, model, modelName,nPara, low, up){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni, the first change
  subOptimalRatio = 0.9 
  if(any(paraNames  == "gamma") || modelName == "BL" ){
    wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9) * subOptimalRatio
  }else{
    wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  tMax = max(tMaxs) # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    nPara = nPara,
                    # real data
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    low = low,
                    up = up,
                    #
                    iti = iti,
                    stepDuration = stepDuration)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(paraNames, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    dplyr::select(-chains) 
  write.table(matrix(unlist(tempt), ncol = length(paraNames) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE) 
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming
  WAIC = waic(log_lik)
  looStat = loo(log_lik)
  save("WAIC", "looStat", file = sprintf("%s_waic.RData", fileName))
  fitSummary <- summary(fit,pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
  
  # detmerine converge
  converge = all(fitSummary[,"Rhat"] < 1.1) & all(fitSummary[, "n_eff"] >100)
  return(converge)  
}


modelFittingCVdb = function(thisTrialData, fileName, paraNames, model, modelName,nPara, low, up){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 100
  
  # determine wIni, the first change
  subOptimalRatio = 0.9 
  if(any(paraNames  == "gamma") || modelName == "BL" ){
    wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9) * subOptimalRatio
  }else{
    wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  }
  
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  tMax = max(tMaxs) # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    nPara = nPara,
                    # real data
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    low = low,
                    up = up,
                    #
                    iti = iti,
                    stepDuration = stepDuration)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  fitSummary <- summary(fit,pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
  
  # detmerine converge
  converge = all(fitSummary[,"Rhat"] < 1.1) & all(fitSummary[, "n_eff"] >100)
  return(converge)  
}