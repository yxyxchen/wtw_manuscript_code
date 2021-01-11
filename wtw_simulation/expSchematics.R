# this script plots delay distributions and reward rates in two environments
expSchematics = function(smallReward, iti){
  # load libararies
  library("tidyverse")
  library(latex2exp)
  source('subFxs/plotThemes.R')
  
  # load experiment parameters 
  load("expParas.RData")
  
  # for display purposes, all variables on the continous time scale
  # are discretized into 0.1 second time bins
  bin = 0.1 # width of a time bin
  smallReward = 0
  time = list(
    HP = seq(bin, delayMaxs[1], by = bin),
    LP = seq(bin, delayMaxs[2], by = bin)
  ) 
  
  # delay CDFS
  ### delays in HP follow unif(0, 20)
  ### delays in LP follow pareto(mu = 0, k =4, sigma = 2), truncated at 40
  HP = 1 / length(time[['HP']]) * (1 : length(time[['HP']]))
  mu = 0; k = 4; sigma = 2; # parameters for the pareto distribution
  LP = 1 - (1 + k * (time[['LP']]- mu) / sigma) ^ (-1 / k)
  LP[length(LP)] = 1 # truncated at 40
  rewardDelayCDFs = list(
    HP = HP,
    LP = LP
  )
  
  # delay PDFs
  HP = diff(c(0, rewardDelayCDFs$HP))
  LP = diff(c(0, rewardDelayCDFs$LP))
  rewardDelayPDFs = list(
    "HP" = HP,
    "LP" = LP
  )
  
  # average waiting durations given different policies
  # Here we assume rewards occur at the middle of each time bin
  HP = cumsum((time[['HP']] - 0.5 * bin) * rewardDelayPDFs$HP) / cumsum(rewardDelayPDFs$HP)
  LP = cumsum((time[['LP']] - 0.5 * bin) * rewardDelayPDFs$LP) / cumsum(rewardDelayPDFs$LP)
  meanRewardDelays = list('HP' = HP, 'LP' = LP)
  
  # rewardRates given different policies
  ## might be different from the values used in expParas.R, 
  ## which are calcuated with a higher temporal resoluation
  HP = (tokenValue * rewardDelayCDFs$HP + smallReward * (1 - rewardDelayCDFs$HP)) /
    ((meanRewardDelays$HP * rewardDelayCDFs$HP + time[['HP']] * (1 - rewardDelayCDFs$HP)) + iti)
  LP = (tokenValue * rewardDelayCDFs$LP + smallReward * (1 - rewardDelayCDFs$LP))/
    ((meanRewardDelays$LP * rewardDelayCDFs$LP + time[['LP']] * (1 - rewardDelayCDFs$LP)) + iti)
  rewardRates = list('HP' = HP, 'LP' = LP)
  
  # optimal raward rates and optimal policies
  optimWaitThresholds = list()
  optimWaitThresholds$HP = time$HP[which.max(HP)]
  optimWaitThresholds$LP = time$LP[which.max(LP)]
  optimRewardRates = list()
  optimRewardRates$HP = max(HP)
  optimRewardRates$LP = max(LP)
  
  # calculate subjective value as a function of elapsed time 
  subjectValues = list()
  for(cIdx in 1 : 2){
    condition = conditions[cIdx]
    delayMax = delayMaxs[cIdx]
    pdf = rewardDelayPDFs[[cIdx]]
    cdf = rewardDelayCDFs[[cIdx]]
    thisTime = time[[cIdx]]
    Rstar = optimRewardRates[[cIdx]]
    ts = seq(0, delayMax, by = 0.1) # elapsed time
    
    # initialize 
    thisSubjectValues = vector(length = length(ts))
    Tstars = vector(length = length(ts))
    # loop over different elapsed time
    for(i in 1 : length(ts)){
      t = ts[i] # this elapsed time
      trctTime = thisTime[thisTime > t]
      
      if(t == delayMax){
        Tstar = t
        gt_max = tokenValue 
      }else{
        # loop over different waiting policies 
        Tstar = t
        gt_max = -100
        for(T in seq(t, delayMax, by = 0.1)){
          trctPDF = pdf[thisTime > t] / sum(pdf[thisTime > t])
          at = tokenValue * sum(trctPDF[trctTime <= T])
          trctPDF[trctTime == T] = trctPDF[trctTime == T]
          bt = sum((trctTime[trctTime <= T] - 0.5 * bin - t) * trctPDF[trctTime <= T]) + 
            + (T - t) * sum(trctPDF[trctTime > T]) 
          gt = at - bt * Rstar
          if(gt > gt_max){
            gt_max = gt
            Tstar = T
          }
        }
      }
      
      thisSubjectValues[i] = gt_max 
      Tstars[i] = Tstar
    }
    subjectValues[[condition]] = thisSubjectValues
  }
  
  # return outputs 
  outputs = list(
    "optimWaitThresholds" = optimWaitThresholds,
    "optimRewardRates" = optimRewardRates,
    "rewardRates" = rewardRates,
    "subjectValues" = subjectValues,
    "time" = time
  )
  return(outputs)
}


  


