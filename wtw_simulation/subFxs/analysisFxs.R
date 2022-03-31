# plots a single subject's trial-by-trial data
trialPlots <- function(thisTrialData) {
  source('./subFxs/plotThemes.R')
  nBlock = length(unique(thisTrialData$blockNum))
  # num of trials in each block
  nTrials = sapply(1:nBlock, function(i) sum(thisTrialData$blockNum == i))
  # num of trials accumulated across blocks
  ac_nTrials = cumsum(nTrials)
  # ensure timeWaited = scheduledWait on rewarded trials
  thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= 0] = scheduledWait[trialEarnings!= 0]})
  # plot
  ## blue : timeWaited on rewarded trials
  ## red : timeWaited on non-rewarded trials
  ## black : scheduledTime on non-rewarded trials
  # p = thisTrialData %>%
  #   ggplot(aes(trialNum, timeWaited,color = factor(trialEarnings))) +
  #   geom_point() + geom_line() + scale_color_manual(values = c("grey", "#fed976")) +
  #   geom_point(data = thisTrialData[thisTrialData$trialEarnings == 0, ],
  #              aes(trialNum, scheduledWait),
  #              color = 'black') + 
  #   geom_vline(xintercept = head(ac_nTrials, nBlock - 1),
  #              color = "grey", linetype = "dashed", size = 2) +
  #   xlab("Trial num") + ylab("Time (s)") + 
  #   myTheme 
  p = thisTrialData %>%
    ggplot(aes(trialNum, timeWaited,color = factor(trialEarnings))) +
    geom_point(size = 4) + geom_line(size = 1.5) + scale_color_manual(values = c("#737373", "#fb6a4a")) +
    geom_point(data = thisTrialData[thisTrialData$trialEarnings == 0, ],
               aes(trialNum, scheduledWait),
               color = 'black', size = 4) +
    xlab("Trial") + ylab("Time (s)") + 
    myTheme 
  print(p)
  return(p)
}

model_based_kmsc = function(wait_minus_quit, tau, tMax){
  nStep = length(wait_minus_quit)
  pWait =  1 /(1  + exp(-wait_minus_quit* tau))
  pSurvival = cumprod(pWait)
  # calc AUC
  kmT = seq(0, nStep - 1)
  kmF = c(1, pSurvival[1: (nStep-1)])
  keepIdx = kmT<=tMax
  kmT <- kmT[keepIdx]
  kmF <- kmF[keepIdx]
  # calculate auc
  auc <- sum(diff(kmT) * head(kmF,-1))
  return(list(kmT=kmT, kmF=kmF, auc=auc))
}
# using kaplan-meier survival analysis to measure:
# average WTW (area under the curve)
# std WTW
kmsc <- function(thisTrialData, tMax, plotKMSC=FALSE, grid) {
  library(survival)
  # ensure timeWaited = scheduledWait on rewarded trials
  thisTrialData = within(thisTrialData, {timeWaited[trialEarnings == tokenValue] = scheduledWait[trialEarnings == tokenValue]})
  # fit a kaplan-meier survival curve
  kmfit = survfit(Surv(thisTrialData$timeWaited, (thisTrialData$trialEarnings ==0), type='right') ~ 1, 
            type='kaplan-meier', conf.type='none', start.time=0, se.fit=FALSE)
  # extract elements of the survival curve object 
  kmT = kmfit$time # time value 
  kmF = kmfit$surv # function value
  # add a point at zero, since "kaplan-meier" starts from the first event
  kmT = c(0, kmT) 
  kmF = c(1, kmF)
  # keep only points up through tMax 
  keepIdx = kmT<=tMax
  kmT <- kmT[keepIdx]
  kmF <- kmF[keepIdx]
  # extend the last value to exactly tMax
  # notice that kmT is not evenly spaced
  kmT <- c(kmT, tMax)
  kmF <- c(kmF, tail(kmF,1))
  # calculate auc
  auc <- sum(diff(kmT) * head(kmF,-1))
  # calculate std WTW
  ## (1-kmF) gives the cdf of WTW
  cdf_WTW = 1 - kmF
  ## by adding 1 at the end, we truncate WTW at tMax
  cdf_WTW = c(cdf_WTW, 1)
  # calcualte the pdf of WTW
  pdf_WTW = diff(cdf_WTW)
  # calculate the std of WTW 
  varWTW = sum((kmT - auc) ^2 * pdf_WTW)
  stdWTW = sqrt(varWTW)
  # plot the k-m survival curve
  if (plotKMSC) {
    p = data.frame(kmT = kmT, kmF = kmF) %>%
      ggplot(aes(kmT, kmF)) + geom_line() + xlab('Delay (s)') +
      ylab('Survival rate') + ylim(c(0,1)) + xlim(c(0,tMax)) +
      ggtitle(sprintf('AUC = %1.1f',auc)) + 
      myTheme
    print(p)
  }
  # resample the survival curve for averaging 
  kmOnGrid = resample(kmF, kmT, kmGrid)
  return(list(kmT=kmT, kmF=kmF, auc=auc, kmOnGrid = kmOnGrid, stdWTW = stdWTW))
}

# calculate willingness to wait (WTW) time-series
wtwTS <- function(thisTrialData, tGrid, wtwCeiling, plotWTW = F) {
  # ensure timeWaited = scheduledWait on rewarded trials
  thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= 0] = scheduledWait[trialEarnings!= 0]})
  
  # initialize the per-trial estimate of WTW
  nTrial = length(thisTrialData$trialEarnings)
  trialWTW = numeric(length = nTrial) 
  
  ## find the longest time waited up through the first quit trial
  ## (or, if there were no quit trials, the longest time waited at all)
  ## that will be the WTW estimate for all trials up to the first quit
  quitIdx = thisTrialData$trialEarnings == 0
  firstQuitTrial = which(quitIdx)[1] 
  if (is.na(firstQuitTrial)) {firstQuitTrial = nTrial} 
  currentTrial = firstQuitTrial
  currentWTW = max(thisTrialData$timeWaited[1 : currentTrial]) 
  trialWTW[1:currentTrial] = currentWTW 
  ## iterate through the remaining trials, updating currentWTW
  ## which is the longest time waited since the recentest quit trial
  if(currentTrial < nTrial){
    for(tIdx in (currentTrial + 1) : nTrial){
      if(quitIdx[tIdx]) {currentWTW = thisTrialData$timeWaited[tIdx]}
      else {currentWTW = max(currentWTW, thisTrialData$timeWaited[tIdx])}
      trialWTW[tIdx] = currentWTW 
    }
  }
  
  # impose a ceiling value, since WTW exceeding some value may be infrequent
  trialWTW = pmin(trialWTW, wtwCeiling)
  
  # convert from per-trial to per-second 
  timeWTW = resample(trialWTW, thisTrialData$sellTime, tGrid)
  
  # plot time WTW
  if(plotWTW){
    p = ggplot(data.frame(tGrid, timeWTW), aes(tGrid, timeWTW)) + geom_line() +
      xlab("Time in block (s)") + ylab("WTW (s)")  + myTheme
    print(p)
  }
  
  # return outputs 
  outputs = list(timeWTW = timeWTW,  trialWTW = trialWTW)
  return(outputs)
}

# resample pair-wise sequences
# inputs:
# ys: y in the original sequence
# xs: x in the original sequence
# Xs: x in the new sequence
# outputs: 
# Ys : y in the new sequence 
resample = function(ys, xs, Xs){
  isBreak = F
  # initialize Ys
  Ys = rep(NA, length = length(Xs))
  for(i in 1 : length(Xs)){
    # for each X in Xs
    X = Xs[i]
    # find the index of cloest x value on the right
    # if closest_right_x_idx exists 
    if(X <= tail(xs,1)) {
      # Y takes the corresonding y value
      closest_right_x_idx = min(which(xs >= X))
      Ys[i] = ys[closest_right_x_idx]
    }else{
      isBreak = T
      lastY = i - 1
      break
    }
  }
  # fill the remaining elements in Ys by the exisiting last element
  if(isBreak){
    Ys[(lastY + 1) : length(Xs)] = Ys[lastY]
  }
  return(Ys)
}

# this function can truncate trials in the simualtion object
# which enables us to zoom in and look and specific trials
truncateTrials = function(data, startTidx, endTidx){
  nVar = length(data)
  varNames = names(data)
  outputs = vector(mode = "list", length = nVar)
  anyMatrix = F
  for(i in 1 : nVar){
    junk = data[[varNames[i]]]
    if(is.matrix(junk)){
      outputs[[i]] = junk[, startTidx:endTidx]
      anyMatrix  = T 
    }else{
      outputs[[i]] = junk[startTidx:endTidx]
    }
  }
  names(outputs) = varNames
  if(!anyMatrix )   outputs = as.data.frame(outputs)
  return(outputs)
}


# correlation plot
# the first col of plotData is x, the second col is y, the third col is the group
plotCorrelation = function(data, dotColor = "black",isRank){
  conditions = c("HP", "LP")
  colnames(data) = c("x", "y", "cond")
  
  # calculate correlations
  corTests = lapply(1:2, function(i) cor.test(data[data$cond == conditions[i], "x"],
                                              data[data$cond == conditions[i], "y"],
                                              method = "spearman"))
  # corTestsPerm = lapply(1:2, function(i) spearman_test(data[data$cond == conditions[i], "x"] ~ data[data$cond == conditions[i], "y"]))
  nus = sapply(1:2, function(i) round(as.numeric(corTests[[i]]$estimate), 3))
  ps = sapply(1:2, function(i) round(corTests[[i]]$p.value, 3))
  # ps = sapply(1:2, function(i) round(as.numeric(pvalue(corTestsPerm[[i]])), 3))
  
  textColors = ifelse(ps < 0.05, "red", "blue")
  textData = data.frame(label = paste(nus, "(p =", ps, ")"),
                        cond= c("HP", "LP"), color = textColors)
  # prepare rank 
  if(isRank){
    plotData = data %>% group_by(cond) %>% mutate(xRank = rank(x), yRank = rank(y))
  }
  
  # plot
  if(isRank){
    p0 = ggplot(plotData, aes(xRank, yRank)) + geom_point(size = 3, color = dotColor, fill = dotColor)
  }else{
    p0 = ggplot(plotData, aes(x, y)) + geom_point(size = 3, color = dotColor, fill = dotColor)
  }
  p = p0  + geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
                      hjust   = -0.2,vjust = -1,color = textColors, size = 5, fontface = 2, color = "#252525") +
    facet_grid(~cond)
  return(p)
} 


getCorrelation = function(data){
  conditions = c("HP", "LP")
  colnames(data) = c("x", "y", "cond")
  
  # calculate correlations
  # since we can't get nu from the later
  corTests = lapply(1:2, function(i) cor.test(data[data$cond == conditions[i], "x"],
                                              data[data$cond == conditions[i], "y"],
                                              method = "kendall") 
  )
  # supposedly, kendall can deal with data with a lot of ties
  # corTestsPerm = lapply(1:2, function(i) spearman_test(data[data$cond == conditions[i], "x"] ~
  #                                                    data[data$cond == conditions[i], "y"]))
  nus = sapply(1:2, function(i) as.numeric(corTests[[i]]$estimate))
  ps = sapply(1:2, function(i) round(corTests[[i]]$p.value, 3))
  # ps = sapply(1:2, function(i) round(pvalue(corTestsPerm [[i]]), 3))
  return(list(nus = nus, ps = ps))
}

getPartCorrelation = function(data){
  library("ppcor")
  
  conditions = c("HP", "LP")
  colnames(data) = c("x", "y", "z", "cond")
  
  # calculate correlations
  # since we can't get nu from the later
  corTests = lapply(1:2, function(i) pcor.test(data[data$cond == conditions[i], "x"],
                                               data[data$cond == conditions[i], "y"],
                                               data[data$cond == conditions[i], "z"],
                                               method = "kendall") 
  )
  # supposedly, kendall can deal with data with a lot of ties
  # corTestsPerm = lapply(1:2, function(i) spearman_test(data[data$cond == conditions[i], "x"] ~
  #                                                    data[data$cond == conditions[i], "y"]))
  nus = sapply(1:2, function(i) as.numeric(corTests[[i]]$estimate))
  ps = sapply(1:2, function(i) round(corTests[[i]]$p.value, 3))
  # ps = sapply(1:2, function(i) round(pvalue(corTestsPerm [[i]]), 3))
  return(list(nus = nus, ps = ps))
}

# integrate stacked data from several blocks 
block2session = function(thisTrialData){
  nBlock = length(unique(thisTrialData$blockNum))
  nTrials = sapply(1:nBlock, function(i) sum(thisTrialData$blockNum == i))
  # accumulated trials 
  ac_nTrials = c(0, cumsum(head(nTrials, -1)))
  # accumulated task durations
  ac_taskTimes = (1:nBlock - 1) * blockSec
  # accumulated totalEarnings sofar
  ac_totalEarnings_s = c(0, thisTrialData$totalEarnings[cumsum(nTrials)[1:(nBlock - 1)]])
  
  # convert within-block variables to accumualting-across-block variables
  thisTrialData = within(thisTrialData,
                         {trialNum = trialNum + rep(ac_nTrials , time = nTrials);
                         sellTime = sellTime + rep(ac_taskTimes, time = nTrials);
                         totalEarnings = totalEarnings + rep(ac_totalEarnings_s, time = nTrials)});
  return(thisTrialData)
}


