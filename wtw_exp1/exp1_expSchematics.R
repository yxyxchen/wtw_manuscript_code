# this script plots delay distributions and reward rates in two environments
expSchematics = function(smallReward, iti, isPlot){
  # load libararies
  library("tidyverse")
  library(latex2exp)
  source('subFxs/plotThemes.R')
  
  # load experiment parameters 
  load("expParas.RData")
  
  # for display purposes, all variables on the continous time scale
  # are discretized into 0.1 second time bins
  bin = 0.1 # width of a time bin
  time = list(
    HP = seq(bin, delayMaxs[1], by = bin),
    LP = seq(bin, delayMaxs[2], by = bin)
  ) 
  stepSec = 1
  
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
  coarseSubjectValues = list()
  Tstars = list()
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
    thisCoarseSubjectValues = vector(length = delayMax / stepSec + 1)
    thisTstars = vector(length = length(ts))
    # loop over different elapsed time
    c_i = 1 # index for thisCoarseSubjectValues
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
      if(abs(t - t %/% stepSec * stepSec) < bin * 0.5){
        thisCoarseSubjectValues[c_i] = gt_max
        c_i = c_i + 1
      }
      thisTstars[i] = Tstar
    }
    subjectValues[[condition]] = thisSubjectValues
    coarseSubjectValues[[condition]] = thisCoarseSubjectValues
    Tstars[[condition]] = thisTstars
  }
  
  if(isPlot){
    # plot PDF
    # the end point in the LP condition is plotted separately 
    plotData =  data.frame(
      pdf = as.numeric(unlist(rewardDelayPDFs)),
      time = as.numeric(unlist(time)),
      condition = c(rep("HP", length(rewardDelayPDFs$HP)), rep("LP", length(rewardDelayPDFs$LP)))
    )
    
    endPointData = plotData[plotData$time == max(time$LP) & plotData$condition == "LP",]

    plotData$pdf[(plotData$time == max(time$LP) & plotData$condition == "LP")] = NA
    figPDF = plotData %>% ggplot(aes(time, pdf)) + geom_line(size = 1.5, aes(color = condition)) +
      facet_grid(~condition)  + 
      geom_point(data = endPointData, aes(time, pdf, color = condition), size = 2, inherit.aes = F) +
      scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)), labels = c("0", max(delayMaxs)/2, max(delayMaxs)),
                         limits = c(0, max(delayMaxs) * 1.1)) + 
      myTheme + xlab('Delay duration (s)') + ylab('Probability density') + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none") + 
      scale_y_continuous(breaks = c(0, 0.05), labels = c(0, 0.05), limits = c(0, 0.08)) +
      ggtitle("Exp.1") 
    
    # plot CDFs 
    ## here we extend the HP CDF to 32s for display purposes
    figCDF = data.frame(CDF = c(0,c(rewardDelayCDFs$HP, rep(1, length(time$LP) - length(time$HP))), 0, rewardDelayCDFs$LP),
               time = c(0, time$LP, 0, time$LP),
               condition =  rep(c('HP', 'LP'), c(length(time$LP) + 1, length(time$LP) + 1))) %>%
      ggplot(aes(time, CDF)) + geom_line(size = 3, aes(color = condition)) + facet_grid(~condition) +
      ylim(c(0,1)) + scale_y_continuous(breaks = c(0,0.5,1)) + 
      scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)), labels = c("0", max(delayMaxs)/2, max(delayMaxs)),
                         limits = c(0, max(delayMaxs) * 1.1)) + 
      myTheme + xlab('Delay duration (s)') + ylab('CDF') + 
      theme(plot.title = element_text(hjust = 0.5, color = themeColor)) + 
      annotate("text", x = max(delayMaxs)/ 2, y = 0.4, label = "+10¢", size = 6) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")
    
    # plot reward rates
    optimData = data.frame(condition = c("HP", "LP"), waitThreshold = as.double(optimWaitThresholds))
    figRewardRate = data.frame(rewardRate = c(rewardRates[[1]], rewardRates[[2]]),
               time = as.numeric(unlist(time)),
               condition = rep(c("HP", "LP"), c(length(time$HP), length(time$LP)))) %>%
      ggplot(aes(time, rewardRate)) +
      geom_line(size = 2, aes(color = condition))  + myTheme + 
      ylab(expression(bold(paste("Reward rate ", Phi['T'], "(¢ s"^"-1", ")")))) + xlab("Give-up time (s)") +
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2), limits = c(0, 1.5)) +
      scale_x_continuous(breaks = c(0, max(delayMaxs)/2, max(delayMaxs))) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none") + facet_grid(~condition) +
      ggtitle("Exp.1") 

    # plot subjective value of waiting 
    delayClock = c(seq(0, delayMaxs[1], by = 0.1), seq(0, delayMaxs[2], by = 0.1)) # time starting from the token onset
    figSV = data.frame(
      value =  c(subjectValues$HP, subjectValues$LP),
      delayClock = delayClock,
      condition = rep(conditions, c(length(subjectValues$HP), length(subjectValues$LP))))%>%
      ggplot(aes(delayClock, value)) +
      geom_line(aes(color = condition), size = 2) +
      myTheme+
      scale_color_manual(values = conditionColors) +
      scale_linetype_manual(values = c(1, 2)) +
      scale_x_continuous(breaks = c(0, max(delayMaxs)/2, max(delayMaxs)))+
      xlab("Trial time (s)") + ylab("Subjective value (¢)")  + 
      theme(legend.position = "none") + facet_grid(~condition)
  }
  
  # return outputs 
  if(isPlot){
    outputs = list(
      "optimWaitThresholds" = optimWaitThresholds,
      "optimRewardRates" = optimRewardRates,
      "rewardRates" = rewardRates,
      "subjectValues" = subjectValues,
      "coarseSubjectValues" = coarseSubjectValues,
      "time" = time,
      "figs" =  list("pdf" = figPDF, "cdf" = figCDF, "rt" = figRewardRate, "sv" = figSV)
    )
  }else{
    outputs = list(
      "optimWaitThresholds" = optimWaitThresholds,
      "optimRewardRates" = optimRewardRates,
      "rewardRates" = rewardRates,
      "subjectValues" = subjectValues,
      "coarseSubjectValues" = coarseSubjectValues,
      "time" = time
    )
  }

  return(outputs)
}


  


