# this script plots delay distributions and reward rates in two environments
expSchematics = function(smallReward, iti, isPlot){
  # small reward 
  smallReward = 0 # get 0 when the agent quits

  # load experiment parameters 
  load("expParas.RData")
  seqHP = rewardDelays$HP
  seqLP = rewardDelays$LP

  ## for display purposes, all variables on continous time
  ## are discretized into 0.1 secondition time bins
  bin = 0.1 # width of a time bin
  time = list(
    HP = seq(bin, delayMaxs[1], by = bin),
    LP = seq(bin, delayMaxs[2], by = bin)
  ) 
  stepSec = 1

  ## delays in both HP and LP are drawn repeatedly from 8 possible values,
  ## In HP, they are uniformly distributed from 0 to 16, 
  ## in LP, they are log-spaced from 0 to 32, therefore forming a heavy tail distribution

  # delay PDFs
  HP = rep(0, length(time$HP))
  LP = rep(0, length(time$LP))
  HP[sapply(1 :8, function(i) which.min(abs(time$HP - seqHP[i])))] = 1 / 8
  LP[sapply(1 :8, function(i) which.min(abs(time$LP - seqLP[i])))] = 1 / 8
  rewardDelayPDFs =list('HP'= HP, 'LP' = LP)

  # delay CDFs
  rewardDelayCDFs = list('HP' = cumsum(rewardDelayPDFs$HP),
                        'LP' = cumsum(rewardDelayPDFs$LP))

  # average waiting durations given different policies
  ## Here we assume rewards occur at the middle of each time bin
  HP = cumsum((time[['HP']]) * rewardDelayPDFs$HP) / cumsum(rewardDelayPDFs$HP)
  LP = cumsum((time[['LP']]) * rewardDelayPDFs$LP) / cumsum(rewardDelayPDFs$LP)
  meanRewardDelays = list('HP' = HP, 'LP' = LP)

  # rewardRates given different policies
  HP = tokenValue * rewardDelayCDFs$HP /
    ((meanRewardDelays$HP * rewardDelayCDFs$HP + time[['HP']] * (1 - rewardDelayCDFs$HP)) + iti)
  HP[is.nan(HP)] = 0
  LP = tokenValue * rewardDelayCDFs$LP /
    ((meanRewardDelays$LP * rewardDelayCDFs$LP + time[['LP']] * (1 - rewardDelayCDFs$LP)) + iti)
  rewardRates = list('HP' = HP, 'LP' = LP)


  # optimal raward rates and optimal policies
  optimWaitThresholds = list()
  optimWaitThresholds$HP = time$HP[which.max(HP)]
  optimWaitThresholds$LP = time$LP[which.max(LP)]
  optimRewardRates = list()
  optimRewardRates$HP = max(HP, na.rm = T)
  optimRewardRates$LP = max(LP, na.rm = T)

  # calculate theoretical value of waiting as a function of elapsed time 
  subjectValues = list()
  coarseSubjectValues = list()
  for(cIdx in 1 : 2){
    condition = conditions[cIdx]
    delayMax = delayMaxs[cIdx]
    pdf = rewardDelayPDFs[[cIdx]]
    cdf = rewardDelayCDFs[[cIdx]]
    thisTime = time[[cIdx]]
    Rstar = optimRewardRates[[cIdx]] # opportunity cost
    ts = seq(0, delayMax, by = 0.1) # elapsed time
    
    # initialize 
    thisSubjectValues = vector(length = length(ts)) # waiting value given the elapsed time value
    thisCoarseSubjectValues = vector(length = delayMax / stepSec + 1)
    Tstars = vector(length = length(ts))  # optimal waiting policy given the elapsed time value
    # loop over different elapsed time values
    c_idx = 1 # index for thisCoarseSubjectValues 
    for(i in 1 : length(ts)){
      t = ts[i] 
      trctTime = thisTime[thisTime > t]
      # given the elapsed time value, t, loop over different waiting policies to find Tstar 
      if(t == delayMax){
        Tstar = t
        gt_max = tokenValue 
      }else{
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
        thisCoarseSubjectValues[c_idx] = gt_max
        c_idx = c_idx + 1
      }
      Tstars[i] = Tstar
    }
    subjectValues[[condition]] = thisSubjectValues
    coarseSubjectValues[[condition]] = thisCoarseSubjectValues
  }

   # plot 
  if(isPlot){
    
    # plot CDFs 
    library('ggplot2')
    source('subFxs/plotThemes.R')
    library("tidyr"); library('dplyr')
    ## here we extend the HP CDF to 32s for display purposes
    figCDF = data.frame(CDF = c(0,c(rewardDelayCDFs$HP, rep(1, length(time$LP) - length(time$HP))), 0, rewardDelayCDFs$LP),
               time = c(0, time$LP, 0, time$LP),
               condition =  rep(c('HP', 'LP'), c(length(time$LP) + 1, length(time$LP) + 1))) %>%
      ggplot(aes(time, CDF)) + geom_line(size = 3, aes(color = condition)) + facet_grid(~condition) +
      ylim(c(0,1)) + scale_y_continuous(breaks = c(0,0.5,1)) + 
      scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)),limits = c(0, max(delayMaxs) * 1.1)) +
      myTheme + xlab('Delay duration (s)') + ylab('CDF')  +
      annotate("text", x = 16, y = 0.4, label =  "+10¢", size = 6) + 
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
      scale_color_manual(values = conditionColors) 
    
    # plot PDFs
    tempt = unlist(rewardDelayPDFs)
    tempt[tempt == 0] = NA
    vlineData = data.frame(
      delay = c(time$HP[rewardDelayPDFs$HP > 0], time$LP[rewardDelayPDFs$LP > 0]),
      condition = c(rep("HP", length(seqHP)), rep("LP", length(seqLP)))
    )
    
    figPDF = data.frame(
      PDF = tempt,
      time = unlist(time),
      condition = c(rep("HP", length(rewardDelayPDFs$HP)), rep("LP", length(rewardDelayPDFs$LP)))
    ) %>% ggplot(aes(time, PDF)) + geom_point(size = 1.5, aes(color = condition)) + facet_grid(~condition) + 
      scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)), labels = c("0", max(delayMaxs)/2, max(delayMaxs)),
                         limits = c(0, max(delayMaxs) * 1.1)) + 
      myTheme + xlab('Delay duration (s)') + ylab('Probability mass') + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none") + scale_y_continuous(breaks = c(0, 0.125), limits = c(0, 0.15), labels = c(0, 0.125)) + 
      geom_segment(data =  vlineData, aes(color = condition, x = delay, xend = delay), y = 0, yend = 1 / 8) +
      ggtitle("Exp.2")
    
    # plot reward rates
    optimData = data.frame(condition = c("HP", "LP"), waitThreshold = as.double(optimWaitThresholds))
    figRewardRate = data.frame(
      time = unlist(time),
      rewardRate = unlist(rewardRates),
      condition = rep(c("HP", "LP"), c(length(time$HP), length(time$LP)))
    ) %>% ggplot(aes(time, rewardRate)) +
      geom_line(size = 2, aes(color = condition)) + facet_grid(~condition) +
      myTheme +
      ylab(expression(bold(paste("Reward rate ", Phi['T'], "(¢ s"^"-1", ")")))) + xlab("Give-up time (s)")  +
      scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2), limits = c(0, 1.3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")  + 
      scale_x_continuous(breaks = c(0, max(delayMaxs) / 2, max(delayMaxs))) +
      ggtitle("Exp.2")
    
    # plot subjective value of waiting 
    figSV = data.frame(
      value =  c(subjectValues$HP, subjectValues$LP),
      t = c(seq(0, delayMaxs[1], by = 0.1), seq(0, delayMaxs[2], by = 0.1)),
      condition = rep(conditions, c(length(subjectValues$HP), length(subjectValues$LP))))%>%
      ggplot(aes(t, value)) +
      geom_line(aes(color = condition), size = 2) +
      myTheme+
      scale_color_manual(values = conditionColors) +
      scale_linetype_manual(values = c(1, 2)) +
      xlab("Elapsed time (s)") + ylab("Subjective value (¢)")  + 
      theme(legend.position = "none") + facet_grid(~condition) + 
      scale_x_continuous(breaks = c(0, max(delayMaxs) / 2, max(delayMaxs)))
  }

  # return outputs 
  if(isPlot){
    outputs = list(
      "optimWaitThresholds" = optimWaitThresholds,
      "optimRewardRates" = optimRewardRates,
      "subjectValues" = subjectValues,
      "coarseSubjectValues" = coarseSubjectValues,
      "figs" =  list("pdf" = figPDF, "cdf" = figCDF, "rt" = figRewardRate, "sv" = figSV)
    )
  }else{
    outputs = list(
      "optimWaitThresholds" = optimWaitThresholds,
      "optimRewardRates" = optimRewardRates,
      "subjectValues" = subjectValues,
      "coarseSubjectValues" = coarseSubjectValues
    )
  }
return(outputs)
}
