# this script plots delay distributions, reward rates and subjective values in two environments
expSchematics = function(smallReward, iti, isPlot){
# small reward 
smallReward = 0 # get 0 when the agent quits
  
# load experiment parameters 
load("expParas.RData")

# for display purposes, all variables on the continous time scale
# are discretized into 0.1 second time bins
bin = 0.1 # width of a time bin
stepSec = 1
time = list(
  HP = seq(bin, delayMaxs[1], by = bin),
  LP = seq(bin, delayMaxs[2], by = bin)
) 

# delay CDFS
## early delay distribution: gamma(k = 2, theta = 2), truncated at 30
## late delay distribution: gamma(k = 6, theta = 2), truncated at 30
fastCDF = pgamma(time$HP,  2, scale = 2, log = FALSE) 
fastCDF[length(fastCDF)] = 1 
slowCDF = pgamma(time$LP,  6, scale = 2, log = FALSE)
slowCDF[length(slowCDF)] = 1 


# delay PDFs
fastPDF = diff(c(0,fastCDF))
slowPDF = diff(c(0,slowCDF))

# average waiting durations given different policies
# Here we assume rewards occur at the middle of each time bin
slowMeanDelay = cumsum((time$HP - 0.5 * bin) * slowPDF) / cumsum(slowPDF)
fastMeanDelay =  cumsum((time$LP - 0.5 * bin) * fastPDF) / cumsum(fastPDF)


# rewardRates given different policies
HPRate = (8 *  slowCDF + (-1) * fastCDF) / 2 / 
  ((slowMeanDelay *  slowCDF + time$HP * (1 - slowCDF)) / 2 +
     (fastMeanDelay *  fastCDF + time$HP * (1 - fastCDF)) / 2 + iti)
LPRate = ((-1) *  slowCDF + (8) * fastCDF) / 2 / 
  ((slowMeanDelay *  slowCDF + time$LP * (1 - slowCDF)) / 2 +
     (fastMeanDelay *  fastCDF + time$LP * (1 - fastCDF)) / 2 + iti)

# optimal raward rates and optimal policies
optimWaitThresholds = list()
optimWaitThresholds$HP = time$HP[which.max(HPRate)]
optimWaitThresholds$LP = time$LP[which.max(LPRate)]
optimRewardRates = list()
optimRewardRates$HP = max(HPRate)
optimRewardRates$LP = max(LPRate)

# calculate the value of waiting as a function of elapsed time 
subjectValues = list()
coarseSubjectValues = list()
cIdx = 1
i = 1
for(cIdx in 1 : 2){
  condition = conditions[cIdx]
  delayMax = delayMaxs[cIdx]
  thisTime = time[[cIdx]]
  Rstar = optimRewardRates[[cIdx]]
  ts = seq(0, delayMax, by = 0.1) # elapsed time
  
  # initialize 
  thisSubjectValues = vector(length = length(ts)) # value of waiting given the elapsed time value
  thisCoarseSubjectValues = vector(length = delayMax / stepSec + 1) 
  Tstars = vector(length = length(ts)) # optimal waiting policy given the elapsed time value
  
  # loop over different elapsed time values
  c_idx = 1 # index for thisCoarseSubjectiveValues
  for(i in 1 : length(ts)){
    t = ts[i] # this elapsed time
    trctTime = thisTime[thisTime > t]
    
    # loop over different waiting policies to find Tstar and gt_max
    if(t == delayMax){
      Tstar = t
      gt_max = max((ifelse(condition == "HP", -1, 8) * tail(fastPDF, 1) +
                  ifelse(condition == "HP", 8, -1)  * tail(slowPDF, 1)) / (tail(fastPDF, 1) + tail(slowPDF, 1)), 0)
    }else{
      Tstar = t
      gt_max = -100
      for(T in seq(t, delayMax, by = 0.1)){
        # normalize the truncated distributions
        trctFastPDF = fastPDF[thisTime > t] / sum(fastPDF[thisTime > t] + slowPDF[thisTime > t])
        trctSlowPDF = slowPDF[thisTime > t] / sum(fastPDF[thisTime > t] + slowPDF[thisTime > t])
        # expected reward
        at = ifelse(condition == "HP", -1, 8) * sum(trctFastPDF[trctTime <= T]) + 
          ifelse(condition == "HP", 8, -1) * sum(trctSlowPDF[trctTime <= T])
        # remaining delay 
        bt = sum((trctTime[trctTime <= T] - 0.5 * bin - t) * trctFastPDF[trctTime <= T]) + 
          sum((trctTime[trctTime <= T] - 0.5 * bin - t) * trctSlowPDF[trctTime <= T])+
          (T - t) * sum(trctFastPDF[trctTime > T] + trctSlowPDF[trctTime > T]) 
        # net value
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

if(isPlot){
  # plot cdf 
  library("tidyverse")
  library("latex2exp")
  source('subFxs/plotThemes.R')
  dir.create("../../figures/wtw_exp3/expSchematics")
  figCDF = data.frame(cdf = c(0, fastCDF, 0, slowCDF),
             time = c(0, time$HP, 0, time$LP),
             cond = rep(c('Early', 'Late'), each = length(fastCDF) + 1)) %>%
    ggplot(aes(time, cdf)) + geom_line(size = 3) + facet_grid(~cond) + 
    myTheme + xlab('Delay duration (s)') + ylab('CDF') + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)),limits = c(0, max(delayMaxs) * 1.1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) + 
    geom_text(data = data.frame(label = c("HP = -1¢", "HP = +8¢"), cond = c("Early", "Late"), x = c(16, 22)),
              aes(label = label, x = x), y = 0.4, size = 5, color = conditionColors[1]) + 
    geom_text(data = data.frame(label = c("LP = +8¢", "LP = -1¢"), cond = c("Early", "Late"), x = c(16, 22)),
              aes(label = label, x = x), y = 0.2, size = 5, color = conditionColors[2])
  ggsave('../../figures/wtw_exp3/expSchematics/CDF.eps', width =4, height = 3)
  # plot pdf
  # the end point is plotted as a separate dot 
  plotData = data.frame(
    pdf = c(fastPDF, slowPDF, fastPDF, slowPDF),
    time = c(time$HP, time$HP, time$LP, time$LP),
    condition = c(rep("HP", 2 * length(time$HP)), rep("LP", 2 * length(time$HP))),
    value = factor(c(rep(min(tokenValue), length(time$HP)),
              rep(max(tokenValue), length(time$HP)),
              rep(max(tokenValue), length(time$HP)),
              rep(min(tokenValue), length(time$HP))))) 
  endPointData = plotData[plotData$time == max(time$HP),]
  plotData$pdf[plotData$time == max(time$HP)] = NA
  
  figPDF = plotData %>% ggplot(aes(time, pdf, color = value)) + geom_line(size = 1) + facet_grid(~condition) +
    myTheme + xlab('Delay duration (s)') + ylab("Probability density") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)),limits = c(0, max(delayMaxs) * 1.1))  +
    scale_color_manual(values = c("grey", "gold")) + 
    geom_point(data = endPointData, aes(time, pdf, color = value), inherit.aes = F, size = 2) +
    scale_y_continuous(breaks = c(0, 0.02), labels = c(0, 0.02), limits = c(0, 0.022)) +
    theme(legend.position = "None") + ggtitle("Exp.3")
  ggsave('../../figures/wtw_exp3/expSchematics/PDF.eps', width =4, height = 3)
  
  # plot reward rates
  figRewardRate = data.frame(rate = c (HPRate, LPRate), 
             time = c(time$HP, time$LP),
             condition = rep(c('HP', 'LP'), each = length(time$HP))) %>%
    ggplot(aes(time, rate)) + geom_line(size = 2, aes(color = condition)) + facet_grid(~cond) +
    myTheme + 
    ylab(expression(bold(paste("Reward rate ", Phi['T'], "(¢ s"^"-1", ")")))) + xlab("Give-up time (s)") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, max(delayMaxs)/2, max(delayMaxs))) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + facet_grid(~condition) + ylim(c(-0.05, 0.5)) + ggtitle("Exp.3")
  ggsave("../../figures/wtw_exp3/expSchematics/reward_rate.eps", width = 4, height = 3)
  
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
    scale_x_continuous(breaks = c(0, max(delayMaxs)/2, max(delayMaxs)))
  ggsave("../../figures/wtw_exp3/expSchematics/subjective.eps", width = 4, height = 3)
}    


# return outputs 
if(isPlot){
  outputs = list(
    "optimWaitThresholds" = optimWaitThresholds,
    "optimRewardRates" = optimRewardRates,
    "time" = time,
    "subjectValues" = subjectValues,
    "fastPDF" = fastPDF,
    "slowPDF" = slowPDF,
    "coarseSubjectValues" = coarseSubjectValues,
    "figs" = list("pdf" = figPDF, "cdf" = figCDF, "rt" = figRewardRate, "sv" = figSV))
}else{
  outputs = list(
    "optimWaitThresholds" = optimWaitThresholds,
    "optimRewardRates" = optimRewardRates,
    "time" = time,
    "subjectValues" = subjectValues,
    "fastPDF" = fastPDF,
    "slowPDF" = slowPDF,
    "coarseSubjectValues" = coarseSubjectValues)
}
}