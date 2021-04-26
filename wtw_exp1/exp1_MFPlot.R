MFPlot = function(){
  library('dplyr')
  library("tidyr")
  library("ggplot2")
  library("ggpubr")
  library("lme4")
  source("subFxs/plotThemes.R")
  source("MFAnalysis.R")
  library("lmerTest")
  library("latex2exp")
  
  # output dir
  dir.create('../../figures/wtw_exp1/MFplot')
  
  # load experiment parameters
  load("expParas.RData")
  
  # normative analysis 
  iti = 2
  source("exp1_expSchematics.R")
  normResults = expSchematics(0, iti, F)
  optimWaitThresholds = normResults$optimWaitThresholds
  
  ##################################################################
  ##           plot WTW timecourses in two environments           ##
  ##################################################################
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  timeWTW_ = MFResults[['timeWTW_']]
  nSub = nrow(sumStats)
  
  # use background color to distinguish used and excluded data 
  yellowData = data.frame(
    xmin = rep(blockSec * (0 : 2), 2), xmax = rep(blockSec * (0 : 2), 2) + blockSec - max(delayMaxs),
    condition = rep(c("HP", "LP"), each = nBlock)
  )
  greyData = data.frame(
    xmin = rep(blockSec * (1 : 3), 2) - max(delayMaxs), xmax = rep(blockSec * (1 : 3), 2),
    condition = rep(c("HP", "LP"), each = nBlock)
  )
  
  # plotData 
  plotData = data.frame(wtw = unlist(timeWTW_),
                        time = rep(tGrid, nSub),
                        condition = rep(sumStats$condition, each = length(tGrid)))
  plotData = plotData %>%
    group_by(condition, time) %>%
    dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                     min = mu - se, max = mu + se) %>% ungroup()
  ## reset at the block onset
  plotData$mu[plotData$time %in% (1:3 * blockSec - 1)] = NA
  plotData$min[plotData$time %in% (1:3 * blockSec - 1)] = NA
  plotData$max[plotData$time %in% (1:3 * blockSec - 1)] = NA
  
  # plot 
  figWTW = plotData %>% ggplot(aes(time, mu, color = condition))  +
    geom_ribbon(aes(ymin=min, ymax=max, fill = condition, color = NA), alpha = 0.5) +
    geom_line(aes(color = condition), size = 1) +
    xlab("Task time (min)") + ylab("Willingness to wait (s)") + 
    myTheme + ylim(0, 20)  +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
    scale_x_continuous(breaks = 0:3 * 420, labels = 0:3 * 7) + 
    theme(legend.position = "none") +
    scale_fill_manual(values = conditionColors) +
    scale_color_manual(values = conditionColors) +
    geom_hline(aes(yintercept = optimWaitThresholds$HP), color = "red", size = 2, linetype = "dashed") +
    geom_hline(aes(yintercept = optimWaitThresholds$LP), color = "red", size = 2, linetype = "dashed")
  
  ggsave("tmp/fig_WTW.pdf", width = 4, height = 4)
  ###################################################################
  ##              compare AUC in the two environments              ##
  ###################################################################
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  sumStats %>% group_by(condition) %>% summarise(median(muWTW),IQR(muWTW))
  wilcox.test(sumStats$muWTW[sumStats$condition == "HP"],
              sumStats$muWTW[sumStats$condition == "LP"])
  
  ################################################################
  ##              plot AUC in the two environments              ##
  ################################################################
  figAUC = sumStats %>% ggplot(aes(condition, muWTW))  +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                       bracket.size = 1, size = 6, label.y = 22) + ylab("AUC (s)") + 
    myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
    scale_y_continuous(breaks = c(0, 12, 24), limits = c(0, 26)) +
    scale_fill_manual(values = conditionColors) +
    theme(legend.position = "none") + xlab("")
  

  ###################################################################
  ##              compare sigma_wtw in the two environments        ##
  ###################################################################
  wilcox.test(sumStats$stdWTW[sumStats$condition == "HP"],
              sumStats$stdWTW[sumStats$condition == "LP"])
  
  ######################################################################
  ##              plot sigma_wtw in the two environments              ##
  ######################################################################
  sumStats %>% group_by(condition) %>% summarise(median(stdWTW))
  sumStats %>% group_by(condition) %>% summarise(median(stdWTW),IQR(stdWTW))
  figSigma = sumStats %>% ggplot(aes(condition, stdWTW))  +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                       bracket.size = 1, size = 6, label.y = 10) + xlab("") +
    ylab(expression(bold(paste(sigma["WTW"], " (s"^2,")"))))+ 
    myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
    scale_y_continuous(breaks = c(0, 5, 10), limits = c(0, 11)) +
    scale_fill_manual(values = conditionColors) +
    theme(legend.position = "none") + xlab("")
  
  ##################################################################
  ##                     plot survival curves                     ##
  ##################################################################
  optim = data.frame(
    t = rep(kmGrid,  2),
    surv = rep(1, length(kmGrid) * 2),
    condition = rep(conditions, each = length(kmGrid)),
    select = rep(1:2, each = length(kmGrid))
  ) 
  optim$surv[optim$condition == "LP" & kmGrid> optimWaitThresholds$LP] = 0 # quit after 2.2 s
  optim$surv[optim$condition == "HP" & kmGrid> optimWaitThresholds$LP] = NA # don't plot after 2.2 s
  optim$select[optim$condition == "HP" & kmGrid <= optimWaitThresholds$LP] = rep(1:2, each = 3) # plot interleaving colors 
  optim$select[optim$condition == "LP" & kmGrid <= optimWaitThresholds$LP] = rep(1:2, each = 3) # plot interleaving colors
  ## stats test
  survCurve_ = MFResults$survCurve_
  plotData = data.frame(survCurve = unlist(survCurve_),
                        time = rep(kmGrid, nSub),
                        condition = rep(sumStats$condition, each = length(kmGrid)))
  
  figCurve = plotData %>%
    group_by(condition, time) %>%
    dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),
                     min = mu- se, max = mu + se) %>%
    ggplot(aes(time, mu, color = condition, fill = condition)) + geom_line() +
    geom_ribbon(aes(time, ymin = min, ymax = max), alpha = 0.5, color = NA) +
    geom_line(data = optim, aes(t, surv, color = condition, linetype = condition, alpha = condition), size = 1.2) +
    geom_line(data = data.frame(t = kmGrid[kmGrid > 2],surv = 1),
              aes(t, surv), color = conditionColors[1], size = 1.2, inherit.aes = F, alpha = 0.8) + 
    scale_fill_manual(values = conditionColors) +
    scale_color_manual(values = conditionColors) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    scale_alpha_manual(values = c(0.8, 1))+
    xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
    theme(legend.position = "none") 
  
  ############# return the output ###############
  outputs = list(
    "wtw" = figWTW,
    "auc" = figAUC,
    "sigma" = figSigma,
    "curve" = figCurve
  )
  return(outputs)
}


