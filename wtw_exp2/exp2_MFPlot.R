MFPlot = function(){
  library("tidyverse")
  library("ggpubr")
  library("latex2exp")
  source("subFxs/plotThemes.R")
  source("exp2_expSchematics.R")
  source('MFAnalysis.R')
  
  # output dir
  dir.create('../../figures/wtw_exp2/MFplot')
  
  # load experiment parameters
  load("expParas.RData")
  
  # normative analysis
  iti = 2
  source("exp2_expSchematics.R")
  normResults = expSchematics(0, iti, F)
  optimWaitThresholds = normResults$optimWaitThresholds
  nCondition = length(conditions)
  
  ######################### plot WTW timecourses in two environments ######################
  MFResults = MFAnalysis(isTrct = F)
  sumStats = MFResults[['sumStats']]
  timeWTW_ = MFResults[['timeWTW_']]
  nSub = nrow(sumStats)
  ## use background color to distinguish used and excluded data 
  colorData = data.frame(
    xmin = c(0, blockSec), xmax = c(blockSec - max(delayMaxs), nBlock * blockSec - max(delayMaxs)),
    condition = c(conditions[2], conditions[1])
  )
  greyData = data.frame(
    xmin = 1:2 * blockSec - max(delayMaxs), xmax = 1:2 * blockSec
  )
  
  plotData = data.frame(wtw = unlist(timeWTW_),
             time = rep(seq(0, blockSec * nBlock -1, by = 1), nSub),
             condition = rep(rep(c("LP", "HP"), each = length(tGrid))), nSub) %>%
    group_by(condition, time) %>%
    dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                     min = mu- se, max = mu + se) %>% ungroup()
  # reset at the block onset
  plotData$mu[plotData$time %in% (1:2 * blockSec - 1)] = NA
  plotData$min[plotData$time %in% (1:2 * blockSec - 1)] = NA
  plotData$max[plotData$time %in% (1:2 * blockSec - 1)] = NA
  
  # plot
  figWTW = plotData %>%
    ggplot(aes(time, mu, color = condition, fill = condition))  +
    geom_rect(data = greyData, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 20,
              fill = "#d9d9d9", inherit.aes = F) + 
    geom_ribbon(aes(ymin=min, ymax=max),  color = NA) +
    geom_line(size = 0.5) +
    xlab("Task time (min)") + ylab("WTW (s)") + 
    scale_y_continuous(breaks = c(0, 8, 16), limits = c(0, 17)) +
    myTheme + scale_x_continuous(breaks = c(0, 600, 1200), labels = c(0, 600, 1200) / 60 ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
    scale_color_manual(values = conditionColors) + scale_fill_manual(values = c("#7fbf7b", "#af8dc3")) +
    theme(legend.position = "none")
  
  ################################################################
  ##              compare AUC in the two environments              ##
  ################################################################
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  sumStats %>% group_by(condition) %>% summarise(median = median(auc),
                                                 IQR(auc))
  wTest = wilcox.test( sumStats[sumStats$condition == "HP", "auc"],
                       sumStats[sumStats$condition == "LP", "auc"],paired = T)
  
  ################################################################
  ##              plot AUC in the two environments              ##
  ################################################################
  figAUC = data.frame(muWTWHP = sumStats$auc[sumStats$condition == 'HP'],
             muWTWLP = sumStats$auc[sumStats$condition == 'LP']) %>%
    ggplot(aes(muWTWLP, muWTWHP)) +
    geom_point(size = 5, shape = 21, stroke =1) +
    geom_abline(slope = 1, intercept = 0) + 
    annotate("text", x = 8, y = 3, label = sprintf('p = %0.3f*', wTest$p.value), size = 6) +
    xlab("LP AUC (s)") + ylab("HP AUC (s)") + 
    myTheme + xlim(c(-1,17)) + ylim(c(-1,17)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
    scale_x_continuous(breaks = c(0, 4, 8, 12, 16), limits = c(0, 16))
  
  
  ################################################################
  ##              compare sigma_wtw in the two environments              ##
  ################################################################
  wTest = wilcox.test( sumStats[sumStats$condition == "HP", "stdWTW"],
                       sumStats[sumStats$condition == "LP", "stdWTW"], paired = T)
  
  ################################################################
  ##              plot sigma_wtw in the two environments              ##
  ################################################################
  sumStats %>% group_by(condition) %>% summarise(median = median(stdWTW),
                                                 IQR = IQR(stdWTW))
  # plot
  figSigma = data.frame(stdWTWHP = sumStats$stdWTW[sumStats$condition == 'HP'],
             stdWTWLP = sumStats$stdWTW[sumStats$condition == 'LP']) %>%
    ggplot(aes(stdWTWLP, stdWTWHP)) +
    geom_point(size = 5, shape = 21, stroke =1) +
    geom_abline(slope = 1, intercept = 0)  +
    xlab(expression(bold(paste("LP ", sigma["WTW"], " (s"^2,")")))) +
    ylab(expression(bold(paste("HP ", sigma["WTW"], " (s"^2,")")))) + 
    myTheme + xlim(c(-1,6)) + ylim(c(-1,6)) +
    annotate("text", x = 5, y = 1, label = sprintf('p = %0.3f*', wTest$p.value))
  
  ################################################################
  ##              plot delta AUC in the two environments        ##
  ################################################################
  LP_end_vs_start = MFResults$sub_auc_[,2] - MFResults$sub_auc_[,1] 
  HP_end_vs_start = MFResults$sub_auc_[,4] - MFResults$sub_auc_[,3] 
  HP_vs_LP = MFResults$sumStats[sumStats$condition == "HP", "auc"] - MFResults$sumStats[sumStats$condition == "LP", "auc"] 
  figDelta = 
    data.frame(
      deltaWTW = c(LP_end_vs_start, HP_end_vs_start, HP_vs_LP),
      contrast = rep(c("LP", "HP", "HP-LP"), each = nrow(sumStats) / nBlock)
    ) %>% ggplot(aes(contrast, deltaWTW)) +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = contrast))  + 
    ylab(TeX('$AUC_{end} - AUC_{start}$')) +
    xlab(" ") + 
    scale_fill_manual(values = conditionColors) +
    myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position =  "none")
  
  wilcox.test(deltaWTW[sumStats$condition == "HP"])
  wilcox.test(deltaWTW[sumStats$condition == "LP"])
  cor.test(deltaWTW[sumStats$condition == "HP"], sumStats$auc[sumStats$condition == "HP"], method=c( "spearman"))
  cor.test(deltaWTW[sumStats$condition == "LP"], sumStats$auc[sumStats$condition == "LP"], method=c( "spearman"))
  
  ###################################################
  ##              plot survival curves            ##
  ##################################################
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
  # stats test
  survCurve_ = MFResults$survCurve_
  plotData = data.frame(survCurve = unlist(survCurve_),
                        time = rep(kmGrid, nSub * nCondition),
                        condition = rep(sumStats$condition, each = length(kmGrid)))
  # plot
  figCurve = plotData %>%
    group_by(condition, time) %>%
    dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),min = mu- se, max = mu + se) %>%
    ggplot(aes(time, mu, color = condition, fill = condition)) +
    geom_ribbon(aes(time, ymin = min, ymax = max), color = NA) +
    geom_line() + 
    geom_line(data = optim, aes(t, surv, color = condition, linetype = condition, alpha = condition), size = 1.2) +
    geom_line(data = data.frame(t = kmGrid[kmGrid > 2],surv = 1),
              aes(t, surv), color = conditionColors[1], size = 1.2, inherit.aes = F, alpha = 0.8) + 
    scale_fill_manual(values = c("#7fbf7b", "#af8dc3")) +
    scale_color_manual(values = conditionColors) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    scale_alpha_manual(values = c(0.8, 1))+
    xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0, 4, 8, 12, 16), limits = c(0, 16))
  
  ############# return the output ###############
  outputs = list(
    "wtw" = figWTW,
    "auc" = figAUC,
    "sigma" = figSigma,
    "curve" = figCurve
  )
  return(outputs)  
}
