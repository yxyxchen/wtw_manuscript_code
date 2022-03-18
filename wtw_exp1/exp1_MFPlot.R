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
  MFResults = MFAnalysis(isTrct = F)
  sumStats = MFResults[['sumStats']]
  timeWTW_ = MFResults[['timeWTW_']]
  nSub = nrow(sumStats)
  
  
  # plotData 
  plotData = data.frame(wtw = unlist(timeWTW_),
                        time = rep(tGrid, nSub),
                        condition = rep(sumStats$condition, each = length(tGrid)))
  plotData = plotData %>%
    group_by(condition, time) %>%
    dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                     ss = sum( (wtw - mean(wtw) )^2 ),
                     min = mu - se, max = mu + se) %>% ungroup()
  ## reset at the block onset
  plotData$mu[plotData$time %in% (1:3 * blockSec - 1)] = NA
  plotData$min[plotData$time %in% (1:3 * blockSec - 1)] = NA
  plotData$max[plotData$time %in% (1:3 * blockSec - 1)] = NA
  
  greyData = data.frame(
    xmin = 1:nBlock * blockSec - max(delayMaxs), xmax = 1:nBlock * blockSec
  )
  
  # plot 
  figWTW = plotData %>% ggplot(aes(time, mu, color = condition)) + 
    geom_rect(data = greyData, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 20,
              fill = "#d9d9d9", inherit.aes = F) +
    geom_ribbon(aes(ymin=min, ymax=max, fill = condition), color = NA) +
    geom_line(aes(color = condition), size = 1) +
    xlab("Task time (min)") + ylab("WTW (s)") + 
    myTheme + ylim(0, 20)  +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
    scale_x_continuous(breaks = 0:3 * 420, labels = 0:3 * 7) + 
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#7fbf7b", "#af8dc3")) +
    scale_color_manual(values = conditionColors)
    # geom_hline(aes(yintercept = optimWaitThresholds$HP), color = "red", size = 2, linetype = "dashed") +
    # geom_hline(aes(yintercept = optimWaitThresholds$LP), color = "red", size = 2, linetype = "dashed")
  
  
  
  ############
  ###########
  sumStats %>% group_by(condition) %>% summarise(median(auc),IQR(auc))
  wilcox.test(sumStats$auc[sumStats$condition == "HP"],
              sumStats$auc[sumStats$condition == "LP"])
  
  
  
  ###################################################################
  ##              compare stats in the two environments              ##
  ###################################################################
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  sumStats$delta = MFResults$sub_auc_[,6] - MFResults$sub_auc_[,1] 
  
  stats_df = sumStats %>% gather("variable", "value", -setdiff(names(sumStats), c("auc", "stdWTW", "delta"))) %>%
    group_by(variable, condition) %>%
    summarise(
      median = round(median(value), 2),
      IQR = round(IQR(value), 2),
      p = round(wilcox.test(value)$p.value, 3)
    )
  
  contrast_df = 
    data.frame(
      variable = c("auc", "delta", "stdWTW"),
      pvalue = c(
        round(wilcox.test( sumStats[sumStats$condition == "HP", "auc"],
                           sumStats[sumStats$condition == "LP", "auc"],paired = F)$p.value, 3),
        round(wilcox.test( sumStats[sumStats$condition == "HP", "delta"],
                           sumStats[sumStats$condition == "LP", "delta"],paired = F)$p.value, 3),
        round(wilcox.test( sumStats[sumStats$condition == "HP", "stdWTW"],
                           sumStats[sumStats$condition == "LP", "stdWTW"],paired = F)$p.value, 3)
      )
    )
  ################################################################
  ##              plot AUC in the two environments              ##
  ################################################################
  figAUC = sumStats %>% ggplot(aes(condition, auc))  +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                       bracket.size = 1, size = 6, label.y = 22,
                       method = "wilcox.test") + 
    ylab("AUC (s)") + myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
    scale_y_continuous(breaks = c(0, 10, 20), limits = c(0, 25)) +
    scale_fill_manual(values = conditionColors) +
    theme(legend.position = "none") + xlab("")
  

  ################################################################
  ##              plot delta AUC in the two environments              ##
  ################################################################
  figDelta = sumStats %>% ggplot(aes(condition, delta)) +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition))  + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                       bracket.size = 1, size = 6, label.y = 9,
                       method = "wilcox.test") +
    ylab(expression(bold(paste(AUC[end] - AUC[start], " (s)")))) +
    xlab(" ") + 
    ylim(c(-10, 12)) +
    scale_fill_manual(values = conditionColors) +
    myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position =  "none") 
  
  ######################################################################
  ##              plot sigma_wtw in the two environments              ##
  ######################################################################
  figSigma = sumStats %>% ggplot(aes(condition, stdWTW))  +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                       bracket.size = 1, size = 6, label.y = 10,
                       method = "wilcox.test") + xlab("") +
    ylab(expression(bold(paste(sigma["WTW"], " (s"^2,")"))))+ 
    myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
    scale_y_continuous(breaks = c(0, 5, 10), limits = c(0, 11)) +
    scale_fill_manual(values = conditionColors) +
    theme(legend.position = "none") + xlab("")

  ##################################################################
  ##                     correlations among task measures        ##
  ##################################################################
  
  pairs(sumStats[, c("auc", "delta", "stdWTW")], condition = sumStats$condition, gap=0,
        lower.panel = my.reg.HP, upper.panel = my.reg.LP, nCmp = 1, lwd = 2) 
  figCorr = recordPlot()
  plot.new()
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
    theme(legend.position = "none") 
  
  ############# return the output ###############
  outputs = list(
    "wtw" = figWTW,
    "auc" = figAUC,
    "sigma" = figSigma,
    "delta" = figDelta,
    "curve" = figCurve,
    "corr" = figCorr,
    "stats_df" = stats_df,
    "contrast_df" = contrast_df
  )
  return(outputs)
}


