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
  nSub = length(unique(sumStats$ID))
  ## use background color to distinguish used and excluded data 
  colorData = data.frame(
    xmin = c(0, blockSec), 
    # xmax = c(blockSec - max(delayMaxs), nBlock * blockSec - max(delayMaxs)),
    xmax = c(blockSec, 2 * blockSec),
    condition = c(conditions[2], conditions[1])
  )
  greyData = data.frame(
    xmin = 1:2 * blockSec - max(delayMaxs), xmax = 1:2 * blockSec
  )
  
  plotData = data.frame(wtw = unlist(timeWTW_),
             time = rep(seq(0, blockSec * nBlock -1, by = 2), nSub),
             condition = rep(rep(c("LP", "HP"), each = length(tGrid))), nSub) %>%
    group_by(condition, time) %>%
    dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                     min = mu- se, max = mu + se) %>% ungroup() 
  # reset at the block onset
  plotData$mu[plotData$time %in% (1:2 * blockSec - 1)] = NA
  plotData$min[plotData$time %in% (1:2 * blockSec - 1)] = NA
  plotData$max[plotData$time %in% (1:2 * blockSec - 1)] = NA
  
  # plot
  # plotData %>%
  #   ggplot(aes(time, mu))  + 
  #   geom_rect(data = colorData, aes(xmin = xmin, xmax = xmax, fill = condition),
  #             ymin = 0, ymax = 20, inherit.aes = F, alpha = 0.7) + 
  #   geom_rect(data = greyData, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 20,
  #             color = "red", inherit.aes = F, fill = NA) +
  #   geom_ribbon(aes(ymin=min, ymax=max),  color = NA, fill = "#878787") +
  #   geom_line(size = 0.5) +
  #   xlab("Task time (min)") + ylab("WTW (s)") + 
  #   scale_y_continuous(breaks = c(0, 8, 16), limits = c(0, 17)) +
  #   myTheme + scale_x_continuous(breaks = c(0, 600, 1200), labels = c(0, 600, 1200) / 60 ) +
  #   theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  #   scale_color_manual(values = conditionColors) +
  #   scale_fill_manual(values = c("#e6f5d0", "#c2a5cf")) +
  #   theme(legend.position = "none")
  
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
  ##              compare stats in the two environments              ##
  ################################################################
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  
  auc_stats = sumStats %>% group_by(condition) %>% summarise(median = round(median(auc), 2),
                                                             IQR = round(IQR(auc), 2),
                                                             pvalue = round(wilcox.test(auc)$p.value, 3)) %>%
    mutate(lower = median - IQR, upper = median + IQR)
  
  auc_contrast_p = round(wilcox.test( sumStats[sumStats$condition == "HP", "auc"],
                                sumStats[sumStats$condition == "LP", "auc"],paired = T)$p.value, 3)
  
  sigma_stats = sumStats %>% group_by(condition) %>% summarise(median = round(median(stdWTW),2),
                                                               IQR = round(IQR(stdWTW), 2),
                                                               pvalue = round(wilcox.test(stdWTW)$p.value, 3)) %>%
    mutate(lower = median - IQR, upper = median + IQR)
  
  sigma_contrast_p = round(wilcox.test( sumStats[sumStats$condition == "HP", "stdWTW"],
                                  sumStats[sumStats$condition == "LP", "stdWTW"], paired = T)$p.value, 3)
  
  delta_contrast_p = round(wilcox.test(MFResults$sub_auc_[,2] - MFResults$sub_auc_[,1],
                      MFResults$sub_auc_[,4] - MFResults$sub_auc_[,3], paired = T)$p.value, 3)
  
  delta_stats = data.frame(
    condition = c("LP", "HP"),
    median = c(round(median(MFResults$sub_auc_[,2] - MFResults$sub_auc_[,1]), 2),
               round(median(MFResults$sub_auc_[,4] - MFResults$sub_auc_[,3]), 2)),
    IQR = c(round(IQR(MFResults$sub_auc_[,2] - MFResults$sub_auc_[,1]), 2),
            round(IQR(MFResults$sub_auc_[,4] - MFResults$sub_auc_[,3]), 2)),
    pvalue = c(round(wilcox.test(MFResults$sub_auc_[,2] - MFResults$sub_auc_[,1])$p.value, 3),
               round(wilcox.test(MFResults$sub_auc_[,4] - MFResults$sub_auc_[,3])$p.value, 3))
  ) %>% 
    mutate(lower = median - IQR, upper = median + IQR)
  
  stats_df = rbind(auc_stats, delta_stats, sigma_stats)
  stats_df['variable'] = rep(c("auc",  "delta", "stdWTW"), each = 2) 
  stats_df = stats_df %>% arrange(variable, condition)
  
  contrast_df = data.frame(
    variable = c("auc", "delta", "stdWTW"),
    pvalue = c(auc_contrast_p, delta_contrast_p, sigma_contrast_p)
  )
  ################################################################
  ##              plot AUC in the two environments              ##
  ################################################################
 figAUC = sumStats %>% ggplot(aes(condition, auc)) +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), 
                       label.x = 1.5,
                       bracket.size = 1, size = 6, label.y = 18, paired = T,
                       method = "wilcox.test") + 
    ylab("AUC (s)") + myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
    scale_y_continuous(breaks = c(0, 8, 16), limits = c(0, 20)) +
    scale_fill_manual(values = conditionColors) +
    theme(legend.position = "none") + xlab("")
  
  ################################################################
  ##              plot sigma_wtw in the two environments              ##
  ################################################################

  # plot
  figSigma = sumStats %>% ggplot(aes(condition, stdWTW)) +
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), 
                       label.x = 1.5,
                       bracket.size = 1, size = 6, label.y = 6.8, paired = T,
                       method = "wilcox.test") + 
    ylab(expression(bold(paste(sigma["WTW"], " (s"^2, ")")))) + myTheme  + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
    scale_y_continuous(breaks = c(0, 3, 6), limits = c(0, 7)) +
    scale_fill_manual(values = conditionColors) +
    theme(legend.position = "none") + xlab("")
  
  ################################################################
  ##              plot delta AUC in the two environments        ##
  ################################################################
  
  df_delta = 
    data.frame(
      delta = c(MFResults$sub_auc_[,2] - MFResults$sub_auc_[,1], MFResults$sub_auc_[,4] - MFResults$sub_auc_[,3]),
      condition = rep(c("LP", "HP"), each = nSub)
    )
  figDelta = df_delta %>% ggplot(aes(condition, delta)) + 
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
    stat_compare_means(comparisons = list(c("HP", "LP")),
                       aes(label = ..p.signif..), 
                       label.x = 1.5,
                       bracket.size = 1, size = 6, label.y = 6.8,
                       paired = T, method = "wilcox.test") + 
    ylab(expression(bold(paste(AUC[end] - AUC[start], " (s)")))) + myTheme + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
    scale_y_continuous(breaks = c(-8, -4, 0, 4, 8), limits = c(-8, 8)) +
    scale_fill_manual(values = conditionColors) +
    theme(legend.position = "none") + xlab("")

  # figDelta = data.frame(
  #   delta_HP =  HP_end_vs_start,
  #   delta_LP = LP_end_vs_start 
  # ) %>% ggplot(aes(delta_LP, delta_HP)) + 
  #   geom_point(size = 5, shape = 21, stroke =1) +
  #   geom_abline(slope = 1, intercept = 0) + 
  #   ggtitle(expression(bold(paste(AUC[end] - AUC[start], " (s)")))) + 
  #   xlab("LP") +
  #   ylab("HP") + 
  #   annotate("text", x = 2, y = -4, label = sprintf('p = %0.3f*', wTest$p.value)) +
  #   xlim(c(-8,5)) + ylim(c(-8,5)) + myTheme +
  #   theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # wilcox.test(deltaWTW[sumStats$condition == "HP"])
  # wilcox.test(deltaWTW[sumStats$condition == "LP"])
  # cor.test(deltaWTW[sumStats$condition == "HP"], sumStats$auc[sumStats$condition == "HP"], method=c( "spearman"))
  # cor.test(deltaWTW[sumStats$condition == "LP"], sumStats$auc[sumStats$condition == "LP"], method=c( "spearman"))
  
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
    "delta" = figDelta,
    "curve" = figCurve,
    "stats_df" = stats_df,
    "contrast_df" = contrast_df
  )
  return(outputs)  
}
