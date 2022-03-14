MFPlot = function(){
  library('dplyr')
  library("tidyr")
  library("ggplot2")
  library("latex2exp")
  library("ggpubr")
  source("subFxs/plotThemes.R")
  source('MFAnalysis.R')
  source("exp3_expSchematics.R")
  
  # normative analysis 
  iti = 2
  normResults = expSchematics(0, iti, F)
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # output dir
  dir.create('../../figures/wtw_exp3/MFplot')
  
  # load experiment parameters
  load("expParas.RData")
  
  ###################################################################
  ##              plot WTW timecourses in the two environments             ##
  ###################################################################
  nBlock = 4
  MFResults = MFAnalysis(isTrct = F)
  blockStats = MFResults[['blockStats']]
  timeWTW_ = MFResults[['timeWTW_']]
  nSub = length(unique(blockStats$id))
  cbal = blockStats$cbal; cbal[cbal == 1] = "HPLP"; cbal[cbal == 2] = "LPHP"
  # define background colors
  greyData = data.frame(
    xmin = (1:nBlock * blockSec - max(delayMaxs)) / 60, xmax = 1:nBlock * blockSec / 60,
    block = rep(c("Block1-2", "Block3-4"), each = 2)
  )
  
  # plot
  plotData = data.frame(wtw = unlist(timeWTW_),
                        time = rep(seq(0, blockSec * nBlock - 1, by = 2) / 60, nSub),
                        condition = rep(blockStats$condition, each = length(tGrid)),
                        blockIdx = rep(rep(1 : nBlock, length(tGrid)), nSub),
                        cbal = rep(ifelse(blockStats$cbal == 1, "HP First", "LP First"), each = length(tGrid))) %>%
    group_by(time, cbal, condition) %>% 
    dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                     min = mu- se, max = mu + se) %>% ungroup()
  # insert NA values to use two colors to draw one line 
  na_df = plotData
  na_df$mu = NA
  na_df$min = NA
  na_df$max = NA
  na_df$condition = ifelse(plotData$condition == "HP", "LP", "HP")
  plotData= rbind(plotData, na_df)

  figWTW = plotData %>%
    ggplot(aes(time, mu, color = condition, fill = condition)) +
      geom_rect(data = greyData, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 20), inherit.aes = F, fill = "#d9d9d9") +
    geom_ribbon(aes(ymin=min, ymax=max), color = NA) +
    geom_line(size = 1) +
    xlab("Task time (min)") + ylab("WTW (s)") + 
    myTheme + facet_grid(cbal~.)  +
    scale_fill_manual(values = c("#7fbf7b", "#af8dc3")) + 
    scale_color_manual(values = conditionColors) +
    theme(legend.position =  "none")

  
  ###################################################################
  ##              compare AUC in the two environments             ##
  ###################################################################
  MFResults = MFAnalysis(isTrct = T)
  # plot the first two blocks and the second two blocks separately 
  blockStats = MFResults[['blockStats']]
  blockStats$block = ifelse(blockStats$blockNum <= 2, "Block1-2", "Block3-4")
  survCurve_ = MFResults$survCurve_
  blockStats %>% group_by(condition, block) %>%
    summarise(median(auc), IQR(auc), median(stdWTW), IQR(stdWTW))
  
  wilcox.test( blockStats[blockStats$condition == "HP" & blockStats$block == "Block1-2", "auc"],
               blockStats[blockStats$condition == "LP" & blockStats$block == "Block1-2", "auc"], paired = T)
  
  wilcox.test( blockStats[blockStats$condition == "HP" & blockStats$block == "Block3-4", "auc"],
               blockStats[blockStats$condition == "LP" & blockStats$block == "Block3-4", "auc"], paired = T)  
  ###################################################################
  ##              plot AUC in the two environments             ##
  ###################################################################
  figAUC = ggdotplot(blockStats, x = "block", y = "auc", fill = "condition") + 
    ggpubr::stat_compare_means(aes(group = condition, label = ..p.signif..),
                               method = "wilcox.test", paired = T, label.y = 31) +
    scale_fill_manual(values = conditionColors) + 
    ylab("AUC (s)") + xlab("") + ylim(c(0, 32)) + 
    theme(legend.position = "None")
    
  ###################################################################
  ##              compare sigma_wtw in the two environments             ##
  ###################################################################
  wilcox.test( blockStats[blockStats$condition == "HP" & blockStats$block == "Block1-2", "stdWTW"],
               blockStats[blockStats$condition == "LP" & blockStats$block == "Block1-2", "stdWTW"], paired = T)
  
  wilcox.test( blockStats[blockStats$condition == "HP" & blockStats$block == "Block3-4", "stdWTW"],
               blockStats[blockStats$condition == "LP" & blockStats$block == "Block3-4", "stdWTW"], paired = T)  
  
  ###################################################################
  ##              plot sigma_wtw in the two environments             ##
  ###################################################################
  figSigma = ggdotplot(blockStats, x = "block", y = "stdWTW", fill = "condition") + 
    ggpubr::stat_compare_means(aes(group = condition, label = ..p.signif..),
                               method = "wilcox.test", paired = T, size = 6, label.y = 16) +
    scale_fill_manual(values = conditionColors) + 
    ylab(expression(bold(paste(sigma["WTW"], " (s"^2, ")")))) +
    xlab("") + ylim(c(0, 18)) + 
    theme(legend.position = "None")

  ###################################################################
  ##              plot adaptation in the two environments             ##
  ###################################################################
  # plot sub_auc first I guess
  figDelta = data.frame(
    delta =  as.vector(t(MFResults$sub_auc_[, c(2, 4, 6, 8)] - MFResults$sub_auc_[, c(1, 3, 5, 7)])),
    condition = blockStats$condition,
    block = blockStats$block,
    id = blockStats$id,
    cbal = blockStats$cbal,
    blockNum = blockStats$blockNum
  ) %>% ggdotplot(x = "block", y = "delta", fill = "condition") + 
    ggpubr::stat_compare_means(aes(group = condition, label = ..p.signif..),
                               method = "wilcox.test", paired = T,  label.y = 13) +
    scale_fill_manual(values = conditionColors) +
    ylab(expression(bold(paste(AUC[end]-AUC[start], " (s)")))) +
    xlab("")  + ylim(c(-15, 15)) +
    theme(legend.position = "None")
    
  ###################################################################
  ##              plot  survival curves            ##
  ###################################################################
  # prepare data
  plotData = data.frame(survCurve = unlist(survCurve_),
                        time = rep(kmGrid, nSub * nCondition * 2),
                        condition = rep(blockStats$condition, each = length(kmGrid)),
                        blockNum = rep(blockStats$blockNum, each = length(kmGrid)),
                        cbal = rep(blockStats$cbal, each = length(kmGrid))) %>%
    mutate(block = ifelse(blockNum <= 2, "Block1-2", "Block3-4")) 
  # optimal strategy
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
  #optim = rbind(optim, optim)
  #optim$block = rep(c("Block1-2", "Block3-4"), each = nrow(optim) / 2)

  figCurve = plotData %>% group_by(condition, time, block) %>%
    dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),
                     min = mu- se, max = mu + se) %>%
    ggplot(aes(time, mu, color = condition, fill = condition)) + facet_grid(.~block) +
    geom_ribbon(aes(time, ymin = min, ymax = max), color = NA)  +
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
  
  # figCurve = plotData %>% filter(block == "Block1-2") %>% 
  #   mutate(cbal = ifelse(cbal == 1, "HPLP", "LPHP")) %>%
  #   group_by(condition, time, cbal) %>%
  #   dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),
  #                    min = mu- se, max = mu + se) %>%
  #   ggplot(aes(time, mu, color = condition, fill = condition)) + 
  #   facet_grid(cbal~.) + 
  #   geom_ribbon(aes(time, ymin = min, ymax = max), color = NA)  +
  #   geom_line() +
  #   geom_line(data = optim, aes(t, surv, color = condition, linetype = condition, alpha = condition), size = 1.2) +
  #   geom_line(data = data.frame(t = kmGrid[kmGrid > 2],surv = 1),
  #             aes(t, surv), color = conditionColors[1], size = 1.2, inherit.aes = F, alpha = 0.8) + 
  #   scale_fill_manual(values = c("#7fbf7b", "#af8dc3")) +
  #   scale_color_manual(values = conditionColors) +
  #   scale_linetype_manual(values = c("solid", "dotted")) +
  #   scale_alpha_manual(values = c(0.8, 1))+
  #   xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
  #   theme(legend.position = "none")
    
  ############# return the output ###############
    outputs = list(
      "wtw" = figWTW,
      "auc" = figAUC,
      "delta" = figDelta,
      "sigma" = figSigma,
      "curve" = figCurve
    )
    return(outputs)
}



