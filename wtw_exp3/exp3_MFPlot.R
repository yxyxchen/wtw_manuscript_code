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
  greenData = data.frame(
    xmin = c(0, blockSec * 2, blockSec, blockSec * 3),
    xmax = c(blockSec, blockSec* nBlock, blockSec * 2, blockSec * 4) - max(delayMaxs),
    cbal = c("HPLP", "HPLP", "LPHP", "LPHP")
  )
  purpleData = data.frame(
    xmin = c(blockSec,  blockSec * 3, 0, blockSec * 2),
    xmax = c(blockSec* 2, blockSec * 4, blockSec, blockSec* 3) - max(delayMaxs),
    cbal = c("HPLP", "HPLP", "LPHP", "LPHP")
  )
  greyData = data.frame(
    xmin = 1:nBlock * blockSec - max(delayMaxs), xmax = 1:nBlock * blockSec
  )
  
  # plot
  ## add NA values
  plotData = data.frame(wtw = unlist(timeWTW_),
                        time = rep(seq(0, blockSec * nBlock - 1, by = 1), nSub),
                        condition = rep(blockStats$condition, each = length(tGrid)),
                        blockIdx = rep(rep(1 : nBlock, length(tGrid)), nSub),
                        cbal = rep(cbal, each = length(tGrid))) %>%
    group_by(time, cbal, condition) %>% 
    dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                     min = mu- se, max = mu + se) %>% ungroup()
  naData = data.frame(
    time = rep((1 : (blockSec* nBlock))- 1, 2),
    condition = rep(factor(c("LP", "HP", "LP", "HP", "HP", "LP", "HP", "LP")), each = blockSec),
    cbal = rep(factor(c("HPLP", "LPHP")), each = blockSec * nBlock),
    mu = rep(NA, blockSec * nBlock * 2),
    se = rep(NA, blockSec * nBlock * 2),
    min = rep(NA, blockSec * nBlock * 2),
    max = rep(NA, blockSec * nBlock * 2)
  )
  newPlotData = rbind(naData, plotData)
  figWTW = newPlotData %>%
    ggplot(aes(time, mu, fill = condition, color = condition)) +
    geom_rect(data = greyData, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 20), inherit.aes = F, fill = "#d9d9d9") +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.5, color = NA) +
    geom_line(size = 1) +
    xlab("Task time (min)") + ylab("WTW (s)") + 
    myTheme +
    facet_grid(.~cbal) +
    scale_x_continuous(breaks = seq(0, blockSec * nBlock, by = blockSec),
                       labels = seq(0, blockSec * nBlock, by = blockSec) / 60) +
    scale_y_continuous(breaks = c(0,  10, 20), limits = c(0, 20))  +
    scale_fill_manual(values = conditionColors) + 
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
    summarise(median(muWTW), IQR(muWTW), median(stdWTW), IQR(stdWTW))
  
  wilcox.test( blockStats[blockStats$condition == "HP" & blockStats$block == "Block1-2", "muWTW"],
               blockStats[blockStats$condition == "LP" & blockStats$block == "Block1-2", "muWTW"], paired = T)
  
  wilcox.test( blockStats[blockStats$condition == "HP" & blockStats$block == "Block3-4", "muWTW"],
               blockStats[blockStats$condition == "LP" & blockStats$block == "Block3-4", "muWTW"], paired = T)  
  ###################################################################
  ##              plot AUC in the two environments             ##
  ###################################################################
  figAUC = data.frame(muWTWHP = blockStats$muWTW[blockStats$condition == 'HP'],
             muWTWLP = blockStats$muWTW[blockStats$condition == 'LP'],
             cbal = blockStats$cbal[blockStats$condition == "HP"],
             block = ifelse(blockStats$blockNum <= 2, "Block1-2", "Block3-4")) %>%
    ggplot(aes(muWTWLP, muWTWHP)) +
    geom_point(size = 5, shape = 21, stroke =1) +
    geom_abline(slope = 1, intercept = 0)  +
    xlab("LP muAUC / (s)") + ylab("HP muAUC / (s)") + 
    myTheme + xlim(c(-1,31)) + ylim(c(-1,31)) + 
    xlab("LP AUC (s)") + ylab("HP AUC (s)") + facet_grid(~block)
  
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
  figSigma = data.frame(stdWTWHP = blockStats$stdWTW[blockStats$condition == 'HP'],
             stdWTWLP = blockStats$stdWTW[blockStats$condition == 'LP'],
             cbal = blockStats$cbal[blockStats$condition == "HP"],
             block = ifelse(blockStats$blockNum <= 2, "Block1-2", "Block3-4")) %>%
    ggplot(aes(stdWTWLP, stdWTWHP)) +
    geom_point(size = 5, shape = 21, stroke =1) +
    geom_abline(slope = 1, intercept = 0)  +
    xlab(expression(bold(paste("LP ", sigma["WTW"], " (s"^2,")")))) +
    ylab(expression(bold(paste("HP ", sigma["WTW"], " (s"^2,")")))) + 
    myTheme + xlim(c(-1,16)) + ylim(c(-1,16)) + facet_grid(~block)
  
  ###################################################################
  ##              plot  survival curves            ##
  ###################################################################
  # prepare data
  plotData = data.frame(survCurve = unlist(survCurve_),
                        time = rep(kmGrid, nSub * nCondition * 2),
                        condition = rep(blockStats$condition, each = length(kmGrid)),
                        blockNum = rep(blockStats$blockNum, each = length(kmGrid))) %>%
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
  optim = rbind(optim, optim)
  optim$block = rep(c("Block1-2", "Block3-4"), each = nrow(optim) / 2)

    ## plot
  figCurve = plotData %>% group_by(condition, time, block) %>%
    dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),
                     min = mu- se, max = mu + se) %>%
    ggplot(aes(time, mu, color = condition, fill = condition)) + geom_line() +
    geom_ribbon(aes(time, ymin = min, ymax = max), alpha = 0.5, color = NA) +
    facet_grid(~block) +
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



