# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName, allData = NULL, MFResults = NULL, repOutputs = NULL){
  set.seed(123)
  # create output directories
  # dir.create("../../figures/wtw_exp1")
  # dir.create("../../figures/wtw_exp1/expModelRep/")
  # dir.create(sprintf("../../figures/wtw_exp1/expModelRep/%s",modelName))
  
  # load experiment parameters
  load('expParas.RData')
  
  # load packages and sub functions 
  library("tidyverse")
  library("latex2exp")
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library("patchwork")
  library(lattice)
  source("exp1_expSchematics.R")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  source("subFxs/modelRep.R")
  
  
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  
  # load empirical data 
  if(is.null(allData)){
    allData = loadAllData()
  }
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id[hdrData$stress == "no_stress"]
  nSub = length(ids)

  # check fit
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp1/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  # compare willingness to wait (WTW) from empirical and replicated data
  ## WTW from empirical data 
  if(is.null(MFResults)){
    source("MFAnalysis.R")
    MFResults = MFAnalysis(isTrct = T)
  }
  sumStats = MFResults[['sumStats']]
  muWTWEmp = sumStats$muWTW
  stdWTWEmp = sumStats$stdWTW
  timeWTW_ = matrix(NA, nrow = 60, ncol = length(tGrid))
  trialWTW_ = MFResults$trialWTW_
  for(i in 1 : 60){
    timeWTW_[i,] = MFResults[['timeWTW_']][[i]]
  }
  
  ## replicate data
  if(is.null(repOutputs)){
    repOutputs =  modelRep(modelName, trialData, ids, nRep, T)
    save(repOutputs, file = sprintf("../../genData/wtw_exp1/expModelRep/%s_trct.RData", modelName))
  }

  ## 
  plotData = data.frame(id = ids, mu =  repOutputs$muWTWRep_mu, std = repOutputs$stdWTWRep_mu,
                        empMu = muWTWEmp, empStd = stdWTWEmp,
                        passCheck = passCheck, 
                        condition = sumStats$condition) %>% filter(passCheck)
  
  ## cacl variance explained by time #
  repTimeWTW_ = repOutputs$timeWTW_
  df = list()
  for(condition in c("HP", "LP")){
    ss_ = apply(timeWTW_[sumStats$condition == condition & passCheck,], FUN = function(x) sum((x - mean(x))^2), MARGIN = 2)
    RL_r2_ = vector(length = length(tGrid))
    baseline_r2_ = vector(length = length(tGrid))
    # auc = sumStats$muWTW[sumStats$condition == condition & passCheck]
    baseline = apply(timeWTW_[sumStats$condition == condition & passCheck,], FUN = mean, MARGIN = 1)
    for(tIdx in 1 : length(tGrid)){
      y = timeWTW_[sumStats$condition == condition & passCheck, tIdx]
      x =  repTimeWTW_[tIdx,sumStats$condition == condition & passCheck]
      RL_r2_[tIdx] = summary(lm(y ~ x))$r.squared
      baseline_r2_[tIdx] = summary(lm(y ~ baseline))$r.squared
    }
    thisdf = data.frame(
      ss = rep(ss_, 2),
      r2 = c(RL_r2_, baseline_r2_),
      predictor = rep(c("RL-replicated WTW", "auc"), each = length(tGrid)), 
      time = rep(tGrid, 2),
      condition = rep(condition, length(tGrid) * 2)
    )
    df[[condition]] = thisdf
  }
  plotdf = rbind(df[['HP']], df[['LP']])
  
  plotdf %>% ggplot(aes(time, r2, color = predictor)) + geom_line() + 
    facet_grid(~condition) + 
    myTheme +
    xlab("Task time (s)") + 
    ylab("Variance explained")
  
  ## Yeah I might as well plot at the p_predicted 
  HP_rep_time_WTW = apply(repTimeWTW_[,sumStats$condition == "LP" & passCheck], mean, MARGIN = 1)
  HP_time_WTW = apply(timeWTW_[sumStats$condition == "LP" & passCheck,], mean, MARGIN = 2)
  data.frame(
    wtw = c(HP_time_WTW, HP_rep_time_WTW),
    time = rep(tGrid, 2),
    type = rep(c("Observed", "Replicated"), each = length(tGrid))
  ) %>% ggplot(aes(time, wtw, color = type)) + geom_line()
  
  HP_rep_time_WTW = apply(repTimeWTW_[,sumStats$condition == "HP" & passCheck], mean, MARGIN = 1)
  HP_time_WTW = apply(timeWTW_[sumStats$condition == "HP" & passCheck,], mean, MARGIN = 2)
  data.frame(
    wtw = c(HP_time_WTW, HP_rep_time_WTW),
    time = rep(tGrid, 2),
    type = rep(c("Observed", "Replicated"), each = length(tGrid))
  ) %>% ggplot(aes(time, wtw, color = type)) + geom_line()
  
  ## calc variance explained by time #, taking off the baseline 
  df = list()
  for(condition in c("HP", "LP")){
    # ss_ = apply(timeWTW_[sumStats$condition == condition & passCheck,], FUN = function(x) sum((x - mean(x))^2), MARGIN = 2)
    RL_r2_ = vector(length = length(tGrid))
    baseline_r2_ = vector(length = length(tGrid))
    # auc = sumStats$muWTW[sumStats$condition == condition & passCheck]
    baseline = apply(timeWTW_[sumStats$condition == condition & passCheck,], FUN = mean, MARGIN = 1)
    rep_baseline = apply(repTimeWTW_[,sumStats$condition == condition & passCheck], FUN = mean, MARGIN = 2)
    for(tIdx in 1 : length(tGrid)){
      y = timeWTW_[sumStats$condition == condition & passCheck, tIdx] - baseline
      x =  repTimeWTW_[tIdx,sumStats$condition == condition & passCheck] - baseline
      RL_r2_[tIdx] = summary(lm(y ~ x))$r.squared
      baseline_r2_[tIdx] = summary(lm(y ~ baseline))$r.squared
    }
    thisdf = data.frame(
      ss = rep(ss_, 2),
      r2 = c(RL_r2_, baseline_r2_),
      predictor = rep(c("RL-replicated WTW", "baseline"), each = length(tGrid)), 
      time = rep(tGrid, 2),
      condition = rep(condition, length(tGrid) * 2)
    )
    df[[condition]] = thisdf
  }
  plotdf = rbind(df[['HP']], df[['LP']])
  plotdf %>% ggplot(aes(time, r2, color = predictor)) + geom_line() + 
    facet_grid(~condition) + 
    myTheme +
    xlab("Task time (s)") + 
    ylab("Variance explained")
  
  ########
  # I need to plot local WTW first 
  rep_local_wtw_ = cbind(repOutputs$timeWTW_[1, ], repOutputs$sub_auc_)
  # rep_local_wtw_ =  repOutputs$sub_auc_
  emp_local_wtw_ = cbind(timeWTW_[,1],  MFResults$sub_auc_)
  # emp_local_wtw_ = cbind(timeWTW_[i,1 ], MFResults$sub_auc_)
  
  data.frame(
    local_wtw = as.vector(emp_local_wtw_),
    time = rep(c(0, seq(1.75, 21, by = 3.5)), each = 60),
    condition = rep(sumStats$condition, 7)
  ) %>% group_by(condition, time) %>%
    summarise(mu = mean(local_wtw), se = sd(local_wtw) / sqrt(length(local_wtw))) %>%
    mutate(min = mu - se, max = mu + se) %>%
    ggplot(aes(time, mu, color = condition)) +
    geom_point() + 
    geom_errorbar(aes(x = time, y = mu, ymin = min, ymax = max))
  
  rep_local_wtw_diff_ = matrix(NA, nrow = 60, ncol =6)
  emp_local_wtw_diff_ = matrix(NA, nrow = 60, ncol =6)
  for(i in 1 : 60){
    rep_local_wtw_diff_[i,] = diff(rep_local_wtw_[i,])
    emp_local_wtw_diff_[i,] =  diff(emp_local_wtw_[i,])
  }
  
  # adaptation levels 
  deltaFig = data.frame(
    emp  = emp_local_wtw_[, 7] - emp_local_wtw_[, 2],
    condition = sumStats$condition,
    rep = rep_local_wtw_[, 7] - rep_local_wtw_[, 2],
    passCheck = passCheck
  ) %>% filter(passCheck) %>% ggplot(aes(emp, rep, color = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) +
    scale_color_manual(values = conditionColors) + myTheme +
    xlab("Observed (s)") + ylab("Model-generated (s)") +
    ggtitle(TeX('$AUC_{end} - AUC_{start}$')) + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = c(-16, -8, 0, 8), limits = c(-18, 9)) +
    scale_y_continuous(breaks = c(-16, -8, 0, 8), limits = c(-18, 9))
 
  tmp = data.frame(
    emp  = emp_local_wtw_[, 7] - emp_local_wtw_[, 2],
    condition = sumStats$condition,
    rep = rep_local_wtw_[, 7] - rep_local_wtw_[, 2],
    passCheck = passCheck
  ) 
  
  summary(lm(emp~rep, tmp[tmp$passCheck & tmp$condition == "HP",]))$r.squared
  summary(lm(emp~rep, tmp[tmp$passCheck & tmp$condition == "LP",]))$r.squared
  
  # replicate these values 
  y = emp_local_wtw_[, 7] - emp_local_wtw_[, 2]
  x = sumStats$muWTW
  # It has nothing to do with
  
  summary(lm(y[sumStats$condition == "LP"] ~ x[sumStats$condition == "LP"]))$r.squared
  
  # a more efficient way 
  
  
  # the biggest individual differences were
  tmp = data.frame(
    auc_diff = as.vector(emp_local_wtw_diff_),
    condition = rep(sumStats$condition, 6),
    time = factor(rep(1:6, each = 60))
  ) %>% group_by(condition, time) %>% summarize(ss = sum((auc_diff - mean(auc_diff))^2)) %>% 
    ungroup() 
  ggplot(tmp, aes(time, ss, fill = condition)) + facet_grid(~condition) + geom_bar(stat = "identity") + 
    scale_fill_manual(values = conditionColors) + 
    myTheme + xlab(TeX("$\\Delta AUC$ component")) +
    ylab("Across-subject variance") +
    theme(legend.position =  "none")
  ggsave(file.path("../figures/cmb","exp1_total_vars.eps"), width = 4, height = 2.5)  

  
  data.frame(
    auc_diff = as.vector(emp_local_wtw_diff_),
    condition = rep(sumStats$condition, 6),
    time = factor(rep(seq(1.75, 19.25, by = 3.5), each = 60))
  ) %>% ggplot(aes(time, auc_diff, color =
                     condition)) + geom_boxplot() +
    scale_fill_manual(values = conditionColors) +
    scale_color_manual(values = conditionColors) +
    myTheme + xlab("Task Time") + ylab(TeX("$\\Delta AUC$")) + theme(legend.position = "none") 
  
  # ggsave(file.path("../figures/cmb","exp1_delta_AUC.eps"), fig1, width = 4, height = 4)
  
  # plot delta AUC 1
  data.frame(
    emp = emp_local_wtw_diff_[passCheck, 1],
    rep = rep_local_wtw_diff_[passCheck, 1],
    condition = sumStats$condition[passCheck]
  ) %>% ggplot(aes(emp, rep, color = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0) +
    myTheme + scale_color_manual(values = conditionColors) +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) + 
    xlab("Observed (s)") + ylab("Model-generated (s)") +
    ggtitle(TeX("First component of $\\Delta AUC$")) + 
    coord_fixed() 
  ggsave(file.path("../figures/cmb","exp1_rep_delta_AUC1.eps"),  width = 4, height = 4)
  
  data.frame(
    emp = emp_local_wtw_diff_[passCheck, 2],
    rep = rep_local_wtw_diff_[passCheck, 2],
    condition = sumStats$condition[passCheck]
  ) %>% ggplot(aes(emp, rep, color = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0) +
    myTheme + scale_color_manual(values = conditionColors) +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) + 
    xlab("Observed (s)") + ylab("Model-generated (s)") +
    ggtitle(TeX("Second component of $\\Delta AUC$")) + theme(aspect.ratio=1) +
    scale_x_continuous(limits = c(-10, 7)) +
    scale_y_continuous(limits = c(-10, 7)) 
  ggsave(file.path("../figures/cmb","exp1_rep_delta_AUC2.eps"),  width = 4, height = 4)  
  
  # plot r explained ratio
  auc = sumStats$muWTW
  r2_ = vector(length = 12)
  auc_r2 = vector(length = 12)
  auc_init_r2 =  vector(length = 12)
  
  for(k in 1 : 6){
    filter = sumStats$condition == "HP" & passCheck
    y = emp_local_wtw_diff_[filter, k]
    x =  rep_local_wtw_diff_[filter, k]
    auc = sumStats$muWTW[filter]
    init_wtw = timeWTW_[filter, 1]
    r2_[k] = summary(lm(y ~ x))$r.squared
    auc_r2[k] = summary(lm(y ~  auc))$r.squared
    auc_init_r2[k] = summary(lm(y ~  auc + init_wtw))$r.squared
    
    filter = sumStats$condition == "LP" & passCheck
    y = emp_local_wtw_diff_[filter, k]
    x =  rep_local_wtw_diff_[filter, k]
    auc = sumStats$muWTW[filter]
    init_wtw = timeWTW_[filter, 1]
    r2_[k + 6] = summary(lm(y ~ x))$r.squared
    auc_r2[k + 6] = summary(lm(y ~  auc))$r.squared
    auc_init_r2[k + 6] = summary(lm(y ~  auc + init_wtw))$r.squared
  }  

  data.frame(
    r2 = c(r2_, auc_r2, auc_init_r2),
    model = rep(c("RL", "auc", "auc_init"), each = 12),
    time = rep(c(1:6, 1:6), 3),
    condition = rep(rep(c("HP", "LP"), each = 6),3)
  ) %>%
    ggplot(aes(time, r2 * 100, color = model)) +
    facet_grid(~condition) + geom_line() + geom_point() + myTheme +
    ylab("Variance explained (%)") + xlab(TeX("$\\Delta AUC$ components")) +
    scale_color_manual(values = c("grey", "#737373",  "#0570b0"))  + 
    scale_x_continuous(breaks = 1:6)
  sggsave(file.path("../figures/cmb","exp1_variance_explained.eps"),  width = 4, height = 2.5) 
  
  
  # I already missed the very first part I assume, so this is not the best way
  data.frame(
    auc = as.vector(emp_local_wtw_),
    condition = rep(sumStats$condition, 7),
    time = factor(rep(1:7, each = 60))
  ) %>% ggplot(aes(time, auc, color = condition)) + 
    geom_boxplot()
  

  # I don't know what I want to do
  
  ##########################
  ##      calc RSME      ##
  ##########################
  summary(lm(plotData$empMu ~ plotData$mu))$r.squared # I can also calculate them separately ...
  
  mu_sqerr = (plotData$empMu - plotData$mu)^2
  std_sqerr = (plotData$empStd - plotData$std)^2
  sqerr_df = data.frame(
    id = plotData$id,
    mu_sqerr = mu_sqerr,
    std_sqerr = std_sqerr)
  
  #################################################################
  ##      compare observed and replicated AUC and sigma_wtw      ##
  #################################################################
  
  # plot to compare average willingess to wait
  aucFig = plotData %>%
  ggplot(aes(empMu, mu)) + 
  geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
  geom_abline(slope = 1, intercept = 0)  + 
  ylab("Model-generated AUC (s)") + xlab("Observed AUC (s)") +
  myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(breaks = c(0, 20), limits = c(-1, 21)) + 
  scale_y_continuous(breaks = c(0, 20), limits = c(-1, 21)) +
  scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + 
    ggtitle(sprintf("%s, N = %d", modelName, sum(passCheck)))
  
  summary(lm(empMu~mu, plotData[plotData$passCheck & plotData$condition == "HP",]))$r.squared
  summary(lm(empMu~mu, plotData[plotData$passCheck & plotData$condition == "LP",]))$r.squared
  
  aucFig = plotData %>%
    ggplot(aes(empMu, mu)) +
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) +
    geom_abline(slope = 1, intercept = 0)  +
    ylab("Model-generated(s)") + xlab("Observed(s)") +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 20), limits = c(-1, 21)) +
    scale_y_continuous(breaks = c(0, 20), limits = c(-1, 21)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") +
    ggtitle('AUC')

  
  
  
  ## plot to compare std willingess to wait
  cipFig = plotData %>%
  ggplot(aes(empStd, std, shape = condition)) + 
  geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
  geom_abline(slope = 1, intercept = 0) +
  ylab(expression(bold(paste("Model-generated ", sigma["WTW"], " (s"^2,")")))) +
    xlab(expression(bold(paste("Observed ", sigma["WTW"], " (s"^2,")")))) +
  myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_x_continuous(breaks = c(0, 10), limits = c(0, 10)) + 
  scale_y_continuous(breaks = c(0, 10), limits = c(0, 10)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none")

  # combind
  rep = aucFig / cipFig +  plot_annotation(title = sprintf("%s, n = %d", modelName, sum(passCheck)),
                                     theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  fileName = sprintf("../../figures/wtw_exp1/expModelRep/%s/rep.eps", modelName) 
  
  ##################################################################
  ##                     example participant                      ##
  ##################################################################
  
  # model generated 
  # repTrialData = repOutputs$repTrialData
  # repNo = repOutputs$repNo
  # sIdx = 1
  # thisRepTrialData = repTrialData[[repNo[1, sIdx]]]
  # thisRepTrialData = data.frame(thisRepTrialData[1:6])
  # modelFig = trialPlots(thisRepTrialData) +  
  #   ggtitle("Model-generated") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none")
  # 
  # # observed
  # thisTrialData = trialData[[ids[sIdx]]]
  # thisTrialData  = thisTrialData %>% filter(trialStartTime <=  blockSec - max(delayMaxs))
  # thisTrialData = block2session(thisTrialData)
  # empFig = trialPlots(thisTrialData) + ggtitle("Observed") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "none")
  # example = modelFig / empFig +  plot_annotation(title = sprintf("%s", modelName),
  #                                      theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))
  # ggsave(sprintf("../../figures/wtw_exp1/expModelRep/%s/example.eps", modelName), example, width = 6, height = 8)
  # 
  
  ################# return figure outputs ###############
  # outputs = list("rep" = rep, "example" = example)
  outputs = list("rep" = rep, "example" = example, "sqerr_df" = sqerr_df)
  return(outputs)
  
}

