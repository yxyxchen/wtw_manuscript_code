

b = kmsc(trialData[['31']], min(delayMaxs), F, grid = seq(0, 20, length.out = 20))$kmOnGrid
c = kmsc(trialData[["17"]], min(delayMaxs), F, grid = seq(0, 20, length.out = 20))$kmOnGrid

figAUC = data.frame(
  ts = rep(seq(0, 20, length.out = 20), 2),
  st = c(b, c),
  id = rep(c("31", "17"), each = 20)
) %>% ggplot(aes(ts, st, color = id)) + geom_line(size = 3) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
  theme(legend.position = 'None', text=element_text(face = "bold", size = 22)) + 
  scale_color_manual(values = c("#d8b365", "#5ab4ac"))
ggsave("~/Downloads/figAUC.png", figAUC, width = 4, height = 4)

###################################
b = c(1, 1, 1, 1, 0.98809524, 
      0.96222349, 0.93238164, 0.86989942, 0.81793093, 0.74389520, 
      0.56346534, 0.45293945, 0.34756812, 0.15798551, 0.07899275, 
      0.02899275, 0.01899275, 0.012, 0.011, 0.011)
c = c(1, 0.985, 0.953, 0.850, 0.820,
      0.810, 0.760, 0.742, 0.641, 0.531,
      0.521, 0.450, 0.430, 0.390, 0.330, 
      0.252, 0.202, 0.153, 0.053, 0.042)
figSigma = data.frame(
  ts = rep(seq(0, 20, length.out = 20), 2),
  st = c(b, c),
  id = rep(c("31", "17"), each = 20)
) %>% ggplot(aes(ts, st, color = id)) + geom_line(size = 3) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
  theme(legend.position = 'None', text=element_text(face = "bold", size = 22)) + 
  scale_color_manual(values = c("#d8b365", "#5ab4ac"))
ggsave("~/Downloads/figSigma.png", figSigma, width = 4, height = 4)

###################################
thisTrialData = trialData[['1']]
b_first = kmsc(thisTrialData[thisTrialData$sellTime < 210 & thisTrialData$blockNum == 1,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid
b_second = kmsc(thisTrialData[thisTrialData$sellTime >= 210 & thisTrialData$blockNum == 3,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid


thisTrialData = trialData[['2']]
c_first = kmsc(thisTrialData[thisTrialData$sellTime < 210 & thisTrialData$blockNum == 1,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid
c_second = kmsc(thisTrialData[thisTrialData$sellTime >= 210 & thisTrialData$blockNum == 3,], min(delayMaxs), F, seq(0, 20, length.out = 20))$kmOnGrid
figDelta = data.frame(
  ts = rep(seq(0, 20, length.out = 20), 4),
  st = c(b_first, b_second, c_first, c_second),
  chunck = rep(rep(c("start", "end"), each = 20), 2),
  id = rep(c("1", "2"), each = 40),
  color = rep(c("1", "2", "3", "4"), each = 20)
) %>% ggplot(aes(ts, st, color = color)) + geom_line(size = 3) + 
  facet_grid(~id) + myTheme +
  scale_color_manual(values = c("#d8b365", "#8c510a", "#5ab4ac", "#01665e")) +
  theme(legend.title = element_blank(), text=element_text(face = "bold", size = 22)) +
  xlab("Elapsed time (s)") + ylab("Survival rate") 
ggsave("~/Downloads/figDelta.eps", figDelta, width = 8, height = 4)
  

