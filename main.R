#-# the main script to run all analyses #-# 

# libraries
library("tidyverse")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(extrafont)
library(patchwork)
library(figpatch)
# extrafont::font_import()


# settings 
exps = c("exp1", "exp2", "exp3")
nExp = length(exps)
pwd = getwd()
wds = c("wtw_exp1", "wtw_exp2", "wtw_exp3")
models = c("QL1", "QL2", "RL1", "RL2", "naive", "omni")
nModel = length(models)
iti = 2 # inter-trial interval in seconds 
smallReward = 0 # initial token value

# create the output fir 
dir.create("../figures")
dir.create("../genData")
dir.create("../figures/cmb")
#################################################################
##                    experiment schematics                    ##
#################################################################
figs_ = vector("list", length = nExp )
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_expSchematics.R", i))
  outs = expSchematics(smallReward, iti, T)
  figs_[[i]] = outs$figs
}
# assemble the figures
setwd(pwd)
expCmb = (figs_[[1]][['pdf']] | figs_[[1]][['rt']]) / (figs_[[2]][['pdf']] | figs_[[2]][['rt']]) /  (figs_[[3]][['pdf']] | figs_[[3]][['rt']]) 
ggsave(file.path("../figures/cmb", "exp.eps"), expCmb, width = 9, height = 9)


#################################################################
##                     model free analysis                     ##
#################################################################
outs_ = vector("list", length = nExp )
for(i in 1 : 3){
  outs = vector("list", length = nExp) # initialize the output 
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_MFPlot.R", i))
  outs = MFPlot()
  outs_[[i]] = outs
}
# plot behavioral results 
setwd(pwd)
figMF12 = ( outs_[[1]][['curve']] | outs_[[1]][['wtw']] | outs_[[1]][['auc']] | outs_[[1]][['delta']] | outs_[[1]][['sigma']]) / ( outs_[[2]][['curve']] | outs_[[2]][['wtw']] | outs_[[2]][['auc']] | outs_[[2]][['delta']] | outs_[[2]][['sigma']]) + plot_annotation(tag_levels = "a")
ggsave(file.path("../figures/cmb","mf12.pdf"), figMF12 , width = 20, height = 8)

figMF3 = ( outs_[[3]][['curve']] | outs_[[3]][['wtw']] ) / ( outs_[[3]][['auc']] | outs_[[3]][['delta']] | outs_[[3]][['sigma']] ) + plot_annotation(tag_levels = "a")
ggsave(file.path("../figures/cmb","mf3.eps"), figMF3 , width = 12, height = 8)

# plot correlations among task measures
pdf("../figures/cmb/exp3_task_corr1.pdf", width = 5, height = 5) 
outs_[[3]][['corr1']]
dev.off()

pdf("../figures/cmb/exp1_task_corr.pdf", width = 5, height = 5) 
outs_[[1]][['corr']]
dev.off()

pdf("../figures/cmb/exp2_task_corr.pdf", width = 5, height = 5) 
outs_[[2]][['corr']]
dev.off()

pdf("../figures/cmb/exp3_task_corr2.pdf", width = 5, height = 5) 
outs_[[3]][['corr2']]
dev.off()

# plot correlations among task measures and self-report measures
source("wtw_exp1/exp1_taskTraitCorr.R")
setwd("wtw_exp1")
figs = taskTraitCorr()
figTaskTrait = figs[['figTraitAUC']] / figs[['figTraitDelta']] / figs[['figTraitSigma']] + plot_annotation(tag_levels = "a")
setwd(pwd)
ggsave(file.path("../figures/cmb","exp1_task_trait.pdf"), figTaskTrait, width = 12, height = 15)

##################################################################
##                 Performance check simulation                 ##
##################################################################
setwd("./wtw_simulation")
load("expParas.RData")
source("performCheck.R")
figs_ = performCheck()
# plot AUC 
setwd(pwd)
figAUC = (figs_[[1]][['learn']] | figs_[[2]]['learn'])
ggsave(file.path("..", "figures", "cmb", "exante_learn_curve.eps"), figAUC, width = 6, height = 3)

# plot the relative value of waiting and the net return of waiting 
normResults = expSchematics(smallReward, iti)
figNetReturn = data.frame(
  net_return = unlist(normResults$subjectValues),
  time = c(seq(0, delayMaxs[1], by = 0.1), seq(0, delayMaxs[2], by = 0.1)),
  condition = rep(c("HP", "LP"), as.numeric(sapply(normResults$subjectValues, length)))
) %>% filter(time < 20) %>%
  ggplot(aes(time, net_return, color = condition)) + geom_point() +
  facet_grid(~condition) + scale_color_manual(values = conditionColors) + 
  xlab("Elapsed time (s)") + ylab("Net return of waiting (¢)") + myTheme + 
  theme(legend.position = "None") + ggtitle("Normative solution")
figRV = (figNetReturn | figs_[[1]][['rv']] | figs_[[2]]['rv']) + plot_annotation(tag_levels = "a")
ggsave(file.path("..", "figures", "cmb", "exante_rv.eps"), figRV, width = 12, height = 4.5)


figSnippet = (plot_spacer() | figs_[[1]][["Gs_"]][1] | figs_[[1]][["Gs_"]][2] | figs_[[1]][["Gs_"]][3]) / (figs_[[1]][["values_"]][[1]] | figs_[[1]][["values_"]][[2]] | figs_[[1]][["values_"]][[3]] | figs_[[1]][["values_"]][[4]]) 
ggsave(file.path("..", "figures", "cmb", "exante_snippet.eps"), figSnippet , width = 12, height = 6)
figPrior = figs_[[1]][["prior"]]
ggsave(file.path("..", "figures", "cmb", "prior.eps"), figPrior, width = 5, height = 5)

#################################################################
##                 Parameter-effect simulation                 ##
#################################################################
setwd("./wtw_simulation")
source("paraEffect.R")
figs_ = paraEffect()
# assemble figures
setwd(pwd)
figPostHoc = figs_[[1]]  
ggsave(file.path("..", "figures", "cmb", "posthoc_para_effect.eps"), figPostHoc , width = 8, height = 4)


###################################################################
##                 Fit models to empirical data               #####
###################################################################
## Warnings: this model fitting process can take hours. 
## To save time, you can directly download the outputs from the link below and move to the next step:
## https://www.dropbox.com/sh/a2yqj3f21fkzj3r/AADly4VA7SMeaBWY5Nkeo83Ga?dl=0
parallel = F # When run on a local PC, you can use set this variabel to true to turn on parallel computing
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_expModelFit.R", i))
  for(modelName in models){
    # fit all participants
    print(sprintf("Model fitting in Exp.%d", i))
    print(sprintf("Model fitting results saved at ../genData/wtw_exp%d/expModelFit/%s", i, modelName))
    sprintf("Stan warning messages saved at stanWarnings/exp_%s.txt", modelName)
    expModelFit(modelName, isFirstFit = T, parallel = parallel)
    
    # check whether the model converge and fit again if necessary 
    print(sprintf("Increase samples to fit participants with disconvergent results in Exp.%d", i))
    sprintf("Stan warning messages saved at stanWarnings/exp_refit_%s.txt", modelName)
    expModelFit(modelName, isFirstFit = F, parallel = parallel)
  }
}

#####################################################################################
##                 Observed vs Model-generated, example participants               ##
#####################################################################################
setwd('./wtw_exp1')
source("exp1_expModelRepInd.R")
figs_ = expModelRepInd()
setwd(pwd)
figRepExample = figs_[['emp']] | figs_[['modelInd']] | figs_[['modelGroup']]
ggsave(file.path("..", "figures", "cmb", "modelRep_example.eps"), figRepExample, width = 14, height = 4)


###############################################################################
##                       qualitative model comparison                       ##
##############################################################################
MFResults_ = list()
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i]))
  source("MFAnalysis.R")
  MFResults = MFAnalysis(isTrct = T)
  MFResults_[[i]] = MFResults
}
# plot ovserved and model-generated AUC and sigma_WTW
outs_ = vector("list", length = nExp )
# sqerr_df_ = vector("list", length = nExp)
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i])) # set the working directory
  source(sprintf("exp%d_expModelRep.R", i)) 
  source("subFxs/loadFxs.R")

  allData = loadAllData() # load all the data 
  outs = vector("list", length = nModel) # initialize the output 
  for(j in 1 : nModel){
    model = models[j]
    output = expModelRep(model, allData, MFResults_[[i]])
    outs[[j]] = output
  }
  outs_[[i]] = outs
}


# plot WTW timecourses
setwd(pwd)
figWTW = (outs_[[1]][[2]][['figWTW']] |outs_[[2]][[2]][['figWTW']] | outs_[[3]][[2]][['figWTW']]) + plot_annotation(tag_levels = "a")
ggsave(file.path("..", "figures", "cmb", "rep_wtw_new.pdf"), figWTW, width = 15, height = 4 )

# combine QL2 results 
setwd(pwd)
figQL2 = (outs_[[1]][[2]]$figStat | outs[[2]]$figStat | outs[[3]]$figStat) + plot_annotation(tag_levels = "a")
ggsave(file.path("..", "figures", "cmb", "figQL2.pdf"), figQL2, width = 12, height = 12)

# look at variance 
cbind(outs_[[1]][[2]]$r2_df, outs_[[2]][[2]]$r2_df, outs_[[3]][[2]]$r2_df)

# combine figures together 
for(i in 1 : nExp){
  outs = outs_[[i]]
  figStats = (outs[[1]]$figStats | outs[[2]]$figStats | outs[[3]]$figStats | outs[[4]]$figStats | outs[[5]]$figStats | outs[[6]]$figStats)
  setwd(pwd)
  ggsave(file.path("..", "figures", "cmb", sprintf("modelRep_exp%d.eps", i)), figStats, width = 4 * 6, height = 3 * 4 )
}

figWTW_ = list()
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_expModelRepGroup.R", i))
  figWTW_[[i]]= expModelRepGroup()
}
setwd(pwd)
figWTW_all = (figWTW_[[1]] | figWTW_[[2]] | figWTW_[[3]])
ggsave(file.path("..", "figures", "cmb", "figWTW_all.pdf"), figWTW_all, width = 15, height = 4)

###############################################################################
##                       quantitative model comparison                       ##
##############################################################################
cmpOuts_ = vector("list", length = nExp)
waic_ave = matrix(NA, 3, 5) 
waic_se = matrix(NA, 3, 5)  
best_fit = matrix(NA, 3, 5) 
for(i in 1 : 3){
  setwd(file.path(pwd, wds[i])) # set the working directory
  source(sprintf("exp%d_expModelCmp.R", i)) 
  cmpOuts_[[i]] = expModelCmp()
  waic_ave[i,] = cmpOuts_[[i]][['full_waic_aves']]
  waic_se[i,] = cmpOuts_[[i]][['full_waic_ses']]
  best_fit[i,] = as.numeric(cmpOuts_[[i]][['bestFit']]$bestFitNum[1:5])
}
round(best_fit / apply(best_fit, MARGIN = 1, sum) * 100)
###################################################################################
##                     Parameter histograms and correlations                     ##
###################################################################################
outs_ = vector("list", length = nExp )
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_expParaAnalysis.R", i))
  outs = expParaAnalysis()
  outs_[[i]] = outs
}
# assemble figures
setwd(pwd)
histPara = (outs_[[1]][["hist"]] / outs_[[2]][['hist']] / outs_[[3]][['hist']]) + plot_layout(heights = c(0.5, 0.5, 0.5))
ggsave(file.path("../figures/cmb", "para_hist.eps"), histPara, width = 6, height = 6)



for(i in 1 : nExp){
  outs = outs_[[i]]
  print(paste(paste0(outs$para_summar$median, " [", outs$para_summar$q1, "-", outs$para_summar$q3, "]"), seq = "&", collapse = ''))
}

 
paste(outs_[[1]]$para_summar$q1, seq = "&", collapse = '')

outs_[[2]]$optim_summary


# figures for correlation analysis (with self-report measures and among parameters) are saved separately for each experiments
# print results for optimism bias 
outs_[[1]]$nuTest # whether nu < 1 in Exp.1
outs_[[1]]$numOptim # number of participants with optimism bias in Exp.1
outs_[[2]]$nuTest # whether nu < 1 in Exp.2
outs_[[2]]$numOptim # number of participants with optimism bias in Exp.2
outs_[[3]]$nuTest # whether nu < 1 in Exp.3
outs_[[3]]$numOptim # number of participants with optimism bias in Exp.3

#################################################################
## parameter recovery, simulation and parameter estimation     ##
#################################################################
parallel = F # When run on a local PC, you can use set this variabel to true to turn on parallel computing
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_simModelFit.R", i))
  source(sprintf("exp%d_simulation.R", i))
  # simulate artificial participants using indvidually fit parameters 
  simulation()
  # fit all participants
  print(sprintf("Model fitting for simulated participants in Exp.%d", i))
  print(sprintf("Model fitting results saved at ../genData/wtw_exp%d/simModelFit/%s", i, "QL2"))
  sprintf("Stan warning messages saved at stanWarnings/sim_%s_%s.txt", "QL2", "QL2")
  simModelFit("QL2", "QL2", isFirstFit = T, parallel = parallel)
  
  # check whether the model converge and fit again if necessary 
  print(sprintf("Increase samples to fit simulated participants with disconvergent results in Exp.%d", i))
  sprintf("Stan warning messages saved at stanWarnings/sim_refit_%s_%s.txt", "QL2", "QL2")
  simModelFit("QL2", "QL2", isFirstFit = F, parallel = parallel)

}

#########################################################
##            parameter recovery, plotting            ##
########################################################
outs_ = vector("list", length = nExp )
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i])) # set the working directory
  source(sprintf("exp%d_paraRecovery.R", i))
  outs = paraRecovery()
  outs_[[i]] = outs
}
# assemble the figures
setwd(pwd)
figParaRecovery = outs_[[1]]$fig / outs_[[2]]$fig / outs_[[3]]$fig
ggsave(file.path("../figures/cmb", "paraRecovery.eps"), figParaRecovery, width = 25, height = 12)
# print the number of included artificial participants for each experiments
outs_[[1]]$nPass
outs_[[2]]$nPass
outs_[[3]]$nPass

########### individual participants ##########
setwd("wtw_exp1")
source("exp1_expModelRepInd.R")
outs = expModelRepInd()
figInd = (outs$emp | outs$modelInd | outs$modelGroup) + plot_annotation(tag_levels = "a")
setwd(pwd)
ggsave(file.path("../figures/cmb", "figInd.eps"), figInd, width = 15, height = 4)

##### combine correlation figures ####3
library("figpatch")
outs = list()
for(i in 1 : nExp){
  this_path = system.file("../figures/cmb", sprintf("exp%d_corr.pdf", i))
  img = fig(this_path)
  outs[[i]] = img
}
pat <- patchwork::wrap_plots(outs[[1]], outs[[2]], outs[[3]])




