#-# the main script to run all analyses #-# 

# libraries
library("tidyverse")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(extrafont)
library(patchwork)
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
figs_ = vector("list", length = nExp )
for(i in 1 : nExp){
  figs = vector("list", length = nExp) # initialize the output 
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_MFPlot.R", i))
  figs = MFPlot()
  figs_[[i]] = figs
}
# assemble the figures
setwd(pwd)
figMF12 = (figs_[[1]][['curve']] | figs_[[1]][['wtw']] | figs_[[1]][['auc']] | figs_[[1]][['sigma']] | figs_[[1]][['delta']]) / \
(figs_[[2]][['curve']] | figs_[[2]][['wtw']] | figs_[[2]][['auc']] | figs_[[2]][['sigma']] | figs_[[2]][['delta']]) / + plot_annotation(tag_levels = "a")
ggsave(file.path("../figures/cmb","mf12.eps"), figMF12 , width = 16, height = 8)
figMF3 = (figs_[[3]][['curve']] | figs_[[3]][['auc']] | figs_[[3]][['sigma']] | figs_[[3]][['wtw']]) + plot_annotation(tag_levels = "a")
ggsave(file.path("../figures/cmb","mf3.eps"), figMF3 , width = 16, height = 8)

##################################################################
##                 Performance check simulation                 ##
##################################################################
setwd("./wtw_simulation")
load("expParas.RData")
source("performCheck.R")
figs_ = performCheck()
# assemble the figures
setwd(pwd)
figAUC = (figs_[[1]][['learn']] | figs_[[2]]['learn'])
ggsave(file.path("..", "figures", "cmb", "exante_learn_curve.eps"), figAUC, width = 6, height = 3)
figRV = (figs_[[1]][['rv']] | figs_[[2]]['rv'] | figs_[[3]]['rv']) + ylim(-2, 12) + ylab("Expected net return")  + plot_annotation(tag_levels = "a")
ggsave(file.path("..", "figures", "cmb", "exante_rv.eps"), figRV, width = 12, height = 3)
figSnippet = (plot_spacer() | figs_[[1]][["Gs_"]][1] | figs_[[1]][["Gs_"]][2] | figs_[[1]][["Gs_"]][3]) / 
  (figs_[[1]][["values_"]][[1]] | figs_[[1]][["values_"]][[2]] | figs_[[1]][["values_"]][[3]] | figs_[[1]][["values_"]][[4]])
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
figPostHoc = (figs_[[1]] + theme(legend.position = "None")) + xlab("Simulation time (min)")
ggsave(file.path("..", "figures", "cmb", "posthoc_para_effect.eps"), figPostHoc , width = 9, height = 4)


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
# plot ovserved and model-generated AUC and sigma_WTW
figs_ = vector("list", length = nExp )
# sqerr_df_ = vector("list", length = nExp)
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i])) # set the working directory
  source(sprintf("exp%d_expModelRep.R", i)) 
  source("subFxs/loadFxs.R")
  source("MFAnalysis.R")
  allData = loadAllData() # load all the data 
  MFResults = MFAnalysis(isTrct = T)
  figs = vector("list", length = nModel) # initialize the output 
  for(j in 1 : nModel){
    model = models[j]
    thisFig = expModelRep(model, allData, MFResults)
    figs[[j]] = thisFig
  }
  figs_[[i]] = figs
}

# assemble the figures 
i = 2
(figs_[[i]][[1]]$figWTW | figs_[[i]][[2]]$ | figs_[[i]][[3]] | figs_[[i]][[4]] | figs_[[i]][[5]] | figs_[[i]][[6]])

figs_[[i]][[j]]$figWTW
# assemble the figures 
setwd(pwd)
for(i in 1 : nExp){
  figRep = (figs_[[i]][[1]] | figs_[[i]][[2]] | figs_[[i]][[3]] | figs_[[i]][[4]] | figs_[[i]][[5]] | figs_[[i]][[6]])
  ggsave(file.path("../figures/cmb", sprintf("exp%s_modelRep.eps", i)), figRep, width = 24, height = 8)
}
for(i in 1 : nExp){
  figRep = (figs_[[i]][[2]] | figs_[[i]][[5]] | figs_[[i]][[6]])
  ggsave(file.path("../figures/cmb", sprintf("exp%s_modelRep_simple.eps", i)), figRep, width = 12, height = 8)
}

###############################################################################
##                       quantitative model comparison                       ##
##############################################################################
cmpOuts_ = vector("list", length = nExp)
waic_ave = matrix(NA, 3, 6) 
waic_se = matrix(NA, 3, 6)  
best_fit_4 = matrix(NA, 3, 4) 
best_fit_6 = matrix(NA, 3, 6) 
for(i in 1 : 3){
  setwd(file.path(pwd, wds[i])) # set the working directory
  source(sprintf("exp%d_expModelCmp.R", i)) 
  cmpOuts_[[i]] = expModelCmp()
  full_waic_aves = full_waic_ %>% apply(2, mean)
  full_waic_ses = full_waic_ %>% apply(2, function(x) sd(x) / sqrt(length(x)))
  waic_ave[i,] = cmpOuts_[[i]][['full_waic_aves']]
  waic_se[i,] = cmpOuts_[[i]][['full_waic_ses']]
  best_fit_4[i,] = as.numeric(cmpOuts_[[i]][['bestFit4']]$bestFitNum[1:4])
  best_fit_6[i,] = as.numeric(cmpOuts_[[i]][['bestFit6']]$bestFitNum[1:6])
}


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



