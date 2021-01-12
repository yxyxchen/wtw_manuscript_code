#-# the main script to run all analyses #-# 

# libraries
library("tidyverse")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(extrafont)
library(patchwork)
extrafont::font_import()
source("subFxs/plotThemes.R")


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
expCmb = (figs_[[1]][['pdf']] | figs_[[2]][['pdf']] | figs_[[3]][['pdf']]) /
  (figs_[[1]][['rt']] | figs_[[2]][['rt']] | figs_[[3]][['rt']]) +
  plot_annotation(tag_levels = 'a')
ggsave(file.path("../figures/cmb", "exp.eps"), expCmb, width = 12, height = 6)


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
figMFCurve = (figs_[[1]][['curve']] / figs_[[2]][['curve']]  / figs_[[3]][['curve']]) + plot_annotation(tag_levels = 'a')
ggsave(file.path("../figures/cmb","mf_curve.pdf"), figMFCurve, width = 4, height = 12)
figMFAUC =  (figs_[[1]][['auc']] / figs_[[2]][['auc']]  / figs_[[3]][['auc']]) + plot_annotation(tag_levels = 'a')
ggsave(file.path("../figures/cmb","mf_auc.eps"), figMFAUC, width = 4, height = 12)
figMFsigma =  (figs_[[1]][['sigma']] / figs_[[2]][['sigma']]  / figs_[[3]][['sigma']]) + plot_annotation(tag_levels = 'a')
ggsave(file.path("../figures/cmb","mf_sigma.eps"), figMFsigma, width = 4, height = 12)
figMFWTW =  (figs_[[1]][['wtw']] / figs_[[2]][['wtw']]  / figs_[[3]][['wtw']]) + plot_annotation(tag_levels = 'a')
ggsave(file.path("../figures/cmb","mf_wtw.pdf"), figMFWTW, width = 6, height = 12)


##################################################################
##                 Performance check simulation                 ##
##################################################################
setwd("./wtw_simulation")
source("performCheck.R")
figs_ = performCheck()
# assemble the figures
setwd(pwd)
figAUC = (figs_[[1]][['learn']] | figs_[[2]]['learn'] | figs_[[3]]['learn'] | figs_[[4]]['learn'] | figs_[[5]]['learn'])
ggsave(file.path("..", "figures", "cmb", "exante_learn_curve.eps"), figAUC, width = 15, height = 3)
figRV = (figs_[[1]][['rv']] | figs_[[2]]['rv'] | figs_[[3]]['rv']) + plot_annotation(tag_levels = "a")
ggsave(file.path("..", "figures", "cmb", "exante_rv.eps"), figRV, width = 15, height = 3)
figSnippet = (plot_spacer() | figs_[[1]][["Gs_"]][1] | figs_[[1]][["Gs_"]][2] | figs_[[1]][["Gs_"]][3]) / 
  (figs_[[1]][["values_"]][[1]] | figs_[[1]][["values_"]][[2]] | figs_[[1]][["values_"]][[3]] | figs_[[1]][["values_"]][[4]])
ggsave(file.path("..", "figures", "cmb", "exante_snippet.eps"), figSnippet , width = 12, height = 6)


#################################################################
##                 Parameter-effect simulation                 ##
#################################################################
setwd("./wtw_simulation")
source("paraEffect.R")
figs_ = paraEffect()
# assemble figures
setwd(pwd)
figPostHoc = (figs_[[1]] + theme(legend.position = "None") |
    figs_[[2]] + theme(legend.position = "None") + theme(legend.position = "None"))
ggsave(file.path("..", "figures", "cmb", "posthoc_para_effect.eps"), figPostHoc , width = 18, height = 4)


##################################################################
##                       model comparison                       ##
##################################################################
# plot ovserved and model-generated AUC and sigma_WTW
figs_ = vector("list", length = nExp )
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i])) # set the working directory
  source(sprintf("exp%d_expModelRep.R", i)) 
  source("subFxs/loadFxs.R")
  source("MFAnalysis.R")
  allData = loadAllData() # load all the data 
  MFResults = MFAnalysis(isTrct = T)
  figs = vector("list", length = nExp) # initialize the output 
  
  for(j in 1 : nModel){ 
    modelName = models[j]
    load(sprintf("../../genData/wtw_exp%d/expModelRep/%s_trct.RData", i, modelName))
    outs = expModelRep(modelName, allData,  MFResults, repOutputs)
    # outs = expModelRep(modelName, allData,  MFResults)
    figs[[j]] = outs[['rep']]
  }
  figs_[[i]] = figs
}
# figures for example participants are saved separately for each experiment
# pie charts and WAIC plots
cmpOuts_ = vector("list", length = nExp)
cmpFigs = vector("list", length = nExp)
comFigFulls = vector("list", length = nExp)
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i])) # set the working directory
  source(sprintf("exp%d_expModelCmp.R", i)) 
  cmpOuts_[[i]] = expModelCmp()
  cmpFigs[[i]] = cmpOuts_[[i]][['4pie']] / cmpOuts_[[i]][['4waic']]
  comFigFulls[[i]] = cmpOuts_[[i]][['6pie']] / cmpOuts_[[i]][['6waic']]
}
# assemble the figures 
setwd(pwd)
for(i in 1 : nExp){
  figLeft = (figs_[[i]][[1]] | figs_[[i]][[2]] | figs_[[i]][[3]] | figs_[[i]][[4]] | figs_[[i]][[5]] | figs_[[i]][[6]])
  figRep = (figLeft + cmpFigs[[i]])
  ggsave(file.path("../figures/cmb", sprintf("exp%s_modelRep.eps", i)), figRep, width = 28, height = 8)
}
figCmp6 = (comFigFulls[[1]] | comFigFulls[[2]] | comFigFulls[[3]]) + plot_annotation(tag_levels = 'a')
ggsave(file.path("../figures/cmb", sprintf("exp%s_modelCmp6.eps", i)), figCmp6, width = 12, height = 8)
# print summary stats 
for(i in 1 : nExp){
  print(sprintf("Exp.%d", i))
  print(cmpOuts_[[i]]$df4)
  print(cmpOuts_[[i]]$df6)
}

#################################################################
##                     Parameter estimates                     ##
#################################################################
outs_ = vector("list", length = nExp )
for(i in 1 : nExp){
  setwd(file.path(pwd, wds[i]))
  source(sprintf("exp%d_expParaAnalysis.R", i))
  outs = expParaAnalysis()
  outs_[[i]] = outs
}
# assemble figures
setwd(pwd)
histPara = (outs_[[1]][["hist"]] / outs_[[2]][['hist']] / outs_[[3]][['hist']]) + plot_layout(heights = c(1, 0.5, 0.5))
ggsave(file.path("../figures/cmb", "para_hist.eps"), histPara, width = 12, height = 8)
# figures for correlation analysis are saved separately for each experiments
# print summary stats
outs_[[1]]$rhoTest
outs_[[1]]$numOptim
outs_[[2]]$rhoTest
outs_[[2]]$numOptim
outs_[[3]]$rhoTest
outs_[[3]]$numOptim

##################################################################
##                      parameter recovery                      ##
##################################################################
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
# print the summary states
outs_[[1]]$nPass
outs_[[2]]$nPass
outs_[[3]]$nPass



