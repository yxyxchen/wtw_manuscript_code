# this is the analysis script for simualtion results

# change the working directory
pwd = getwd()
setwd("./wtw_simulation")

# libraries 
source("simExAnte.R")
source("expSchematics.R")
source(file.path("subFxs", "simPostHoc.R"))
source(file.path("subFxs", "plotThemes.R"))
source(file.path("subFxs", "helpFxs.R"))
source(file.path("subFxs", "taskFxs.R"))
library(ggpubr)
library("latex2exp")
library(patchwork)
# load expParas
load("expParas.RData")

#################################################################
##                      ExAnte Simulation                      ##
#################################################################
modelName_ = c("QL1", "RL1", "omni", "naive", "naive")
modelLabel_ = c("QL1", "RL1", "omni", "naive75", "naive1")

# it is good to use different seqs acoss simulations yet keep them the same across models
nSim = 10
set.seed(123)
delays_ = vector(mode = "list", length = nSim)
for(i in 1 : nSim){
  delays_[[i]] = list(
    HP = replicate(120 * 60  / 2, drawSample("HP")),
    LP = replicate(120 * 60  / 2, drawSample("LP"))
  )
}
# simulation parameters 
paras_ = list(
  list("HP" = c(0.05, 5, 0.85, 6), "LP" = c(0.05, 5, 0.85, 6)), # "QL2"
  list("HP" = c(0.05, 5, 10, 0.04), "LP" = c(0.05, 5, 10, 0.04)), # "RL2"
  list("HP" = c(5), "LP" = c(5)), # omni
  list("HP" = c(0.75), "LP" = c(0.75)), # naive p = 0.75
  list("HP" = c(1), "LP" = c(1)) # naive p = 0.75
)

modelName = modelName_[1]
modelLabel = modelLabel_[1]
paras = paras_[[1]]

# simulate and plot 
figs_ = vector(mode = "list", length = length(modelName_))
for(i in 1 : length(modelName_)){
  outs = simExAnte(modelName_[i], modelLabel_[i], paras_[[i]], delays_)
  figs_[[i]] = outs
}
setwd(pwd)
figAUC = (figs_[[1]][['learn']] | figs_[[2]]['learn'] | figs_[[3]]['learn'] | figs_[[4]]['learn'] | figs_[[5]]['learn'])
ggsave(file.path("..", "..", "figures", "cmb", "exante_learn_curve.eps"), figAUC, width = 15, height = 3)
figRV = (figs_[[1]][['rv']] | figs_[[2]]['rv'] | figs_[[3]]['rv'])
ggsave(file.path("..", "..", "figures", "cmb", "exante_rv.eps"), figRV, width = 15, height = 3)
figSnippet = (plot_spacer() | figs_[[1]][["Gs_"]][1] | figs_[[1]][["Gs_"]][2] | figs_[[1]][["Gs_"]][3]) / 
                (figs_[[1]][["values_"]][[1]] | figs_[[1]][["values_"]][[2]] | figs_[[1]][["values_"]][[3]] | figs_[[1]][["values_"]][[4]])
ggsave(file.path("..", "figures", "cmb", "exante_snippet.eps"), figSnippet , width = 12, height = 6)

#################################################################
##                      Parameter effects                      ##
#################################################################
# simulate with different parameter combinations
modelName_ = c("QL2", "RL2")
paraLabels_ = list(
  c("alpha", "rho", "tau", "gamma", "eta"),
  c("alpha", "rho", "tau", "eta", "beta")
  )
nCut = 4
paraSamples_ = list(
  list("HP" = data.frame(
    alpha = seq(0.01, 0.05, length.out = nCut),
    rho = seq(0.03, 3, length.out = nCut),
    tau = exp(seq(log(0.5), log(8), length.out = nCut)),
    gamma = seq(0.85, 0.95, length.out = nCut),
    prior = seq(0, 1, length.out = nCut)
    ), "LP" = data.frame(
    alpha = seq(0.01, 0.05, length.out = nCut),
    rho = seq(0.03, 3, length.out = nCut),
    tau = exp(seq(log(0.5), log(8), length.out = nCut)),
    gamma = seq(0.85, 0.95, length.out = nCut),
    prior = seq(2, 6, length.out = nCut)
    )),
  list("HP" = data.frame(
    alpha = seq(0.01, 0.05, length.out = nCut),
    rho = seq(0.03, 3, length.out = nCut),
    tau = exp(seq(log(0.5), log(8), length.out = nCut)),
    prior = seq(0, 1, length.out = nCut),
    beta = seq(0.0105, 0.025, length.out = nCut)
    ), "LP" = data.frame(
    alpha = seq(0.01, 0.05, length.out = nCut),
    rho = seq(0.03, 3, length.out = nCut),
    tau = exp(seq(log(0.5), log(8), length.out = nCut)),
    prior = seq(2, 6, length.out = nCut),
    beta = seq(0.0105, 0.025, length.out = nCut)
    ))
  )

for(i in 1 : length(modelName_)){
  simPostHoc(modelName_[i], paraLabels_[[i]], paraSamples_[[i]])
}

# ###################### robustness analysis ##################
# # simulate with different itis and non-zero quit values
# itis = c(2, 4, 8)
# smallReward.mag = 1
# smallReward.probs= c(0, 0.4, 0.8)
# 
# # normative analysis for different itis and non-zero quit values
# rewardRate_ = list(); condition_ = list(); iti_ = list(); prob_ = list(); t_= list()
# count = 1
# for(iti in itis){
#   for(prob in smallReward.probs){
#     normResults = expSchematics(prob * smallReward.mag, iti, F)
#     rewardRate_[[count]] = c(normResults$rewardRates$HP, normResults$rewardRates$LP)
#     condition_[[count]] = rep(conditions, sapply(normResults$rewardRates, length))
#     iti_[[count]] = rep(iti, length(condition_[[count]]))
#     prob_[[count]] = rep(prob, length(condition_[[count]]))
#     t_[[count]] = unlist(normResults$time)
#     count = count + 1
#   }
# }
# iti.labs = c("ITI = 2 s", "ITI = 4 s", "ITI = 8 s")
# names(iti.labs) = c("2", "4", "8")
# data.frame(
#   rewardRate = unlist(rewardRate_),
#   condition = unlist(condition_),
#   iti = unlist(iti_),
#   prob = unlist(prob_),
#   t = unlist(t_)
# ) %>%
#   ggplot(aes(t, rewardRate, color = as.factor(prob))) + geom_line() + 
#   facet_grid(condition ~ iti, labeller = labeller(iti = iti.labs)) +
#   myTheme + ylab(expression(bold("Reward rate (¢ s"^"-1"*")"))) + xlab("Waiting policy (s)") +
#   scale_colour_manual(values = c("#bdbdbd", "#737373", "black")) +
#   labs(color = TeX(c("p_{1¢}"))) + 
#   scale_x_continuous(breaks = c(0, 20, 40))
# ggsave("figures/simulation/extended/rewardRate.png", width = 7, height = 3)
# 
# # simulation 
# p_ = list(length = 9)
# count = 1
# for(iti in itis){
#   for(prob in smallReward.probs){
#     p_[[count]] = simExtended(smallReward.mag, prob, iti)
#     count = count + 1
#   }
# }
# d = ggarrange(p_[[1]] + rremove("legend") + rremove("ylab") + rremove("xlab"), 
#               p_[[2]] + rremove("legend") + rremove("ylab") + rremove("xlab"),
#               p_[[3]] + rremove("legend") + rremove("ylab") + rremove("xlab"),
#               p_[[4]] + rremove("legend")  + rremove("ylab") + rremove("xlab"),
#               p_[[5]] + rremove("legend") + rremove("ylab") + rremove("xlab"),
#               p_[[6]] + rremove("legend") + rremove("ylab") + rremove("xlab"),
#               p_[[7]] + rremove("legend") + rremove("ylab") + rremove("xlab"),
#               p_[[8]] + rremove("legend") +  rremove("ylab") + rremove("xlab"),
#               p_[[9]] + rremove("legend") + rremove("ylab") + rremove("xlab"),
#               ncol= 3, nrow = 3)
# ggsave("figures/simulation/extended/paraEffect.png", d, width = 10, height = 10)
# 
# ########################### opportunity cost irrelevant version ###################
# simOneOff()





