figs_ = vector("list", length = nExp )
sqerr_df_ = vector("list", length = nExp)
for(expIdx in 1 : nExp){
  setwd(file.path(pwd, wds[expIdx])) # set the working directory
  source(sprintf("exp%d_expModelRep.R", expIdx)) 
  source("subFxs/loadFxs.R")
  source("MFAnalysis.R")
  allData = loadAllData() # load all the data 
  MFResults = MFAnalysis(isTrct = T)
  figs = vector("list", length = nExp) # initialize the output 
  dfs = vector("list", length = nExp) # initialize the output 
  for(mIdx in 1 : nModel){
    modelName = models[mIdx]
    if(file.exists(sprintf("../../genData/wtw_exp%d/expModelRep/%s_trct.RData", expIdx, modelName))){
      load(sprintf("../../genData/wtw_exp%d/expModelRep/%s_trct.RData", expIdx, modelName))
      outs = expModelRep(modelName, allData,  MFResults, repOutputs)
    }else{
      outs = expModelRep(modelName, allData,  MFResults)
    }
    # record generated figures 
    figs[[mIdx]] = outs$rep
    dfs[[mIdx]] = outs$sqerr_df
    # merge sqerr tables across models
    # this_sqerr = outs$sqerr_df
    # if(expIdx == 1){
    #   names(this_sqerr) = c("id", paste("mu_sqerr", modelName, sep = "_"), paste("std_sqerr", modelName, sep = "_"))
    # }else{
    #   names(this_sqerr) = c("id", "condition", paste("mu_sqerr", modelName, sep = "_"), paste("std_sqerr", modelName, sep = "_"))
    # }
    # 
    # if(mIdx == 1){
    #   sqerr_df = this_sqerr
    # }else{
    #   if(expIdx == 1){
    #     sqerr_df = merge(sqerr_df,  this_sqerr, by = "id")
    #   }else{
    #     sqerr_df = merge(sqerr_df,  this_sqerr, by = c("id", "condition"))
    #   }
    # }
  }
  figs_[[expIdx]] = figs
  sqerr_df_[[expIdx]] = dfs
}

# 
plotdf = data.frame()
for(expIdx in 1 : nExp){
  included_ids = Reduce(intersect, lapply(sqerr_df_[[expIdx]], FUN = function(x) x$id))
  # plotdf = data.frame()
  for(mIdx in 1 : nModel){
    modelName = models[mIdx]
    this_df = sqerr_df_[[expIdx]][[mIdx]]
    this_df['model'] = modelName
    this_df['exp'] = exps[expIdx]
    if(expIdx == 1){
      this_df['condition'] = NaN
    }
    plotdf  = rbind(plotdf, this_df[this_df$id %in% included_ids,])
  }
}

plotdf$model = factor(plotdf$model, levels = models)
plotdf %>% ggplot(aes(model, sqrt(mu_sqerr))) + geom_boxplot() + facet_grid(exp~., scales = "free") + myTheme +
  xlab("") + ylab("Square root error in predicting AUC")

plotdf %>% ggplot(aes(model, sqrt(std_sqerr))) + geom_boxplot() + facet_grid(exp~., scales = "free") +
  xlab("") + ylab("Square root error in predicting sigma_WTW") + myTheme

