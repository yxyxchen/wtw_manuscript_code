paraRecovery = function(){
  # libraries and scripts
  library("ggplot2")
  library("dplyr")
  library("tidyr")
  source("subFxs/helpFxs.R")
  source("subFxs/loadFxs.R")
  source("subFxs/plotThemes.R")
  library("latex2exp")
  library("patchwork")
  
  # load model names
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id            
  nSub = length(ids) 
  
  # check fit
  modelName = "QL2"
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  expPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp2/expModelFit/%s", modelName))
  expPassCheck = checkFit(paraNames, expPara)
  expPara = expPara[expPassCheck, ]
  simPara = loadExpPara(paraNames, sprintf("../../genData/wtw_exp2/simModelFit/%s/%s", modelName, modelName)) 
  passCheck = checkFit(paraNames, simPara)
  nPass = sum(passCheck)
  
  
  # plot 
  plotData = expPara[,1:nPara] %>% gather(key = "para", value = "truth")
  tempt = hdrData$condition
  plotData$condition = rep(tempt[expPassCheck], nPara)
  plotData$passCheck = rep(passCheck, nPara)
  addData = simPara[,1:nPara] %>% gather(key = "para", value = "recovered")
  plotData$recovered = addData$recovered
  paraLabels = c("$\\alpha$", "$\\rho$", "$\\tau$", "$\\gamma$", "$\\eta$")
  paraHats = c("$\\hat{\\alpha}$", "$\\hat{\\rho}$", "$\\hat{\\tau}$", "$\\hat{\\gamma}$", "$\\hat{\\eta}$")
  figs = vector(mode = "list", length = nPara) # initialize the outputs 
  rankTaus = vector(length = nPara)
  ps = rankTaus = vector(length = nPara)
  for(i in 1 : nPara){
    paraName = paraNames[i]
    paraLabel = paraLabels[i]
    paraHat = paraHats[i]
    plotData = data.frame(
      passCheck = passCheck,
      truth = expPara[[paraName]],
      recovered = simPara[[paraName]]
    ) %>% filter(passCheck) 
    corRes = cor.test(plotData$truth, plotData$recovered, method = "kendall")
    rankTaus[[i]] = corRes$estimate
    ps[[i]] = corRes$p.value
    lowLim = min(min(plotData$recovered), min(plotData$truth))
    upLim = max(max(plotData$recovered), max(plotData$truth)) 
    ran = upLim - lowLim
    lowLim = lowLim - ran * 0.05
    upLim = upLim + ran * 0.05
    figs[[i]] = plotData %>% ggplot(aes(truth, recovered)) + geom_point(color = "#636363") +
      myTheme + xlab(TeX(paraLabel)) +
      ylab(TeX(paraHat)) + geom_abline(xintercept = 0, slope = 1, linetype = "dashed") + 
      theme(legend.position = "none") +
      annotate("text", label = paste("list(r ==", round(rankTaus[[i]], 2), ", p < 0.001 ", ")"),
               x= upLim - 0.4 * ran, y =max(plotData$recovered), parse = T, size = 5) +
      coord_fixed(ratio = 1) +
      xlim(c(lowLim, upLim)) +
      ylim(c(lowLim, upLim)) 
    # ggsave(file.path("..", "figures", "exp2", "paraRecovery", sprintf("%s.eps", paraName)), figs[[i]], width = 4, height = 4)
  }
  
  fig = (figs[[1]] | figs[[2]] | figs[[3]] | figs[[4]] | figs[[5]])
  
  ########## return outputs #######
  outputs = list(
    "p" = ps,
    "rankTau" = rankTaus,
    "fig" = fig,
    "nPass" = nPass
  )
}


   
