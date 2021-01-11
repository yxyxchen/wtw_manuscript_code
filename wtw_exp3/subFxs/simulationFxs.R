# define functions
simulate = function(modelName, nRep, paraTable, scheduledWaitList, cond){
  # create the output file
  dir.create("../../genData/wtw_exp3/simulation")
  dir.create(sprintf("../../genData/wtw_exp3/simulation/%s", modelName))
  
  # choose modelFun
  repModelFun = getRepModelFun(modelName)
  
  # determine paraComb and nSeq
  paraComb = getParaComb(paraTable)
  nComb = length(paraComb) / length(paraTable)
  nSeq = length(scheduledWaitList)
  simNo = array(t(seq(1 : (nComb * nSeq * nRep))), dim = c(nRep, nSeq, nComb)) 
  
  # initialize outputs
  trialData = vector(length = nComb * nSeq * nRep, mode ='list')
  
  # loop
  for(cbIdx in 1 : nComb){
    para = as.double(paraComb[cbIdx,])
    for(sqIdx in 1: nSeq){
      seq = scheduledWaitList[[sqIdx]]
      # using dim = c(nRep, nSeq, nComb), we get continous idex for the repetitions with the same seq and the same comb
      trialData[simNo[, sqIdx,cbIdx]] = lapply(1 : nRep, function(rpIdx) repModelFun(para, cond, seq))
    }
  }
  return(trialData)
}
