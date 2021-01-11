loadAllData = function() {
  # loads hdrData and trialData
  
  # outputs:
  # hdrData: exp info for each participant, like ID, condition, questionnaire measurements
  # trialData: a length = nSub list, each element containing trial-wise data for each partcipant. 
  
  # each elemant of trialData is formatted similarly to this example:
  # blockNum : [130x1 int]
  # trialNum : [130x1 int]
  # trialStartTime : [130x1 num] # when participant start a trial
  # sellTime : [130x1 num] # when participants sell a token
  # nKeyPresses : [130x1 int]
  # scheduledWait : [130x1 num] # delay durations for rewards
  # rewardTime : [130x1 num] # actual delay durations (constrainted by the screen update freq), NA for non-rewarded trials
  # timeWaited : [130x1 num] # persistence durations, namely sellTime - trialStartTime
  # trialEarnings : [130x1 int] trial-wise payments, either 10 or 0
  # totalEarnings : [130x1 int] cumulative payments
  
  # load experiment paras
  load('expParas.RData')
  
  # load hdrData
  hdrData = read.csv(file.path('data', 'hdrData.csv'), comment = "#")
  hdrData$stress = ifelse( hdrData$stress == 'stress', 'stress', 'no_stress')
  hdrData$condition = ifelse(hdrData$condition == 1, 'HP', 'LP')
  hdrData$id = as.character(hdrData$id)
  
  # load trialData
  trialData = list()
  nSub = nrow(hdrData)
  trialDataNames = c('blockNum', 'trialNum', 'trialStartTime', 'nKeyPresses', 'scheduledWait',
                      'rewardTime', 'timeWaited', 'sellTime', 'trialEarnings','totalEarnings')
  # loop over subjects
  for (sIdx in 1:nSub) {
    id = hdrData$id[sIdx]
    condition = hdrData$condition[hdrData$id == id]
    stress = hdrData$stress[hdrData$id == id]
    thisTrialData = vector(nBlock, mode = 'list')
    # loop over blocks
    for (bkIdx in 1 : nBlock){
      blockFile = list.files(path="data", pattern=(sprintf('wtw_stress_SDG%s_bk%d_1.txt',id, bkIdx)))
      if (length(blockFile) != 1) {
        cat('Could not identify a single data file for subject',thisID,' block', bkIdx, '\n')
        browser()
      }
      thisTrialData[[bkIdx]] = read.csv(file.path("data", blockFile), header = F, col.names = trialDataNames)
      thisTrialData[[bkIdx]]$blockNum = rep(bkIdx, nrow(thisTrialData[[bkIdx]]))
      thisTrialData[[bkIdx]]$condition = rep(condition, nrow(thisTrialData[[bkIdx]]))
    }
    trialData[[id]] = bind_rows(thisTrialData)
  } # end of loop over subjects
  # return outputs
  outputs = list(hdrData=hdrData, trialData=trialData)
  return(outputs)
} 


loadExpPara = function(paraNames, dirName){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 4 + 1)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    idIndexs = str_locate(fileName, "s[0-9]+")
    idList[i] = substr(fileName, idIndexs[1]+1, idIndexs[2])
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
    expPara[i, nE * 4 + 1] = junk[1,11]
  }

  # transfer expPara to data.frame
  expPara = data.frame(expPara, stringsAsFactors = F)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"), 'nDt')
  expPara$id = idList # ensure the levels are consistent, usually not that problematic though
  return(expPara)
}


