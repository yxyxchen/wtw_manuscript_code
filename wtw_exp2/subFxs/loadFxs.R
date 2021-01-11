loadAllData = function() {
  # loads hdrData and trialData
  
  # outputs:
  # hdrData: exp info for each participant
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
  
  library("gtools")

  # get the filename for each participant 
  fileNames = list.files(path= "data", pattern=('wtw-timing-fixed_[0-9]{3}_[0-9].txt'))
  fileNames = mixedsort(sort(fileNames))
  nSub = length(fileNames)
  if(any(duplicated(fileNames))){
    print("duplicated files!")
    break
  }else{
    sprintf("load data for %d participants", nSub) 
  }
  
  # initialize output variables
  hdrData = matrix(ncol = 1, nrow = nSub) # hdrData only contains ids 
  colnames(hdrData) = "id"
  trialData = list()
  trialDataNames = c('blockNum', 'trialNum', 'trialStartTime', 'nKeyPresses', 'scheduledWait',
               'rewardTime', 'timeWaited', 'sellTime', 'trialEarnings','totalEarnings')
  
  # loop over participants
  for(i in 1 : nSub){
    fileName = fileNames[i]
    id = substr(fileName, 18,20)
    hdrData[i,1] = id
    thisTrialData = read.csv(sprintf("%s/%s", "data", fileName), col.names = trialDataNames, header = F)
    thisTrialData = within(thisTrialData, {condition = ifelse(blockNum == 1, "LP", "HP")})
    trialData[[id]] = thisTrialData
  }
 
  # return outputs
  hdrData = as.data.frame(hdrData, stringsAsFactors = F)
  outputs= list(
    hdrData=hdrData,
    trialData=trialData)
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
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat")), "nDt")
  expPara$id = idList # ensure the levels are consistent, usually not that problematic though
  return(expPara)
}

# I also need to load 2.5% and 97.5%
loadCVPara = function(paraNames, dirName, pattern){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern= pattern)
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 6)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    sIndexs = str_locate(fileName, "s[0-9]+")
    s = substr(fileName, sIndexs[1]+1, sIndexs[2])
    fIndexs = str_locate(fileName, "f[0-9]+")
    f = substr(fileName, fIndexs[1]+1, fIndexs[2])
    idList[i] = sprintf("s%s_f%s", s, f)
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
    expPara[i, (4*nE + 1) : (5 * nE)] = junk[,4]
    expPara[i, (5 * nE + 1) : (6 * nE)] = junk[,8]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"),
                        paste0(junk, "2.5"),paste0(junk, "97.5"))
  expPara$id = idList
  return(expPara)
}

