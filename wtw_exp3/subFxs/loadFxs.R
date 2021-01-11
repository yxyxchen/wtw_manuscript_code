loadAllData = function() {
  
  library("gtools")
  # load hdrData fileNames
  fileNames = list.files(path= "data", pattern=('*_hdr.txt'))
  fileNames = mixedsort(sort(fileNames))
  nFile = length(fileNames)
  if(any(duplicated(fileNames))){
    print("duplicated hdrData files!")
    break
  }else{
    sprintf("load %d hdrData files", nFile) 
  }
  # load hdrData
  id = vector(length = nFile)
  cbal = vector(length = nFile)
  for(i in 1 : nFile){
    fileName = fileNames[i]
    thisId = substr(fileName, 12,14)
    junk = read.csv(sprintf("%s/%s", "data", fileName))
    id[i] = thisId
    cbal[i] = as.double(substr(junk[[1]][1], 7, 7))
  }
  hdrData = data.frame(id = id, cbal = cbal, stringsAsFactors = F)
  
  # define column names 
  colNames = c('blockNum', 'trialNum', 'trialStartTime', 'nKeyPresses', 'scheduledWait',
               'rewardTime', 'timeWaited', 'sellTime', 'trialEarnings','totalEarnings')
  
  # get trialData fileNames
  trialFileNames = list.files(path= "data", pattern=('wtw-work-7_[0-9]{3}_[0-9].txt'))
  trialFileNames = mixedsort(sort(trialFileNames))
  nFile = length(trialFileNames)
  if(any(duplicated(trialFileNames))){
    print("duplicated trialData files!")
    break
  }else{
    sprintf("load %d trialData files", nFile) 
  }
  
  # load trialData
  trialData = list()
  for(i in 1 : nFile){
    fileName = trialFileNames[i]
    thisId = substr(fileName, 12,14)
    junk = read.csv(sprintf("%s/%s", "data", fileName), col.names = colNames, header = T)
    if(hdrData$cbal[which(hdrData$id == thisId)] == 1){
      condition = ifelse(junk$blockNum %% 2 == 1, "HP", "LP")
    }else{
      condition = ifelse(junk$blockNum %% 2 == 1, "LP", "HP")
    }
    junk$condition = condition
    trialData[[thisId]] = junk
  }
  
  outputData = list(hdrData=hdrData, trialData=trialData)
  return(outputData)
  
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
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"),
                        "nDt")
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




