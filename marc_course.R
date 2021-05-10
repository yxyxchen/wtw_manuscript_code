
# summary statistics 
MFResults$sumStats %>% filter(ID %in% c("461", "478", "531"))

# trialplots
load("expParas.RData")

# load exp data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$id
nSub = length(ids)                    # n
cat('Analyzing data for',nSub,'subjects.\n')

# mkdir 
dir.create("tmp")

# plot 
for(sIdx in 1 : 3){
  ID = ids[sIdx]
  thisTrialData = trialData[[ID]]
  trialPlots(block2session(thisTrialData)) +
    theme(legend.position = "None") +
    geom_vline(xintercept = max(thisTrialData$trialNum[thisTrialData$condition == "LP"]))
  ggsave(sprintf("tmp/sub%d.png", sIdx), width = 5, height = 4) 
  
  # 
  wtwTS(block2session(thisTrialData), tGrid = seq(0, blockSec*2-1, by = 1), min(delayMaxs), T) 
  ggsave(sprintf("tmp/subwtw%d.png", sIdx), width = 5, height = 4) 
  
  
}


