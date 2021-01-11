conditions = c("HP", "LP")
delayMaxs = c(16, 32) # max trial durations in secs
nBlock = 2
blockMin = 10 # block duration in mins
blockSec = blockMin * 60 # block duration in secs
tokenValue = 10 # value of the matured token
# possible values of reward dlays
rewardDelays = list(HP = seq(2, 16, by = 2),
                    LP = c(0.2759766, 0.7589358, 1.6041142, 3.0831765,
                             5.6715355, 10.2011638, 18.1280133, 32)) 
# analyses parameters
tGrid = seq(0, blockSec-1, by = 1) # time grid for wtw time courses
kmGrid = seq(0, min(delayMaxs), by = 0.1) # time grid for Kaplan-Meier survival curves
save("conditions" = conditions,
     "delayMaxs" = delayMaxs,
     "blockMin" = blockMin,
     "blockSec" = blockSec,
     "nBlock" = nBlock,
     "tokenValue" = tokenValue,
     "rewardDelays" = rewardDelays,
     'tGrid' = tGrid,
     'kmGrid' = kmGrid,
     file = "expParas.RData")