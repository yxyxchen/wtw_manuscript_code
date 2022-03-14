conditions = c("HP", "LP")
nCondition = length(conditions)
delayMaxs = c(30, 30) # max trial durations in secs
blockMin = 10 # block duration in mins
blockSec = blockMin * 60 # block duration in secs
tokenValue = c(-1, 8)
# analyses parameters
tGrid = seq(0, blockSec - 1, by = 2) # time grid for wtw time courses
kmGrid = seq(0, min(delayMaxs), by = 0.1) # time grid for Kaplan-Meier survival curves
save("conditions" = conditions,
     "nCondition" = nCondition,
     "delayMaxs" = delayMaxs,
     "blockMin" = blockMin,
     "blockSec" = blockSec,
     "tokenValue" = tokenValue,
     'tGrid' = tGrid,
     'kmGrid' = kmGrid,
     file = "expParas.RData")