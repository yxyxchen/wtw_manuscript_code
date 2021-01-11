# experiment design paramters
conditions = c("HP", "LP")
delayMaxs = c(20, 40) # max delay durations in secs
nBlock = 3
blockMin = 7 # block duration in mins
blockSec = blockMin * 60 # block duration in secs
tokenValue = 10 # value of the matured token

# 
pareto = list()
pareto[['k']] = 3
pareto[['mu']] = 0
pareto[['sigma']] = 1.5

# analyses parameters
tGrid = seq(0, blockSec * nBlock - 1, by = 1) # time grid for wtw time courses, open interval 
kmGrid = seq(0, min(delayMaxs) - 0.2, by = 0.2) # time grid for Kaplan-Meier survival curves

# save 
save("conditions" = conditions,
     "delayMaxs" = delayMaxs,
     "blockMin" = blockMin,
     "blockSec" = blockSec,
     "nBlock" = nBlock,
     "tokenValue" = tokenValue,
     "pareto" = pareto,
     'tGrid' = tGrid,
     'kmGrid' = kmGrid,
     file = "expParas.RData")

