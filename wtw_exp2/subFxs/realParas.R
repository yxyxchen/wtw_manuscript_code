# monte
nPara = 4
paraNames = c('phi', 'tau', 'gamma', "QwaitIni")
nValue = 3
nComb = nValue ^ nPara
realParas = matrix(NA, nValue^nPara, nPara)
realParas[,1] = rep(seq(0.01, 0.05, length.out = 5), nValue^(nPara - 1)) # phi
realParas[,2] = rep(rep(seq(10, 22, length.out = 5), each = nValue), nValue^(nPara - 2)) # tau
realParas[,3] = rep(seq(0.8, 0.98, length.out = 5), each = nValue^2)
nRep = 5
simIdx = t(matrix(1 : (nComb * nRep), nRep, nComb))
save("realParas", "nComb", "nValue", "nPara", "paraNames", "simIdx", "nRep", file = "../../genData/wtw_exp2/simulation/monte/simParas.RData")

paraValues = data.frame(phi = c())
