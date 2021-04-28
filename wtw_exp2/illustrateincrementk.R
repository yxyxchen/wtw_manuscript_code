source("normalization.R")

# set up some s values
# how many?
nums <- 1000 
# place to put them
Ss <- array(0,nums)
# smallest s value
mins <- .001
Ss[1] <- mins
# This is c as in 1+c in a bunch of papers
c <- .01
# this makes the values of s logarithmically spaced.
for(i in 2:nums){
	Ss[i] <- (1+c) * Ss[i-1]
}

# this is difference between each index on a log scale.  It's the same for
# every index 
deltalogs <- log(1+c)

myk <- 2
myn <- 6

# make some analytic ftildes.  Suppose you have a delta function tau units
# in the past  Here tau is 10 (i happen to know that fits into these svalues
# but that's the only reason.
analyticftilde2 <-  analyticftilde(tau=10,svals=Ss,k=myk)
analyticftilde12 <-  analyticftilde(tau=10,svals=Ss,k=(myk*myn))

# This is analytic  ftilde raised to the kfactor power and scaled as a
# function of s^{kfactor - 1} 
exponentk12 <- multiplyk(analyticftilde2,svals=Ss,n=myn,k=myk) 

# center on peaks
peakanalyticindex <- which(analyticftilde12==max(analyticftilde12))
peakexpindex <- which(exponentk12==max(exponentk12))


# I should also be able to shift indices by this amount and find the same
# value,  note that there's noise if log(1+c) is not an integer.
shiftsindex <- log(myn)/deltalogs


print(paste("peak of analyticftilde12 at", peakanalyticindex))
print(paste("peak of exponentk12 at", peakexpindex))
print(paste("difference should be log(n)/log(1+c) = ", shiftsindex))
print(paste("difference is", peakanalyticindex - peakexpindex))
# if the notes are correct (and there's no bug in this file), exponentk12 should be
# equal to analyticftilde12, but at different indices.  

layout(1:3)
par(mar=c(4,4,1,1))
par(tcl=0.2)

# these should be translated versions of one another
plot(-log(Ss) +log(Ss[peakexpindex]), exponentk12/max(exponentk12), type="l",
	 lwd=3, 
	xlab=" delta log tau from peak", 
   ylab = "Analytic and exponentiated ftilde"	
	 )
lines(-log(Ss) +log(Ss[peakanalyticindex]), analyticftilde12/max(analyticftilde12), col="blue")
#lines(-log(Ss), analyticftilde2/max(analyticftilde2), col="red")



# how does exponentk12 scale with tau?
# make some ftildes with different taus
analyticftildetau1 <-  analyticftilde(tau=2,svals=Ss,k=myk)
analyticftildetau10 <-  analyticftilde(tau=4,svals=Ss,k=myk)
analyticftildetau100 <-  analyticftilde(tau=8,svals=Ss,k=myk)

analyticftildetau1n <-  analyticftilde(tau=2,svals=Ss,k=myk*myn)
analyticftildetau10n <-  analyticftilde(tau=4,svals=Ss,k=myk*myn)
analyticftildetau100n <-  analyticftilde(tau=8,svals=Ss,k=myk*myn)

# exponentiate these
exponenttau1 <- multiplyk(analyticftildetau1,svals=Ss,n=myn,k=myk) 
exponenttau10 <- multiplyk(analyticftildetau10,svals=Ss,n=myn,k=myk) 
exponenttau100 <- multiplyk(analyticftildetau100,svals=Ss,n=myn,k=myk) 

# here's what we want (albeit translated)
plot(-log(Ss), analyticftildetau1n, type="l", xlab="log tau", 
	 ylab=paste("analytic ftilde with k=",myk*myn))
lines(-log(Ss), analyticftildetau10n)
lines(-log(Ss), analyticftildetau100n)



# here's what we got
plot(-log(Ss), exponenttau1, type="l", xlab="log tau", 
	 ylab=paste("ftilde with k=",myk,", n=",myn))
lines(-log(Ss), exponenttau10)
lines(-log(Ss), exponenttau100)

print(paste("max exponenttau1 is", max(exponenttau1)))
print(paste("max exponenttau10 is", max(exponenttau10)))
print(paste("max exponenttau100 is", max(exponenttau100)))

print(paste("max analyticftildetau1n is", max(analyticftildetau1)))
print(paste("max analyticftildetau10n is", max(analyticftildetau10)))
print(paste("max analyticftildetau100n is", max(analyticftildetau100)))


# if the calculation is correct, these indices should line up like so:
#plot(
#	 analyticftilde8[floor(shiftsindex) :		nums],
#	 exponentk8		[		1			:(nums - floor(shiftsindex))],
#	 xlab="Analytic ftilde_8",
#	 ylab="ftilde_2^6 with indices shifted"
#)
# and we should get a straight line in that plot


