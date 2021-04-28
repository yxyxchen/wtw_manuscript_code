# This file gives GetLk
source("Lk.R")

# This file has some simple functions for normalization
source("normalization.R")

# this file has hw 4 as a function and returns M after a list of three things
# decided to just make my own F's to test this out.
#source("simpleexp.R")


# if you don't want to make the figures set this to FALSE
makefigs=T

# k is the order of the Post approximation
myk <- 8 
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
for(i in 2:nums){
	Ss[i] <- (1+c) * Ss[i-1]
}
# this makes the values of s logarithmically spaced.

# for plotting things with the \Lk
taustars <- myk/Ss

# now we've got our values, lets check on them
print(paste("maximum value of s is", max(Ss)))
print(paste("maximum value of taustar is", max(taustars)))
print(paste("minimum value of s is", min(Ss)))
print(paste("minimum value of taustar is", min(taustars)))



# for if'n we want to do the 2nd deriv + divisive normalization
logtaus <- -log(Ss)
# make a deltas for the 2nd deriv + divisive normalization
deltalogs <- array(Ss[2] - Ss[1], length(Ss) + 1)

# mask to get rid of edge effects.
Lkmask <- array(1, length(Ss))
Lkmask[1:(2*myk-1)] <- 0
Lkmask[length(Ss):(length(Ss) - (2*myk-1))] <- 0


# this should be the Post approximation with k (and edge effects removed)
MyLk <- diag(Lkmask) %*% GetLk(svals=Ss, power=myk, mults=T)

# i want to use this to get the second derivative with 
SecondDerivLogs <- diag(Lkmask) %*% 
	GetLk(svals = logtaus, power=2, deltas=deltalogs, mults=F)



#############################################
#			 MAKE SOME Fs to invert			# 
#############################################

# Notice that the Ss are the same, but \tau varies
Ft100 <- exp(-Ss * 100) 
Ft10 <- exp(-Ss * 10) 
Ft1 <- exp(-Ss * 1) 


#############################################
#		 Done making Fs to invert			# 
#############################################

# Make some ftildes
ftilde100 <- MyLk %*% Ft100
ftilde10 <- MyLk %*% Ft10
ftilde1 <- MyLk %*% Ft1

# How does changing k affect these functions?  How would individual \taustars
# behave through time?  That is, if you pick one cell with one taustar, how
# does it behave as \tau goes from zero to infinity?

# alternate forms of approximating the inverse transform
# this is just the second derivative with respect to -log s
ftildelog100 <- clip(SecondDerivLogs %*% Ft100)
ftildelog10 <- clip(SecondDerivLogs %*% Ft10)
ftildelog1 <- clip(SecondDerivLogs %*% Ft1)

if(makefigs){
		# make some pictures of them
		layout(cbind(1:2,3:4))
		par(mar=c(4,4,1,1))


		plot(taustars,ftilde100, type="l", xlab="taustar", xlim=c(0,600))
		lines(taustars,ftilde10)
		lines(taustars,ftilde1)


		plot(log(taustars),ftilde100, type="l", xlab="log taustar")
		lines(log(taustars),ftilde10)
		lines(log(taustars),ftilde1)

		plot(1/Ss,ftildelog100, type="l", xlab="taustar", xlim=c(0,200))
		lines(1/Ss,ftildelog10)
		lines(1/Ss,ftildelog1)


		plot(logtaus,ftildelog100, type="l", xlab="log taustar")
		lines(logtaus,ftildelog10)
		lines(logtaus,ftildelog1)
}


# this is the second derivative run through divisive normalization
# How does this depend on the exponent?  on sigma?
divnormftilde <- divisivenorm(ftildelog10, sigma=0, exponent=2)

# You can also run through a softmax.  How does the width of the output depend on beta?
softmaxfitlde <- softmax(ftildelog10, betascale=1)
