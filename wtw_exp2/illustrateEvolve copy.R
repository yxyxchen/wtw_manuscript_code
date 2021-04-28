# Ensure that Lk is working ok
source("Lk.R")
# need this for analytic ftilde
source("normalization.R")
source("evolveF.R")

# if you've recently run this, you can save some time by setting this to FALSE
makemats <- FALSE
# set this to FALSE if you don't want to see the pictures.
makepix <- TRUE

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

if(makemats){
		# this should be the Post approximation with k (and edge effects removed)
		MyLk2 <- GetLk(svals=Ss, power=myk, mults=T)
		MyLk8 <- GetLk(svals=Ss, power=myk*4, mults=T)
}

#############################################
#			 MAKE SOME Fs to invert			# 
#############################################


# some vectors for items
A <- c(1,0,0) # start of the trial
B <- c(0,1,0) # the reward 
C <- c(0,0,1) 

# the times at which to present things
prestimes <- c(0,10,20)

# the items presented at each of these times
fvecs <- cbind(A, B, C)

# initialize an F(t).
Ft <- array(0, c(length(Ss),3))

# evolve it over the list
Ft <- EvolveF(PrevF = Ft, svals=Ss, dt=0, fin=fvecs[,1])
for(i in 2:length(prestimes)){
	deltat <- prestimes[i] - prestimes[i-1]
	Ft <- EvolveF(PrevF = Ft, svals=Ss, dt=deltat, fin=fvecs[,i])
}

# Let's evolve it a smidge so we don't get numerical errors on the C column:
Ft <- EvolveF(PrevF = Ft, svals=Ss, dt=10, fin = c(0,0,0))

# Now let's make some ftildes and see what happens.
ftilde2 <- MyLk2 %*% Ft
ftilde8 <- MyLk8 %*% Ft

if(makepix){
		layout(1:2)
		par(mar=c(4,4,1,1))
		par(tcl=0.2)

		# these are the three stimulus dimensions.  A is red, B is green, C is
		# blue
		plot(-log(Ss), Ft[,1], type="l",
			 lwd=1,  col=rgb(1,0,0,0.5),
			xlab=" log tau", 
		   ylab = "F(s)"	
			 )
		lines(-log(Ss),Ft[,2] , col=rgb(0,1,0,0.5))
		lines(-log(Ss),Ft[,3] , col=rgb(0,0,1,0.5))

		# same for ftilde8
		plot(-log(Ss), ftilde8[,3], type="l",
			 lwd=1,  col=rgb(0,0,1,0.5),
			xlab=" log tau", 
		   ylab = "ftilde, k=8"	
			 )
		lines(-log(Ss),ftilde8[,2] , col=rgb(0,1,0,0.5))
		lines(-log(Ss),ftilde8[,1] , col=rgb(1,0,0,0.5))
}
