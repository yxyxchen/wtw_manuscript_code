# Ensure that Lk is working ok
source("Lk.R")
# need this for analytic ftilde
source("normalization.R")

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


# make some analytic ftildes.  Suppose you have a delta function tau units
# in the past 
analyticftilde2 <-  analyticftilde(tau=10,svals=Ss,k=myk)
analyticftilde8 <-  analyticftilde(tau=10,svals=Ss,k=(myk*4))


if(makemats){
		# this should be the Post approximation with k (and edge effects removed)
		MyLk2 <- GetLk(svals=Ss, power=myk, mults=T)
		MyLk8 <- GetLk(svals=Ss, power=myk*4, mults=T)
}

#############################################
#			 MAKE SOME Fs to invert			# 
#############################################

# Notice that the Ss are the same, but \tau varies
#Ft100 <- exp(-Ss * 100) 
Ft10 <- exp(-Ss * 10) 
#Ft1 <- exp(-Ss * 1) 



# Make some ftildes.  One for each of the k's.  These should resemble analytic
# ftildes
ftilde102 <- MyLk2 %*% Ft10
ftilde108 <- MyLk8 %*% Ft10

if(makepix){
		par(mar=c(4,4,1,1))
		par(tcl=0.2)

		layout(rbind(1:2,3:4))

		# these should be translated versions of one another
		plot(-log(Ss), analyticftilde2, type="l",
			 lwd=3, 
			xlab=" log tau", 
		   ylab = "ftilde, k=2"	
			 )
		lines(-log(Ss), ftilde102, col="blue")

		# these should be translated versions of one another
		plot(-log(Ss), analyticftilde8, type="l",
			 lwd=3, 
			xlab=" log tau", 
		   ylab = "ftilde, k=8"	
			 )
		lines(-log(Ss), ftilde108, col="blue")

		plot( analyticftilde2, ftilde102, type="p", pch=19, col=rgb(0.5,0.5,0.5,0.5),
			 xlab="Analytic", ylab="Computed") 
		abline(0,1)

		plot( analyticftilde8, ftilde108, type="p", pch=19, col=rgb(0.5,0.5,0.5,0.5),
			 xlab="Analytic", ylab="Computed") 
		abline(0,1)
}
