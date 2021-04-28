
# here's a function to evolve F(t).  Takes in an F set of s values, a time since
# the last update and an input vector.  Checks the dimensionality of
# everything

EvolveF <- function(PrevF, svals, dt=1, fin){
	numfs <- length(fin)
	if(numfs != length(PrevF[1,])){
		stop(paste("EvolveF -> fin has dim", numfs, "but PrefF has dim",
				   length(PrevF[1,])))
	}
	nums <- length(Ss)
	if(nums != length(PrevF[,1])){
		stop(paste("EvolveF -> svals has dim", nums, "but PrefF has dim",
				   length(PrevF[,1])))
	}
	evolvefac <- exp(-dt*svals)
	retval <- diag(evolvefac) %*% PrevF
	retval <- retval + array(1,nums) %*% t(fin)
	return(retval)
}

