# Taken from JNeuro '14 code.
# This generates \L^{-1}_k .  Power is k.   I'm asking for svals and deltas
# separately because i want the last s value to make sense.   if mults is
# true, it multiplies each row by s^k to make it \Lk
# mask eliminates edges. 
# exponent appears as s^{k+1 + exponent}.  Set exponent to 0 to get \Lk.  Set
# it to 1 to have all the cells have the same peak.  Other values are not
# recommended. 
GetLk <- function(svals=array(1,100), deltas=NULL, power=2, mults=TRUE,
				  mask=T, exponent=0){
	# how many svals do you have
	size <- length(svals)
	# make sure s is sorted
	svals <- sort(svals)
	# This is the matrix that discretely approximates the first derivative
	Delta <- array(0, c(size,size))

	# if you don't specify the deltas, the function figures it out from
	# the s values 
	if(is.null(deltas)){
		deltas <- array(0,size+1)
	    for(i in 2:(size)){
			deltas[i] = Ss[i] - Ss[i-1]
		}
		# try to minimize the edge effect
		deltas[1] <- deltas[2]
	}
	deltas[size+1] <- deltas[size]
	print(paste("length of deltas", length(deltas)))
	#diag(retval) <- 1/(size:1) * 1/(size:1)
	for(i in 1:(size)){
		# referring to Eq 18 in the JMachLearnR paper,
		# s_o - s_{-1} = deltas[i]
		# s_1 - s_o = deltas[i+1]
		# s_1 - s_{-1} = deltas[i] + deltas[i+1]
		s2 <- deltas[i] + deltas[i+1]
		if(i < size){
			rat1 <- deltas[i]/(deltas[i+1]*s2)
			Delta[i,i+1] <- rat1
			Delta[i,i] <- Delta[i,i] - rat1
		}
		if(i > 1){
			rat2 <- deltas[i+1]/(deltas[i] * s2)
			Delta[i,i-1] <- -rat2 
			Delta[i,i] <- Delta[i,i] + rat2
		}
	}
#	print(paste("rat1", rat1,"   rat2", rat2))
	retval <- Delta
#	return(retval)
	# raise Delta^power
	if (power > 1){
		print(paste("Raising matrix Delta to power",power))
		for(i in 2:power){
			print(".")
			retval <- retval %*% Delta
		}
	}
	if(mults){
			# at the end, make it be Linvk by multiplying each row by s^{k+1}.
			for(i in 1:size){
				retval[i,] <- retval[i,] * ((svals[i])^(power+1-exponent))
			}
	}
	
	Lkmask <- array(1, length(svals))
	Lkmask[1:(2*power-1)] <- 0
	Lkmask[length(svals):(length(svals) - (2*power-1))] <- 0

	if(mask){
	 retval <- diag(Lkmask) %*% retval	
	}



	return(retval)
}
