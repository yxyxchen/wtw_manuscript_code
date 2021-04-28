# these functions take in an input and return an output of the same size with
# some form of normalization

# want to return only positive values
clip <- function(x){
	pos <- x > 0
	return(pos * x)
}

# Based on Eq 5 in the Carrandini and Heeger 2011 paper
divisivenorm <- function(input, sigma=0, exponent = 2){
	num <- input^exponent
	denom <- sum(num) + sigma
	return(num/denom)	
}

# generic softmax function
softmax <- function(input, betascale=1){
	beta <- betascale/max(input)
	num <- exp(beta * input)
	denom <- sum(num) 
	return(num/denom)	
}

# here's a function to ``increment the k'' of ftilde by exponentiation.
# compare to the notes from April 17, 2021
multiplyk <- function(ftildek, svals, n=6, k=2){
	prefactor <- n^{n * k + 1}
	sexp <- 1-n*(k+1)
	return((ftildek^n)*(svals^(sexp)))
}

# here's a function for impulse response of ftilde.  I'm ignoring Ck
analyticftilde <- function(tau=1,svals,k){
	stauterm <- svals*tau
	return((svals*(stauterm^k) * exp(-stauterm)))
}

# here's a function to ``increment the k'' of ftilde by exponentiation.
# compare to the notes from April 17, 2021
# i'm ignoring Ck
multiplyk <- function(ftildek, svals, n=6, k=2){
	prefactor <- n^(n * k + 1)
	sexp <- 1-n
#	\scalefunction = \frac{C_{nk}}{C_k^n}\ n^{nk+1}\  s^{1-n(k+1)} 
	Gofs <- prefactor * (svals^(sexp))
	return((ftildek^n)*(Gofs))
}



