###3
seqAppend = function(seq,nUnique){
  library('stringr')
  #
  avail = 1:nUnique; # potential elements
  # % only balance 1st-order transitions (length 2)
  # % this introduces a weak temporal autocorrelation structure
  # % (autocorrelation from balancing 0th-order transitions is much larger)
  subLength = 2; # length of subsequence to balance (hardcoded)
  candidates = avail;
  nCand = length(candidates);
  cost = rep(0,nCand);
  if(length(seq)>=subLength){
    for(i in 1:nCand){  # evaluate each current candidate
    provisSeq = c(seq, candidates[i]); # the full sequence that candidate i would create
    newSubseq = provisSeq[(length(provisSeq)-subLength+1):length(provisSeq)] # the relevant subsequence
    seqStr = toString(seq)
    seqStr = gsub(", ", "", seqStr)
    newSeqStr = toString(newSubseq)
    newSeqStr = gsub(", ", "", newSeqStr)
    cost[i] = str_count(seqStr, newSeqStr) # number of subsequence repetitions
    } #provided this subsequence length exists
  
  # only keep candidates with the fewest subsequence repetitions
  candidates = candidates[cost==min(cost)];
  
  }
  # choose randomly among the candidates
  nextItem = candidates[sample(length(candidates), 1)];
  
  # append it
  newSeq = c(seq, nextItem)
  
  # return
  outputs = list("seq" = newSeq, "nextItem" = nextItem)
  return(outputs)
}

drawSample = function(distrib,seq){
  # % generates a sample from next quantile of the designated distribution
  # % timing parameters are specified within this function
  # 
  # % identify the distribution
  n = 8
  dist = matrix(NA, 2, n)
  if(distrib == 'unif16'){
    dist[1,] = seq(2, 16, by = 2);
    dist[2,] = rep(1, ncol(dist));
  }else{
    m = 32; # maximum delay
    fac = 1.75;
    d1 = log(m/(fac^n - 1));
    d2 = log(m + m/(fac^n - 1));
    tempt = exp(seq(d1,d2,length.out = n+1));
    dist[1, ] = round(tempt[2:length(tempt)] - tempt[1], 3); # subtract 1st to normalize relative to zero.
    dist[2, ] = rep(1,  ncol(dist));
  }
  
  # for any of the real distributions
  # find the index of the next sequence element
  nUnique = sum(dist[2,]); # number of elements to randomize
    seqOutputs = seqAppend(seq,nUnique)


  nextItem = seqOutputs[['nextItem']]
  seq = seqOutputs[['seq']]
  
  # % map the selected randomization element (stored in nextItem) to a
  # % delay value, allowing for distributions in which delay values have
  # % unequal probabilities.
  idxCol = cumsum(dist[2,]);
  sampIdx = min(which(nextItem<=idxCol));
  delay = dist[1,sampIdx];
  
  # 
  outputs = list('delay' = delay, 'seq' = seq, 'dist' = dist)
  return(outputs)
}


