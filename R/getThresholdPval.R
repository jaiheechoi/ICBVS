getThresholdPval <- function(p, pvals, causal, T1E){
  check <- round(p/2)
  numNotCausal <- p - length(causal)
  
  
  t1eRates <- rep(NA, check)
  
  for(i in 1:check){
    selected <- order(pvals, decreasing = F)[1:i]
    correct <- selected[which(selected %in% causal)]
    notCorrect <- selected[which(!(selected %in% causal))]
    t1eRates[i] <- length(notCorrect)/(numNotCausal)
  }
  
  threshold <- sort(pvals, decreasing = F)[max(which(t1eRates < T1E))]
  
  return(threshold)
}
