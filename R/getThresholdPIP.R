#' Type I error threshold for estimated posterior inclusion probabilities. 
#' 
#' @param q Total number of SNPs being fine-mapped.  
#' @param pip A 1 x q vector of the estimated posterior inclusion probabilities. 
#' @param causal True causal variants.
#' @param T1E Desired Type 1 Error rate. 
#' @returns A single value indicating the threshold for beta estimate. 

getThresholdPIP <- function(p, pip, causal, T1E){
  check <- round(p/2)
  numNotCausal <- p - length(causal)
  
  
  t1eRates <- rep(NA, check)
  
  for(i in 1:check){
    selected <- order(pip, decreasing = T)[1:i]
    correct <- selected[which(selected %in% causal)]
    notCorrect <- selected[which(!(selected %in% causal))]
    t1eRates[i] <- length(notCorrect)/(numNotCausal)
  }
  
  threshold <- sort(pip, decreasing = T)[max(which(t1eRates < T1E))]
  
  return(threshold)
}