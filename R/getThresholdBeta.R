#' Type I error threshold for estimated genetic effects. 
#' 
#' @param q Total number of SNPs being fine-mapped.  
#' @param estBeta Estimated effect sizes.
#' @param causal True causal variants.
#' @param T1E Desired Type 1 Error rate. 
#' @returns A single value indicating the threshold for beta estimate. 

getThresholdBeta <- function(q, estBeta, causal, T1E){
  check <- round(p/2)
  t1eRates <- rep(NA, check)
  numNotCausal <- p - length(causal)  
  for(i in 1:check){
    selected <- order(estBeta, decreasing = T)[1:i]
    correct <- selected[which(selected %in% causal)]
    notCorrect <- selected[which(!(selected %in% causal))]
    t1eRates[i] <- length(notCorrect)/(numNotCausal)
  }
  
  # threshold <- sort(abs(estBeta), decreasing = T)[max(which(t1eRates < T1E))]
  threshold <- sort(estBeta, decreasing = T)[max(which(t1eRates < T1E))]
  
  return(threshold)
}