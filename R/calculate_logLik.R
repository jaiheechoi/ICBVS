# calculate the log likelihood of the data - give it the linear predictors
calculate_logLik <- function(tposInd, obsInd, alphaLeft, alphaRight, etaPart, betaPart) {
  # the survival terms
  # the design matrix should be the non-SNPs, then SNPs, then spline terms
  HL <- as.numeric(exp(alphaLeft + etaPart + betaPart))
  HR <- as.numeric(exp(alphaRight + etaPart + betaPart))
  SL <- ifelse(tposInd == 0, 1, exp(-HL))
  SR <- ifelse(obsInd == 0, 0, exp(-HR))
  diffS <- SL - SR
 
  # sometimes the log hazard will be decreasing - in that case just return -Inf for ll
  hazCheck <- which(diffS < 0)
  if (length(hazCheck) > 0) {
    return(-Inf)
  }
  
  # check for SL = SR
  equalCheck <- which(diffS == 0)
  if (length(equalCheck) > 0) {
    diffS[equalCheck] <- min(diffS[-equalCheck])
  }
  
  return(sum(log(diffS)))
}
