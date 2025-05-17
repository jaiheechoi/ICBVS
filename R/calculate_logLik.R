#' Calculate log-likelihood for acceptance ratio in MCMC scheme.
#' 
#' @param tposInd A n x 1 vector indicating whether the event was observed before last follow-up.
#' @param obsInd A n x 1 vector indicating whether the event was observed after follow-up started.
#' @param alphaLeft A n x p matrix of the left time point fixed coefficients times the covariates.
#' @param alphaRight  A n x p matrix of the right time point fixed coefficients times the covariates.
#' @param etaPart  A n x (nKnots + 2) matrix of the spline coefficients times cubic spline values.
#' @param betaPart A n x q matrix of the estimated genetic coefficients times the SNP data.
#' @returns A scalar value summing the log likelihood of the left and right time points. 

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