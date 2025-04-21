# run and record
run_finemap_spikeslabsparse <- function(xMat, q, leftTimes, rightTimes, nKnots=1, tposInd, obsInd, sigmaAprop, sigmaEprop, sigmaBprop,
                                  sigmaAprior, sigmaEprior, sigmaBprior, B=10000, quant_r=NULL, addProb, removeProb,
                                  seed=0, burnIn=0.2, alphaPriorMean=NULL, etaPriorMean = NULL, checkpoint=FALSE) {
  
  # initialize
  #set.seed(0)
  alphaCurr <- rep(1, nKnots + 2)
  etaCurr <- rep(0, q)
  betaCurr <- rep(0.1, ncol(xMat) - q)
  gammaCurr <- rbinom(n=ncol(xMat) - q, size=1, prob=0.2)
  betaCurr <- betaCurr * gammaCurr
  #piCurr <- exp(runif(n=1, min=log(1 / (ncol(xMat) - q)), max=0))
  piCurr = .2
  # make design matrix
  dmats <- make_IC_dmat(xMat = xMat, lt=leftTimes, rt=rightTimes, nKnots=nKnots, tpos_ind=tposInd,
                        obs_ind = obsInd, quant_r=quant_r)
  
  # for saving
  betaChain <- matrix(data=NA, nrow=B, ncol=length(betaCurr))
  etaChain <- matrix(data=NA, nrow=B, ncol=length(etaCurr))
  alphaChain <- matrix(data=NA, nrow=B, ncol=length(alphaCurr))
  gammaChain <- matrix(data=NA, nrow=B, ncol=length(betaCurr))
  piChain <- matrix(data=NA, nrow=B, ncol=1)
  accMat <- matrix(data=NA, nrow=B, ncol=nKnots + 3 + q + 2*length(betaCurr))
  
  # loop
  for (mc_it in 1:B) {
    # run MCMC one interation
    mcmcOut <- mcmc_step_sparse(leftDmat = dmats$left_dmat, rightDmat =  dmats$right_dmat, tposInd = tposInd, 
                                obsInd = obsInd, alphaCurr = alphaCurr, etaCurr = etaCurr, betaCurr = betaCurr, 
                                sigmaAprop = sigmaAprop, sigmaEprop = sigmaEprop, sigmaBprop = sigmaBprop,
                                sigmaAprior = sigmaAprior, sigmaEprior = sigmaEprior, sigmaBprior = sigmaBprior,
                                gammaCurr = gammaCurr, piCurr = piCurr, removeProb=removeProb, addProb=addProb,
                                alphaPriorMean = alphaPriorMean, etaPriorMean = etaPriorMean)
    #leftDmat = dmats$left_dmat; rightDmat =  dmats$right_dmat
    
    # update current
    betaCurr <- mcmcOut$betaCurr
    etaCurr = mcmcOut$etaCurr
    alphaCurr <- mcmcOut$alphaCurr
    gammaCurr <- mcmcOut$gammaCurr
    piCurr <- mcmcOut$piCurr
    
    # record
    betaChain[mc_it, ] <- betaCurr
    etaChain[mc_it, ] <- etaCurr
    alphaChain[mc_it, ] <- alphaCurr
    gammaChain[mc_it, ] <- gammaCurr
    piChain[mc_it, ] <- piCurr
    accMat[mc_it, ] <- mcmcOut$accVec
    
    if (checkpoint) {if(mc_it%%100 == 0) {cat(mc_it, '\n')}}
  }
  
  # remove burn in
  burn_idx <- round(burnIn * B) + 1
  return(list( betaChain=betaChain[burn_idx:B, ], alphaChain=alphaChain[burn_idx:B, ], 
               etaChain=etaChain[burn_idx:B, ],
               gammaChain=gammaChain[burn_idx:B, ], piChain = piChain[burn_idx:B, ],
               accMat = accMat[burn_idx:B, ]) )
}


