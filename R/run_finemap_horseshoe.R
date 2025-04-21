# run and record
run_finemap_horseshoe <- function(xMat, p, q, leftTimes, rightTimes, nKnots=1, tposInd, obsInd, sigmaAprop, sigmaEprop, sigmaBprop,
                          sigmaAprior, sigmaEprior, sigmaBprior, A, B, subPct, quant_r=NULL, addProb, removeProb,
                          seed=0, burnIn=0.2, alphaPriorMean=NULL, etaPriorMean = NULL, checkpoint=FALSE) {
  
  # initialize
  #set.seed(0)
  alphaCurr <- rep(1, nKnots + 2)
  etaCurr <- rep(0, q)
  betaCurr <- rnorm(p, mean = 0, sd = 1)
  gammaCurr <- rgamma(n = p, shape = 1, rate = 1)
  zetaCurr <- rgamma(n = p, shape = 1, rate = 1)
  
  # make design matrix
  dmats <- make_IC_dmat(xMat = xMat, lt=leftTimes, rt=rightTimes, nKnots=nKnots, tpos_ind=tposInd,
                        obs_ind = obsInd, quant_r=quant_r)
  
  # for saving
  betaChain <- matrix(data=NA, nrow=B, ncol=length(betaCurr))
  etaChain <- matrix(data=NA, nrow=B, ncol=length(etaCurr))
  alphaChain <- matrix(data=NA, nrow=B, ncol=length(alphaCurr))
  zetaChain <- matrix(data=NA, nrow=B, ncol=length(betaCurr))
  gammaChain <- matrix(data=NA, nrow=B, ncol=length(betaCurr))
  accMat <- matrix(data=NA, nrow=B, ncol=(2 + length(etaCurr) + length(alphaCurr) + length(betaCurr) ))
  
  # loop
  for (mc_it in 1:B) {
    # run MCMC one interation
    #print(mc_it)
    mcmcOut <- mcmc_step_hs(leftDmat = dmats$left_dmat, rightDmat =  dmats$right_dmat, tposInd = tposInd, 
                             obsInd = obsInd, zetaCurr = zetaCurr, alphaCurr = alphaCurr, etaCurr = etaCurr, betaCurr = betaCurr, 
                             sigmaAprop = sigmaAprop, sigmaEprop = sigmaEprop, sigmaBprop = sigmaBprop,
                             sigmaAprior = sigmaAprior, sigmaEprior = sigmaEprior, sigmaBprior = sigmaBprior,
                             gammaCurr = gammaCurr, removeProb=removeProb, addProb=addProb,
                             alphaPriorMean = alphaPriorMean, etaPriorMean = etaPriorMean, A = A, subPct = subPct)
    #leftDmat = dmats$left_dmat; rightDmat =  dmats$right_dmat
    
    # update current
    betaCurr <- mcmcOut$betaCurr
    etaCurr = mcmcOut$etaCurr
    alphaCurr <- mcmcOut$alphaCurr
    gammaCurr <- mcmcOut$gammaCurr
    zetaCurr <- mcmcOut$zetaCurr
    
    # record
    betaChain[mc_it, ] <- betaCurr
    etaChain[mc_it, ] <- etaCurr
    alphaChain[mc_it, ] <- alphaCurr
    gammaChain[mc_it, ] <- gammaCurr
    zetaChain[mc_it, ] <- zetaCurr
    accMat[mc_it, ] <- mcmcOut$accVec
    
    if (checkpoint) {if(mc_it%%100 == 0) {cat(mc_it, '\n')}}
  }
  
  # remove burn in
  burn_idx <- round(burnIn * B) + 1
  return(list( betaChain=betaChain[burn_idx:B, ], alphaChain=alphaChain[burn_idx:B, ], 
               etaChain=etaChain[burn_idx:B, ],
               zetaChain = zetaChain[burn_idx:B, ],
               gammaChain = gammaChain[burn_idx:B,],
               accMat = accMat[burn_idx:B, ]) )
}


