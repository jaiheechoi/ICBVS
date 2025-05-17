#' Function for MCMC for interval-censored horseshoe.
#'
#' @param xMat A n x p matrix of covariate data.
#' @param q Total number of SNPs.
#' @param rightTimes A n x 1 vector of observed left times of the event of interest.
#' @param leftTimes A n x 1 vector of observed right times of the event of interest.
#' @param nKnots Number of internal knots in cubic spline for estimating baseline hazard function.
#' @param tposInd A n x 1 vector indicating whether the event was observed before last follow-up.
#' @param obsInd A n x 1 vector indicating whether the event was observed after follow-up started.
#' @param sigmaAprop Hyperparameter for the standard deviation of the proposal distribution for the spline coefficients.
#' @param sigmaBprop Hyperparameter for the standard deviation of the proposal distribution for the genetic effect coefficients.
#' @param sigmaEprop Hyperparameter for the standard deviation of the proposal distribution for the fixed effect coefficients.
#' @param sigmaAprior Hyperparameter for the standard deviation of the prior distribution for the spline coefficients.
#' @param sigmaBprior Hyperparameter for the standard deviation of the prior distribution for the genetic effect coefficients.
#' @param sigmaEprior Hyperparameter for the standard deviation of the prior distribution for the fixed effect coefficients.
#' @param A Hyperparameter in rate parameter in gamma prior. Default set to 1.
#' @param B Number of MCMC iterations (excluding burn-in). (Default set to 10000 iterations)
#' @param subPct A single value from 0 to 1 indicating the percentage of genetic variants to update each iteration of the MCMC. Default set to 0.5.
#' @param quant_r Quantiles of time to use in constructing the spline.
#' @param addProb Probability of adding another causal variant in selection scheme.
#' @param removeProb Probability of removing another causal variant in selection scheme.
#' @param seed Random seed. (Default set to 0)
#' @param burnIn Percentage of first beginning iterations removed from MCMC. (Default set to 0.2)
#' @param alphaPriorMean A nKnot + 1 vector for original spline coefficient estimates. (Default set to NULL)
#' @param A p x 1 vector for original fixed effect coefficient estimates. (Default set to NULL)
#' @param checkpoint Boolean indicating whether to print results every 100 iterations. (Default set to FALSE)
#' @return A list with the elements:
#' \item{betaChain}{A (B*burnIn) x q matrix containing all genetic effect coefficient estimates from all MCMC iterations, excluding burn-in.}
#' \item{alphaChain}{A (B*burnIn) x (nKnot + 2) matrix containing all spline coefficient estimates from all MCMC iterations, excluding burn-in.}
#' \item{etaChain}{A (B*burnIn) x p matrix containing all fixed effect coefficient estimates from all MCMC iterations, excluding burn-in.}
#' \item{gammaChain}{A (B*burnIn) x q matrix containing all inclusion indicator from all MCMC iterations, excluding burn-in.}
#' \item{piChain}{A (B*burnIn) x 1 vector containing all sparsity parameter estimates from all MCMC iterations, excluding burn-in.}
#' \item{accMat}{A (B*burnIn) x (q*2 + (nKnot + 2) + p + 1) matrix containing indicators for whether the proposal was accepted (= 1) or not (=0).}
#' @examples
#' run_finemap_horseshoe(xMat = cbind(xMat, gMat), p = 2, q = 100, leftTimes = lt, rightTimes = rt, nKnots = nKnots,
#' tposInd = tposInd, obsInd = obsInd,
#' sigmaAprior = sqrt(2), sigmaEprior = sqrt(2),
#' sigmaBprior = sqrt(2), A = 1, B=100000, subPct = 0.5,
#' sigmaAprop = 1, sigmaEprop = 1, sigmaBprop = 1,
#' addProb = 0.4, removeProb = 0.4, seed=0, burnIn=0.2,
#' alphaPriorMean=rep(1, 3), etaPriorMean=rep(0, 2), checkpoint=TRUE)


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
