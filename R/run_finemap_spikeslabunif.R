#' Function for MCMC for interval-censored spike and slab using uniform setting.
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
#' @param B Number of MCMC iterations (excluding burn-in). (Default set to 10000 iterations)
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
#' run_finemap_spikeslabunif(xMat = cbind(xMat, gMat), q = 2, leftTimes = lt, rightTimes = rt, nKnots = nKnots,
#'   tposInd = tposInd, obsInd = obsInd,
#'   sigmaAprior = sigmaAprior, sigmaEprior = sigmaEprior,
#'   sigmaBprior = sigmaBprior, B=B,
#'   sigmaAprop = sigmaAprop, sigmaEprop = sigmaEprop, sigmaBprop = sigmaBprop,
#'   addProb = addProb, removeProb = removeProb, seed=0, burnIn=0.2,
#'   alphaPriorMean=rep(1, nKnots + 2), etaPriorMean=rep(0, ncol(xMat)), checkpoint=TRUE)

# run and record
run_finemap_spikeslabunif<- function(xMat, q, leftTimes, rightTimes, nKnots=1, tposInd, obsInd, sigmaAprop, sigmaEprop, sigmaBprop,
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
    mcmcOut <- mcmc_step_unif(leftDmat = dmats$left_dmat, rightDmat =  dmats$right_dmat, tposInd = tposInd,
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
