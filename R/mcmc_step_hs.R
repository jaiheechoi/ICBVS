#' One step of the MCMC scheme using the horseshoe prior specification.
#'
#' @param leftDmat Left design matrix
#' @param rightDmat Right design matrix
#' @param tposInd A n x 1 vector indicating whether the event was observed before last follow-up.
#' @param obsInd A n x 1 vector indicating whether the event was observed after follow-up started.
#' @param zetaCurr Current zeta values.
#' @param alphaCurr Current alpha values.
#' @param etaCurr Current eta values.
#' @param betaCurr Current beta values.
#' @param sigmaAprop Hyperparameter for the standard deviation of the proposal distribution for the spline coefficients.
#' @param sigmaBprop Hyperparameter for the standard deviation of the proposal distribution for the genetic effect coefficients.
#' @param sigmaEprop Hyperparameter for the standard deviation of the proposal distribution for the fixed effect coefficients.
#' @param sigmaAprior Hyperparameter for the standard deviation of the prior distribution for the spline coefficients.
#' @param sigmaBprior Hyperparameter for the standard deviation of the prior distribution for the genetic effect coefficients.
#' @param sigmaEprior Hyperparameter for the standard deviation of the prior distribution for the fixed effect coefficients.
#' @param gammaCurr Current gamma values.
#' @param addProb Probability of adding another causal variant in selection scheme.
#' @param removeProb Probability of removing another causal variant in selection scheme.
#' @param alphaPriorMean A 1 x (nKnot + 2) vector of the spline coefficients.
#' @param etaPriorMean A 1 x p vector of the fixed effects.
#' @param A Hyperparameter in rate parameter in gamma prior. Default set to 1.
#' @param subPct A single value from 0 to 1 indicating the percentage of genetic variants to update each iteration of the MCMC. Default set to 0.5.
#' @returns A n x q matrix of genetic data.
#' @return A list with the elements:
#' \item{betaCurr}{Matrix of current beta values.}
#' \item{alphaCurr}{Matrix of current alpha values.}
#' \item{etaCurr}{Matrix of current eta values.}
#' \item{gammaCurr}{Matrix of current gamma value}
#' \item{accMat}{A (B*burnIn) x (q*2 + (nKnot + 2) + p + 1) matrix containing indicators for whether the proposal was accepted (= 1) or not (=0).}

mcmc_step_hs<- function(leftDmat, rightDmat, tposInd, obsInd, zetaCurr, alphaCurr, etaCurr, betaCurr,
                        sigmaAprop, sigmaEprop, sigmaBprop, sigmaAprior, sigmaEprior, sigmaBprior,
                        gammaCurr, removeProb, addProb,alphaPriorMean, etaPriorMean, A = 1, subPct = 0.5) {

  # acceptance rates
  # holds zeta, gamma, eta, alpha, beta
  accVec <- rep(0, 2 + length(etaCurr) + length(alphaCurr) + length(betaCurr) )# last is beta

  # linear predictors - to speed things up
  p <- length(betaCurr)
  q <- length(etaCurr)
  if (q> 0) {
    etaPartCurr <- as.numeric(leftDmat[, 1:q] %*% etaCurr)
  } else {
    etaPartCurr <- rep(0, nrow(leftDmat))
  }
  betaPartCurr <- as.numeric(leftDmat[, (q+1):(q+length(betaCurr))] %*% betaCurr)
  alphaLeftCurr <- as.numeric(leftDmat[, (q+length(betaCurr)+1):ncol(leftDmat)] %*% alphaCurr)
  alphaRightCurr <- as.numeric(rightDmat[, (q+length(betaCurr)+1):ncol(rightDmat)] %*% alphaCurr)


  # propose alpha
  for (alpha_it in 1:length(alphaCurr)) {
    # propose alpha from normal centered at current value
    alphaProp <- alphaCurr
    alphaProp[alpha_it] <- rnorm(n=1, mean=alphaCurr[alpha_it], sd=sigmaAprop)

    # get the acceptance probability
    alphaLeftProp <- as.numeric(leftDmat[, (q+length(betaCurr)+1):ncol(leftDmat)] %*% alphaProp)
    alphaRightProp <- as.numeric(rightDmat[, (q+length(betaCurr)+1):ncol(rightDmat)] %*% alphaProp)
    alphaNum <- dnorm(x=alphaProp[alpha_it], mean=alphaPriorMean[alpha_it], sd=sigmaAprior, log=TRUE) +
      calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftProp,
                       alphaRight = alphaRightProp, etaPart = etaPartCurr, betaPart = betaPartCurr)
    alphaDenom <- dnorm(x=alphaCurr[alpha_it], mean=alphaPriorMean[alpha_it], sd=sigmaAprior, log=TRUE) +
      calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                       alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartCurr)
    # acceptance probabilites for alpha
    if(alphaNum == -Inf & alphaDenom == -Inf){
      accProb = 0
    } else{
      accProb <- min(1, exp(alphaNum - alphaDenom))}
    # accProb <- min(1, exp(alphaNum - alphaDenom))
    # roll
    draw <- runif(n=1)
    if (draw <= accProb) {
      alphaCurr <- alphaProp
      accVec[alpha_it + 2 + length(etaCurr)] <- 1
    }
  } # done proposing alpha
  # update alphaLeftCurr and alphaRightCurr
  alphaLeftCurr <- as.numeric(leftDmat[, (q+length(betaCurr)+1):ncol(leftDmat)] %*% alphaCurr)
  alphaRightCurr <- as.numeric(rightDmat[, (q+length(betaCurr)+1):ncol(rightDmat)] %*% alphaCurr)

  # propose eta
  if (length(etaCurr) > 0) {
    for (eta_it in 1:length(etaCurr)) {
      # propose eta from normal centered at existing eta
      etaProp <- etaCurr
      etaProp[eta_it] <- rnorm(n=1, mean=etaCurr[eta_it], sd=sigmaEprop)

      # get the acceptance probability
      etaPartProp <- as.numeric(leftDmat[, 1:q] %*% etaProp)
      etaNum <- dnorm(x=etaProp[eta_it], mean=etaPriorMean[eta_it], sd=sigmaEprior, log=TRUE) +
        calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                         alphaRight = alphaRightCurr, etaPart = etaPartProp, betaPart = betaPartCurr)
      etaDenom <- dnorm(x=etaCurr[eta_it], mean=etaPriorMean[eta_it], sd=sigmaEprior, log=TRUE) +
        calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                         alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartCurr)
      # acceptance probabilites for eta
      accProb <- min(1, exp(etaNum - etaDenom))
      # roll
      draw <- runif(n=1)
      if (draw <= accProb) {
        etaCurr <- etaProp
        accVec[2 + eta_it] <- 1
      }
    }
  } # done eta


  # sample eta
  zetaProp = rgamma(n = p, shape = 1,
                    rate = betaCurr^2/2 + gammaCurr)
  zetaCurr <- zetaProp

  # sample gamma
  gammaProp = rgamma(n = p, shape = 1, rate = zetaCurr + 1/A^2)
  gammaCurr <- gammaProp
  accVec[1:2] <- 1

  subTot <- round(length(betaCurr)*subPct)
  betaSubIdx <- sort(sample(1:length(betaCurr), subTot))

  for (beta_it in 1:length(betaSubIdx)) {
    betaIdx <- betaSubIdx[beta_it]
    betaProp <- betaCurr
    betaProp[betaIdx] <- rnorm(n=1, mean=betaCurr[betaIdx], sd=sigmaBprop)


    # get the acceptance probability
    betaPartProp <- leftDmat[, (q+1):(q+length(betaCurr))] %*% betaProp

    # previously I forget this part
    priorTerm <- dnorm(x=betaProp[betaIdx], mean=0, sd=(zetaCurr[betaIdx])^(-1), log=TRUE)-
      dnorm(x=betaCurr[betaIdx], mean=0, sd=(zetaCurr[betaIdx])^(-1), log=TRUE)

    # full acceptance probability for gamma + beta change
    logLikTerm <- calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                                   alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartProp) -
      calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                       alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartCurr)
    # acceptance probability for gamma + beta change
    accProb <- min(1, exp(priorTerm + logLikTerm))
    # roll
    draw <- runif(n=1)
    if (draw <= accProb) {
      betaCurr <- betaProp
      accVec[length(alphaCurr) + length(etaCurr) + 2 + betaIdx] <- 1
    }
  }


  # return
  return(list(betaCurr=betaCurr, alphaCurr=alphaCurr, etaCurr=etaCurr, zetaCurr = zetaCurr, gammaCurr=gammaCurr, accVec=accVec))
}

