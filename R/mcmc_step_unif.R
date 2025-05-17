#' One step of the MCMC scheme using SS-Unif prior specification.
#'
#' @param leftDmat The design matrix for the left times.
#' @param rightDmat The design matrix for the right times.
#' @param tposInd A n x 1 vector indicating whether the event was observed before last follow-up.
#' @param obsInd A n x 1 vector indicating whether the event was observed after follow-up started.
#' @param piCurr A vector containing of the current pi values.
#' @param alphaCurr A matrix of current alpha values.
#' @param etaCurr A matrix of the current eta values.
#' @param betaCurr A matrix of the current beta values.
#' @param sigmaAprop Hyperparameter for the standard deviation of the proposal distribution for the spline coefficients.
#' @param sigmaBprop Hyperparameter for the standard deviation of the proposal distribution for the genetic effect coefficients.
#' @param sigmaEprop Hyperparameter for the standard deviation of the proposal distribution for the fixed effect coefficients.
#' @param sigmaAprior Hyperparameter for the standard deviation of the prior distribution for the spline coefficients.
#' @param sigmaBprior Hyperparameter for the standard deviation of the prior distribution for the genetic effect coefficients.
#' @param sigmaEprior Hyperparameter for the standard deviation of the prior distribution for the fixed effect coefficients.
#' @param gammaCurr A matrix of the current gamma values.
#' @param addProb Probability of adding another causal variant in selection scheme.
#' @param removeProb Probability of removing another causal variant in selection scheme.
#' @param alphaPriorMean A 1 x (nKnot + 2) vector of the spline coefficients.
#' @param etaPriorMean A 1 x p vector of the fixed effects.
#' @returns A n x q matrix of genetic data.
#' @return A list with the elements:
#' \item{betaCurr}{Matrix of current beta values.}
#' \item{alphaCurr}{Matrix of current alpha values.}
#' \item{etaCurr}{Matrix of current eta values.}
#' \item{gammaCurr}{Matrix of current gamma value}
#' \item{accMat}{A (B*burnIn) x (q*2 + (nKnot + 2) + p + 1) matrix containing indicators for whether the proposal was accepted (= 1) or not (=0).}

mcmc_step_unif <- function(leftDmat, rightDmat, tposInd, obsInd, piCurr, alphaCurr, etaCurr, betaCurr,
                           sigmaAprop, sigmaEprop, sigmaBprop, sigmaAprior, sigmaEprior, sigmaBprior,
                           gammaCurr, removeProb, addProb,alphaPriorMean, etaPriorMean) {

  # acceptance rates
  # holds pi, eta, gamma, beta, alpha
  accVec <- rep(0, 1 + length(etaCurr) + 2 * length(betaCurr) + length(alphaCurr))

  # linear predictors - to speed things up
  p <- length(gammaCurr)
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
    accProb <- min(1, exp(alphaNum - alphaDenom))
    # roll
    draw <- runif(n=1)
    if (draw <= accProb) {
      alphaCurr <- alphaProp
      accVec[alpha_it + 1 + q + 2*p] <- 1
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
        accVec[1 + eta_it] <- 1
      }
    }
  } # done eta
  # update etaPartCurr
  etaPartCurr <- as.numeric(leftDmat[, 1:q] %*% etaCurr)

  # propose gamma and beta
  zeroIdx <- which(gammaCurr == 0)
  oneIdx <- which(gammaCurr == 1)
  gammaProp <- gammaCurr
  betaProp <- betaCurr
  gammaRoll <- runif(n=1, min=0, max=1)

  # propose gamma
  currIn <- sum(gammaCurr)
  # subtract
  if (currIn == p | (gammaRoll <= removeProb & currIn > 0)) {
    # if full, must subtract
    changeIdx <- oneIdx[ceiling(runif(n=1, min=0, max=length(oneIdx)))]
    gammaProp[changeIdx] <- 0
    betaSub <- betaCurr[changeIdx]
    betaProp[changeIdx] <- 0

    # multiply this term by the ratio of likelihoods to get the acceptance probability
    if (currIn == 1) {
      cSub <- (1 / p) / removeProb
    } else if (currIn == length(gammaCurr)) {
      cSub <- addProb / (1 / p)
    } else {
      cSub <- (addProb / (p - currIn + 1)) / (removeProb / currIn)
    }

    # propose pi - acc/rej
    if(currIn == 0){
      piProp <- rbeta(n=1, shape1=1, shape2=length(gammaCurr))
    } else{
      piProp <- rbeta(n=1, shape1=currIn, shape2=length(gammaCurr) - currIn + 1)
    }
    # multiply by the part that accounts for the new beta
    # I used to have mean=betaPriorMean and sd=sigmaBprior here, but that's just wrong because this
    # is the proposal part
    propTerm <- log(cSub) + dnorm(betaSub, mean=0, sd=sigmaBprop, log=TRUE)
    # previously I forgot this prior term
    priorTerm <- dnorm(x=0, mean=0, sd=sigmaBprior, log=T) +
      log(currIn)  - log(p - currIn + 1) -
      ((currIn - 1) * log(piProp) ) + ((p - currIn + 1) * log(1 - piProp)) -
      (currIn * log(piCurr)) + ((p - currIn)*log(1 - piCurr)) +
      dnorm(x=betaSub, mean=0, sd=sigmaBprior, log=T) +
      log((1/piProp)*(1/(log(1) - log(1/10000))))  -
      log((1/piCurr)*(1/(log(1) - log(1/10000))))
  } else if (currIn == 0 | (gammaRoll > removeProb & gammaRoll <= removeProb + addProb & currIn < p) ) {
    # if nothing, must add
    changeIdx <- zeroIdx[ceiling(runif(n=1, min=0, max=length(zeroIdx)))]
    gammaProp[changeIdx] <- 1
    betaAdd <- rnorm(n=1, mean=0, sd=sigmaBprop)
    betaProp[changeIdx] <- betaAdd

    # multiply this term by the ratio of likelihoods to get the acceptance probability
    if (currIn == 0) {
      cAdd <- removeProb / (1 / p)
    } else if (currIn == p - 1) {
      cAdd <- (1/p) / addProb
    } else {
      cAdd <- (removeProb / (currIn + 1)) / (addProb / (p - currIn))
    }

    if(currIn == 0){
      piProp <- rbeta(n=1, shape1=1, shape2=length(gammaCurr))
    } else{
      piProp <- rbeta(n=1, shape1=currIn, shape2=length(gammaCurr) - currIn + 1)
    }
    # I used to have mean=betaPriorMean and sd=sigmaBprior here, but that's just wrong because this
    # is the proposal part
    propTerm <- log(cAdd) -  dnorm(betaAdd, mean=0, sd=sigmaBprop, log=TRUE)
    # previously I forgot this prior term
    priorTerm <- dnorm(x=betaAdd, mean=0, sd=sigmaBprior, log=T) +
      log(p - currIn) -  log(currIn + 1) -
      dnorm(x=0, mean=0, sd=sigmaBprior, log=T) +
      ((currIn + 1)* log(piProp)) + ((p - currIn - 1)* log(1 - piProp)) -
      (currIn* log(piCurr)) + ((p - currIn)* log(1 - piCurr)) +
      log((1/piProp)*(1/(log(1) - log(1/10000)))) -
      log((1/piCurr)*(1/(log(1) - log(1/10000))))
  } else {
    # swap
    addIdx <- zeroIdx[ceiling(runif(n=1, min=0, max=length(zeroIdx)))]
    subIdx <- oneIdx[ceiling(runif(n=1, min=0, max=length(oneIdx)))]
    changeIdx <- c(addIdx, subIdx)
    gammaProp[addIdx] <- 1
    gammaProp[subIdx] <- 0
    betaAdd <- rnorm(n=1, mean=0, sd=sigmaBprop)
    betaProp[addIdx] <- betaAdd
    betaSub <- betaCurr[subIdx]
    betaProp[subIdx] <- 0
    piProp = piCurr
    # I used to have mean=betaPriorMean and sd=sigmaBprior here, but that's just wrong because this
    # is the proposal part
    propTerm <- dnorm(betaSub, mean=0, sd=sigmaBprop, log=TRUE) - dnorm(betaAdd, mean=0, sd=sigmaBprop, log=TRUE)
    # previously I forget this part
    priorTerm <- dnorm(betaAdd, mean=0, sd=sigmaBprior, log=TRUE) - dnorm(betaSub, mean=0, sd=sigmaBprior, log=TRUE)
  }




  # full acceptance probability for gamma + beta change
  betaPartProp <- leftDmat[, (q+1):(q+length(betaCurr))] %*% (betaProp * gammaProp)
  logLikTerm <- calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                                 alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartProp) -
    calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                     alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartCurr)
  # acceptance probability for gamma + beta change
  accProb <- min(1, exp(propTerm + priorTerm + logLikTerm))

  # roll
  draw <- runif(n=1)
  if (draw <= accProb) {
    gammaCurr <- gammaProp
    betaCurr <- betaProp
    piCurr = piProp
    accVec[1] <- 1
    accVec[length(etaCurr) + 1 + changeIdx] <- 1
    accVec[length(etaCurr) + length(betaCurr) + 1 + changeIdx] <- 1
  }
  #print(piCurr)
  #print(gammaCurr)


  # loop through the other beta terms
  otherBetaIdx <- which(betaCurr != 0)
  if (length(otherBetaIdx) > 0) {
    otherBetaIdx <- otherBetaIdx[which(!(otherBetaIdx %in% changeIdx))]
  }
  if (length(otherBetaIdx) > 0) {
    for (beta_it in 1:length(otherBetaIdx)) {
      tempIdx <- otherBetaIdx[beta_it]
      betaProp <- betaCurr
      betaProp[tempIdx] <- rnorm(n=1, mean=betaCurr[tempIdx], sd=sigmaBprop)

      # get the acceptance probability
      betaPartProp <- leftDmat[, (q+1):(q+length(betaCurr))] %*% (betaProp * gammaCurr)
      betaNum <- dnorm(x=betaProp[tempIdx], mean=0, sd=sigmaBprior, log=TRUE) +
        calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                         alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartProp)
      betaDenom <- dnorm(x=betaCurr[tempIdx], mean=0, sd=sigmaBprior, log=TRUE) +
        calculate_logLik(tposInd = tposInd, obsInd = obsInd, alphaLeft = alphaLeftCurr,
                         alphaRight = alphaRightCurr, etaPart = etaPartCurr, betaPart = betaPartCurr)
      # acceptance probabilites for beta
      accProb <- min(1, exp(betaNum - betaDenom))
      # roll
      draw <- runif(n=1)
      if (draw <= accProb) {
        betaCurr <- betaProp
        accVec[length(alphaCurr) + length(etaCurr) + length(betaCurr) + 1 + tempIdx] <- 1
      }
    }
  } # done proposing beta

  # return
  return(list(betaCurr=betaCurr, alphaCurr=alphaCurr, etaCurr=etaCurr, piCurr=piCurr, gammaCurr=gammaCurr, accVec=accVec))
}
