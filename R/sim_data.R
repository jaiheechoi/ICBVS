#' Simulate sample interval-censored data.
#'
#' @param seed Random seed.
#' @param n Total number of observations.
#' @param q Total number of SNPs.
#' @param bnFunInv The inverse of the baseline hazard function.
#' @param obsTimes Vector of the intended observation times.
#' @param windowHalf The amount of time before or after the intended obsTimes that a visit might take place.
#' @param trueEta A vector of true effect sizes for the fixed effects.
#' @param trueBeta A vector of true effect sizes for the random effects.
#' @param gMat A n x q matrix of the genetic data.
#' @returns A list containing the outcome times, covariate data, and genetic data.

sim_data <- function(seed=0, n, q, bhFunInv, obsTimes, windowHalf,
                     trueEta=NULL, trueBeta, gMat) {
  # usually you'd call it with something like
  # seed=0; n=3000; p=10; mafVec = rep(0.3, p); numCausal=3; bhFunInv = function(x) {x};
  # obsTimes = 1:5; windowHalf = 0.1; trueEta = c(-1, 0); trueBeta = rep(0, p);
  # trueBeta[c(2, 7, 8)] <- c(-1.5, 1.5, 1.5)

  #set.seed(seed)
  # this part is the non-SNP covariates forced in
  if (!is.null(trueEta)) {
    xMat <- matrix(data=rnorm(n=length(trueEta) * n), nrow=n)
    xEff <- xMat %*% trueEta
  } else {
    xEff <- rep(0, n)
  }

  # the SNPs
  #gMat <- matrix(data=rbinom(n=n*p, size=2, prob=mafVec), nrow=n, ncol=p)
  linearPred <- xEff + gMat %*% trueBeta

  # generate otucome
  outcomeDat <- ICSKAT::gen_IC_data(bhFunInv = bhFunInv, obsTimes = obsTimes, windowHalf = windowHalf,
                                    probMiss = 0, etaVec = linearPred)

  # return
  return(list(outcomeDat = outcomeDat, xMat = xMat, gMat = gMat))
}
