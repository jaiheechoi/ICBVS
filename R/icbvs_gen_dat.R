#' Generating interval censored data.
#' 
#' @param n Total number of observations.  
#' @param gMat A n x q matrix of the genetic data. 
#' @param noCausal Number of causal genetic variants. 
#' @param rho Common pairwise correlation parameter.
#' @param ES Effect size of a causal genetic variant. 
#' @param bhFunInv The inverse of the baseline hazard function.
#' @param obsTimes Vector of the intended observation times.
#' @param windowHalf The amount of time before or after the intended obsTimes that a visit might take place.
#' @param trueEta A vector of true effect sizes for the fixed effects. 
#' @return A list with the elements:
#' \item{outcomeDat}{List containing the simulated left and right times.}
#' \item{trueBeta}{A 1 x q vector of the true genetic effects used to simulate the data.}
#' \item{xMat}{A n x p matrix of the simulated covariate data.}

icbvs_gen_dat <- function(n, gMat, noCausal, rho, ES, bhFunInv = function(x){x}, obsTimes = 1:5, windowHalf = 0.1, trueEta = c(-1, 0)){
  
  p = ncol(gMat)
  
  # generate vector of true effects
  trueBeta <- rep(0, p)	
  trueBeta[sample(c(1:p), noCausal)] <- ES
  #trueBeta[which.max(colSums(gMat))] <- ES
  
  
  dat <- sim_data_new(seed=0, n=n, p=length(trueBeta), bhFunInv = bhFunInv,
                      obsTimes = obsTimes, windowHalf = windowHalf, trueEta=trueEta, trueBeta = trueBeta,
                      gMat = gMat)
  
  # save the new data
  outcomeDat <- dat$outcomeDat
  
  return(list(outcomeDat = outcomeDat, trueBeta = trueBeta, xMat = dat$xMat))}