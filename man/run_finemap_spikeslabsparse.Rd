% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_finemap_spikeslabsparse.R
\name{run_finemap_spikeslabsparse}
\alias{run_finemap_spikeslabsparse}
\title{Function for MCMC for interval-censored spike and slab using sparse setting.}
\usage{
run_finemap_spikeslabsparse(
  xMat,
  q,
  leftTimes,
  rightTimes,
  nKnots = 1,
  tposInd,
  obsInd,
  sigmaAprop,
  sigmaEprop,
  sigmaBprop,
  sigmaAprior,
  sigmaEprior,
  sigmaBprior,
  B = 10000,
  quant_r = NULL,
  addProb,
  removeProb,
  seed = 0,
  burnIn = 0.2,
  alphaPriorMean = NULL,
  etaPriorMean = NULL,
  checkpoint = FALSE
)
}
\arguments{
\item{xMat}{A n x p matrix of covariate data.}

\item{q}{Total number of SNPs.}

\item{leftTimes}{A n x 1 vector of observed right times of the event of interest.}

\item{rightTimes}{A n x 1 vector of observed left times of the event of interest.}

\item{nKnots}{Number of internal knots in cubic spline for estimating baseline hazard function.}

\item{tposInd}{A n x 1 vector indicating whether the event was observed before last follow-up.}

\item{obsInd}{A n x 1 vector indicating whether the event was observed after follow-up started.}

\item{sigmaAprop}{Hyperparameter for the standard deviation of the proposal distribution for the spline coefficients.}

\item{sigmaEprop}{Hyperparameter for the standard deviation of the proposal distribution for the fixed effect coefficients.}

\item{sigmaBprop}{Hyperparameter for the standard deviation of the proposal distribution for the genetic effect coefficients.}

\item{sigmaAprior}{Hyperparameter for the standard deviation of the prior distribution for the spline coefficients.}

\item{sigmaEprior}{Hyperparameter for the standard deviation of the prior distribution for the fixed effect coefficients.}

\item{sigmaBprior}{Hyperparameter for the standard deviation of the prior distribution for the genetic effect coefficients.}

\item{B}{Number of MCMC iterations (excluding burn-in). (Default set to 10000 iterations)}

\item{quant_r}{Quantiles of time to use in constructing the spline.}

\item{addProb}{Probability of adding another causal variant in selection scheme.}

\item{removeProb}{Probability of removing another causal variant in selection scheme.}

\item{seed}{Random seed. (Default set to 0)}

\item{burnIn}{Percentage of first beginning iterations removed from MCMC. (Default set to 0.2)}

\item{alphaPriorMean}{A nKnot + 1 vector for original spline coefficient estimates. (Default set to NULL)}

\item{checkpoint}{Boolean indicating whether to print results every 100 iterations. (Default set to FALSE)}

\item{A}{p x 1 vector for original fixed effect coefficient estimates. (Default set to NULL)}
}
\value{
A list with the elements:
\item{betaChain}{A (B\emph{burnIn) x q matrix containing all genetic effect coefficient estimates from all MCMC iterations, excluding burn-in.}
\item{alphaChain}{A (B}burnIn) x (nKnot + 2) matrix containing all spline coefficient estimates from all MCMC iterations, excluding burn-in.}
\item{etaChain}{A (B\emph{burnIn) x p matrix containing all fixed effect coefficient estimates from all MCMC iterations, excluding burn-in.}
\item{gammaChain}{A (B}burnIn) x q matrix containing all inclusion indicator from all MCMC iterations, excluding burn-in.}
\item{piChain}{A (B\emph{burnIn) x 1 vector containing all sparsity parameter estimates from all MCMC iterations, excluding burn-in.}
\item{accMat}{A (B}burnIn) x (q*2 + (nKnot + 2) + p + 1) matrix containing indicators for whether the proposal was accepted (= 1) or not (=0).}
}
\description{
Function for MCMC for interval-censored spike and slab using sparse setting.
}
\examples{
run_finemap_spikeslabsparse(xMat = cbind(xMat, gMat), q = 2, leftTimes = lt, rightTimes = rt, nKnots = nKnots,
 tposInd = tposInd, obsInd = obsInd,
sigmaAprior = sigmaAprior, sigmaEprior = sigmaEprior,
sigmaBprior = sigmaBprior, B=B,
sigmaAprop = sigmaAprop, sigmaEprop = sigmaEprop, sigmaBprop = sigmaBprop,
addProb = addProb, removeProb = removeProb, seed=0, burnIn=0.2,
alphaPriorMean=rep(1, nKnots + 2), etaPriorMean=rep(0, ncol(xMat)), checkpoint=TRUE)
}
