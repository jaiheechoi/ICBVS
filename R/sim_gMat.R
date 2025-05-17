#' Generate simulated genetic data set
#'
#' @param n Total number of observations.
#' @param q Total number of SNPs.
#' @param rho Common pairwise correlation parameter.
#' @param maxMAF Maximum minor allele frequency for mean parameter vector.
#' @returns A n x q matrix of genetic data.
#' @examples
#' sim_gMat(1000, 100, 0.1, .05)
sim_gMat <- function(n,q,rho, mafVec){

  # Construct a binary correlation matrix
  cmat <- matrix(rho, nrow = q, ncol = q)
  diag(cmat) <- 1

  # create mean parameter vector using from Unif(0, maxMAF) distribution
  meanparam <- runif(q, .01, maxMAF)

  # create genetic data from two binary RVs
  x <- bindata::rmvbin(n, margprob = meanparam, bincorr = cmat)  + bindata::rmvbin(n, margprob = meanparam, bincorr = cmat)
  return(x)
}
