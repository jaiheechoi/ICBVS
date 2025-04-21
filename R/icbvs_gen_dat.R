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


