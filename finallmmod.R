function(formulat, dataFrame=NULL) {
  formulat <- as.formula(formulat)
  m3 <- as.matrix(model.frame(terms(formulat), dataFrame))
  nn <- ncol(m3)
  Ymat <- m3[,1]
  Xmat <- m3[,2:nn]
  m4 <- as.data.frame(m3)
  
  Xmatrix <- as.matrix(Xmat)
  ymatrix <- as.matrix(Ymat)
  xmatrix <- as.matrix(cbind(1, Xmatrix))
  
  betas <- as.matrix(solve(t(xmatrix)%*%xmatrix)%*%t(xmatrix)%*%ymatrix)
  
  yhat <- xmatrix%*%betas
  
  resid <- as.matrix((as.numeric(ymatrix))-yhat)
  
  
  n <- nrow(m3) #number of observations (rows)
  k <- ncol(m3) #of variables (columns)
  bb <- nrow(betas)
  VCV <- (1/(n-k))*as.numeric(t(resid)%*%resid)*solve(t(xmatrix)%*%xmatrix)
  SE <- sqrt(diag(VCV))
  p_val <-2*pt(abs(betas/SE), df=n-k,lower.tail = FALSE)
  
  names <- colnames(m4[2:nn])
  
  #3commneting out, but alternative name <-deparse(substitute(head(m3[2:nn])))
  betaout <- as.data.frame(cbind(c("Intercept", names), betas, SE, p_val))
  names(betaout) <- c("coefficients", "estimate", "standard Error", "P_Value")
  return(betaout)}


