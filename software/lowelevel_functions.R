
# binary model ------------------------------------------------------------

#usage: binaryTrait ~ b1*cov1 + b2*cov2 + ... + bX*snp
#output: Chisq;p-value 
mb = sample(c(NA,2),size = 100,replace = T)
snp = sample(c(0,1,2),size = 100,replace = T)
covar = matrix(rnorm(200),nrow = 100,dimnames = list(NULL,c("cov1","cov2")))

quantitative_model = function (snp, mb, covar){
 
  y = as.integer(!is.na(mb))
  X = cbind(1,covar,snp)
  
  sigmoid <- function(z)
  {
    g <- 1/(1+exp(-z))
    return(g)
  }
  
}



# quantitative model ------------------------------------------------------------


#usage: binaryTrait ~ b1*cov1 + b2*cov2 + ... + bX*snp
#estimate betas and betas' SD
quantitative_model = function (snp, mb, covar){
  
  y = mb[!is.na(mb)]
  X = cbind(1,covar,snp)[!is.na(mb)]
  Nres = length(y) - ncol(X)

  betaHat = solve(t(X) %*% X) %*% t(X) %*% y
  MSE = sum((y - X %*% betaHat)^2)/Nres
  betaSD =  sqrt(diag(MSE * solve(t(X)%*% X)))
  beta_snp = betaHat[nrow(betaHat)]
  betaSD_snp = betaSD[nrow(betaHat)]
  q= (beta_snp]/betaSD[length(betaSD)]) ^ 2
  pvalue = pf(q,1,Nres,lower.tail = FALSE)

  return(c(length(y),q,beta_snp,betaSD_snp))
}
# the question if PF function. it's a distribution function for F discribution. In R it calls embeded C function, we need to find the same one in Java

