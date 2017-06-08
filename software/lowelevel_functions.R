

# quantitative model ------------------------------------------------------------
# snp is a vector
# mb is a vector
# covar is a matrix of covariates 

quantitative_model = function (snp, mb, covar){
  y = mb[!is.na(mb)]
  X = cbind(1,covar,snp)[!is.na(mb),]
  Nres = length(y) - ncol(X)
  reSNP = round(X[,ncol(X)])
  reMAF = min( 
    (2*sum(reSNP == 0) + sum(reSNP == 1))/(2*length(reSNP)),
    (2*sum(reSNP == 2) + sum(reSNP == 1))/(2*length(reSNP))
  )

  betaHat = try(solve(t(X) %*% X) %*% t(X) %*% y,silent = T)
  if(class(betaHat) == "try-error") {
    return(c(length(y),reMAF,NA,NA,NA,NA))
  } else {
     MSE = sum((y - X %*% betaHat)^2)/Nres
  betaSD =  sqrt(diag(MSE * solve(t(X)%*% X)))
  beta_snp = betaHat[nrow(betaHat)]
  betaSD_snp = betaSD[nrow(betaHat)]
  q= (beta_snp/betaSD[length(betaSD)])^2
  pvalue = pf(q,1,Nres,lower.tail = FALSE)
  return(N=c(length(y),reMAF=reMAF,beta=beta_snp,betaSD=betaSD_snp,F = q,pvalue = pvalue))
  }
}
# the question if PF function. it's a distribution function for F discribution. In R it calls embeded C function, we need to find the same one in Java

# run analysis per bacteria for all snps (assuming in genotype file snps are in rows)
#bac is vector of bacteria
#genotypes is matrix with snps in rows

run_per_bac_parallel = function(genotypes,bac,covar,cluster){
  res = parApply(cl=cluster,genotypes,1,function(x){quantitative_model(x,bac,covar)})
  res = t(res)
  colnames(res) = c("Nnz","MAF","beta","betaSD","Fstat","pvalue")
  res
}


run_per_bac_single = function(genotypes,bac,covar){
  res = apply(genotypes,1,function(x){quantitative_model(x,bac,covar)})
  res = t(res)
  colnames(res) = c("Nnz","MAF","beta","betaSD","Fstat","pvalue")
  res
}



LRT.glm.fit <- function(glm1,glm0){
  df.null <- glm0$df.residual
  df.mod <- glm1$df.residual
  dev.null <- glm0$deviance
  dev.mod <- glm1$deviance
  dev.diff <- glm0$deviance - glm1$deviance
  p.value <- 1 - pchisq(dev.null - dev.mod, df.null - df.mod)
  output <- c(round(df.null), round(df.mod), dev.null, dev.mod, p.value)
  names(output) <- c("df.null", "df.mod", "dev.null", "dev.mod", "p.value")
  output
}

bin_model = function(snp,mb,covar){
  y = as.integer(!is.na(mb))
  X = cbind(1,covar,snp)
  
  glm1 = speedglm.wfit(y,X,family = binomial())
  return(summary(glm1)$coef[4,])
}

run_spearman = function(snp,mb,covar){
  y = mb[!is.na(mb)]
  covar  = covar[!is.na(mb),]
  snp = snp[!is.na(mb)]
  reSNP = round(snp)
  reMAF = min( 
    (2*sum(reSNP == 0) + sum(reSNP == 1))/(2*length(reSNP)),
    (2*sum(reSNP == 2) + sum(reSNP == 1))/(2*length(reSNP))
  )
  yres = resid(lm(y ~ covar))
  c1 = cor.test(yres,snp,method = "spearman")
  return(c(N=length(y),reMAF = reMAF,cor=c1$estimate,pvalue=c1$p.value))
}

