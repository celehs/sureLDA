# sureLDA.R: Contains sureLDA function. See Ahuja et al. (2019), JAMIA (in revision) for details.
# Author: Yuri Ahuja
# Last Updated: 9/19/2019

# install.packages(c("Rcpp","RcppArmadillo","Matrix","flexmix","stats"), repos = "http://cran.us.r-project.org")
# library(Rcpp)
# library(RcppArmadillo)
# library(Matrix)
# library(flexmix)
# library(stats)

# sourceCpp("../src/sureLDA.cpp")
# source("../R/FUN_PheNorm_Publish_ZH.R")

# INPUT:
# X = nPatients x nFeatures matrix of feature counts
# weights = nPhenotypes x nFeatures matrix of phenotype-specific feature weights
# ICD = nPatients x nPhenotypes matrix of main ICD surrogate counts
# NLP = nPatients x nPhenotypes matrix of main NLP surrogate counts
# HU = nPatients-dimensional hospital utilization vector
# filter = nPatients x nPhenotypes binary matrix indicating filter-positives
# nEmpty = Number of 'empty' topics to include in LDA step
# alpha, beta = LDA hyperparameters
# burnin = number of burnin Gibbs iterations; ITER = number of subsequent iterations for inference

# OUTPUT:
# scores = nPatients x nPhenotypes matrix of sureLDA phenotype scores
# probs = nPatients x nPhenotypes matrix of posterior probabilities (post-clustering of sureLDA scores)

#' ...
#' @param X nPatients x nFeatures matrix of feature counts
#' @param weight nPhenotypes x nFeatures matrix of phenotype-specific feature weights
#' @param ICD nPatients x nPhenotypes matrix of main ICD surrogate counts
#' @param NLP nPatients x nPhenotypes matrix of main NLP surrogate counts
#' @param HU nPatients-dimensional hospital utilization vector
#' @param filter nPatients x nPhenotypes binary matrix indicating filter-positives
#' @param nEmpty Number of 'empty' topics to include in LDA step
#' @param alpha LDA hyperparameters
#' @param beta LDA hyperparameters
#' @param burnin number of burnin Gibbs iterations
#' @param ITER number of subsequent iterations for inference
#' @export
sureLDA <- function(X,weight,ICD,NLP,HU,filter,nEmpty=20,alpha=1,beta=1,burnin=50,ITER=150){
  knowndiseases = ncol(ICD)
  D = knowndiseases + nEmpty
  W = ncol(X)
  N = nrow(X)
  
  
  ## PheNorm (Step 1) ##
  print("Starting PheNorm")
  prior <- sapply(1:knowndiseases, function(i){
    print(paste("On disease",i))
    
    mat = Matrix(data=cbind(log(ICD[,i]+1), log(NLP[,i]+1), log(ICD[,i]+NLP[,i]+1)), sparse=TRUE)
    note = Matrix(data=log(HU+1), sparse=TRUE)
    keep = which(filter[,i]==1)
    data = cbind(keep,note[keep],mat[keep,])
    
    fit.phenorm = PheNorm.Prob(c(3:ncol(data)), 2, data, nm.X=NULL, corrupt.rate=0.3, train.size=10000)
    score = rep(0,dim(mat)[1])
    score[keep] = as.vector(fit.phenorm$probs)
    score
  })
  
  
  ## Guided LDA (Step 2) ##
  Add_probs = matrix(0,ncol=(D-knowndiseases),nrow=N)
  prior = t(cbind(prior,Add_probs)) ##MAP_initial_probs is a matrix of N rows, 10
  weight = t(cbind(weight,matrix(1,W,nEmpty)))
  
  xx=data.frame("V1"=rep(1:N,rep(W,N)),"variable"=rep(1:W,N),"value"=as.vector(t(as.matrix(X))))
  xx = xx[xx$value>0,]
  d = rep(xx$V1,xx$value) - 1
  w = rep(xx$variable,xx$value) - 1
  z = rep(0,length(d))
  
  res1 = lda_rcpp(d,w,z,weight,prior,alpha,beta,D,knowndiseases,burnin,ITER)[1:knowndiseases,]
  res2 = lda_rcpp(d,w,z,weight,prior,alpha,beta,D,knowndiseases,burnin,ITER)[1:knowndiseases,]
  res3 = lda_rcpp(d,w,z,weight,prior,alpha,beta,D,knowndiseases,burnin,ITER)[1:knowndiseases,]
  res4 = lda_rcpp(d,w,z,weight,prior,alpha,beta,D,knowndiseases,burnin,ITER)[1:knowndiseases,]
  
  LDA_Ndk_predicted = t(res1 + res2 + res3 + res4) / (4*ITER)

  
  ## Clustering of surrogates with sureLDA score (Step 3) ##
  print("Starting final clustering step")
  posterior <- sapply(1:knowndiseases, function(i){
    print(paste("On disease",i))
    
    mat = Matrix(data=cbind(log(ICD[,i]+1), log(NLP[,i]+1), log(LDA_Ndk_predicted[,i]+1)), sparse=TRUE)
    note = Matrix(data=log(HU+1), sparse=TRUE)
    keep = which(filter[,i]==1)
    data = cbind(keep,note[keep],mat[keep,])
    
    fit.phenorm = PheNorm.Prob(c(3:ncol(data)), 2, data, nm.X=NULL, corrupt.rate=0.3, train.size=10000)
    score = rep(0,dim(mat)[1])
    score[keep] = as.vector(fit.phenorm$probs)
    score
  })
  
  
  return(list("scores"=LDA_Ndk_predicted, "probs"=posterior))
}

## dat: all data columns need to be log-transformed and need column names; ##
## nm.logS.ori is the name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
## nm.utl: is the name of healthcare utlization (e.g. note count, encounter_num etc)
## nm.X: additional features other than the main ICD and NLP
PheNorm.Prob = function(nm.logS.ori,nm.utl,dat, nm.X=NULL,corrupt.rate=0.3,train.size=10000){
  dat = as.matrix(dat)
  S.ori = dat[,nm.logS.ori,drop=F]; utl = dat[,nm.utl]
  a.hat = apply(S.ori, 2, function(S){findMagicNumber(S,utl)$coef})
  S.norm = S.ori - VTM(a.hat,nrow(dat))*utl
  if(!is.null(nm.X)){
    X = as.matrix(dat[,nm.X])
    SX.norm = cbind(S.norm,X,utl)
    id = sample(1:nrow(dat), train.size, replace=T)
    SX.norm.corrupt = apply(SX.norm[id,],2,function(x){ifelse(rbinom(length(id),1,corrupt.rate),mean(x),x)})
    b.all = apply(S.norm, 2, function(ss){lm(ss[id]~SX.norm.corrupt-1)$coef})
    b.all[is.na(b.all)] = 0
    S.norm = as.matrix(SX.norm)%*%b.all
    b.all = b.all[-dim(b.all)[1],]
  }
  else{
    b.all = NULL
  }
  if(length(nm.logS.ori)>1){
    postprob = apply(S.norm,2,function(x){fit = normalmixEM2comp2(x, lambda=0.5, mu=quantile(x,probs=c(1/3,2/3)), sigsqrd=sd(S.norm)/2);fit$posterior[,2]})
    list("probs"=rowMeans(postprob,na.rm = T), "betas"=b.all)
    
  }else{
    fit = normalmixEM2comp2(unlist(S.norm), lambda=0.5, mu=quantile(S.norm,probs=c(1/3,2/3)), sigsqrd=sd(S.norm)/2)
    list("probs"=fit$posterior[,2], "betas"=b.all)
  }
}

## dat: all data columns need to be log-transformed and need column names; ##
## nm.logS.ori is the name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
## nm.utl: is the name of healthcare utlization (e.g. note count, encounter_num etc)
## nm.X: additional features other than the main ICD and NLP
PheNorm = function(nm.logS.ori,nm.utl,dat, nm.X=NULL,corrupt.rate=0.3,train.size=100000){
  dat = as.matrix(dat)
  S.ori = dat[,nm.logS.ori,drop=F]; utl = dat[,nm.utl]
  a.hat = apply(as.matrix(S.ori), 2, function(S){findMagicNumber(S,utl)$coef})
  S.norm = S.ori - VTM(a.hat,nrow(dat))*utl
  if(!is.null(nm.X)){
    X = as.matrix(dat[,nm.X])
    SX.norm = cbind(S.norm,X,utl)
    id = sample(1:nrow(dat), train.size, replace=T)
    SX.norm.corrupt = apply(SX.norm[id,],2,function(x){ifelse(rbinom(length(id),1,corrupt.rate),mean(x),x)})
    b.all = apply(S.norm, 2, function(ss){lm(ss[id]~SX.norm.corrupt-1)$coef})
    b.all[is.na(b.all)] = 0
    S.norm = as.matrix(SX.norm)%*%b.all
    b.all = b.all[-dim(b.all)[1],]
  }
  else{
    b.all = NULL
  }
  if(length(nm.logS.ori)>1){
    postprob = apply(S.norm,2,function(x){fit = normalmixEM2comp2(x, lambda=0.5, mu=quantile(x,probs=c(1/3,2/3)), sigsqrd=1);fit$posterior[,2]})
    keep = apply(postprob,1,function(x){if(sum(x>0.5)>=2) x[which(x<0.5)]=NA else x[which(x>0.5)]=NA; x})
    keep = as.matrix(1*(keep>=0)); if(nrow(keep)!=nrow(dat)){keep=t(keep)}
    list("scores"=rowMeans(S.norm*keep,na.rm = T), "betas"=b.all)    
  }else{
    list("scores"=unlist(S.norm), "betas"=b.all)
  }
}

normalmixEM2comp2 <- function (x, lambda, mu, sigsqrd, eps = 1e-08, maxit = 1000, verb = FALSE) {
  arbvar <- (length(sigsqrd) == 2)
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigsqrd1 <- sigsqrd[1]
  sigsqrd2 <- sigsqrd[arbvar + 1]
  mx <- mean(x)
  const <- length(x) * 0.918938533204673
  dl <- 1 + eps
  iter <- 0
  ll <- rep(0, maxit + 1)
  a1 <- (x - mu1)^2
  b1 <- (lambda/sqrt(sigsqrd1)) * exp(-a1/2/sigsqrd1)
  a2 <- (x - mu2)^2
  b2 <- ((1 - lambda)/sqrt(sigsqrd2)) * exp(-a2/2/sigsqrd2)
  l <- sum(log(b1 + b2))
  while (dl > eps && iter < maxit) {
    iter <- iter + 1
    ll[iter] <- l
    postprobs <- b1/(b1 + b2)
    lambda <- mean(postprobs)
    mu1 <- mean(postprobs * x)/lambda
    mu2 <- (mx - lambda * mu1)/(1 - lambda)
    if (arbvar) {
      sigsqrd1 <- mean(postprobs * a1)/lambda
      sigsqrd2 <- mean((1 - postprobs) * a2)/(1 - lambda)
    }
    else {
      sigsqrd1 <- sigsqrd2 <- mean(postprobs * a1 + (1 - 
                                                       postprobs) * a2)
    }
    a1 <- (x - mu1)^2
    b1 <- (lambda/sqrt(sigsqrd1)) * exp(-a1/2/sigsqrd1)
    a2 <- (x - mu2)^2
    b2 <- ((1 - lambda)/sqrt(sigsqrd2)) * exp(-a2/2/sigsqrd2)
    oldl <- l
    l <- sum(log(b1 + b2))
    dl <- l - oldl
    if (verb) {
      cat("iteration =", iter, " log-lik diff =", dl, " log-lik =", 
          l - const, "\n")
    }
  }
  ##cat("number of iterations=", iter, "\n")
  iter <- iter + 1
  ll[iter] <- l
  postprobs <- cbind(postprobs, 1 - postprobs)
  colnames(postprobs) <- c(paste("comp", ".", 1:2, sep = ""))
  out <- list(x = x, lambda = c(lambda, 1 - lambda), mu = c(mu1, 
                                                            mu2), sigma = sqrt(c(sigsqrd1, sigsqrd2)[1:(1 + arbvar)]), 
              loglik = l - const, posterior = postprobs, all.loglik = ll[1:iter] - 
                const, restarts = 0, ft = "normalmixEM")
  class(out) <- "mixEM"
  out
}

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


findMagicNumber = function(surrogate, log_note_count, n.boot=10) {
  a.values = rep(1,n.boot)
  err.values = rep(Inf,n.boot)
  for(k in 1:n.boot) {
    idx = sample(1:length(log_note_count), replace=n.boot>1) # bootstrap to make the process more robust
    coefs = seq(0,1.2,0.05)
    err = rep(0, length(coefs))
    err.prev = Inf
    for (j in 1:length(coefs)) {
      score = surrogate[idx] - coefs[j] * log_note_count[idx]
      fit = normalmixEM2comp2(score, lambda=0.5, mu=quantile(score,probs=c(1/3,2/3)), sigsqrd=1)
      # the following approximates \int|Fn(x)-Fmix(x)|dx
      Fn = ecdf(score) # empirical cdf
      Fmix = function(x) {fit$lambda[1]*pnorm(x,fit$mu[1],fit$sigma) + fit$lambda[2]*pnorm(x,fit$mu[2],fit$sigma)}
      id.small = which.min(fit$mu); id.large = which.max(fit$mu)
      fit.range = seq(fit$mu[id.small]-5*fit$sigma, fit$mu[id.large]+5*fit$sigma, length.out=1000)
      x.step = (10*fit$sigma + fit$mu[id.large] - fit$mu[id.small])/length(fit.range)
      fit.err = function(x){abs(Fn(x)-Fmix(x))}
      err[j] = sum(fit.err(fit.range))*x.step
      if (err[j] > err.prev)
        break
      else
        err.prev = err[j]
    }
    a.values[k] = coefs[j-1]
    err.values[k] = err[j-1]
  }
  return(list("coef"=mean(a.values), "error"=mean(err.values)))
}

