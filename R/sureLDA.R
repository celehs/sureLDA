# sureLDA.R: Contains sureLDA function. See Ahuja et al. (2020), JAMIA for details.
# Author: Yuri Ahuja
# Last Updated: 9/12/2020

library(Rcpp)
library(RcppArmadillo)
require(Matrix)
require(flexmix)
require(stats)
library(foreach)
library(doParallel)
library(glmnet)

# source("../sureLDA/PheNorm.R")
# source("../sureLDA/MAP.R")
# sourceCpp("../sureLDA/lda_rcpp.cpp")


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
# prior = 'PheNorm', 'MAP', or nPatients x nPhenotypes matrix of prior probabilities
# weight = 'beta', 'uniform', or nPhenotypes x nFeatures matrix of feature weights
# labels (optional) = NA when missing, non-NA for observed

# OUTPUT:
# scores = nPatients x nPhenotypes matrix of sureLDA phenotype scores
# probs = nPatients x nPhenotypes matrix of posterior probabilities (post-clustering of sureLDA scores)

logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}

sureLDA <- function(X,ICD,NLP,HU,filter,prior='PheNorm',weight='beta',nEmpty=20,alpha=100,beta=100,burnin=50,ITER=150,phi=NULL,labeled=NULL){
	knowndiseases = ncol(ICD)
	D = knowndiseases + nEmpty
	W = ncol(X)
	N = nrow(X)
	
	
	## PheNorm/MAP (Step 1) ##
	
	if (typeof(prior) != 'character'){
	  print('Prior supplied')
	}
	else if (prior == 'PheNorm' & (typeof(weight) != 'character' | weight == 'uniform')){
	  print("Starting PheNorm")
	  prior <- sapply(1:knowndiseases, function(i){
	    print(paste("Predicting disease",i))
	    
	    mat = Matrix(data=cbind(log(ICD[,i]+1), log(NLP[,i]+1), log(ICD[,i]+NLP[,i]+1)), sparse=TRUE)
	    note = Matrix(data=log(HU+1), sparse=TRUE)
	    filterpos = which(filter[,i]==1)
	    data = cbind(filterpos,note[filterpos],mat[filterpos,],log(X[filterpos,]+1))
	    
	    fit.phenorm = PheNorm.Prob(3:5, 2, data, nm.X=6:ncol(data), corrupt.rate=0.3, train.size=10000)
	    score = rep(0,N)
	    score[filterpos] = as.vector(fit.phenorm$probs)
	    score
	  })
	}
	else if (prior == 'PheNorm' & weight == 'beta'){
	  print("Starting PheNorm")
	  prior <- matrix(nrow=N,ncol=knowndiseases)
	  weight <- matrix(nrow=knowndiseases,ncol=W)
	  for (i in 1:knowndiseases){
	    print(paste("Predicting disease",i))
	    
	    mat = Matrix(data=cbind(log(ICD[,i]+1), log(NLP[,i]+1), log(ICD[,i]+NLP[,i]+1)), sparse=TRUE)
	    note = Matrix(data=log(HU+1), sparse=TRUE)
	    filterpos = which(filter[,i]==1)
	    data = cbind(filterpos,note[filterpos],mat[filterpos,],log(X[filterpos,]+1))
	    
	    fit.phenorm = PheNorm.Prob(3:5, 2, data, nm.X=6:ncol(data), corrupt.rate=0.3, train.size=10000)
	    score = rep(0,N)
	    score[filterpos] = as.vector(fit.phenorm$probs)
	    prior[,i] <- score
	    weight[i,] <- as.vector(fit.phenorm$betas[,3])
	  }
	  print("Finishing PheNorm")
	}
	else if (prior == 'MAP'){
	  print("Starting MAP")
	  prior <- sapply(1:knowndiseases, function(i){
	    print(paste("Predicting prior",i))
	    
	    filterpos = which(filter[,i]==1)
	    mat = Matrix(data=cbind(ICD[filterpos,i],NLP[filterpos,i]), sparse=TRUE)
	    colnames(mat) = c('ICD','NLP')
	    note = Matrix(HU[filterpos], sparse=TRUE)
	    
	    score = rep(0,N)
	    score[filterpos] = as.vector(MAP(mat=mat, note=note)$scores)
	    score
	  })
	  
	  if (anyNA(prior)){
	    stop('MAP output has NAs')
	  }
	  
	  if (typeof(weight)=='character' & weight == 'beta'){
	    weight <- t(sapply(1:knowndiseases, function(i){
	      print(paste("Predicting weight",i))
	      filterpos = which(filter[,i]==1)
	      
	      SX.norm.corrupt <- apply(X[filterpos,],2,function(x){
	        ifelse(rbinom(length(filterpos),1,0.3), mean(x), x)
	      })
	      prior_bounded <- pmax(1e-5,pmin(1-1e-5,prior[filterpos,i]))
	      logit_prior <- logit(prior_bounded)
	      reg.weights <- prior_bounded * (1-prior_bounded)
	      
	      coef(glmnet(SX.norm.corrupt,logit_prior,weights=reg.weights,intercept=FALSE), s=0)[-1]
	    }))
	  }
	  print("Finishing MAP")
	}
	
	if (typeof(weight) == 'character' & weight == 'uniform'){
	  weight = matrix(100,nrow=knowndiseases,ncol=W)
	}
	
	weight[weight<0] <- 0
	weight <- round(100*weight/mean(weight))
	weight <- rbind(weight,matrix(alpha,nrow=nEmpty,ncol=W))
	
	if (!is.null(labeled)){
		prior[!is.na(labeled)] = labeled[!is.na(labeled)]	
	}

	
	
	## Guided LDA (Step 2) ##
	if (is.null(phi)){
	  print('Starting Guided LDA Step')
	  
	  Add_probs = matrix(0,ncol=(D-knowndiseases),nrow=N)
	  priorLDA = t(cbind(prior,Add_probs)) ##MAP_initial_probs is a matrix of N rows, 10
	  
	  xx=data.frame("V1"=rep(1:N,rep(W,N)),"variable"=rep(1:W,N),"value"=as.vector(t(as.matrix(X))))
	  xx = xx[xx$value>0,]
	  d = rep(xx$V1,xx$value) - 1
	  w = rep(xx$variable,xx$value) - 1
	  z = rep(0,length(d))
	  
	  print('Starting Gibbs Sampling')
	  
	  res = foreach(it=1:3) %do% {
	    print(paste('On iteration',it))
	    lda_rcpp(d,w,z,weight,priorLDA,alpha,beta,D,knowndiseases,burnin,ITER)[1:knowndiseases,]
	  }
	  
	  print('Finishing Gibbs Sampling')
	  
	  if (knowndiseases == 1){
	    LDA_Ndk_predicted = as.matrix(res[[1]] + res[[2]] + res[[3]]) / (3*alpha*ITER)
	  }
	  else{
	    LDA_Ndk_predicted = t(res[[1]] + res[[2]] + res[[3]]) / (3*alpha*ITER)
	  }
	}
	else{
	  print("Inferring theta given provided phi")
	  
	  LDA_Ndk_predicted <- foreach(i=1:N, .combine=rbind) %dopar% {
	    prior_i <- c(prior[i,],rep(0,nEmpty)); prior_i <- prior_i/sum(prior_i)
	    post_i <- t(prior_i * phi); post_i <- post_i / rowSums(post_i)
	    z_i <- c((t(post_i)*weight) %*% X[i,])
	    old <- rep(0,D)
	    while (any(z_i-old >= 0.1)){
	      old <- z_i
	      prior_i <- c(prior[i,],rep(1,nEmpty)) + z_i; prior_i <- prior_i/sum(prior_i)
	      post_i <- t(prior_i * phi); post_i <- post_i / rowSums(post_i)
	      z_i <- c((t(post_i)*weight) %*% X[i,])
	    }
	    z_i
	  }
	}

	
	## Clustering of surrogates with sureLDA score (Step 3) ##
	print("Starting final clustering step")
	posterior <- sapply(1:knowndiseases, function(i){
	  print(paste("Predicting posterior",i))
	  
		mat = Matrix(data=log(LDA_Ndk_predicted[,i]+1), sparse=TRUE)
		note = Matrix(data=log(HU+1), sparse=TRUE)
		keep = which(filter[,i]==1)
		data = cbind(keep,note[keep],mat[keep,])
		
		fit.phenorm = PheNorm.Prob(c(3:ncol(data)), 2, data, nm.X=NULL, corrupt.rate=0.3, train.size=10000)
		score = rep(0,dim(mat)[1])
		score[keep] = as.vector(fit.phenorm$probs)
		score
	})
	
	
	return(list("scores"=LDA_Ndk_predicted, "probs"=posterior, "ensemble"=(prior+posterior)/2,
	            "prior"=prior, "weights"=weight[1:knowndiseases,]))
}

