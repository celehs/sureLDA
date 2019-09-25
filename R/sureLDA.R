# sureLDA.R: Contains sureLDA function. See Ahuja et al. (2019), JAMIA (in revision) for details.
# Author: Yuri Ahuja
# Last Updated: 9/19/2019

# install.packages(c("Rcpp","RcppArmadillo","Matrix","flexmix","stats"), repos = "http://cran.us.r-project.org")
library(Rcpp)
library(RcppArmadillo)
require(Matrix)
require(flexmix)
require(stats)

sourceCpp("../src/lda_rcpp.cpp")
source("../R/FUN_PheNorm_Publish_ZH.R")

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
