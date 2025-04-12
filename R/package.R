#' sureLDA: A Novel Multi-Disease Automated Phenotyping Method for the Electronic Health Record
#' 
#' Surrogate-guided ensemble Latent Dirichlet Allocation (sureLDA) is a label-free multidimensional phenotyping method. It first uses the PheNorm algorithm to initialize probabilities based on two surrogate features for each target disease, and then leverages these probabilities to guide the LDA topic model to generate phenotype-specific topics. Finally, it combines phenotype-feature counts with surrogates via clustering ensemble to yield final phenotype probabilities. 
#' 
#' @name sureLDA-package
#' @keywords package
#' @useDynLib sureLDA
#' @import Matrix Rcpp RcppArmadillo glmnet
#' @importFrom stats ecdf
#' @importFrom stats lm
#' @importFrom stats pnorm
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats sd
NULL
