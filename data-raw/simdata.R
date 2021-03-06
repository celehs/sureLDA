library(bindata) ##For Multiple Y

Simulate<-function(N=10000,K=10,M=50){
  ##Generate Paramters
  Alpha <- rep(0.3,K)
  Beta <- rep(0.25,K)
  Lambda <- matrix(rep(0.2,M*K-2*K))
  
  mu_positive <- rbind(
    t(sapply(1:5,function(o) c(4.5,2.1,seq(1.9,1.7,length=M-2)))),
    t(sapply(6:K,function(o) c(4.5,2.1,seq(1.9,1.7,length=M-2)))) )
  mu_negative <- rbind(
    t(sapply(1:5,function(o) c(0.51,0.81,seq(0.85,0.82,length=M-2)))),
    t(sapply(6:K,function(o) c(0.51,0.81,seq(0.85,0.81,length=M-2)))) )
  
  ##Prior Probability of Y
  Pi <- c(0.4,0.2,0.15,0.13,0.1,0.03,0.03,0.03,0.03,0.03)
  
  ##Correlation of Multiple Y
  sigma <- matrix(0, K, K)
  sigma[1,2] <- 0.5
  sigma[2,1] <- 0.5
  diag(sigma)<-1
  ##Generate Y
  Y <- rmvbin(N, margprob=Pi, sigma=sigma)
  Y_table = apply(Y,1,function(x){t(2^(0:(K-1)))%*%t(t(x))})
  table(Y_table)
  length(table(Y_table))
  
  ##Generate X
  H <- rpois(N,lambda = 2)
  random <- t(apply(Y,1,function(y){
    t(sapply(1:K, function(i){
      yi = y[i]
      if(yi==1) rgamma(M, scale=mu_positive[i,], shape=1)
      else rgamma(M, scale=mu_negative[i,], shape=1)
    }))
  }))
  
  logX <- cbind(log(H+1)%*%t(Alpha),log(H+1)%*%t(Beta),log(H+1)%*%t(Lambda)) + log(random)
  X <- floor(exp(logX))
  ICD <- X[,1:K] 
  X <- cbind(H,X)
  colnames(X) <- c('H',paste('ICD',1:K,sep=''),paste('NLP',1:K,sep=''),paste('V',1:(M*K-2*K),sep=''))
  ##Filter is defined as ICD >= 1
  filter <- I(ICD>=1)*1
  
  i = 1
  for(i in 1:K){
    print(table(filter[,i],Y[,i],dnn=c(paste('filter',i,sep=''),paste('Diseases',i,sep='')))/N)
  }
  list(X=X,Y=Y,filter=filter)
}

set.seed(123)

N=500  ##Number of patients
K=10     ##Number of diseases
M=10     ##Number of features for each diseases, two of them are ICD and NLP.
Sim <- Simulate(N,K,M)
X <- Sim$X[,-1]
Y <- Sim$Y
filter <-Sim$filter 
diseases <- seq(K)
HU <- Sim$X[,1]
ICD <- Sim$X[,2:11]
NLP <- Sim$X[,12:21]
nPatients = N

simdata <- list(Y = Y, X = X, ICD = ICD, NLP = NLP, HU = HU, filter = filter)
str(simdata)

usethis::use_data(simdata, overwrite = TRUE)
