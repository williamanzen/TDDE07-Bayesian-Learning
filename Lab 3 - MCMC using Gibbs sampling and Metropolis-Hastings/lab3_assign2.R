ebay <- read.table("eBayNumberOfBidderData.dat", header=TRUE)

# a) Obtain the maximum likelihood estimator of Beta in the Poisson regression model
# for the eBay data [Hint: glm.R, don't forget that glm() adds its own intercept
# so don't input the covariate Const]. Which covariates are significant
ebayNew <- ebay[,-2]
ebayGLM <- glm(nBids~., family=poisson, data=ebayNew)

summary(ebayGLM)

## The covariates that are significant can be seen in the summary
## They are VerifyID, Sealed, MajBlem, LogBook and MinBidShare.

## b) Let's now do a Bayesian analysis of the Poisson regression.
## let the prior be Beta ~ N(0, 100*(t(X)%*%X)^-1) where X is the
## n * p covariate matrix. This is a commonly used prior which is called 
## Zellner's g-prior. Assume first that the posterior density is approximately
## multivariate normal:
## Beta|y ~ N( BetaTilde, Inv J(BetaTilde))
## where BetaTilde is the posterior mode and J(BetaTilde) is the negative Hessian
# att the posterior mode. BetaTilde and J(BetaTilde) can be obtained by numerical 
## optimization (optim.R) exactly like you already did for the logistic
## regression in Lab2 (but with the log posterior function replaced by the corresponding 
## one for the Poisson model, which you have to code up.).

library(mvtnorm)

#Setting up the prior
X <- as.matrix(ebay[,-1])
y <- as.vector(ebay[,1])
numberPara <- dim(X)[2]
mu <- rep(0,numberPara)
Sigma <- 100*solve(t(X)%*%X)
covNames <- names(ebay)[2:length(names(ebay))];


LogLikePoisson <- function(beta, y, x, mu, sigma){
  n<-length(y)
  xBeta <- beta%*%t(x)
  logLik <- sum(-log(factorial(y))+xBeta*y-exp(xBeta))
  prior = dmvnorm(beta, mean=mu, sigma=sigma, log=TRUE)
  return(logLik + prior)
}

set.seed(12345)
initVal <- rnorm(dim(X)[2]); 

##Optimizing loglikelihood with second derivatives (BFGS)
optimRes <-optim(initVal,LogLikePoisson,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

postMode <- optimRes$par
postCovariance <- -solve(optimRes$hessian) # Posterior covariance is -Inv(hessian) because we multiplied loglikelihood with -1 earlier
names(postMode) <- covNames # Naming the coefficient by covariates
approxStd <- sqrt(diag(postCovariance)) # Computing approximate standard deviations.
names(approxStd) <- covNames # Naming the coefficient by covariates
print('The posterior mode is:')
print(postMode)
print('The approximate posterior standard deviation is:')
print(approxStd)

# Approximation through optimization of the BetaVector gives us the posterior mode BetaTilde and the Hessian J(BetaTilde)

# (c) Now, let's simulate from the actual posterior of beta using the Metropolis algorithm and
# compare with the approximate results in b). Program a general
# function that uses the Metropolis algorithm to generate random draws from an
# arbitrary posterior density. In order to show that it is a general function for
# any model, I will denote the vector of model parameters by theta. Let the proposal
# density be the multivariate normal density mentioned in Lecture 8 (random walk Metropolis):
# thetaP | theta(i-1) ~ N(theta(i-1),c*Sigma)
# where Sigma = J(BetaTilde)^-1 obtained in b). The value c is a tuning parameter and
# should be an input to your Metropolis function. The user of your Metropolis function 
# should be able to supply her own posterior density function, not
# necessarily for the Poisson regression, and still be able to use your Metropolis
# function. This is not so straightforward, unless you have come across function
# objects in R and the triple dot (...) wildcard argument. I have posted a note
# (HowToCodeRWM.pdf) on the course web page that describes how to do this
# in R
# Now, use your new Metropolis function to sample from the posterior of beta in
# the Poisson regression for the eBay dataset. Assess MCMC convergence by
# graphical method


RWMSampler <- function(previousVal, c, hessian, MyFunction, ...){
  proposalVal <- rmvnorm(1,mean=previousVal, sigma=c*hessian)
  alpha=min(1, exp(MyFunction(proposalVal, ...)-MyFunction(previousVal,...))) ## subtraction because of log p(beta|.)
  
  u <- runif(1)
  if(u<alpha){
    return(proposalVal)
  } else {
    return(previousVal)
  }
}

## Setting up the RWM
nDraws <- 5000
sampledBetas <- matrix(0,nDraws,ncol(X))
sampledBetas[1,] <- initVal

set.seed(12345)
for(i in 1:nDraws){
  if(i<nDraws){
  sampledBetas[i+1,]<-RWMSampler(sampledBetas[i,],c=0.5,postCovariance, LogLikePoisson, y, X, mu, Sigma)
  }
}

avgAlpha <- dim(sampledBetas[!duplicated(sampledBetas),])[1]/dim(sampledBetas)[1]

## We compute the average acceptance rate for the RWM to 32.6% which is close to the recommended acceptance rate between 25-30%

iter=seq(1,nDraws,1)
for(i in 1:9){
  plot(iter, sampledBetas[,i], type="l", main=paste("Convergence plot for covariate", covNames[i]))
}

## As we see in the plots, each covariate converges to a value


# (d) Use the MCMC draws from c) to simulate from the predictive distribution of
# the number of bidders in a new auction with the characteristics below. Plot
# the predictive distribution. What is the probability of no bidders in this new
# auction

set.seed(12345)
xVector <- c(1,1,1,1,0,0,0,1,0.5) ## Observation given
trueBeta <- sampledBetas[1001:5000,]

predictiveDist <- exp(t(xVector%*%t(trueBeta)))
hist(predictiveDist)

predictedVals <- rpois(5000, predictiveDist)
hist(predictedVals, main="simulated predictive distribution")

probZero = sum(predictedVals==0)/length(predictedVals)

## The probability of no bidders is calculated to 36.02% given the specified characteristics.


