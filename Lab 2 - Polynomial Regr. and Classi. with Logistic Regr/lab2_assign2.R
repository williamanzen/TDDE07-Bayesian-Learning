# (a) Consider the logistic regression Pr(y = 1|x) = exp (xT b) / (1 + exp (xT b))
# where y is the binary variable with y = 1 if the woman works and y = 0 if she
# does not. x is a 8-dimensional vector containing the eight features (including
# a one for the constant term that models the intercept)
# The goal is to approximate the posterior distribution of the 8-dim parameter
# vector b with a multivariate normal distribution b|y, X ~ N(btild, Jy(btild)^???1)
# where btild is the posterior mode and J(btild) = ???d2ln p(b|y)/dbdb^T |b=btild is the
# observed Hessian evaluated at the posterior mode. Note that d2ln p(b|y)/dbdb^T is an 8x8 matrix with
# second derivatives on the diagonal and cross-derivatives d2ln p(b|y)/dbidbj
# on the offdiagonal. It is actually not hard to compute this derivative by hand, but don't
# worry, we will let the computer do it numerically for you. Now, both btild and J(btild)
# are computed by the optim function in R. See my code 
# zip where I have coded everything up for the spam prediction example (it also
# does probit regression, but that is not needed here). I want you to implement
# you own version of this. You can use my code as a template, but I want you
# to write your own file so that you understand every line of your code. Don't
# just copy my code. Use the prior b ~ N (0, T^2*I), with T = 10.
# Your report should include your code as well as numerical values for btild and Jy(btild)^-1
# for the WomenWork data. Compute an approximate 95% credible interval for the variable Nsmallchild.
# Would you say that this feature is an important determinant of the probability that a women works?
  # Hint To verify that your results are reasonable, you can compare to you get by
   # estimating the parameters using maximum likelihood:
# glmModel <- glm(Work ~ 0 + ., data = WomenWork, family = binomial).

#####USER INPUT#####
tau=10


#Reading data
data <- read.table("WomenWork.dat",header=TRUE)

# library(mvtnorm)
#convert data to vectors and matrixes
y <- as.vector(data[,1])
x <- as.matrix(data[,2:9])
covNames <- names(data)[2:length(names(data))];




numberPara <- dim(x)[2]
#Setting up the prior
mu <- rep(0,numberPara)
Sigma <- tau^2*diag(numberPara)

logLikelihood <- function(betaVec,y,x,mu,sigma){
  nPara <- length(betaVec)
  xBeta <- x%*%betaVec
  #defining loglikelihood
  logLik <- sum(y*xBeta-log(1+exp(xBeta)))
  #defining prior
  logPrior <- dmvnorm(betaVec, matrix(0,nPara,1), sigma, log=TRUE)
  return(logLik + logPrior)
  
}
set.seed(12345)
initVal <- rnorm(dim(x)[2]);

## Optimizing loglikelihood with second derivatives (BFGS)
optimRes <-optim(initVal,logLikelihood,gr=NULL,y,x,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

postMode <- optimRes$par
postCovariance <- -solve(optimRes$hessian) # Posterior covariance is -Inv(hessian) because we multiplied loglikelihood with -1 earlier
names(postMode) <- covNames # Naming the coefficient by covariates
approxStd <- sqrt(diag(postCovariance)) # Computing approximate standard deviations.
names(approxStd) <- covNames # Naming the coefficient by covariates
print('The posterior mode is:')
print(postMode)
print('The approximate posterior standard deviation is:')
print(approxStd)

## Compute distribution for feature NSmallChild
nSmallChildMode <- postMode["NSmallChild"]
nSmallChildStd <- approxStd["NSmallChild"]
set.seed(12345)
nSCdist <- rnorm(10000,mean=nSmallChildMode,sd=nSmallChildStd)
sortedNSCdist <- sort(nSCdist)
set.seed(12345)
credintNSC <- c(sortedNSCdist[(0.025*length(sortedNSCdist)+1)],sortedNSCdist[(0.975*length(sortedNSCdist))])
credIntQuant <- qnorm(c(0.025,0.975),mean=nSmallChildMode , sd = nSmallChildStd)
print(paste("The lower bound of the credible interval is:",round(credintNSC[1],6)," and the upper bound is:",round(credintNSC[2],6)))

glmModel <- glm(Work ~ 0+., data=data, family=binomial)
glmModel$coefficients
# Since the variable has a credible interval which is negative, the feature contributes negatively towards the response variable
# making it go towards 0 i.e. women not working. This seems reasonable since a mother with a small child is likely to be home
# taking care of the child instead of working. This is confirmed by seeing that the glmModel which estimates through ML
# also confirms a negative coefficient, which indicates that our code is reasonable.

# (b) Write a function that simulates from the predictive distribution of the response 
# variable in a logistic regression. Use your normal approximation from 2(a). 
# Use that function to simulate and plot the predictive distribution for the Work 
# variable for a 40 year old woman, with two children (3 and 9 years old), 8 years 
# of education, 10 years of experience. and a husband with an income of 10. 
# Hints: The R package mvtnorm will again be handy. Remember my discussion 
# on how Bayesian prediction can be done by simulation. 

logPred <- function(pred){
  return(exp(pred)/(1+exp(pred)))
}


predLogReg <- function(obs,mean,sigma,nDraw){
  
  bDraws <- rmvnorm(n=nDraw, mean=mean, sigma=sigma)
  predY = bDraws %*% obs
  return(logPred(predY))
}

woman=c(1, 10, 8, 10, (10/10)^2, 40, 1, 1)
predictions <- predLogReg(woman,postMode,-solve(optimRes$hessian),10000)
logDistr=c()
for( i in predictions){
  logDistr=c(logDistr,rbinom(1,1,i))
}
barplot(table(logDistr),xlab="probabiltiy", ylab="density", main="Density of predicted probabilites")

# We see in the plot that the mode for the density function is around 0.2 which tells us that there is approx. 20% probability
# that the woman is working.

# (c) Now, consider 10 women which all have the same features as the woman in 2(b).
# Rewrite your function and plot the predictive distribution for the number of
# women, out of these 10, that are working. [Hint: Which distribution can be
# described as a sum of Bernoulli random variables?]


predProb <- function(obs,mean,sigma,nDraw, n){
  
  multiplePred = c()
  for(i in 1:nDraw){
    bDraws <- predLogReg(obs,mean,sigma,1)
    multiplePred=c(multiplePred, rbinom(1,n,bDraws))
  }
  
  barplot(table(multiplePred),main=paste("Distribution for predictions made on ", n," women"))
  
}

predProb(woman, postMode, -solve(optimRes$hessian),10000,10)

## 
