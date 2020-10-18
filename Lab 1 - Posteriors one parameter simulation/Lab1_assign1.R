#1. Bernoulli ... again.
#Let y1, ..., yn|theta ~ Bern(theta), and assume that you have obtained a sample with s = 5
#successes in n = 20 trials. Assume a Beta(alpha0, beta0) prior for theta and let alpha0 = beta0 = 2.
#
#(a) Draw random numbers from the posterior theta|y ~ Beta(alpha0 + s, beta0 + f), y =
# y1..yn and verify graphically that the posterior mean and standard deviation converges to 
# the true values as the number of random draws grows
# large.
set.seed(12345)
a0 = 2
b0 = 2
n = 20
s = 5
f = n-s
a1 = a0+s
b1 = b0+f


posteriorMean = (a0+s)/(a0+s+b0+f)


calcStdDev=function(alpha,beta){
  return(sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1))))
}

stdDev = calcStdDev(a1,b1)


drawBetaVal <- function(n, alpha, beta){
  return(rbeta(n,alpha,beta))
}

calcMSE <- function(n, mean, data){
  return(sqrt(1/(n-1)*sum((data-mean)^2)))
}

calcMean <- function(alpha,beta){
  return(alpha/(alpha+beta))  
}
  
nVector=c(seq(1,10000,5))
meanVector=c()
stdVector=c()
for(i in nVector){
  set.seed(12345)
  betaValues <- drawBetaVal(i, a1,b1)
  meanVector = c(meanVector, mean(betaValues))
  stdVector=c(stdVector, calcMSE(i,mean(betaValues),betaValues))
}

plot(nVector,meanVector,main="Mean converges to 0.29")
plot(nVector,stdVector, main="stdDev converges to 0.09")
print(stdDev)
print(posteriorMean)
#OK!

#Use simulation (nDraws = 10000) to compute the posterior probability
# Pr(theta > 0.3|y) and compare with the exact value [Hint: pbeta()].
set.seed(12345)
trueProb <- 1-pbeta(0.3,a1,b1)
set.seed(12345)
draws=rbeta(10000,a1,b1)
probHat = sum(draws>0.3)
prob = probHat/10000

phi=log(draws/(1-draws))
hist(phi, breaks=20, main="Distribution of the log-odds")

