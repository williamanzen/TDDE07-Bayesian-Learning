# The dataset TempLinkoping.txt contains daily average temperatures (in Celcius degrees)
# at Malmslatt, Linkoping over the course of the year 2018. The response
# variable is temp and the covariate is time =  (the number of days since beginning of year)/365

# The task is to perform a Bayesian analysis of a quadratic regression
# temp = beta0 + beta1 * time + beta2 * time2 + eps, eps ~iid N(0, sigmasq).

# a) Determining the prior distribution of the model parameters. Use the conjugate
# prior for the linear regression model. Your task is to set the prior hyperparameters my0, omega0, v0 and sigmasq0
# to sensible values. Start with my0 = (-10, 100, -100)T
# omega0 = 0.01 * I3, v0 = 4 and sigmasq0 = 1  Check if this prior agrees with your prior
# opinions by simulating draws from the joint prior of all parameters and for
# every draw compute the regression curve. This gives a collection of regression
# curves, one for each draw from the prior. Do the collection of curves look reasonable?
# If not, change the prior hyperparameters until the collection of prior
# regression curves agrees with your prior beliefs about the regression curve.
# Hint: the  R package mvtnorm will be handy. And use your Inv-chisqsimulator from Lab 1]

temp = read.table("TempLinkoping.txt", header=TRUE)

# install.packages("mvtnorm")
  library(mvtnorm)
beta0 = -10
beta1 = 100
beta2 = -100
my0 = c(beta0, beta1, beta2)
v0=365
sigmasq0=0.5
omega0=0.5*diag(3)
invomega0= solve(omega0)

calcRegr <- function(b,rows,x){
  return(b[rows,1]+b[rows,2]*x+b[rows,3]*x^2)
}

drawBeta <- function(my, sigmasq,invomega){
  return(rmvnorm(1,mean=my,sigma=sigmasq*invomega))
}


nDraws = 1000
drawX=rchisq(nDraws,df=v0)
sigmasq=(v0)*sigmasq0/drawX
bmatrix <- matrix(0, nDraws, 3)

plot.new()
plot.window(xlim=c(0,1),ylim=c(-50,50))
axis(side=1)
axis(side=2)
title(main="simulated temps for different times")
set.seed(12345)

for(i in 1:nDraws ){
  bmatrix[i,]=drawBeta(my0,sigmasq[i],invomega0)
  lines(temp$time, bmatrix[i,1]+bmatrix[i,2]*temp$time+bmatrix[i,3]*temp$time^2,col=rgb(0,0,0,0.2))
}

#We changed beta0 to 0 in order to lift the collection of prior regression curves to a slightly higher value as we
#seemed fit. It now looks reasonable as it is low (around 0) in the beginning of the year, then the temperatures
# rise during the summer and goes back down again to around 0 during the winter.

# (b) Write a program that simulates from the joint posterior distribution of beta0,
# beta1,beta2 and sigma2 .Plot the marginal posteriors for each parameter as a histogram.
# Also produce another figure with a scatter plot of the temperature data and
# overlay a curve for the posterior median of the regression function f(time) =
# beta0+beta1 *time+beta2 *time^2
#  computed for every value of time. Also overlay curves
# for the lower 2.5% and upper 97.5% posterior credible interval for f(time).
# That is, compute the 95% equal tail posterior probability intervals for every
# value of time and then connect the lower and upper limits of the interval by
# curves. Does the interval bands contain most of the data points? Should they?

#Defining parameters for the posterior distribution
vn = v0+length(temp$temp)
X=cbind(1,temp$time,temp$time^2)
Y=temp$temp
betaHat <- solve(t(X)%*%X)%*%t(X)%*%Y
omegan <- t(X)%*%X+omega0

myn=solve(t(X)%*%X+omega0)%*%(t(X)%*%X%*%betaHat+omega0%*%my0)
sigmasqn=(v0*sigmasq0+(t(Y)%*%Y+t(my0)%*%omega0%*%my0-t(myn)%*%omegan%*%myn))/vn

#simulate the joint posterior

sigmasqPost=(vn)*c(sigmasqn)/drawX
bmatrixPost <- matrix(0, nDraws, 3)

# plot.new()
# plot.window(xlim=c(0,1),ylim=c(-50,50))
# axis(side=1)
# axis(side=2)

invomegan = solve(omegan)



for(i in 1:nDraws){
  bmatrixPost[i,]=drawBeta(myn,sigmasqPost[i],invomegan)
}

hist(bmatrixPost[,1],breaks=100,main="b0 postdist")
hist(bmatrixPost[,2],breaks=100, main="b1 postdist")
hist(bmatrixPost[,3],breaks=100,main="b2 postidst")
hist(sigmasqPost)
plot(temp$time,temp$temp,col="red",main="Plot of the temp data over times")
responsePostTemp=matrix(0,nDraws,length(temp$time))

for(i in 1:nDraws){
  betaTemp=sapply(temp$time,calcRegr,b=bmatrixPost,rows=i)
  responsePostTemp[i,]=betaTemp
}

responsePost=c()
credInterval=matrix(0,length(temp$time),2)
for(i in 1:length(temp$time)){
  sortedTemp=sort(responsePostTemp[,i])
  responsePost=c(responsePost,(sortedTemp[500]+sortedTemp[501])/2)
  credInterval[i,]=quantile(sortedTemp,probs=c(0.025,0.975))
}

CredIntervalVec <- apply(X=responsePostTemp,MARGIN = 1,FUN=function(x) quantile(x,c(0.025,0.975)))

lines(temp$time,responsePost)
lines(temp$time,credInterval[,2], col="green")
lines(temp$time,credInterval[,1], col="green")

#The interval bands contain most of the data points except around 4% of them. They should contain most of the 
# data points as longs as the model is able to capture the reality. Since the model is quite simple
# (quadratic polynomial) it can be difficult to capture the curve of temperature over the year. And therefore we miss
# a couple of data points

# (c) It is of interest to locate the time with the highest expected temperature (that
# is, the time where f(time) is maximal). Let's call this value xtilde. Use the
# simulations in b) to simulate from the posterior distribution of xtilde. [Hint: the
# regression curve is a quadratic. You can find a simple formula for xtilde given b0, b1 and b2.]
maxtempTime<-c()
calcMax <- function(b,rows){
  return(-b[rows,2]/(2*b[rows,3]))
}
for(i in 1:nDraws){
 maxtempTime<-c(maxtempTime, calcMax(bmatrixPost,i))
}

hist(maxtempTime,breaks=1000,xlim=c(0,1))

#The time with the expected highest temperature is approx 0.55 which corresponds
# to the middle/end of June which is reasonable to be the hottest period of the year.

# (d) Say now that you want to estimate a polynomial model of order 7, but you
# suspect that higher order terms may not be needed, and you worry about overfitting.
# Suggest a suitable prior that mitigates this potential problem. You do
# not need to compute the posterior, just write down your prior. [Hint: the task
# is to specify my0 and omega0 in a smart way.]

#a suitable prior for this task would be to set my0 to 0 in order to obtain increased shrinkage
# for coefficients close to zero. You also want to set omega0 to Lambda*I,
# this would mean for larger values of lambda beta values
# would be close to zero since the spread of the distrubtion of beta values would decrease
# In this case, when worrying about overfitting, it is a good idea to use a larger value for Lambda,
# in order to penalize more => less overfitting
