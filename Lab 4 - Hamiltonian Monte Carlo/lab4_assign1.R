# remove.packages("rstan")
# if (file.exists(".RData")) file.remove(".RData")
#  install.packages("rstan")
 
# Write a function in R that simulates data from the AR(1)-process
# X(t) = mu + phi(X(t-1)-mu)+epsilon(t), epsilon(t) ~ N(0,sigma2)
# for given values of mu, phi and sigma2  Start the process at x1 = mu and then simulate
# values for x(t)
# for t = 2,3,..., T and return the vector x(1:T) containing all time
# points. Use mu = 10, sigma2 = 2 and T = 200 and look at some different realizations
# (simulations) of x(1:T) for values of mu between ???1 and 1 (this is the interval of
#  phi where the AR(1)process is stable). Include a plot of at least one realization
# in the report. What effect does the value of phi have on x(1:T) ?
  
#Parameters
t=200
sigma2=2
mu=10
phiVector=c(seq(-0.9,0.9,0.1))

ARprocess <- function(sigma2, mu, phi, t){
  xVector <- c(mu)
  for(i in 2:t){
    epsilon=rnorm(1,0,sqrt(sigma2))
    xVector <- c(xVector,mu+phi*(xVector[i-1]-mu)+epsilon)
  }
  return(xVector)
}

set.seed(12345)
xMatrix <- matrix(0,t,length(phiVector))
xMatrix[1,]=10
count=0
for(phi in phiVector){
  count=count+1
  xMatrix[,count]=ARprocess(sigma2,mu,phi,t)
}

for(i in 1:length(phiVector)){
  plot(1:t,xMatrix[,i],type="l",main=paste("Plot of realization of AR-process"),xlab="Iteration",
       sub=paste("Phi:",phiVector[i]),ylab="Value", col="grey")
}


## With phi-values below 0, x(t) goes in the opposite direction of x(t-1) which makes the value
## oscillate faster around mu. As phi -> 1, the correlation between x(t) and x(t-1) increases, 
## causing the process to oscillate slower around mu.

## Use your function from a) to simulate two AR(1)-processes x(1:t) with
## phi=0.3 and y(1:t) with phi =0.95. Now treat your simulated vectors as synthetic
## data, and treat the values of mu, phi and sigmasq as unknown and estimate
## them using MCMC. Implement Stan-code that samples from the posterior of the three
# parameters, using suitable non-informative priors of your choice. Hint: Look
# at the time-series models examples in the Stan user's guide/reference manual,
# and note the different parameterization used here.

# i. Report the posterior mean, 95% credible intervals and the number of effective posterior 
# samples for the three inferred parameters for each of the
# simulated AR(1)-process. Are you able to estimate the true values?

Vector03=rep(0,t)
Vector095=rep(0,t)
set.seed(12345)
Vector03=ARprocess(sigma2,mu,0.3,t)
set.seed(12345)
Vector095=ARprocess(sigma2,mu,0.95,t)

StanModel ='
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
}

model {
  for (n in 2:N)
    y[n] ~ normal(mu + phi * (y[n-1]-mu), sigma);
}'

library(rstan)

fit095 = stan(model_code=StanModel,
              data=list(N=t, y=Vector095))

fit03 = stan(model_code= StanModel,
             data=list(N=t, y=Vector03))

print(fit095)
print(fit03)


##For the simulated data with phi=0.3 we receive a posterior mean of mu = 10.06, phi=0.39, sigma=1.36
##For the simulated data with phi=0.95 we recieve a posterior mean of mu = 10.61, phi=0.98, sigma=1.44
##The 95% credible intervals were for phi=0.3: 9.72<mu<10.39 , 0.26<phi<0.52 , 1.23<sigma<1.52
##The 95% credible intervals were for phi=0.95: -31.62<mu<

##Are you able to estimate the true values? It is possible to estimate the true values for the dataset which were
##simulated with phi=0.3 but not for the dataset simulated with phi=0.95. This is due to the oscillation around the
##true value of mu being much slower for phi=0.95 compared to phi=0.3

burnin=1000
niter=2000
postDraws03 <- extract(fit03)
postDraws095 <- extract(fit095)

par(mfrow=c(1,1))

plot(postDraws03$mu[1000:2000],type="l",ylab="mu",main="Traceplot for mu(03)")

# Do traceplots of the first chain for dataset X
par(mfrow = c(1,1))
plot(postDraws03$mu[1000:2000], postDraws03$phi[1000:2000],ylab="phi",xlab="mu",main="Traceplot for dataset phi03")
# Do automatic traceplots of all chains
traceplot(fit03)
# Bivariate posterior plots
pairs(fit03)
# Do traceplots of the first chain for dataset Y
par(mfrow = c(1,1))
plot(postDraws095$mu[1000:2000],postDraws095$phi[1000:2000],ylab="phi",xlab="mu",main="Traceplot for dataset 095")
# Do automatic traceplots of all chains
traceplot(fit095)
# Bivariate posterior plots
pairs(fit095)



## The convergence of the samplers are quite different. For the sample which used phi =0.3 we see a convergence around mu=10.3
## and phi = 0.3 while the sample which used phi=0.95 doesn't converge at all. This correlates with the fact that the credible
## intervals for the non-converging sample was very large. We can see for posterior distribution of the sample with phi = 0.95
## that for lower values of phi the distribution centers around a value between 10-20 while for higher values of phi around 1
## it spreads widely with a range from -40 to 80. This agrees with the fact that the posterior distribtuion with a phi = 0.3
## centers around a value of 10

# (c) The data campy.dat contain the number of cases of campylobacter infections
# in the north of the province Quebec (Canada) in four week intervals from
# January 1990 to the end of October 2000. It has 13 observations per year and
# 140 observations in total. Assume that the number of infections c(t) at each
# time point follows an independent Poisson distribution when conditioned on a
# latent AR(1)-process x(t), that is 
# c(t)|x(t) ~ Poisson(exp(x(t))),
# where x(t) is an AR(1)-process as in a). Implement and estimate the model in
# Stan, using suitable priors of your choice. Produce a plot that contains both the
# data and the posterior mean and 95% credible intervals for the latent intensity
# theta(t)=exp(x(t)) over time. [Hint: should x(t) be seen as data or parameters?]

campy <- read.table("campy.dat", header=TRUE)
T=dim(campy)[1]

StanPois ='
data {
  int<lower=0> T;
  int c[T];
}

parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  vector[T] x;
}
model {
phi ~ uniform(-1,1);
  for (t in 2:T)
    x[t] ~ normal(mu + phi * (x[t-1]-mu), sigma);
    
  for(t in 1:T)
  c[t] ~ poisson(exp(x[t]));

}

generated quantities{
vector[T] post_mean;
post_mean = exp(x);
}'

fitPois <- stan(model_code=StanPois,
                data=list(T=T,c=campy$c))
print(fitPois)
summary(fitPois)

poisMeanList <- fitPois@.MISC$summary$msd
postMean <- poisMeanList[grep("post_mean", rownames(poisMeanList)),]

plot(campy$c, col="blue", ylab="No. of Infections")
points(postMean[,1], col="black", type="l")

quantiles=fitPois@.MISC$summary$quan

quantilesPostMean=quantiles[grep("post_mean", rownames(quantiles)),]
credIntervalPostMean=matrix(0,dim(quantilesPostMean)[1], 2)
credIntervalPostMean[,1]=quantilesPostMean[,1]
credIntervalPostMean[,2]=quantilesPostMean[,ncol(quantilesPostMean)]

lines(credIntervalPostMean[,1], col="gray", lty=21)
lines(credIntervalPostMean[,2], col="gray", lty=21)
title(main="Plot of data vs approximated posterior")
legend("topleft", box.lty=1,legend=c("Data","Posterior mean","95% Cred. Interval"),col=c("blue","black","gray"),
       pch=c(1,NaN,NaN), lwd=c(NaN,1,1), lty=c(NaN,1,21))


## As seen in the plot the posterior mean follows the data accurateley. Almost all of the datapoints
## are inside the credible intervals which aren't that wide which indicates that the approximated posterir
## resembles the reality shown by the data well.

# (d) Now, assume that we have a prior belief that the true underlying intensity theta(t)
# varies more smoothly than the data suggests. Change the prior for sigmasq
# so that
# it becomes informative about that the AR(1)-process increments epsilon(t) should be
# small. Re-estimate the model using Stan with the new prior and produce the
# same plot as in c). Has the posterior for theta(t) changed?
  
StanPoisPrior ='
data {
  int<lower=0> T;
  int c[T];
}

parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  vector[T] x;
}
model {
phi ~ uniform(-1,1);
sigma ~ scaled_inv_chi_square(140,0.15);
  for (t in 2:T)
    x[t] ~ normal(mu + phi * (x[t-1]-mu), sigma);
    
  for(t in 1:T)
  c[t] ~ poisson(exp(x[t]));

}

generated quantities{
vector[T] post_mean;
post_mean = exp(x);
}'


fitPoisPrior <- stan(model_code=StanPoisPrior,
                data=list(T=T,c=campy$c))
print(fitPoisPrior)
summary(fitPoisPrior)

poisMeanListPrior <- fitPoisPrior@.MISC$summary$msd
postMeanPrior <- poisMeanListPrior[grep("post_mean", rownames(poisMeanListPrior)),]

plot(campy$c, col="blue", ylab="No. of Infections")
points(postMeanPrior[,1], col="black", type="l")

quantilesPrior=fitPoisPrior@.MISC$summary$quan

quantilesPostMeanPrior=quantilesPrior[grep("post_mean", rownames(quantilesPrior)),]
credIntervalPostMeanPrior=matrix(0,dim(quantilesPostMeanPrior)[1], 2)
credIntervalPostMeanPrior[,1]=quantilesPostMeanPrior[,1]
credIntervalPostMeanPrior[,2]=quantilesPostMeanPrior[,ncol(quantilesPostMeanPrior)]

lines(credIntervalPostMeanPrior[,1], col="gray", lty=21)
lines(credIntervalPostMeanPrior[,2], col="gray", lty=21)
title(main="Plot of data vs approximated posterior with prior")
legend("topleft", box.lty=1,legend=c("Data","Posterior mean","95% Cred. Interval"),col=c("blue","black","gray"),
       pch=c(1,NaN,NaN), lwd=c(NaN,1,1), lty=c(NaN,1,21))


## With a small specified prior of sigma the posterior mean varies less, and moves more smoothly. The consequence of this
## is that more datapoints lie outside of the credible interval, which suggests that the approximated postieror cannot capture
## the reality described by the data as accurately as before. However, by doing this one can avoid overfitting when the
## model is applied to a new dataset.

