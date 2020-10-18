#Assignment 2
#Assume that you have asked 10 randomly selected persons about their monthly income 
#(in thousands Swedish Krona) and obtained the following ten observations: 44,
#25, 45, 52, 30, 63, 19, 50, 34 and 67. A common model for non-negative continuous
#variables is the log-normal distribution. The log-normal distribution log N (mu, sigma^2)
# has density function
# p(y|mu,sigma^2)=log-normal distribution density function (see assignment)
#for y > 0, mu > 0 and sigma^2 > 0. The log-normal distribution is related to the
#normal distribution as follows: if y ~log N (mu, sigma^2) then log y ~ N (mu, sigma^2). 
#Let y1..yn|mu, sigma^2 iid~ log N (mu, sigma^2), where mu = 3.7 is assumed to be known but sigma^2
# is unknown with non-informative prior p(sigma^2) proportional to 1/sigma^2. The posterior for sigma^2
#is the inv - chi2(n, tao^2) distribution, where
# tao^2= sum(log(yi-mu)^2/n

#(a) Simulate 10,000 draws from the posterior of sigma^2(assuming mu = 3.7) 
#and compare with the theoretical Inv - chi2(n, tao^2) posterior distribution.


calcMSE <- function(n, mean, data){
  return((1/(n-1)*sum((data-mean)^2)))
}

x <- c(44,25,45,52,30,63,19,50,34,67)
n=length(x)

#function for tao^2 calculation
calcTao<-function(data,mu,n){
  return(sum((log(data)-mu)^2)/n)
}

taosq<-calcTao(x,3.7,n)
drawX=rchisq(10000,n)
sigmasq=(n)*taosq/drawX

calcGamma <- function(x, df, taosq){
first = ((taosq*df/2)^(df/2))/gamma(df/2)
second = (exp((-df*taosq)/(2*x)))/(x^(1+df/2))
return(first*second)
}


plot(density(sigmasq),main="distribution simulated in black vs theoretical in red")
xvals=seq(0.001, 3, 0.001)
set.seed(12345)
lines(xvals, calcGamma(xvals,10,taosq), col="red")

#as seen in the plot the simulated follows the theoretical in red with good precision

# (b) The most common measure of income inequality is the Gini coefficient, G,
# where 0 <= G <= 1. G = 0 means a completely equal income distribution,
# whereas G = 1 means complete income inequality. See Wikipedia for more
# information. It can be shown that G = 2phi(sigma/sqrt(2))-1 when incomes follow a
# log N (my, sigma2) distribution. phi(z) is the cumulative distribution function (CDF)
#   for the standard normal distribution with mean zero and unit variance. Use
# the posterior draws in a) to compute the posterior distribution of the Gini
# coefficient G for the current data set

G <- 2* pnorm(sqrt(sigmasq/2),0,1)-1
plot(density(G))

# as seen in the plot the gini coefficient most probably lies somewhere around 0.2 or a little higher.

# (c) Use the posterior draws from b) to compute a 90% equal tail credible interval
# for G. A 90% equal tail interval (a, b) cuts off 5% percent of the posterior
# probability mass to the left of a, and 5% to the right of b. Also, do a kernel
# density estimate of the posterior of G using the density function in R with
# default settings, and use that kernel density estimate to compute a 90% Highest
# Posterior Density interval for G. Compare the two interval

sortedG = sort(G)
sortedG =sortedG[(0.05*length(G)+1):(0.95*length(G))]
credibleIntervalG = c(min(sortedG),max(sortedG))
abline(v = credibleIntervalG[1])
abline(v = credibleIntervalG[2])

densityG <- density(G)
plot(densityG)
df.densityG <- data.frame(x=densityG$x, y=densityG$y)
df.densityG <- df.densityG[order(-df.densityG$y),]
index = dim(df.densityG)[1]
df.densityG <- data.frame(x= df.densityG$x, y=df.densityG$y, t=cumsum(df.densityG$y)/sum(df.densityG$y))
rows = which(df.densityG$t<0.90)
df.densityG <- df.densityG[rows,]
credibleDensityIntervalG = c(min(df.densityG$x),max(df.densityG$x))
print(credibleDensityIntervalG)
abline(v = credibleDensityIntervalG[1])
abline(v = credibleDensityIntervalG[2])
# cumDensity = cumDensity[cumDensity>0.05]
# cumDensity = cumDensity[cumDensity<0.95]
# cumDensity
