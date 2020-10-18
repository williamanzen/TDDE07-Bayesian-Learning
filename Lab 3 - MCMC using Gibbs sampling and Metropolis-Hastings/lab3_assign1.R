## Assignment 1:
# (a) Normal model
# Assume the daily precipitation {y1, ..., yn} are independent normally distributed,
# y1, ..., yn|my, sigma2 ~ N (my, sigma2) where both my and sigma2 are unknown. 
#Let my ~ N (my0, tao20) independently of sigma2 ~ Inv-X2(v0, sigma20).
# i. Implement (code!) a Gibbs sampler that simulates from the joint posterior
# p(my, sigma2|y1, ..., yn). The full conditional posteriors are given on the slides
# from Lecture 7


data <- read.table("rainfall.dat")

#Setup
v0 <- 1
tausq0 <- 100
meandata <- mean(data$V1)
sigmasq0 <- 1
n= length(data$V1)
nDraws <- 5000
mu0=14.79
vn = v0+n
calcTau <- function(n, tau, sigmasq){
  return(1/(n/sigmasq+1/tau))
}

calcMy <- function(sigmasq, n, my0,mean, tausq){
  
  return(calcW(n, sigmasq, tausq)*mean+((1-calcW(n,sigmasq,tausq))*my0))
}

calcW <- function(n, sigmasq, tausq){
  return((n/sigmasq)/((n/sigmasq)+(1/tausq)))
}

calcSigmaHat <- function(v0, sigmasq0, data, my, n){
  return((v0*sigmasq0+sum((data-my)^2))/(v0+n))
}


#Gibbs sampling
myVector <- c()
sigmaVector <- c(sigmasq0)
tauVector <- c(tausq0)
for( i in 1:nDraws){
  myVector <- c(myVector, rnorm(1, calcMy(sigmaVector[i],n,mu0,meandata, tausq0), calcTau(n,tausq0,sigmaVector[i])))
  if(i<nDraws){
  drawX <- rchisq(1,vn)
  sigmaVector <- c(sigmaVector, vn*calcSigmaHat(v0,sigmasq0,data$V1,myVector[i],n)/drawX)
  }
}

plot(myVector,sigmaVector)

# ii. Analyze the daily precipitation using your Gibbs sampler in (a)-i. Evaluate
# the convergence of the Gibbs sampler by suitable graphical methods, for
# example by plotting the trajectories of the sampled Markov chains.


iter=seq(1001,5000,1)
plot(iter, myVector[1001:5000], type="l", main="Convergence plot for my")
plot(iter, sigmaVector[1001:5000], type="l", main="convergence plot for sigma")



# (b) Mixture normal model.
# Let us now instead assume that the daily precipitation {y1, ..., yn} follow an iid
# two-component mixture of normals model:
# p(yi|my, sigmasq, pi) = pi*N(yi|my1, sigmasq1) + (1 ??? pi)N (yi|sigmasq, sigmasq2),
# where my = (my1, my2) and sigmasq = (sigmasq1, sigmasq2).
# Use the Gibbs sampling data augmentation algorithm in NormalMixtureGibbs.R
# (available under Lecture 7 on the course page) to analyze the daily precipitation data.
# Set the prior hyperparameters suitably. Evaluate the convergence of the sampler.

##########    BEGIN USER INPUT #################
x <- as.matrix(data$V1)

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 1*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(14.79,nComp) # Prior mean of mu
tau2Prior <- rep(100,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(1,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws


# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
# Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(data$V1)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)


# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))

wVec <- c()
paramMatrix <- matrix(0, nIter,4) 


for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    wVec <- c(wVec, wPrior)
    
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
    paramMatrix[k,j]=mu[j]
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
    paramMatrix[k,j+2]=sigma2[j]
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topleft", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    
  }
  
}
plot(200:1000,wVec[200:1000],type="l")
plot(200:1000,paramMatrix[200:nrow(paramMatrix),1], type="l")
plot(200:1000,paramMatrix[200:nrow(paramMatrix),2], type="l")
plot(200:1000,paramMatrix[200:nrow(paramMatrix),3], type="l")
plot(200:1000,paramMatrix[200:nrow(paramMatrix),4], type="l")

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

#########################    Helper functions    ##############################################

## with prior alpha = 1, mu = 14.79, tau2=100, nu0=1, the final fitted density follows the histogram quite accurately with
# a mixture distribution with a mode slightly below 20. 

# (c) Graphical comparison.
# Plot the following densities in one figure: 1) a histogram or kernel density
# estimate of the data. 2) Normal density N (yi|my, sigma2) in (a); 3) Mixture of
# normals density p(yi|mean, sigma2, pi) in (b). Base your plots on the mean over all
# posterior draws.

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(myVector), sd = mean(sqrt(sigmaVector))), type = "l", lwd = 2, col = "blue")

## As seen in the plot, the blue line, representing the gibbs sampled distibribution is very similar to the true distribution 
# derived from the data in the previous plot (also the blue line). 


