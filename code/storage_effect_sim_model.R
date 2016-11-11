##  Two Species Annual Plant Population Model: simulates the dynamics of
##  the seedbanks of two competing annual plants. The plants coxist by the
##  storage effect. This is just to show two time series (one with coexistence,
##  one with exclusion) and to show that reducing environmental variability
##  causes exclusion and a reduction in ecosystem stability, here measured by
##  total community abundance.
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: November 11, 2016


rm(list=ls(all.names = TRUE))



####
#### Parameters ----------------------------------------------------------------
####

sigE   <- c(0.05,0.5)
rho    <- c(0.5,0,-0.5)
s      <- c(0.5, 0.5)
alpha  <- c(1,1)
lambda <- c(101,99)
nTime  <- 2000



####
####  Libraries ----------------------------------------------------------------
####
library(mvtnorm)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(zoo)



####
####  Annual Plant Model Functions ---------------------------------------------
####

### Update seedbank function
updateN <- function(g, s, alpha, lambda, lastN){
  newN <- numeric(2)
  for(i in 1:2){
    newN[i] <- lastN[i]*s[i]*(1-g[i]) + ((lambda[i]*g[i]*lastN[i]) / (1 + (alpha[i]*g[i]*lastN[i] + alpha[-i]*g[-i]*lastN[-i])))
    if(newN < 1) { newN <- 0 }
  }
  return(newN)
}

### Get germination fractions function
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e      <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g      <- exp(e) / (1+exp(e))
  return(g)
}



####
#### Simulations ---------------------------------------------------------------
####

par(mfrow=c(1,3), mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1), las=1)

### Low variability // competitive exclusion
gSeries <- getG(sigE[1], rho[2], nTime)
N       <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,]   <- c(10,10)
for(t in 2:nTime){
  N[t,] <- updateN(g=gSeries[t,], s, alpha, lambda, lastN=N[t-1,])
}
matplot(c(1:nTime), N, type="l", lwd=1, lty=1, bty="n", xlab="Time", 
        ylab="Population Size (N)", main="Competitive Exclusion")
totpop1 <- rowSums(N)

### High variability // coexistence
gSeries <- getG(sigE[2], rho[2], nTime)
N       <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,]   <- c(10,10)
for(t in 2:nTime){
  N[t,] <- updateN(g=gSeries[t,], s, alpha, lambda, lastN=N[t-1,])
}
matplot(c(1:nTime), N, type="l", lwd=1, lty=1, bty="n", xlab="Time", 
        ylab="Population Size (N)", main="Species Coexistence")
totpop2 <- rowSums(N)

### Plot the rolling CV
mycv <- function(x) {sd(x) / mean(x)}
rolling_cv <- data.frame(cv1 = rollapply(totpop1, width=10, FUN=mycv, fill=NA),
                         cv2 = rollapply(totpop2, width=10, FUN=mycv, fill=NA),
                         iteration = 1:length(totpop1))
matplot(rolling_cv[,c(1:2)], type="l", col=c("grey45","steelblue"), bty="n")



