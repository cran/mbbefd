### R code from vignette source 'mbbefd.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
  options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
	set.seed(123)


###################################################
### code chunk number 2: load:ec
###################################################
library(mbbefd)
ecunif(0:4/4)
ecbeta(0:4/4, 3, 2)
eecf(rbeta(100, 3, 2))


###################################################
### code chunk number 3: oinfl
###################################################
doibeta(0:4/4, 3, 2, 1/2)
poibeta(0:4/4, 3, 2, 1/2)


###################################################
### code chunk number 4: example1
###################################################
net<-mbbefdExposure(x=1/2, a=0.2,b=0.04)*40000
ceded<-40000-net


###################################################
### code chunk number 5: example1b
###################################################
expectedLoss<-1/dG(x=0,a=0.2,b=0.04)*40000
expectedLoss


###################################################
### code chunk number 6: survivalPlot
###################################################
x<-seq(from=0,to=1,by=.1)
smbbefd<-function(x,a,b) 1-pmbbefd(q=x,a=a,b=b)
y<-sapply(x, smbbefd,a=.2,b=.04)

plot(x, y, type="l",lwd=2,col="steelblue",xlab="normalized loss",
     ylab="survival function",main="Survival function for MBBEFD at a=0.2, b=0.04",ylim=c(0,1))


###################################################
### code chunk number 7: example1c
###################################################
pTotalLoss<-1-pmbbefd(q=1,a=0.2,b=0.04)
pTotalLoss


###################################################
### code chunk number 8: example1d
###################################################
qmbbefd(p=0.6,a=0.2,b=0.04)


###################################################
### code chunk number 9: example1e
###################################################
100*(1-pmbbefd(q=0.8,a=0.2,b=0.04))


###################################################
### code chunk number 10: example1e2
###################################################
simulatedLosses<-rmbbefd(n=10000,a=0.2,b=0.04)
mean(simulatedLosses)
sum(simulatedLosses==1)/length(simulatedLosses)


###################################################
### code chunk number 11: distrPlot
###################################################
hist(simulatedLosses, main="Simulated MBBEFD distribution",probability=TRUE,col="steelblue")


###################################################
### code chunk number 12: example1f
###################################################
integrate(dmbbefd,lower=0, upper=1, a=0.2, b=0.04)


###################################################
### code chunk number 13: fitting
###################################################
#get data
data1<-rmbbefd(n=1000,a = .2,b=.04)
data(loss, package = "copula")
data2<-pmin(1,pmax(0,loss$loss/loss$limit)) #capping loss data to lim

#functions used to initialize the parameters
#using one iteration of Method of Moments

#method of moments

giveFunction2Minimize<-function(mu,g) {
  out = function(b) (mu - (log(g*b)*(1 - b))/( log(b)*(1 - g*b)) )^2
  return(out)
}

giveFunction2Integrate<-function(b,g) {
  out = function(x) x^2*dmbbefd(x,b=b,g=g)
  return(out)
}

giveInits<-function(x) {
  m0<-mean(x)
  m2<-mean(x^2)
  
  #p<=1/g
  
  p0=m2 #m2 upper limit of p0
  g=1/p0
  
  #equate 1rst moment to get the mean
  myMin<-giveFunction2Minimize(mu=m0,g=g)
  b<-nlm(f=myMin,p=.1)$estimate
  
  #return a
  a=(g-1)*b/(1-g*b)
  out<-list(a=a, b=b)
  return(out)
}

###fitting process

library(fitdistrplus)
#using close starting points
est1<-fitdist(data=data1,distr = "mbbefd",method = "mle",start=list(a=.9,b=.14))
est1
#using estimated starting points
inits2<-giveInits(x=data2)
est2<-fitdist(data=data2,distr = "mbbefd",method = "mle",start=inits2)
est1


