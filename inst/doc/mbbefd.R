### R code from vignette source 'mbbefd.Rnw'

###################################################
### code chunk number 1: setup
###################################################
	options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
	set.seed(123)


###################################################
### code chunk number 2: load
###################################################
library(mbbefd)


###################################################
### code chunk number 3: figureBackoffice
###################################################
x=seq(0,1,by=0.01)
y=sapply(x,"mbbefdExposure",a=0.2,b=0.04)


###################################################
### code chunk number 4: drateplot
###################################################
plot(x, y, type="l",lwd=2,col="steelblue",xlab="normalized loss",ylab="exposure curve",main="MBBEFD Exposure Curve for a=0.2, b=0.04")


###################################################
### code chunk number 5: example1
###################################################
net<-mbbefdExposure(x=1/2, a=0.2,b=0.04)*40000
ceded<-40000-net


###################################################
### code chunk number 6: example1b
###################################################
expectedLoss<-1/dG(x=0,a=0.2,b=0.04)*40000
expectedLoss


###################################################
### code chunk number 7: example1c
###################################################
pTotalLoss<-1-pmbbefd(x=1,a=0.2,b=0.04)
pTotalLoss


###################################################
### code chunk number 8: example1d
###################################################
qmbbefd(p=0.6,a=0.2,b=0.04)


###################################################
### code chunk number 9: example1e
###################################################
100*(1-pmbbefd(x=0.8,a=0.2,b=0.04))


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


