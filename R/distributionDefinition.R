#exposure curve

.G<-function(x,a,b,g)
{
  if(missing(a)) a=((g-1)*b)/(1-g*b)
  out<-(log(a+b^x)-log(a+1))/(log(a+b)-log(a+1))
  return(out)
}

dG<-function(x,a,b,g)
{
  if(missing(a)) a<-((g-1)*b)/(1-g*b)
  out<-( (log(b)*b^x) /(a+b^x) )/(log(a+b)-log(a+1))
  return(out)
}

.Sx<-function(x,a,b,g)
{
  if(missing(a)) a<-((g-1)*b)/(1-g*b)
  out<-dG(x=x,a=a,b=b)/dG(x=0,a=a,b=b)
  return(out)
}

mbbefdExposure<-function(x, a, b,g)
{
  if(missing(a)) a<-((g-1)*b)/(1-g*b)
  if(x>1|x<0) stop("Error! x should be between 0 and 1")
  out<-.G(x=x, a=a, b=b)
  return(out)
}

####################################
#classical functions
#distirbution function
pmbbefd<-function(x,a,b,g)
{
  if(missing(a)) a<-((g-1)*b)/(1-g*b)
  if(x>1|x<0) stop("Error! x should lie between 0 and 1")
  out<-1-.Sx(x=x,a=a,b=b)
  return(out)
}
#inverse distribution funtcion (quantile)
qmbbefd<-function(p,a,b,g)
{
  if(missing(a)) a<-((g-1)*b)/(1-g*b)
  if(p>1|p<0) stop("Error! p should lie between 0 and 1")
  out<-log((a*(1-p))/(a+p))/log(b)
  #p > 1 - prob total loss would make to exceed 1. So it will return one.
  out<-ifelse(out>1, 1,out)
  return(out)
}

#density functio

dmbbefd<-function(x,a,b,g)
{
  if(missing(a)) a<-((g-1)*b)/(1-g*b)
  if(x>1|x<0) stop("Error! x should lie between 0 and 1")
  out<- -((a+1)*a*log(b)*b^x)/((a+b^x)^2)
  return(out)
}


#random generation function



.f4Random<-function(x,a,b) min(1,(log((a*(1-x))/(a+x)))/(log(b)))
# 
# b=12.648
# g=4.22069

rmbbefd<-function(n,a,b,g)
{
  if(missing(a)) a<-((g-1)*b)/(1-g*b)
  out<-numeric(n)
  u<-runif(n=n,min=0,max=1)
 out<-sapply(u,.f4Random,a=a,b=b)
  return(out)
}


####Swiss Re curves####

swissRe<-function(c)
{
  out<-numeric(2)
  b <- exp(3.1 - 0.15*c*(1+c))
  g <-exp(c*(0.78 + 0.12*c))
  out<-c(b,g)
  names(out)<-c("b","g")
  return(out)
}