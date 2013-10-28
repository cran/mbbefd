\name{mbbefd-package}
\alias{mbbefd-package}
\alias{mbbefd}
\docType{package}
\title{
Package to Model MBBEFD distribution
}
\description{
This package models MBBEFD distribution as well as Exposure Curves for Reinsurance exposure rating.
}
\details{
\tabular{ll}{
Package: \tab mbbefd\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-10-28\cr
License: \tab GPL 2\cr
}
}
\author{
Giorgio Spedicato
Maintainer:  Giorgio Spedicato <spedicato_giorgio@yahoo.it>
}
\references{
BERNEGGER, STEFAN. THE SWISS RE EXPOSURE CURVES AND THE MBBEFD DISTRIBUTION CLASS. Astin Bulletin (1997): 99.
}
\keyword{package}
\examples{
#sample variates
losses<-rmbbefd(n=10,a=0.2,b=0.04)
}
