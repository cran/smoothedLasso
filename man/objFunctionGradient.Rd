% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothedLasso.r
\name{objFunctionGradient}
\alias{objFunctionGradient}
\title{Auxiliary function which computes the (non-smooth) gradient of an L1 penalized regression operator.}
\usage{
objFunctionGradient(betavector, w, du, dv, dw)
}
\arguments{
\item{betavector}{The vector of regression coefficients.}

\item{w}{The function encoding the dependence structure among the regression coefficients.}

\item{du}{The derivative (gradient) of the objective of the regression operator.}

\item{dv}{The derivative (gradient) of the penalty of the regression operator.}

\item{dw}{The derivative (Jacobian matrix) of the function encoding the dependence structure among the regression coefficients.}
}
\value{
The value of the gradient for the input \eqn{betavector}.
}
\description{
Auxiliary function which computes the (non-smooth) gradient of an L1 penalized regression operator.
}
\examples{
library(smoothedLasso)
n <- 100
p <- 500
betavector <- runif(p)
X <- matrix(runif(n*p),nrow=n,ncol=p)
y <- X \%*\% betavector
lambda <- 1
temp <- standardLasso(X,y,lambda)
print(objFunctionGradient(betavector,temp$w,temp$du,temp$dv,temp$dw))

}
\references{
Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
}
