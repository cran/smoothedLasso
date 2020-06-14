# **************************************************************************************************************************************
# Methodology for smoothing the LASSO objective function using Nesterov smoothing.
# 
# The main functions provided in the package are:
# 
# 1) The smoothed LASSO objective function 'objFunctionSmooth' and its gradient 'objFunctionSmoothGradient'.
# 2) The minimization procedure 'solveSmoothedLASSO' returning the regression estimates for the smoothed LASSO.
# 3) The progressive smoothing approach 'solveSmoothedLASSOSequence' returning the regression estimates for the smoothed LASSO.
# 
# **************************************************************************************************************************************





# ************************************************************************
# Auxiliary function needed for the Nesterov smoothing framework
# 
# All the functions in this section are not exported.
# ************************************************************************

# Definition of a one dimensional PWA (piecewise affine function).
# Since each PWA is assumed to go to the origin (as does the absolute value), each PWA is parameterized by (p1,0) with only the slope p1 given in input vector 'p'.
# That is, 'p' is the vector of all p1 for as many PWAs as the length of 'p'.
# The function evaluates all PWAs at all one dimensional points given in the vector 'x'.
# Each row returned by the function PWA contains the evaluation of all pieces of the PWA parameterized by 'p' at each point 'x'.
PWA <- function(p,x) {
	outer(x,p)
}

# Fast version of the Michelot (1986) projection algorithm. The input 'x' can be a matrix where each row is projected separately.
michelot <- function(x) {
	d <- ncol(x)
	u <- t(apply(x,1, function(z) sort(z,decreasing=T) ))
	v <- t((1-apply(u,1,cumsum))/(1:d))
	r <- apply(u+v,1, function(z) max(which(z>0)) )
	l <- v[cbind(1:nrow(x),r)]
	return(pmax(x+l,0))
}

# Nesterov with entropy and squared error prox smoothing.
# The function smoothes several functions simultaneously.
# The input 'fx' can be a matrix, where each row contains the values per PWA. The input 'mu' is the Nesterov smoothing parameter.
nesterov <- function(fx,mu,entropy=F) {
	if(entropy) {
		# analytic expression of Nesterov smoothing obtained with entropy prox function
		return( mu*log( rowMeans(exp(fx/mu)) ))
	}
	else {
		# squared error prox function
		m <- ncol(fx)
		w <- michelot(fx/mu-1/m)
		return( rowSums(w*fx) - mu*rowSums((w-1/m)**2)/2 )
	}
}

# Nesterov's smoothed PWA. The function smoothes several functions simultaneously.
# The function inputs are the vector 'p' of slopes for the one dimensional PWAs, the vector of one dimensional points 'x', the Nesterov smoothing parameter 'mu',
# and a boolean flag 'entropy' to switch between the entropy prox function and the squared error prox function.
nesterovPWA <- function(p,x,mu,entropy=F) {
	nesterov(PWA(p,x),mu,entropy)
}

# Gradient of Nesterov smoothing with respect to 'x' for the entropy prox function.
# The input 'p' is a vector of one dimensional slopes defining the PWAs, 'x' is a vector of one dimensional points, and 'mu' is the Nesterov smoothing parameter.
nesterovArgumentGradient_entropy <- function(p,x,mu) {
	fx <- PWA(p,x)
	r <- t(exp(fx/mu))*p
	colSums(r)/rowSums(exp(fx/mu))
}

# Gradient of Nesterov smoothing with respect to 'x' for the squared error prox function.
# The input 'p' is a vector of one dimensional slopes defining the PWAs, 'x' is a vector of one dimensional points, and 'mu' is the Nesterov smoothing parameter.
nesterovArgumentGradient_sq <- function(p,x,mu) {
	fx <- PWA(p,x)
	# Michelot projection
	m <- ncol(fx)
	w <- michelot(fx/mu-1/m)
	# gradient
	colSums(t(w)*p)
}

# Gradient of Nesterov's smoothed PWA.
# The input 'p' is a vector of one dimensional slopes defining the PWAs, 'x' is a vector of one dimensional points, and 'mu' is the Nesterov smoothing parameter.
# The boolean flag 'entropy' switches between the entropy prox function and the squared error prox function.
nesterovArgumentGradient <- function(p,x,mu,entropy=F) {
	if(entropy) nesterovArgumentGradient_entropy(p,x,mu)
	else nesterovArgumentGradient_sq(p,x,mu)
}





# ************************************************************************
# Definition of the LASSO objective function
# ************************************************************************

#' Auxiliary function defining the LASSO objective function.
#' 
#' @param beta The \eqn{p}-vector of coefficients.
#' @param X The data matrix of dimensions \eqn{n \times p}.
#' @param y The \eqn{n}-vector of responses.
#' @param lambda The LASSO regularization parameter.
#' 
#' @return The value of the LASSO objective function.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' beta <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% beta
#' lambda <- 1
#' print(objFunction(beta,X,y,lambda))
#' 
#' @export
objFunction <- function(beta,X,y,lambda) {
	n <- nrow(X)
	1/n*sum((y - X %*% beta)**2) + lambda*sum(abs(beta))
}

#' Auxiliary function which computes the (non-smooth) gradient of the LASSO objective function with respect to \eqn{\beta}.
#' 
#' @param beta The \eqn{p}-vector of coefficients.
#' @param X The data matrix of dimensions \eqn{n \times p}.
#' @param y The \eqn{n}-vector of responses.
#' @param lambda The LASSO regularization parameter.
#' 
#' @return The value of the gradient of the LASSO objective function at \eqn{\beta}.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' beta <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% beta
#' lambda <- 1
#' print(objFunctionGradient(beta,X,y,lambda))
#' 
#' @export
objFunctionGradient <- function(beta,X,y,lambda) {
	n <- nrow(X)
	as.numeric(-2/n * t(y - X %*% beta) %*% X) + lambda*sign(beta)
}





# ************************************************************************
# Smoothed LASSO objective function and its gradient
# ************************************************************************

#' Auxiliary function defining the smoothed LASSO objective function.
#' 
#' @param beta The \eqn{p}-vector of coefficients.
#' @param X The data matrix of dimensions \eqn{n \times p}.
#' @param y The \eqn{n}-vector of responses.
#' @param lambda The LASSO regularization parameter.
#' @param mu The Nesterov smoothing parameter.
#' @param entropy A boolean switch to select the entropy prox function (default) or the squared error prox function.
#' 
#' @return The value of the smoothed LASSO objective function.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' beta <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% beta
#' lambda <- 1
#' print(objFunctionSmooth(beta,X,y,lambda,mu=0.1))
#' 
#' @export
objFunctionSmooth <- function(beta,X,y,lambda,mu,entropy=T) {
	n <- nrow(X)
	# two lines through the origin with slopes -1 and +1 define the PWA for the absolute value
	p <- c(-1,+1)
	# smooth each component of beta separately and sum them up to obtain the L1 norm of beta
	1/n*sum((y - X %*% beta)**2) + lambda*sum(nesterovPWA(p,beta,mu,entropy))
}

#' Auxiliary function which computes the gradient of the smoothed LASSO objective function with respect to \eqn{\beta}.
#' 
#' @param beta The \eqn{p}-vector of coefficients.
#' @param X The data matrix of dimensions \eqn{n \times p}.
#' @param y The \eqn{n}-vector of responses.
#' @param lambda The LASSO regularization parameter.
#' @param mu The Nesterov smoothing parameter.
#' @param entropy A boolean switch to select the entropy prox function (default) or the squared error prox function.
#' 
#' @return The value of the gradient of the LASSO objective function at \eqn{\beta}.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' beta <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% beta
#' lambda <- 1
#' print(objFunctionSmoothGradient(beta,X,y,lambda,mu=0.1))
#' 
#' @export
objFunctionSmoothGradient <- function(beta,X,y,lambda,mu,entropy=T) {
	n <- nrow(X)
	# two lines through the origin with slopes -1 and +1 define the PWA for the absolute value
	p <- c(-1,+1)
	as.numeric(-2/n * t(y - X %*% beta) %*% X) + lambda*nesterovArgumentGradient(p,beta,mu,entropy)
}





# ************************************************************************
# Linear regression with the unsmoothed and smoothed LASSO
# 
# Additionally, the progressive smoothing procedure is defined here.
# ************************************************************************

#' Minimize the unsmoothed LASSO objective function with respect to \eqn{\beta}.
#' Three options are available: BFGS with analytical gradient (\eqn{method=0}), BFGS with numerical gradient (\eqn{method=1}), and simulated annealing which is gradient free (\eqn{method=2}).
#' The default is \eqn{method=0}.
#' 
#' @param X The data matrix of dimensions \eqn{n \times p}.
#' @param y The \eqn{n}-vector of responses.
#' @param lambda The LASSO regularization parameter.
#' @param method The method used for minimization: BFGS with analytical gradient (\eqn{method=0}), BFGS with numerical gradient (\eqn{method=1}), and simulated annealing which is gradient free (\eqn{method=2}). The default is \eqn{method=0}.
#' 
#' @return The LASSO estimator \eqn{\beta}.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' beta <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% beta
#' lambda <- 1
#' print(solveUnsmoothedLASSO(X,y,lambda))
#' 
#' @export
solveUnsmoothedLASSO <- function(X,y,lambda,method=0) {
	p <- ncol(X)
	# analytical gradient
	if(method==0) {
		optim(	par = runif(p),
				fn = function(z) objFunction(z,X,y,lambda),
				gr = function(z) objFunctionGradient(z,X,y,lambda),
				method = "BFGS")$par
	}
	# numerical gradient
	else if(method==1) {
		optim(	par = runif(p),
				fn = function(z) objFunction(z,X,y,lambda),
				method = "BFGS")$par
	}
	# gradient free
	else if(method==2) {
		optim(	par = runif(p),
				fn = function(z) objFunction(z,X,y,lambda),
				method = "SANN")$par
	}
}

#' Minimize the smoothed LASSO objective function with respect to \eqn{\beta} using BFGS.
#' 
#' @param X The data matrix of dimensions \eqn{n \times p}.
#' @param y The \eqn{n}-vector of responses.
#' @param lambda The LASSO regularization parameter.
#' @param mu The Nesterov smoothing parameter.
#' @param entropy A boolean switch to select the entropy prox function (default) or the squared error prox function.
#' 
#' @return The LASSO estimator \eqn{\beta}.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' beta <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% beta
#' lambda <- 1
#' print(solveSmoothedLASSO(X,y,lambda,mu=0.1))
#' 
#' @export
solveSmoothedLASSO <- function(X,y,lambda,mu,entropy=T) {
	p <- ncol(X)
	optim(	par = runif(p),
			fn = function(z) objFunctionSmooth(z,X,y,lambda,mu,entropy),
			gr = function(z) objFunctionSmoothGradient(z,X,y,lambda,mu,entropy),
			method = "BFGS")$par
}

#' Minimize the smoothed LASSO objective function with respect to \eqn{\beta} using the progressive smoothing algorithm.
#' 
#' @param X The data matrix of dimensions \eqn{n \times p}.
#' @param y The \eqn{n}-vector of responses.
#' @param lambda The LASSO regularization parameter.
#' @param muSeq The sequence of Nesterov smoothing parameters. The default is \eqn{2^{-n}} for \eqn{n \in \{0,\ldots,5\}}.
#' @param entropy A boolean switch to select the entropy prox function (default) or the squared error prox function.
#' 
#' @return The LASSO estimator \eqn{\beta}.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' library(smoothedLasso)
#' require(Matrix)
#' n <- 100
#' p <- 500
#' beta <- runif(p)
#' X <- Matrix(sample(0:1,size=n*p,replace=TRUE),nrow=n,ncol=p,sparse=TRUE)
#' y <- X %*% beta
#' lambda <- 1
#' print(solveSmoothedLASSOSequence(X,y,lambda))
#' 
#' @export
solveSmoothedLASSOSequence <- function(X,y,lambda,muSeq=2**seq(3,-6),entropy=T) {
	p <- ncol(X)
	res <- runif(p)
	for(mu in muSeq) {
		res <- optim(	par = res,
						fn = function(z) objFunctionSmooth(z,X,y,lambda,mu,entropy),
						gr = function(z) objFunctionSmoothGradient(z,X,y,lambda,mu,entropy),
						method = "BFGS")$par
	}
	return(res)
}
