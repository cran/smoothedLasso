# **************************************************************************************************************************************
# Methodology for smoothing L1 penalized regression operators using Nesterov smoothing.
# 
# The main functions provided in the package are:
# 
# 1) The functions defining popular L1 penalized regression operators, for instance 'standardLasso'.
# 2) The functions 'objFunction' and 'objFunctionGradient' defining the unsmoothed objective function of an L1 penalized regression operator and its gradient, as well as 'objFunctionSmooth' and 'objFunctionSmoothGradient' defining the smoothed objective function of an L1 penalized regression operator.
# 3) The minimization procedure 'solveFunction' returning the regression estimates for an unsmoothed or smoothed regression operator.
# 4) The progressive smoothing approach 'solveSmoothedSequence' returning the regression estimates for a smoothed regression operator.
# **************************************************************************************************************************************





# *****************************************************************************
# Auxiliary function needed for the Nesterov smoothing framework
# 
# All the functions in this section are not exported.
# *****************************************************************************

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
	u <- t(apply(x,1, function(z) sort(z,decreasing=TRUE) ))
	v <- t((1-apply(u,1,cumsum))/(1:d))
	r <- apply(u+v,1, function(z) max(which(z>0)) )
	l <- v[cbind(1:nrow(x),r)]
	return(pmax(x+l,0))
}

# Nesterov with entropy and squared error prox smoothing.
# The function smoothes several functions simultaneously.
# The input 'fx' can be a matrix, where each row contains the values per PWA. The input 'mu' is the Nesterov smoothing parameter.
nesterov <- function(fx,mu,entropy=TRUE) {
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
nesterovPWA <- function(p,x,mu,entropy=TRUE) {
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
nesterovArgumentGradient <- function(p,x,mu,entropy=TRUE) {
	if(entropy) nesterovArgumentGradient_entropy(p,x,mu)
	else nesterovArgumentGradient_sq(p,x,mu)
}

# Parametrization of a positive definite matrix of dimensions d*d via its Cholesky decomposition: x should be of length d*(d+1)/2 and will be put on the lower triangle of the Cholesky matrix
toCholesky <- function(d,x) {
	Theta <- matrix(0,d,d)
	Theta[lower.tri(Theta,diag=TRUE)] <- x
	Theta %*% t(Theta)
}





# *****************************************************************************
# Definition of L1 penalized regression operators
# 
# All functions return a list with the following 3+3 functions defining
# the regression operator and its derivative: u,v,w,du,dv,dw
# *****************************************************************************

#' Auxiliary function which returns the objective, penalty, and dependence structure among regression coefficients of the Lasso.
#' 
#' @param X The design matrix.
#' @param y The response vector.
#' @param lambda The Lasso regularization parameter.
#' 
#' @return A list with six functions, precisely the objective \eqn{u}, penalty \eqn{v}, and dependence structure \eqn{w}, as well as their derivatives \eqn{du}, \eqn{dv}, and \eqn{dw}.
#' 
#' @importFrom Rdpack reprompt
#' @references Tibshirani, R. (1996). Regression Shrinkage and Selection Via the Lasso. J Roy Stat Soc B Met, 58(1):267-288.
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' lambda <- 1
#' temp <- standardLasso(X,y,lambda)
#' 
#' @export
standardLasso <- function(X,y,lambda) {
	n <- nrow(X)
	p <- ncol(X)
	u <- function(beta) 1/n * sum((y - X %*% beta)**2)
	v <- function(z) lambda * sum(z)
	w <- function(beta) beta
	du <- function(beta) -2/n * as.vector(t(y - X %*% beta) %*% X)
	dv <- function(z) rep(lambda,p)
	dw <- function(beta) diag(p)
	list(u=u,v=v,w=w,du=du,dv=dv,dw=dw)
}

#' Auxiliary function which returns the objective, penalty, and dependence structure among regression coefficients of the elastic net.
#' 
#' @param X The design matrix.
#' @param y The response vector.
#' @param alpha The regularization parameter of the elastic net.
#' 
#' @return A list with six functions, precisely the objective \eqn{u}, penalty \eqn{v}, and dependence structure \eqn{w}, as well as their derivatives \eqn{du}, \eqn{dv}, and \eqn{dw}.
#' 
#' @importFrom Rdpack reprompt
#' @references Zou, H. and Hastie, T. (2005). Regularization and variable selection via the elastic net. J Roy Stat Soc B Met, 67(2):301-320.
#' @references Friedman, J., Hastie, T., Tibshirani, R., Narasimhan, B., Tay, K., Simon, N., and Qian, J. (2020). glmnet: Lasso and Elastic-Net Regularized Generalized Linear Models. R-package version 4.0.
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' alpha <- 0.5
#' temp <- elasticNet(X,y,alpha)
#' 
#' @export
elasticNet <- function(X,y,alpha) {
	n <- nrow(X)
	p <- ncol(X)
	u <- function(beta) 1/(2*n) * sum((y - X %*% beta)**2) + 1/2*(1-alpha) * sum(beta**2)
	v <- function(z) alpha * sum(z)
	w <- function(beta) beta
	du <- function(beta) -1/n * as.vector(t(y - X %*% beta) %*% X) + (1-alpha)*beta
	dv <- function(z) rep(alpha,p)
	dw <- function(beta) diag(p)
	list(u=u,v=v,w=w,du=du,dv=dv,dw=dw)
}

#' Auxiliary function which returns the objective, penalty, and dependence structure among regression coefficients of the fused Lasso.
#' 
#' @param X The design matrix.
#' @param y The response vector.
#' @param E The adjacency matrix which encodes with a one in position \eqn{(i,j)} the presence of an edge between variables \eqn{i} and \eqn{j}. Note that only the upper triangle of \eqn{E} is read.
#' @param lambda The first regularization parameter of the fused Lasso.
#' @param gamma The second regularization parameter of the fused Lasso.
#' 
#' @return A list with six functions, precisely the objective \eqn{u}, penalty \eqn{v}, and dependence structure \eqn{w}, as well as their derivatives \eqn{du}, \eqn{dv}, and \eqn{dw}.
#' 
#' @importFrom Rdpack reprompt
#' @references Tibshirani, R., Saunders, M., Rosset, S., Zhu, J., and Knight, K. (2005). Sparsity and Smoothness via the Fused Lasso. J Roy Stat Soc B Met, 67(1):91-108.
#' @references Arnold, T.B. and Tibshirani, R.J. (2020). genlasso: Path Algorithm for Generalized Lasso Problems. R package version 1.5.
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' E <- matrix(sample(c(TRUE,FALSE),p*p,replace=TRUE),p)
#' lambda <- 1
#' gamma <- 0.5
#' temp <- fusedLasso(X,y,E,lambda,gamma)
#' 
#' @export
fusedLasso <- function(X,y,E,lambda,gamma) {
	n <- nrow(X)
	p <- ncol(X)
	E <- upper.tri(E,diag=TRUE)
	indices <- unname(which(E,arr.ind=TRUE))
	m <- p + sum(E)
	u <- function(beta) 1/2 * sum((y - X %*% beta)**2)
	v <- function(z) lambda * sum(z[1:p]) + lambda*gamma * sum(z[(p+1):m])
	w <- function(beta) c(beta, beta[indices[,1]]-beta[indices[,2]])
	du <- function(beta) -1 * as.vector(t(y - X %*% beta) %*% X)
	dv <- function(z) c(rep(lambda,p), rep(lambda*gamma,m-p))
	temp <- matrix(0,m,p)
	for(j in 1:p) temp[j,j] <- 1
	for(j in (p+1):m) {
		temp[j,indices[j-p,1]] <- 1
		temp[j,indices[j-p,2]] <- -1
	}
	dw <- function(beta) temp
	list(u=u,v=v,w=w,du=du,dv=dv,dw=dw)
}

#' Auxiliary function which returns the objective, penalty, and dependence structure among regression coefficients of the graphical Lasso.
#' 
#' @param S The sample covariance matrix.
#' @param lambda The regularization parameter of the graphical Lasso.
#' 
#' @return A list with three functions, precisely the objective \eqn{u}, penalty \eqn{v}, and dependence structure \eqn{w}. Not all derivatives are available in closed form, and thus computing the numerical derivative of the entire objective function is recommended.
#' 
#' @importFrom Rdpack reprompt
#' @references Friedman, J., Hastie, T., and Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9(3):432-441.
#' @references Friedman, J., Hastie, T., and Tibshirani, R. (2019). glasso: Graphical Lasso: Estimation of Gaussian Graphical Models. R package version 1.11.
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' p <- 30
#' S <- matrix(rWishart(1,p,diag(p)),p)
#' lambda <- 1
#' temp <- graphicalLasso(S,lambda)
#' 
#' @export
graphicalLasso <- function(S,lambda) {
	p <- ncol(S)
	u <- function(beta) {
		Theta <- toCholesky(p,beta)
		sum(diag(S %*% Theta)) - as.numeric(determinant(Theta,logarithm=TRUE)$modulus)
	}
	v <- function(z) lambda * sum(z)
	w <- function(beta) {
		toCholesky(p,beta)
	}
	list(u=u,v=v,w=w)
}

#' Auxiliary function which returns the objective, penalty, and dependence structure among regression coefficients of the Lasso for polygenic risk scores (prs).
#' 
#' @param X The design matrix.
#' @param y The response vector.
#' @param s The shrinkage parameter used to regularize the design matrix.
#' @param lambda The regularization parameter of the prs Lasso.
#' 
#' @return A list with six functions, precisely the objective \eqn{u}, penalty \eqn{v}, and dependence structure \eqn{w}, as well as their derivatives \eqn{du}, \eqn{dv}, and \eqn{dw}.
#' 
#' @importFrom Rdpack reprompt
#' @references Mak, T.S., Porsch, R.M., Choi, S.W., Zhou, X., and Sham, P.C. (2017). Polygenic scores via penalized regression on summary statistics. Genet Epidemiol, 41(6):469-480.
#' @references Mak, T.S. and Porsch, R.M. (2020). lassosum: LASSO with summary statistics and a reference panel. R package version 0.4.5.
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' s <- 0.5
#' lambda <- 1
#' temp <- prsLasso(X,y,s,lambda)
#' 
#' @export
prsLasso <- function(X,y,s,lambda) {
	p <- ncol(X)
	r <- as.vector(t(X) %*% y)
	u <- function(beta) (1-s)*sum((X %*% beta)**2) - 2*sum(beta*r) + s*sum(beta**2)
	v <- function(z) 2 * lambda * sum(z)
	w <- function(beta) beta
	du <- function(beta) (1-s)*2*as.vector((t(X) %*% X) %*% beta) - 2*r + 2*s*beta
	dv <- function(z) rep(2*lambda,p)
	dw <- function(beta) diag(p)
	list(u=u,v=v,w=w,du=du,dv=dv,dw=dw)
}





# *****************************************************************************
# Definition of the unsmoothed objective function and its (non-smooth) gradient
# *****************************************************************************

#' Auxiliary function to define the objective function of an L1 penalized regression operator.
#' 
#' @param betavector The vector of regression coefficients.
#' @param u The function encoding the objective of the regression operator.
#' @param v The function encoding the penalty of the regression operator.
#' @param w The function encoding the dependence structure among the regression coefficients.
#' 
#' @return The value of the L1 penalized regression operator for the input \eqn{betavector}.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' lambda <- 1
#' temp <- standardLasso(X,y,lambda)
#' print(objFunction(betavector,temp$u,temp$v,temp$w))
#' 
#' @export
objFunction <- function(betavector,u,v,w) {
	u(betavector) + v(abs(w(betavector)))
}

#' Auxiliary function which computes the (non-smooth) gradient of an L1 penalized regression operator.
#' 
#' @param betavector The vector of regression coefficients.
#' @param w The function encoding the dependence structure among the regression coefficients.
#' @param du The derivative (gradient) of the objective of the regression operator.
#' @param dv The derivative (gradient) of the penalty of the regression operator.
#' @param dw The derivative (Jacobian matrix) of the function encoding the dependence structure among the regression coefficients.
#' 
#' @return The value of the gradient for the input \eqn{betavector}.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' lambda <- 1
#' temp <- standardLasso(X,y,lambda)
#' print(objFunctionGradient(betavector,temp$w,temp$du,temp$dv,temp$dw))
#' 
#' @export
objFunctionGradient <- function(betavector,w,du,dv,dw) {
	temp <- abs(w(betavector))
	du(betavector) + as.vector( t(dw(betavector)) %*% (sign(w(betavector))*dv(temp)) )
}





# *****************************************************************************
# Definition of the smoothed objective function and its gradient
# *****************************************************************************

#' Auxiliary function to define the objective function of the smoothed L1 penalized regression operator.
#' 
#' @param betavector The vector of regression coefficients.
#' @param u The function encoding the objective of the regression operator.
#' @param v The function encoding the penalty of the regression operator.
#' @param w The function encoding the dependence structure among the regression coefficients.
#' @param mu The Nesterov smoothing parameter.
#' @param entropy A boolean switch to select the entropy prox function (default) or the squared error prox function.
#' 
#' @return The value of the smoothed regression operator for the input \eqn{betavector}.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' lambda <- 1
#' temp <- standardLasso(X,y,lambda)
#' print(objFunctionSmooth(betavector,temp$u,temp$v,temp$w,mu=0.1))
#' 
#' @export
objFunctionSmooth <- function(betavector,u,v,w,mu,entropy=TRUE) {
	p <- c(-1,+1)
	temp <- nesterovPWA(p,w(betavector),mu,entropy)
	u(betavector) + v(temp)
}

#' Auxiliary function which computes the gradient of the smoothed L1 penalized regression operator.
#' 
#' @param betavector The vector of regression coefficients.
#' @param w The function encoding the dependence structure among the regression coefficients.
#' @param du The derivative (gradient) of the objective of the regression operator.
#' @param dv The derivative (gradient) of the penalty of the regression operator.
#' @param dw The derivative (Jacobian matrix) of the function encoding the dependence structure among the regression coefficients.
#' @param mu The Nesterov smoothing parameter.
#' @param entropy A boolean switch to select the entropy prox function (default) or the squared error prox function.
#' 
#' @return The value of the gradient for the input \eqn{betavector}.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' lambda <- 1
#' temp <- standardLasso(X,y,lambda)
#' print(objFunctionSmoothGradient(betavector,temp$w,temp$du,temp$dv,temp$dw,mu=0.1))
#' 
#' @export
objFunctionSmoothGradient <- function(betavector,w,du,dv,dw,mu,entropy=TRUE) {
	p <- c(-1,+1)
	temp <- nesterovArgumentGradient(p,w(betavector),mu,entropy)
	du(betavector) + as.vector( t(dw(betavector)) %*% (temp*dv(temp)) )
}





# *****************************************************************************
# Minimization of the unsmoothed and smoothed objective functions
# 
# Additionally, the progressive smoothing procedure is defined here.
# *****************************************************************************

#' Minimize the objective function of an unsmoothed or smoothed regression operator with respect to \eqn{betavector} using BFGS.
#' 
#' @param p The dimension of the unknown parameters (regression coefficients).
#' @param obj The objective function of the regression operator as a function of \eqn{betavector}.
#' @param objgrad The gradient function of the regression operator as a function of \eqn{betavector}.
#' 
#' @return The estimator \eqn{betavector} (minimizer) of the regression operator.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' lambda <- 1
#' temp <- standardLasso(X,y,lambda)
#' obj <- function(z) objFunctionSmooth(z,temp$u,temp$v,temp$w,mu=0.1)
#' objgrad <- function(z) objFunctionSmoothGradient(z,temp$w,temp$du,temp$dv,temp$dw,mu=0.1)
#' print(minimizeFunction(p,obj,objgrad))
#' 
#' @export
minimizeFunction <- function(p,obj,objgrad) {
	optim(par=runif(p), fn=obj, gr=objgrad, method="BFGS")$par
}

#' Minimize the objective function of a smoothed regression operator with respect to \eqn{betavector} using the progressive smoothing algorithm.
#' 
#' @param p The dimension of the unknown parameters (regression coefficients).
#' @param obj The objective function of the regression operator. Note that in the case of the progressive smoothing algorithm, the objective function must be a function of both \eqn{betavector} and \eqn{mu}.
#' @param objgrad The gradient function of the regression operator. Note that in the case of the progressive smoothing algorithm, the gradient must be a function of both \eqn{betavector} and \eqn{mu}.
#' @param muSeq The sequence of Nesterov smoothing parameters. The default is \eqn{2^{-n}} for \eqn{n \in \{-3,\ldots,6\}}.
#' 
#' @return The estimator \eqn{betavector} (minimizer) of the regression operator.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Lutz, S., Laha, N., and Lange, C. (2020). A framework to efficiently smooth L1 penalties for linear regression. bioRxiv:2020.09.17.301788.
#' 
#' @examples
#' library(smoothedLasso)
#' n <- 100
#' p <- 500
#' betavector <- runif(p)
#' X <- matrix(runif(n*p),nrow=n,ncol=p)
#' y <- X %*% betavector
#' lambda <- 1
#' temp <- standardLasso(X,y,lambda)
#' obj <- function(z,m) objFunctionSmooth(z,temp$u,temp$v,temp$w,mu=m)
#' objgrad <- function(z,m) objFunctionSmoothGradient(z,temp$w,temp$du,temp$dv,temp$dw,mu=m)
#' print(minimizeSmoothedSequence(p,obj,objgrad))
#' 
#' @export
minimizeSmoothedSequence <- function(p,obj,objgrad,muSeq=2**seq(3,-6)) {
	res <- runif(p)
	for(mu in muSeq) {
		res <- optim(par=res, fn=function(z) obj(z,mu), gr=function(z) objgrad(z,mu), method="BFGS")$par
	}
	return(res)
}
