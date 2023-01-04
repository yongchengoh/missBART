#' Sampling once from a MVN(mu, Q^{-1}) distribution (Keefe)
#'
#' @param mu mean vector
#' @param Q precision matrix
#'
#' @return One sample from a MVN(mu, Q^{-1}) distribution
#' @export
#'
#' @examples rMVN(c(0,0), diag(1,2))
rMVN <- function(mu, Q) {
  p  <- NCOL(Q)
  z  <- stats::rnorm(p)
  if(p == 1) mu + z/sqrt(Q) else mu + backsolve(PD_chol(Q), z)
}

#keefe: some versions below use this function to draw from MVN(bQ^{-1},Q^{-1}) a little faster than MVN(mu=bQ^{-1},Q^{-1})
# note that b being a matrix is allowed as per rMVN (removes loops in update_y_miss_old_matrix!)
#' Sampling from MVN(bQ^{-1},Q^{-1})
#'
#' @param b a vector
#' @param Q precision matrix
#' @param is_chol logical
#'
#' @return One sample from a MVN(bQ^{-1},Q^{-1}) distribution
#' @export
#'
#' @examples rMVNc(c(0,0), diag(1,2))
rMVNc   <- function(b, Q, is_chol = FALSE) {
  p     <- NCOL(Q)
  Z     <- if(is.matrix(drop(b))) matrix(stats::rnorm(p * ncol(b)), nrow=p) else stats::rnorm(p)
  if(p  == 1) {
    U   <- drop(if(isTRUE(is_chol)) Q else sqrt(Q))
    drop((b/U + Z)/U)
  } else     {
    U   <- if(isTRUE(is_chol)) Q else chol(Q)
    backsolve(U, backsolve(U, b, transpose=TRUE, k=p) + Z, k=p)
  }
}

#' Sampling multiple times from a multivariate normal distribution
#'
#' @param mean_mat matrix containing n mean vectors
#' @param precision precision matrix
#'
#' @return a matrix where each row corresponds to a draw from a MVN distribution with mean equal to the rows of \code{mean_mat} and precision \code{precision}
#' @export
#'
#' @examples multi_rMVN(matrix(0, ncol=2, nrow=10), diag(1,2))
multi_rMVN = function(mean_mat, precision){
  n = nrow(mean_mat)
  p = ncol(mean_mat)
  Y = matrix(nrow = n, ncol = p)
  for(i in 1:n){
    Y[i,] = rMVN(mean_mat[i,], precision)
  }
  return(Y)
}

PD_chol <- function(x, ...) tryCatch(chol(x, ...), error=function(e) {
  chol(Matrix::nearPD(x, base.matrix=TRUE)$mat, ...)
})

#' Metropolis-Hastings sampler
#'
#' @param p1 scalar or vector
#' @param l1 scalar or vector
#' @param p2 scalar or vector
#' @param l2 scalar or vector
#' @param q1 scalar or vector
#' @param q2 scalar or vector
#'
#' @return scalar or vector for accepting/rejecting proposed values
#' @export
#'
#' @examples
MH <- function(p1, l1, p2, l2, q1 = 0, q2 = 0){
  ratio = (p2 + l2 + q2) - (p1 + l1 + q1)
  accept = ratio >= 0 || - stats::rexp(1L) <= ratio
  return(accept)
}

#' Sampling from a matrix normal distribution MV(mu, U, V)
#' @param mu mean matrix (real nxp matrix)
#' @param U scale matrix (positive-definite real nxn matrix)
#' @param V scale matrix (positive-definite real pxp matrix)
#'
#' @return one sample from a MV(mu, U, V) distribution
#' @export
#'
#' @examples matrnorm(matrix(0, nrow=5, ncol=2), diag(1,5), diag(1,2))
matrnorm <- function(mu, U, V) {
  r <- NCOL(U)
  p <- NCOL(V)
  X <- matrix(stats::rnorm(r * p), nrow=r)
  X <- mu + crossprod(PD_chol(U), X) %*% (if(p == 1) sqrt(V) else PD_chol(V))
  return(X)
}

#' Sampling mu at tree terminal node level
#' @param node_partial_res matrix of partial residuals from trees
#' @param kappa hyperparameter
#' @param omega precision matrix
#'
#' @return a single mu (scalar or vector) at terminal node level
#' @export
#'
#' @examples
sim_mu_node = function(node_partial_res, kappa, omega){
  n = nrow(node_partial_res)
  p = ncol(node_partial_res)
  Q = n*omega + diag(kappa, p)
  b = omega %*% colSums(node_partial_res)
  return(rMVNc(b, Q))
}

#' MCMC sample of mu for all terminal nodes of a BART tree
#'
#' @param change_points vector
#' @param partial_res matrix
#' @param kappa hyperparameter
#' @param omega precision matrix
#'
#' @return mu for all terminal nodes of a tree
#' @export
#'
#' @examples
sim_mu <- function(change_points, partial_res, kappa, omega) {
  return(matrix((sapply(split.data.frame(partial_res, change_points), sim_mu_node, kappa = kappa, omega = omega)), nrow = length(unique(change_points)), byrow = TRUE))
}

#' MCMC sample of BART precision matrix
#'
#' @param y data
#' @param y_hat sum of BART outputs
#' @param alpha hyperparameter
#' @param Vinv hyperparameter
#'
#' @return one MCMC sample for BART precision matrix
#' @export
#'
#' @examples
sim_omega = function(y, y_hat, alpha, Vinv){
  n = nrow(y)
  df = alpha + n
  scale = crossprod(y-y_hat) + Vinv
  return(stats::rWishart(1, df, scale)[,,1])
}

#' MCMC sample of kappa
#'
#' @param mu list of mu's from all BART trees in the current iteration
#' @param a hyperparameter
#' @param b hyperparameter
#'
#' @return one MCMC sample for kappa
#' @export
#'
#' @examples
sim_kappa = function(mu, a, b){ #tree_mu[[i]]
  mu_mat = Reduce(rbind, mu)
  p = ncol(mu_mat)
  n_mu = nrow(mu_mat)
  shape = n_mu*p/2
  rate = b + sum(apply(mu_mat, 1, crossprod))/2
  return(stats::rgamma(1, shape, rate))
}

###--------------------------------- PROBIT UPDATES ---------------------------------###
#' MCMC sample of the n latent variables in the probit model (z)
#'
#' @param Y matrix of probit covariates
#' @param m matrix of missing data indicators
#' @param B matrix of probit parameters
#' @param R correlation matrix
#'
#' @return one MCMC sample of z
#' @export
#'
#' @examples
update_z = function(Y, m, B, R){
  n = nrow(Y)
  mu = Y %*% B
  U = diag(1, n)
  z = matrnorm(mu, U, R)
  z[intersect(which(z<0), which(m==1))] = 0
  z[intersect(which(z>=0), which(m==0))] = 0
  return(z)
}

#' MCMC sample of B, the matrix of probit parameters
#'
#' @param x regression covaraites
#' @param y data
#' @param z latent variables
#' @param tau_b hyperparameter
#' @param include_x logical
#' @param include_y logical
#' @param center logical
#'
#' @return one MCMC sample of B
#' @export
#'
#' @examples
update_B = function(x, y, z, tau_b = 0.01, include_x = TRUE, include_y = TRUE, center = TRUE){
  if(center){
    # x = scale(x, center = TRUE, scale = FALSE)
    y = scale(y, center = TRUE, scale = FALSE)
  }
  n = nrow(y)
  Y = matrix(1,nrow=n, ncol=1)
  if(include_x & !include_y){
    Y = cbind(Y, x)
  } else if (!include_x & include_y){
    Y = cbind(Y, y)
  } else if (include_x & include_y){
    Y = cbind(Y, x, y)
  }
  k = ncol(Y)
  p = ncol(y)
  V = diag(p)
  U = chol2inv(PD_chol(diag(tau_b, k) + t(Y)%*%Y))
  M = t(t(z) %*% Y %*% U)
  return(matrnorm(M, U, V))
}

#' MCMC sample of W, the latent variable in probit regression introduced by \cite{Talhouk, A., Doucet, A., & Murphy, K. (2012). Efficient Bayesian inference for multivariate probit models with sparse inverse correlation matrices. Journal of Computational and Graphical Statistics, 21(3), 739-757.}
#'
#' @param R correlation matrix
#' @param z probit latent variables
#'
#' @return one MCMC sample of W
#' @export
#'
#' @examples
#' @importFrom extraDistr "rinvgamma"
update_W = function(R, z){
  p = ncol(R)
  R_inv = chol2inv(PD_chol(R))
  shape = (p+1)/2
  scale_vec = diag(R_inv)/2
  d2 = rinvgamma(p, shape, scale_vec)
  D = diag(sqrt(d2),p)
  W = z %*% D
  return(W)
}

#' MCMC sample of \eqn{\Sigma} and \eqn{\gamma}, the latent variables in probit regression introduced by \cite{Talhouk, A., Doucet, A., & Murphy, K. (2012). Efficient Bayesian inference for multivariate probit models with sparse inverse correlation matrices. Journal of Computational and Graphical Statistics, 21(3), 739-757.}
#'
#' @param p scalar
#' @param Y matrix
#' @param Psi matrix
#' @param W matrix
#'
#' @return one MCMC sample of \eqn{\Sigma} and \eqn{\gamma}
#' @export
#'
#' @examples
#' @importFrom CholWishart "rInvWishart"
update_sigma_gamma = function(p, Y, Psi, W){
  n = nrow(Y)
  r = ncol(Y)
  Psi_inv = chol2inv(PD_chol(Psi))
  Xi_inv = t(Y) %*% Y + Psi_inv
  Xi = chol2inv(PD_chol(Xi_inv))
  M = Xi %*% t(Y) %*% W
  df = 2 + n
  df = df - p + 1
  scale = t(W)%*%W + diag(1,p) - t(M) %*% Xi_inv %*% M
  sim_Sigma = rInvWishart(1, df, scale)[,,1]
  sim_gamma = matrnorm(M, Xi, sim_Sigma)
  return(list(Sigma = sim_Sigma, gamma = sim_gamma))
}

#' MCMC sample of R (correlation matrix) and B (probit parameters)
#'
#' @param Sigma matrix
#' @param gamma matrix
#'
#' @return one MCMC sample of R and B
#' @export
#'
#' @examples
update_RB = function(Sigma, gamma){
  if(is.matrix(Sigma)){
    D_inv = diag(1/sqrt(diag(Sigma)))
  } else {
    D_inv = 1/sqrt(Sigma)
  }
  # D_inv = chol2inv(PD_chol(D))
  R = D_inv %*% Sigma %*% D_inv
  B = gamma %*% D_inv
  return(list(R = R, B = B, D = solve(D_inv)))
}

#' MCMC sample of \eqn{\Psi} from \cite{Talhouk, A., Doucet, A., & Murphy, K. (2012). Efficient Bayesian inference for multivariate probit models with sparse inverse correlation matrices. Journal of Computational and Graphical Statistics, 21(3), 739-757.}
#'
#' @param nu degrees of freedom
#' @param S scale matrix
#' @param B matrix of probit parameters
#' @param R correlation matrix
#'
#' @return one MCMC sample of \eqn{\Psi}
#' @export
#'
#' @examples
#' @importFrom CholWishart "rInvWishart"
update_Psi = function(nu, S, B, R){
  r = nrow(B)
  df = nu + r
  Sigma = S + B %*% chol2inv(PD_chol(R)) %*% t(B)
  return(rInvWishart(1, df, Sigma)[,,1])
}

#' MCMC sample of missing data
#'
#' @param x regression covariates
#' @param y_hat sum of BART outputs
#' @param m matrix of missing data indicator
#' @param y data
#' @param z probit regression latent variables
#' @param B matrix of probit parameters
#' @param R correlation matrix
#' @param omega precision matrix
#' @param include_x logical
#' @param include_y logical
#'
#' @return a matrix of imputed values
#' @export
#'
#' @examples
update_y_miss_reg <- function(x, y_hat, m, y, z, B, R, omega, include_x = TRUE, include_y = TRUE) {
  M = apply(m, 2, as.logical)
  n <- nrow(y_hat)
  p <- ncol(y_hat)
  q <- ncol(x)
  if(include_x && !include_y) {
    A   <- B[2:(q+1),, drop=FALSE]
    B1  <- matrix(0, nrow=p, ncol=p)
  } else if(!include_x && include_y) {
    A   <- matrix(0, nrow=q, ncol=p)
    B1  <- B[-1,, drop=FALSE]
  } else if(include_x  && include_y) {
    A   <- B[2:(q+1),, drop=FALSE]
    B1  <- B[(q+2):(p+q+1),, drop=FALSE]
  } else {
    A   <- matrix(0, nrow=q, ncol=p)
    B1  <- matrix(0, nrow=p, ncol=p)
  }
  b1    <- t(B[1,,drop=FALSE])
  prec  <- omega + if(p == 1) tcrossprod(B1)/R else B1 %*% solve(R, t(B1))
  Q     <- if(p == 1) sqrt(prec) else PD_chol(prec)
  Oyhat <- tcrossprod(omega, y_hat)
  zXA   <- if(include_x) z - x %*% A else z
  BRinv <- B1 %*% chol2inv(PD_chol(R))
  b     <- Oyhat + tcrossprod(BRinv, sweep(zXA, 2, b1, FUN="-", check.margin=FALSE))
  Y     <- rMVNc(b, Q, is_chol=TRUE)
  Y     <- if(p == 1) matrix(Y, ncol=p) else t(Y)
  return(Y)
}
