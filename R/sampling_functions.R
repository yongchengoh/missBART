#' Sampling once from a MVN(mu, Q^{-1}) distribution (Keefe)
#'
#' @param mu a mean vector
#' @param Q precision matrix
#'
#' @return One sample from a MVN(mu, Q^{-1}) distribution
#' @export
#'
#' @examples rMVN(mu, Q)
rMVN <- function(mu, Q) {
  p  <- NCOL(Q)
  z  <- rnorm(p)
  if(p == 1) mu + z/sqrt(Q) else mu + backsolve(PD_chol(Q), z)
}

#keefe: some versions below use this function to draw from MVN(bQ^{-1},Q^{-1}) a little faster than MVN(mu=bQ^{-1},Q^{-1})
# note that b being a matrix is allowed as per rMVN (removes loops in update_y_miss_old_matrix!)
#' Sampling from MVN(bQ^{-1},Q^{-1})
#'
#' @param b
#' @param Q
#' @param is_chol
#'
#' @return One sample from MVN(bQ^{-1},Q^{-1})
#' @export
#'
#' @examples rMVNc(b, Q)
rMVNc   <- function(b, Q, is_chol = FALSE) {
  p     <- NCOL(Q)
  Z     <- if(is.matrix(drop(b))) matrix(rnorm(p * ncol(b)), nrow=p) else rnorm(p)
  if(p  == 1) {
    U   <- drop(if(isTRUE(is_chol)) Q else sqrt(Q))
    drop((b/U + Z)/U)
  } else     {
    U   <- if(isTRUE(is_chol)) Q else chol(Q)
    backsolve(U, backsolve(U, b, transpose=TRUE, k=p) + Z, k=p)
  }
}

# Sampling multiple times from a MVN distribution
#' Sampling multiple times from a MVN distribution
#'
#' @param n
#' @param mean_mat mean matrix
#' @param precision precision
#'
#' @return a matrix
#' @export
#'
#' @examples
multi_rMVN = function(n, mean_mat, precision){
  p = ncol(mean_mat)
  Y = matrix(nrow = n, ncol = p)
  for(i in 1:n){
    Y[i,] = rMVN(mean_mat[i,], precision)
  }
  return(Y)
}

# Sampling once from a multivariate t-distribution (Keefe's code)
#' Sampling once from a multivariate t-distribution
#'
#' @param mu mean vector
#' @param nu scalar
#' @param Sigma matrix
#'
#' @return a vector
#' @export
#'
#' @examples
rMVt  <- function(mu, nu, Sigma) {
  p <- ncol(Sigma)
  Z <- rnorm(p)
  W <- sqrt(nu/rchisq(1, nu))
  if(p == 1) {
    U <- sqrt(Sigma)
    mu + W * U * Z
  } else {
    U <- PD_chol(Sigma)
    mu + t(W * U %*% Z)
  }
}

PD_chol <- function(x, ...) tryCatch(chol(x, ...), error=function(e) {
  chol(Matrix::nearPD(x, base.matrix=TRUE)$mat, ...)
})

# Metropolis Hastings
#' Title
#'
#' @param p1
#' @param l1
#' @param p2
#' @param l2
#' @param q1
#' @param q2
#'
#' @return
#' @export
#'
#' @examples
MH <- function(p1, l1, p2, l2, q1 = 0, q2 = 0){
  ratio = (p2 + l2 + q2) - (p1 + l1 + q1)
  accept = ratio >= 0 || - stats::rexp(1L) <= ratio
  return(accept)
}

#' Title
#'
#' @param mu
#' @param U
#' @param V
#'
#' @return
#' @export
#'
#' @examples
matrnorm <- function(mu, U, V) {
  r <- NCOL(U)
  p <- NCOL(V)
  X <- matrix(rnorm(r * p), nrow=r)
  X <- mu + crossprod(PD_chol(U), X) %*% (if(p == 1) sqrt(V) else PD_chol(V))
  return(X)
}

#' Title
#'
#' @param node_partial_res
#' @param kappa
#' @param omega
#'
#' @return
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

#' Title
#'
#' @param change_points
#' @param partial_res
#' @param kappa
#' @param omega
#'
#' @return
#' @export
#'
#' @examples
sim_mu <- function(change_points, partial_res, kappa, omega) {
  return(matrix((sapply(split.data.frame(partial_res, change_points), sim_mu_node, kappa = kappa, omega = omega)), nrow = length(unique(change_points)), byrow = TRUE))
}

#' Title
#'
#' @param y
#' @param y_hat
#' @param alpha
#' @param Vinv
#'
#' @return
#' @export
#'
#' @examples
sim_omega = function(y, y_hat, alpha, Vinv){
  n = nrow(y)
  df = alpha + n
  scale = crossprod(y-y_hat) + Vinv
  return(rWishart(1, df, scale)[,,1])
}

#' Title
#'
#' @param mu
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
sim_kappa = function(mu, a, b){ #tree_mu[[i]]
  mu_mat = Reduce(rbind, mu)
  p = ncol(mu_mat)
  n_mu = nrow(mu_mat)
  shape = n_mu*p/2
  rate = b + sum(apply(mu_mat, 1, crossprod))/2
  return(rgamma(1, shape, rate))
}

#' Title
#'
#' @param Y
#' @param m
#' @param B
#' @param R
#'
#' @return
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

###--------------------------------- PROBIT UPDATES ---------------------------------###
#' Title
#'
#' @param x
#' @param y
#' @param z
#' @param sigma
#' @param tau_b
#' @param include_x
#' @param include_y
#' @param center
#'
#' @return
#' @export
#'
#' @examples
update_B = function(x, y, z, sigma, tau_b, include_x = TRUE, include_y = TRUE, center = TRUE){
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
  # precision = chol2inv(PD_chol(V %x% U))
  # return(matrix(rMVN(as.vector(M), precision), nrow=k))
  return(matrnorm(M, U, V))
}

#' Title
#'
#' @param R
#' @param z
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom CholWishart "rinvgamma"
update_W = function(R, z){
  ## Update W = ZD. We need to sample D first, then compute W.
  p = ncol(R)
  R_inv = chol2inv(PD_chol(R))
  shape = (p+1)/2
  scale_vec = diag(R_inv)/2
  d2 = rinvgamma(p, shape, scale_vec) #1/sqrt(d_inv2)
  D = diag(sqrt(d2),p)
  W = z %*% D
  return(W)
}

#' Title
#'
#' @param p
#' @param Y
#' @param Psi
#' @param W
#' @param center
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom CholWishart "rInvWishart"
update_sigma_gamma = function(p, Y, Psi, W, center = FALSE){
  n = nrow(Y)
  # if(center){
  #   x = scale(x, center = TRUE, scale = FALSE)
  #   y = scale(y, center = TRUE, scale = FALSE)
  # }
  r = ncol(Y)
  Psi_inv = chol2inv(PD_chol(Psi))
  Xi_inv = t(Y) %*% Y + Psi_inv
  Xi = chol2inv(PD_chol(Xi_inv))
  M = Xi %*% t(Y) %*% W
  df = 2 + n
  ### Testing
  df = df - p + 1
  scale = t(W)%*%W + diag(1,p) - t(M) %*% Xi_inv %*% M
  # Sigma_inv = rWishart(1, df, chol2inv(PD_chol(scale)))[,,1]
  sim_Sigma = rInvWishart(1, df, scale)[,,1] #chol2inv(PD_chol(Sigma_inv))
  # Xi = chol2inv(PD_chol(Xi_inv))
  # Changing Xi %x% sim_Sigma to sim_Sigma %x% Xi
  sim_gamma = matrnorm(M, Xi, sim_Sigma)
  # Q = chol2inv(PD_chol(sim_Sigma %x% Xi))
  # sim_gamma = matrix(rMVN(mu = as.vector(M), Q = Q), nrow = r, ncol = p)
  return(list(Sigma = sim_Sigma, gamma = sim_gamma))
}

#' Title
#'
#' @param Sigma
#' @param gamma
#'
#' @return
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

#' Title
#'
#' @param nu
#' @param S
#' @param B
#' @param R
#'
#' @return
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

#' Title
#'
#' @param x
#' @param y_hat
#' @param m
#' @param y
#' @param z
#' @param B
#' @param R
#' @param omega
#' @param include_x
#' @param include_y
#'
#' @return
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
