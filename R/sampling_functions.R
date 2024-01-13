#' Sampling once from a MVN(mu, Q^{-1}) distribution (Keefe)
#'
#' @param mu mean vector
#' @param Q precision matrix
#'
#' @return One sample from a MVN(mu, Q^{-1}) distribution
#' @export
#'
#' @examples rMVN(c(0,0), diag(2))
rMVN <- function(mu, Q) {
  p  <- NCOL(Q)
  z  <- stats::rnorm(p)
  if(p == 1) mu + z/sqrt(Q) else mu + backsolve(PD_chol(Q), z, k=p)
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
#' @examples rMVNc(c(0,0), diag(2))
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
#' @examples multi_rMVN(matrix(0, ncol=2, nrow=10), diag(2))
multi_rMVN <- function(mean_mat, precision) {
  n <- nrow(mean_mat)
  p <- ncol(mean_mat)
  Y <- matrix(nrow = n, ncol = p)
  for(i in seq_len(n)) {
    Y[i,] <- rMVNmu0(precision)
  }
  return(Y + mean_mat)
}

#' Sampling from a matrix normal distribution MV(mu, U, V)
#' @param mu mean matrix (real nxp matrix)
#' @param U scale matrix (positive-definite real nxn matrix)
#' @param V scale matrix (positive-definite real pxp matrix)
#'
#' @return one sample from a MV(mu, U, V) distribution
#' @export
#'
#' @examples matrnorm(matrix(0, nrow=5, ncol=2), diag(5), diag(2))
matrnorm <- function(mu, U, V) {
  r <- NCOL(U)
  p <- NCOL(V)
  X <- matrix(stats::rnorm(r * p), nrow=r)
  X <- mu + crossprod(PD_chol(U), X) %*% (if(p == 1) sqrt(V) else PD_chol(V))
  return(X)
}
