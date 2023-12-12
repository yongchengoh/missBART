#' BART prior parameters
#'
#' @param mu0 hyperparameter
#' @param kappa hyperparameter
#' @param alpha hyperparameter
#' @param V hyperparameter
#'
#' @return list of hyperparameters
#' @export
#'
#' @examples hypers_list(mu0 = 5, kappa = 2, alpha = 1, V = 2)
hypers_list <- function(mu0 = 0, kappa = NULL, alpha = NULL, V = NULL, df = 10, q = 0.75) {
  return(list(mu0 = mu0, kappa = kappa, alpha = alpha, V = V, df = df, q = q))
}

#' Tree prior parameters
#'
#' @param prior_alpha hyperparameter
#' @param prior_beta hyperparameter
#' @param min_node minimum number of observations that should fall into a single terminal node.
#' @param max_attempt maximum number of attempts to find a suitable tree
#'
#' @return list of hyperparameters
#' @export
#'
#' @examples tree_list(prior_alpha = 0.95, prior_beta = 3)
tree_list <- function(prior_alpha = 0.95, prior_beta = 2, min_node = 1, max_attempt = 1) {
  if(prior_alpha >= 1 || prior_alpha <= 0) stop("prior_alpha must be between (0,1)")
  if(prior_beta < 0) stop("prior_beta must be positive")
  return(list(prior_alpha = prior_alpha, prior_beta = prior_beta, min_node = min_node, max_attempt = max_attempt))
}
