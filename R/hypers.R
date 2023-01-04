#' Title
#'
#' @param mu0 hyperparameter
#' @param kappa hyperparameter
#' @param alpha hyperparameter
#' @param V hyperparameter
#'
#' @return list of hyperparameters
#' @export
#'
#' @examples
hypers_list <- function(mu0 = 0, kappa = 1, alpha = 1, V = 1) {
  return(list(mu0 = mu0, kappa = kappa, alpha = alpha, V = V))
}

#' Title
#'
#' @param prior_alpha hyperparameter
#' @param prior_beta hyperparameter
#' @param min_node minimum number of observations that should fall into a single terminal node.
#' @param max_attempt maximum number of attempts to find a suitable tree
#'
#' @return list of hyperparameters
#' @export
#'
#' @examples
tree_list <- function(prior_alpha = 0.95, prior_beta = 2, min_node = 5, max_attempt = 100) {
  return(list(prior_alpha = prior_alpha, prior_beta = prior_beta, min_node = min_node, max_attempt = max_attempt))
}
