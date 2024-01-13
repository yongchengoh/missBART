#' Title
#'
#' @param x covariates
#' @param y response
#' @param x_predict out-of-sample covariates
#' @param n_trees number of trees
#' @param burn burn-in samples
#' @param iters post burn-in samples
#' @param thin thinning
#' @param predict whether or not to make out-of-sample predictions
#' @param tree_prior_params prior parameters for BART trees
#' @param hypers prior parameters for BART parameters
#' @param scale logical. Whether to scale data to range (-0.5, 0.5).
#' @param show_progress logical.
#' @param progress_every integer value stating how often to update the progress bar.
#' @param true_trees_data true trees for BART component
#' @param ... Catches unused arguments
#'
#' @return An object of class \code{"BART"}.
#' @export
#'
#' @examples
#' # data <- sim_data_friedman(n = 100, p = 2)
#' # bart_out <- mvBART(data$x, data$y, n_trees = 90, burn = 500,
#' #                    iters = 1000, thin = 2, predict = FALSE)
mvBART <- function(x, y, x_predict = NA, n_trees = 100, burn = 1000, iters = 1000, thin = 2, predict = TRUE, tree_prior_params = tree_list(...), hypers = hypers_list(...),
                  scale = TRUE, show_progress = TRUE, progress_every = 10, true_trees_data = NA, ...) {

  if(is.null(x_predict)) predict <- FALSE
  y <- as.matrix(y)
  x <- as.matrix(x)

  for(l in 1:ncol(x)) {
    if(any(is.na(x[,l]))) {
      x <- cbind(x, 1 - as.integer(is.na(x[,l])))
      if(predict) {
        x_predict <- cbind(x_predict, 1 - as.integer(is.na(x_predict[,l])))
        colnames(x_predict) <- colnames(x)
      }
    }
  }

  min_y <- apply(y, 2, min, na.rm = TRUE)
  max_y <- apply(y, 2, max, na.rm = TRUE)
  if(scale) {
    y <- t(apply(sweep(y, 2, min_y), 1, function(x) x/(max_y - min_y))) - 0.5
    if(nrow(y) == 1) y <- t(y)
  }

  #####-------------------- GET PARAMETERS --------------------#####
  p <- ncol(y) # No. of y variables
  q <- ncol(x)
  x_vars <- seq_len(q)
  n <- nrow(y)
  thinned <- iters
  total_iters <- burn + thin * iters

  #####-------------------- GET BART PRIOR PARAMETERS --------------------#####
  mu0 <- rep(hypers$mu0, p)
  kappa <- 4 * (stats::qnorm(0.9))^2 * n_trees

  nu <- hypers$df
  qchi <- stats::qchisq(1 - hypers$q, nu)
  sigest <- rep(0, p)
  for(i in seq_len(p)) {
    sigest[i] <- summary(stats::lm(y[,i]~x))$sigma
  }
  lambda <- (sigest^2) * qchi/nu
  # print(paste("nu =", nu))
  # print(paste("lambda =", lambda))
  if(p == 1) {
    print(paste("sd_ols", (sigest * (max_y - min_y))^2))
    print(paste("E(sd_original_scale) =", (1/sqrt(1/lambda/(max_y - min_y)^2))^2))
  } else {
    print(paste("sd_ols", (sigest  *(max_y - min_y))^2))
    print(paste("E(sd_original_scale) =", (1/sqrt(1/lambda/(max_y - min_y)^2))^2))

    alpha <- ifelse(is.null(hypers$alpha), nu, hypers$alpha)
    if(is.null(hypers$V)) {
      V <- diag(1/(lambda * alpha), p)
      Vinv <- diag(lambda * alpha,  p)
    } else {
      V <- hypers$V
      Vinv <- solve(V)
    }
  }
  # for(j in seq_len(p)) {
  #   curve(dgamma(x, shape = nu/2, rate = nu * lambda[j]/2), from = 0, to = 100)
  # }

  #####-------------------- GET TREE PRIOR PARAMETERS --------------------#####
  prior_alpha <- tree_prior_params$prior_alpha
  prior_beta <- tree_prior_params$prior_beta
  min_node <- max(tree_prior_params$min_node, p + 1)
  max_attempt <- tree_prior_params$max_attempt

  #####-------------------- CREATE STORAGE FOR BART --------------------#####
  accepted_trees <- lapply(vector(mode = "list", length = n_trees), as.list)
  change_id <- vector(mode = "list", length = n_trees)

  tree_mu <- vector(mode = "list", length = thinned)
  tree_phi <- vector(mode = "list", length = n_trees)
  omega_post <- vector(mode = "list", length = 0)

  tree_prior <- lapply(vector(mode = "list", length = n_trees), as.list) # tree prior for accepted trees
  tree_likely <- lapply(vector(mode = "list", length = n_trees), as.list) # likelihood for accepted trees
  tree_accept <- lapply(vector(mode = "list", length = n_trees), as.list) # accept/reject status for all trees: tree_accept[[j]][[i]]
  tree_moves <- NULL

  partial_res <- matrix(nrow = n, ncol = p)

  #####-------------------- SET INITIAL VALUES FOR BART --------------------#####
  df <- data.frame(matrix(ncol = 8, nrow = 1))
  colnames(df) <- c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction", "NA_direction")
  df[1,] <- c(0,0,1,0,1,0,0,NA)

  accepted_trees <- lapply(seq_len(n_trees), function(x) accepted_trees[[x]] = df)
  change_id <- lapply(seq_len(n_trees), function(x) change_id[[x]] = rep(1, n))

  tree_mu[[1]] <- sapply(seq_len(n_trees), function(x) list(rMVN(mu = matrix(mu0, nrow=p), Q = diag(kappa, p))))
  tree_phi <- lapply(seq_len(n_trees), function(x) rep(tree_mu[[1]][[x]], n))

  tree_prior <- lapply(seq_len(n_trees), function(x) tree_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  tree_accept <- lapply(seq_len(n_trees), function(x) tree_accept[x][[1]] = TRUE)

  y_post <- vector(mode = "list", length = 0)
  y_pred <- vector(mode = "list", length = 0)

  new_omega <- diag(stats::rgamma(p, shape = nu/2, rate = nu * lambda/2), p)

  #####----- OUT-OF-SAMPLE PREDICTIONS -----#####
  if(predict) {
    predict_change_id <- vector(mode = "list", length = n_trees)
    tree_phi_pred <- vector(mode = "list", length = n_trees)
    n_new <- nrow(x_predict)
    new_y_post <- vector(mode = "list", length = 0)
  }

  #####----- PROGRESS BAR -----#####
  if(show_progress) {
    progress <- utils::txtProgressBar(min = 1, max = burn + iters * thin, style = 3, width = 60, title = 'Running rBART...')
  }

  #####-------------------- BART LOOP --------------------#####
  for(i in seq_len(total_iters)) {

    if(show_progress) {
      if(i == 1 || i %% progress_every == 0) {
        utils::setTxtProgressBar(progress, i)
      }
    }

    tree_mu[[i]] <- vector(mode = "list", n_trees)

    for(j in seq_len(n_trees)) {

      ###----- Compute partial residuals -----###
      if(n_trees == 1) {
        partial_res <- y
      } else {
        partial_res <- y - Reduce("+", tree_phi[-j])
      }

      ###----- Set likelihood of first stump of each tree -----###
      if(length(tree_likely[[j]]) == 0) tree_likely[[j]] <- log_marginal_likelihood(node_partial_res = y, kappa = kappa, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      df <- accepted_trees[[j]]
      new_tree <- propose_tree(df, x, min_node, max_attempt, i) # Propose new tree for tree j
      new_df <- new_tree$new_df
      change_points <- new_tree$change_points # Get change points for new tree
      decent_tree <- new_tree$decent_tree
      # new_df <- true_trees_data[[j]]
      # change_points <- get_change_points(new_df, x)
      # decent_tree <- TRUE
      # accept <- TRUE

      ###----- Metropolis-Hastings for accepting/rejecting proposed tree -----###
      if(decent_tree) {
        p2 <- tree_priors(new_df = new_df, prior_alpha, prior_beta)
        l2 <- sum(sapply(split.data.frame(partial_res, change_points),
                         log_marginal_likelihood, kappa = kappa, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        p1 <- tree_prior[[j]]

        l1 <- sum(sapply(split.data.frame(partial_res, change_id[[j]]),
                         log_marginal_likelihood, kappa = kappa, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)) #tree_likely[[j]]
        ratio <- (p2 + l2) - (p1 + l1)
        accept <- ratio >= 0 || - stats::rexp(1L) <= ratio
        tree_accept[[j]][[i]] <- accept
      } else {
        tree_accept[[j]][[i]] <- accept <- FALSE
      }

      ###----- Storing updated tree if accepted -----###
      if(accept) {
        accepted_trees[[j]] <- new_df
        tree_prior[[j]] <- p2
        tree_likely[[j]] <- l2
        change_id[[j]] <- change_points
      }

      ###----- Tree node updates -----###
      ##--Sample mu and compute phi, then store
      L_mu <- sim_mu(change_id[[j]], partial_res, kappa, new_omega)
      L_phi <- L_mu[change_id[[j]],, drop=FALSE]
      tree_mu[[i]][[j]] <- L_mu
      tree_phi[[j]] <- L_phi # Used to calculate partial res.

      ###----- Out-of-Sample Predictions -----###
      if(predict) {
        predict_change_id[[j]] <- get_change_points(accepted_trees[[j]], rbind(x, x_predict))[-seq_len(n)] # assigns each x to a terminal node
        tree_phi_pred[[j]] <- L_mu[predict_change_id[[j]],, drop=FALSE] # mu_j's for new x
      }

    } # End of j iterations

    ###----- BART updates -----###
    #--Get BART predictions
    y_hat <- Reduce("+", tree_phi)

    #--Sample data precision and kappa
    if(p == 1) {
      new_omega <- sim_omega(y = y, y_hat = y_hat, nu = nu, lambda = lambda)
    } else {
      new_omega <- sim_omega(y = y, y_hat = y_hat, alpha = alpha, Vinv = Vinv)
    }
    # kappa <- sim_kappa(tree_mu[[i]], kappa_a, kappa_b)

    ###----- Out-of-Sample Predictions -----###
    if(predict) {
      new_pred_mean <- Reduce("+", tree_phi_pred)
      new_predictions <- multi_rMVN(new_pred_mean, new_omega)
    }

    ###----- Store posterior samples after burn-in, accounting for thinning -----###
    if(i > burn && i %% thin == 0){
      y_post <- append(y_post, list(y_hat))
      omega_post <- append(omega_post, list(new_omega))
      if(predict) new_y_post <- append(new_y_post, list(new_predictions))

      pred <- multi_rMVN(y_hat, new_omega)
      y_pred <- append(y_pred, list(pred))
    }
  } # End of i iterations

  if(!predict) new_y_post <- NA

  mvBART_out <- structure(list(y_post = y_post, omega_post = omega_post,
                               new_y_post = new_y_post, accepted_trees = accepted_trees,
                               burn = burn, iters = iters, thin = thin,
                               max_y = max_y, min_y = min_y,
                               y_pred = y_pred, y = y, tree_mu = tree_mu),
                          class = "BART")
  return(mvBART_out)
}
