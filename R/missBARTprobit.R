#' BART and probit regression for data with non-ignorable missingness in the response.
#'
#' @param x BART covariates
#' @param y data
#' @param x_predict out-of-sample covariates. If not specificied, the default is set to NA and no out-of-sample predictions will be made.
#' @param n_trees number of BART trees. Default is set to 90.
#' @param burn burn-in samples. Default is set to 1000.
#' @param iters post-burn-in samples (after thinning). Default is set to 1000.
#' @param thin thinning. Default is set to 3.
#' @param predict whether or not to make out-of-sample predictions
#' @param tree_prior_params prior parameters for BART trees
#' @param hypers prior parameters for BART parameters
#' @param scale logical. Whether to scale data to range (-0.5, 0.5).
#' @param include_x logical. Include x in probit model?
#' @param include_y logical. Include y in probit model?
#' @param show_progress logical.
#' @param progress_every integer value stating how often to update the progress bar.
#' @param ... Catches unused arguments
#'
#' @return a list containing BART predictions and imputed values
#' @export
#'
#' @examples
#' x <- matrix(runif(6), ncol = 2)
#' y <- matrix(runif(6), ncol = 2) %*% matrix(rnorm(4), ncol=2)
#' missBARTprobit(x, y, n_trees = 2, burn = 2, iters = 2, thin = 1, scale = FALSE)
missBARTprobit = function(x, y, x_predict = c(), n_trees = 150, burn = 1000, iters = 1000, thin = 3, predict = TRUE, tree_prior_params = tree_list(), hypers = hypers_list(),
                           scale = TRUE, include_x = TRUE, include_y = TRUE,
                          show_progress = TRUE, progress_every = 10, pdp_range = c(-0.5, 0.5), make_pdp = FALSE, mice_impute = FALSE,
                          true_trees_data = NA, true_change_points = NA, ...) {


  if(is.null(x_predict)) predict = FALSE
  y = as.matrix(y)
  x = as.matrix(x)

  for(l in 1:ncol(x)){
    if(any(is.na(x[,l]))){
      x = cbind(x, 1-as.integer(is.na(x[,l])))
      if(predict) {
        x_predict = cbind(x_predict, 1-as.integer(is.na(x_predict[,l])))
        colnames(x_predict) = colnames(x)
      }
    }
  }

  min_y = apply(y, 2, min, na.rm = TRUE) #min(y, na.rm = TRUE)
  max_y = apply(y, 2, max, na.rm = TRUE) #max(y, na.rm = TRUE)
  if(scale){
    y = t(apply(sweep(y, 2, min_y), 1, function(x) x/(max_y-min_y))) - 0.5
    if(nrow(y)==1) y = t(y)
  }

  #####-------------------- GET PARAMETERS --------------------#####
  p = ncol(y) # No. of y variables
  q = ncol(x)
  r = 1 + p*include_y + q*include_x
  n = nrow(y)
  thinned = iters
  total_iters = burn + thin*iters

  missing_index = which(is.na(y))
  obs_index = which(!is.na(y))

  m = matrix(1, nrow=n, ncol=p)
  m[is.na(y)] = 0

  #####-------------------- GET BART PRIOR PARAMETERS --------------------#####
  mu0 <- rep(hypers$mu0, p)
  kappa = 16*n_trees
  # alpha <- p + 1 #max(5, p + hypers$alpha)
  # sample_t = 1/(apply(y, 2, sd, na.rm=TRUE))^2
  # V = -diag(sample_t, p)/(1.28*sqrt(2*max(4, alpha))-max(4, alpha)) #diag(1/(apply(y, 2, sd, na.rm=TRUE))^2/alpha, p)
  # V = diag(n_trees, p) #diag(1/(apply(y, 2, sd, na.rm=TRUE))^2/alpha, p)
  # Vinv = solve(V)
  # kappa_a = 16
  # kappa_b = 1/n_trees
  nu = hypers$df
  # lambda = (sqrt(2/nu)*qnorm(1-hypers$q) + 1)/sample_t
  qchi = qchisq(1-hypers$q, nu)
  sigest = rep(0, p)
  for(i in 1:p){
    sigest[i] = summary(lm(y[,i]~x))$sigma
  }
  sigest = summary(lm(y~x))$sigma
  lambda = (sigest^2)*qchi/nu
  print(lambda)

  Psi = rInvWishart(1, r+1, diag(1,r))[,,1]

  #####-------------------- GET TREE PRIOR PARAMETERS --------------------#####
  prior_alpha <- tree_prior_params$prior_alpha
  prior_beta <- tree_prior_params$prior_beta
  min_node <- max(tree_prior_params$min_node, p+1)
  max_attempt <- tree_prior_params$max_attempt

  #####-------------------- CREATE STORAGE FOR BART --------------------#####
  accepted_trees <- lapply(vector(mode = "list", length = n_trees), as.list)
  change_id <- vector(mode = "list", length = n_trees)

  reg_mu <- vector(mode = "list", length = thinned)
  tree_phi <- vector(mode = "list", length = n_trees)
  omega_post <- vector(mode = "list", length = 0)
  omega_post2 <- vector(mode = "list", length = 0)
  omega_post3 <- vector(mode = "list", length = 0)

  tree_prior <- lapply(vector(mode = "list", length = n_trees), as.list) # tree prior for accepted trees
  tree_likely <- lapply(vector(mode = "list", length = n_trees), as.list) # likelihood for accepted trees
  tree_accept <- lapply(vector(mode = "list", length = n_trees), as.list) # accept/reject status for all trees: tree_accept[[j]][[i]]
  tree_moves = NULL

  partial_res <- matrix(nrow = n, ncol = p)

  #####-------------------- SET INITIAL VALUES FOR BART --------------------#####
  df <- data.frame(matrix(ncol = 8, nrow = 1))
  colnames(df) = c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction", "NA_direction")
  df[1,] <- c(0,0,1,0,1,0,0,NA)

  accepted_trees <- lapply(seq_len(n_trees), function(x) accepted_trees[[x]] = df)
  change_id <- lapply(seq_len(n_trees), function(x) change_id[[x]] = rep(1, n))

  reg_mu[[1]] = sapply(seq_len(n_trees), function(x) list(rMVN(mu = matrix(mu0, nrow=p), Q = kappa*diag(p))))
  tree_phi = lapply(seq_len(n_trees), function(x) rep(reg_mu[[1]][[x]], n))

  tree_prior <- lapply(seq_len(n_trees), function(x) tree_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  # tree_accept <- lapply(seq_len(n_trees), function(x) tree_accept[x][[1]] = TRUE)

  #####-------------------- CREATE STORAGE AND SET INITIAL VALUES FOR PROBIT MODEL --------------------#####
  B_post = vector(mode = "list", length = 0)
  R_post = vector(mode = "list", length = 0)
  y_post = vector(mode = "list", length = 0)
  y_pred = vector(mode = "list", length = 0)
  y_impute = vector(mode = "list", length = 0)
  z_post = vector(mode = "list", length = 0)

  new_B = matrix(0, ncol=p, nrow=r) #matrix(rnorm(r*p, mean = 0, sd = 1/sqrt(tau_b)), ncol = p)
  new_R = diag(1, p)
  D = diag(1, p)

  if(mice_impute){
    imputed = mice::complete(mice::mice(cbind(y, x), print = FALSE))[,1:p]
    y[missing_index] = imputed[missing_index]
  } else {
    y[missing_index] = 0
  }
  z = matrix(rep(1, n*p), nrow=n, ncol=p)
  z[missing_index] = -1
  Y = probit_predictors(x, y, include_x = include_x, include_y = include_y, intercept = TRUE)

  # new_omega = sim_omega(y = y, y_hat = Reduce("+", tree_phi), alpha = alpha, Vinv = Vinv)
  new_omega = diag(rgamma(p, shape = nu/2, rate = nu*lambda/2), p)

  #####----- OUT-OF-SAMPLE PREDICTIONS -----#####
  if(predict){
    predict_change_id = vector(mode = "list", length = n_trees)
    tree_phi_pred = vector(mode = "list", length = n_trees)
    n_new = nrow(x_predict)
    new_y_post = vector(mode = "list", length = 0)
  }

  if(p==1 && make_pdp){
    pdp_out = vector(mode = "list", length = total_iters)
    if(include_x){
      pdp_list = pdp_param_mat_list(x = x, intercept = TRUE, y_range = pdp_range, include_x = include_x, n = n)
    } else {

    }
  }

  #####----- PROGRESS BAR -----#####
  if(show_progress) {
    progress = utils::txtProgressBar(min = 1, max = burn+iters*thin, style = 3, width = 60, title = 'Running rBART...')
  }

  #####-------------------- BART LOOP --------------------#####
  for(i in 1:total_iters) {

    if(show_progress){
      if(i==1 || i %% progress_every == 0){
        utils::setTxtProgressBar(progress, i)
      }
    }

    reg_mu[[i]] = vector(mode = "list", n_trees)

    for(j in seq_len(n_trees)) {

      ###----- Compute partial residuals -----###
      if (n_trees==1) {
        partial_res = y
      } else {
        partial_res = y - Reduce("+", tree_phi[-j])
      }

      ###----- Set likelihood of first stump of each tree -----###
      if(length(tree_likely[[j]]) == 0) tree_likely[[j]] = log_marginal_likelihood(node_partial_res = y, kappa = kappa, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      # df <- accepted_trees[[j]]
      # new_tree <- propose_tree(df, x, min_node, max_attempt, i) # Propose new tree for tree j
      # new_df <- new_tree$new_df
      # change_points <- new_tree$change_points # Get change points for new tree
      # decent_tree <- new_tree$decent_tree
      new_df = true_trees_data[[j]]
      change_points = get_change_points(new_df, x)
      decent_tree = TRUE
      accept = TRUE

      ###----- Metropolis-Hastings for accepting/rejecting proposed tree -----###
      # accept = FALSE
      # if(decent_tree) {
      #   p2 = tree_priors(nodes = row.names(new_df), parents = unique(new_df$parent), depth = new_df$depth, prior_alpha, prior_beta)
      #   l2 = sum(sapply(split.data.frame(partial_res, change_points), log_marginal_likelihood, kappa = kappa, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha))
      #   p1 = tree_prior[[j]]
      #   l1 = sum(sapply(split.data.frame(partial_res, change_id[[j]]), log_marginal_likelihood, kappa = kappa, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)) #tree_likely[[j]]
      #   ratio = (p2 + l2) - (p1 + l1)
      #   accept = ratio >= 0 || - stats::rexp(1L) <= ratio
      # }

      ###----- Storing updated tree if accepted -----###
      if(accept) {
        accepted_trees[[j]] <- new_df
        # tree_prior[[j]] <- p2
        # tree_likely[[j]] <- l2
        change_id[[j]] <- change_points
      }

      ###----- Tree node updates -----###
      ##--Sample mu and compute phi, then store
      L_mu <- sim_mu(change_id[[j]], partial_res, kappa, new_omega)
      L_phi <- L_mu[change_id[[j]],, drop=FALSE]
      reg_mu[[i]][[j]] <- L_mu
      tree_phi[[j]] <- L_phi # Used to calculate partial res.

      ###----- Out-of-Sample Predictions -----###
      if(predict){
        # predict_change_id[[j]] = get_change_points(accepted_trees[[j]], rbind(x, x_predict))[-c(1:n)] # assigns each x to a terminal node
        predict_change_id[[j]] = get_change_points(accepted_trees[[j]], t(cbind(t(x), t(x_predict))))[-c(1:n)]
        tree_phi_pred[[j]] = L_mu[predict_change_id[[j]],, drop=FALSE] # mu_j's for new x
      }

    } # End of j iterations

    ###----- BART updates -----###
    #--Get BART predictions
    y_hat = Reduce("+", tree_phi)

    #--Sample missing values
    if(include_y){
      y_miss = update_y_miss_reg(x = x, y_hat = y_hat, m = m, z = z, B = new_B, R = new_R, omega = new_omega, include_x = include_x)
    } else {
      y_miss = multi_rMVN(mean_mat = y_hat, precision = omega)[missing_index]
    }
    # y_miss = update_y_miss_reg(x = x, y_hat = y_hat, m = m, y = y, z = z, B = new_B, R = new_R, omega = new_omega, include_x = include_x, include_y = include_y)[missing_index]

    y[missing_index] = y_miss
    Y = probit_predictors(x, y, include_x = include_x, include_y = include_y, intercept = TRUE)

    #--Sample data precision and kappa
    # new_omega = sim_omega(y = y, y_hat = y_hat, alpha = alpha, Vinv = Vinv)
    new_omega = sim_omega(y = y, y_hat = y_hat, nu = nu, lambda = lambda)
    # kappa = sim_kappa(reg_mu[[i]], kappa_a, kappa_b)

    ###----- Probit updates -----###
    z = update_z(Y, m, new_B, new_R)

    if(p==1){
      new_R = diag(1, p)
      new_B = update_B(y, z, Y, tau_b = 0.01)
    } else {
      W = update_W(R = new_R, z = z)
      Sigma_gamma = update_sigma_gamma(p, Y, Psi, W)
      Sigma = Sigma_gamma$Sigma
      gamma = Sigma_gamma$gamma
      RB = update_RB(Sigma = Sigma, gamma = gamma)
      new_R = RB$R
      new_B = RB$B
      D = RB$D

      #--Sample Psi--#
      Psi = update_Psi(r+1, diag(1,r), new_B, new_R)
    }

    ###----- Out-of-Sample Predictions -----###
    if(predict){
      new_pred_mean = Reduce("+", tree_phi_pred)
      new_predictions = multi_rMVN(new_pred_mean, new_omega)
    }

    ###----- Store posterior samples after burn-in, accounting for thinning -----###
    if(i > burn && i%%thin == 0){
      if(scale){
        y_post = append(y_post, list(unscale(y_hat, min_y, max_y)))
        y_impute = append(y_impute, list(unscale(y[missing_index], min_y, max_y)))
        omega_post = append(omega_post, list(1/sqrt(new_omega/(max_y - min_y)^2)))
      } else {
        y_post = append(y_post, list((y_hat)))
        y_impute = append(y_impute, list(y[missing_index]))
        omega_post = append(omega_post, list(1/sqrt(new_omega)))
      }
      R_post = append(R_post, list(new_R))
      B_post = append(B_post, list(new_B))
      if(predict) new_y_post = append(new_y_post, list(unscale(new_predictions, min_y, max_y)))
      z_post = append(z_post, list(z))

      # pred = multi_rMVN(y_hat, new_omega)
      # pred[missing_index] = y[missing_index]
      # y_pred = append(y_pred, list(pred))
    }

    if(p==1 && make_pdp){
      pdp_out[[i]] = pnorm(Reduce(cbind, lapply(pdp_list, function(x) x %*% new_B)))
    }
  } # End of i iterations

  if(!predict) new_y_post = NA
  if(!make_pdp || p>1) pdp_out = NA

  return(list(y_post = y_post, omega_post = omega_post, R_post = R_post, B_post = B_post,
              new_y_post = new_y_post, accepted_trees = accepted_trees, y_impute = y_impute,
              burn = burn, iters = iters, thin = thin,
              reg_mu = reg_mu,
              max_y = max_y, min_y = min_y,
              y_pred = y_pred,
              z_post = z_post,
              pdp_out = pdp_out))
}

