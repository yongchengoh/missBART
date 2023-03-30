#' Title
#'
#' @param x covariates
#' @param y response
#' @param x_predict out-of-sample covariates. If not specificied, the default is set to NA and no out-of-sample predictions will be made.
#' @param n_reg_trees number of BART trees
#' @param n_class_trees number of probit BART trees
#' @param burn burn-in samples
#' @param iters post-burn-in samples
#' @param thin thinning
#' @param predict make out-of-sample predictions?
#' @param tree_prior_params prior parameters for BART trees
#' @param hypers prior parameters for BART parameters
#' @param scale scale data?
#' @param include_x Include x in probit model?
#' @param include_y Include y in probit model?
#' @param show_progress logical
#' @param progress_every integer value stating how often to update the progress bar.
#' @param ... Catches unused arguments
#'
#' @return
#' @export
#'
#' @examples
missBART2 = function(x, y, x_predict = NA, n_reg_trees = 20, n_class_trees = 20, burn = 100, iters = 100, thin = 2, predict = TRUE,
                    tree_prior_params = tree_list(), hypers = hypers_list(),
                    scale = TRUE, include_x = TRUE, include_y = TRUE, show_progress = TRUE, progress_every = 10, ...) {


  y = as.matrix(y)
  x = as.matrix(x)
  for(l in 1:ncol(x)){
    if(any(is.na(x[,l]))){
      x = cbind(x, 1-as.integer(is.na(x[,l])))
    }
  }

  if(is.na(x_predict)) predict = FALSE

  missing_index = which(is.na(y))
  obs_index = which(!is.na(y))
  miss_row = apply(y, 1, function(x) any(is.na(x)))

  if(scale){
    min_y = apply(y, 2, min, na.rm = TRUE) #min(y, na.rm = TRUE)
    max_y = apply(y, 2, max, na.rm = TRUE) #max(y, na.rm = TRUE)
    y = t(apply(sweep(y, 2, min_y), 1, function(x) x/(max_y-min_y))) - 0.5
    if(nrow(y)==1) y = t(y)
    if(predict) y_predict = t(apply(sweep(y_predict, 2, min_y), 1, function(x) x/(max_y-min_y))) - 0.5
  }

  #####-------------------- GET PARAMETERS --------------------#####
  p = ncol(y) # No. of y variables
  q = ncol(x) # No. of x variables
  r = p*include_y + q*include_x
  Y_vars = seq_len(r)
  n = nrow(y)
  thinned = iters
  total_iters = burn + thin*iters

  #####-------------------- GET BART PRIOR PARAMETERS --------------------#####
  mu0 = rep(hypers$mu0, p)
  alpha = max(5, p + hypers$alpha)
  V = diag(hypers$V, p) #diag(1/apply(y, 2, sd, na.rm=TRUE), p)
  Vinv = solve(V)
  # Psi = diag(1, r)
  # kappa_a_reg = kappa_a_class = 16
  # kappa_b_reg = 1/n_reg_trees
  # kappa_b_class = 1/n_class_trees

  #####-------------------- GET TREE PRIOR PARAMETERS --------------------#####
  prior_alpha = tree_prior_params$prior_alpha
  prior_beta = tree_prior_params$prior_beta
  min_node = max(tree_prior_params$min_node, p+1)
  max_attempt = tree_prior_params$max_attempt

  #####-------------------- CREATE STORAGE FOR REGRESSION BART --------------------#####
  accepted_reg_trees = lapply(vector(mode = "list", length = n_reg_trees), as.list)
  reg_change_id = vector(mode = "list", length = n_reg_trees)

  reg_mu = vector(mode = "list", length = thinned)
  reg_phi = vector(mode = "list", length = n_reg_trees)
  omega_post = vector(mode = "list", length = 0)

  reg_prior = lapply(vector(mode = "list", length = n_reg_trees), as.list) # tree prior for accepted trees
  reg_likely = lapply(vector(mode = "list", length = n_reg_trees), as.list) # likelihood for accepted trees
  reg_accept = lapply(vector(mode = "list", length = n_reg_trees), as.list) # accept/reject status for all trees: reg_accept[[j]][[i]]
  reg_moves = NULL

  #####-------------------- CREATE STORAGE FOR PROBIT BART --------------------#####
  m = matrix(1, nrow=n, ncol=p)
  m[is.na(y)] = 0

  R_post = vector(mode = "list", length = 0)
  y_post = vector(mode = "list", length = 0)
  y_pred = vector(mode = "list", length = 0)

  accepted_class_trees = lapply(vector(mode = "list", length = n_class_trees), as.list)
  class_change_id = vector(mode = "list", length = n_class_trees)

  class_mu = vector(mode = "list", length = thinned)
  class_phi = vector(mode = "list", length = n_class_trees)
  omega_post = vector(mode = "list", length = 0)

  class_prior = lapply(vector(mode = "list", length = n_class_trees), as.list) # tree prior for accepted trees
  class_likely = lapply(vector(mode = "list", length = n_class_trees), as.list) # likelihood for accepted trees
  class_accept = lapply(vector(mode = "list", length = n_class_trees), as.list) # accept/reject status for all trees: class_accept[[j]][[i]]
  class_moves = NULL

  y_miss_accept = matrix(nrow = total_iters, ncol = length(missing_index)) #vector(mode = "list", length = total_iters) # Metropolis-Hastings acceptance/rejection for missing y's

  partial_res_y = matrix(nrow = n, ncol = p)
  partial_res_z = matrix(nrow = n, ncol = p)

  partial_reg = vector(mode = "list", length = n_reg_trees)

  #####-------------------- SET INITIAL VALUES FOR REGRESSION BART --------------------#####
  df = data.frame(matrix(ncol = 8, nrow = 1))
  colnames(df) = c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction", "NA_direction")
  df[1,] <- c(0,0,1,0,1,0,0,NA)

  accepted_reg_trees = lapply(seq_len(n_reg_trees), function(x) accepted_reg_trees[[x]] = df)
  reg_change_id = lapply(seq_len(n_reg_trees), function(x) reg_change_id[[x]] = rep(1, n))

  kappa_reg = 16*n_reg_trees
  kappa_class = 16*n_class_trees
  reg_mu[[1]] = sapply(seq_len(n_reg_trees), function(x) list(rMVN(mu = matrix(0, nrow=p), Q = kappa_reg*diag(p))))
  reg_phi = lapply(seq_len(n_reg_trees), function(x) rep(reg_mu[[1]][[x]], n))

  reg_prior = lapply(seq_len(n_reg_trees), function(x) reg_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  reg_accept = lapply(seq_len(n_reg_trees), function(x) reg_accept[x][[1]] = TRUE)

  #####-------------------- SET INITIAL VALUES FOR PROBIT MODEL --------------------#####
  accepted_class_trees = lapply(seq_len(n_class_trees), function(x) accepted_class_trees[[x]] = df)
  class_change_id = lapply(seq_len(n_class_trees), function(x) class_change_id[[x]] = rep(1, n))

  # new_omega = diag(1,p) #drop(rWishart(1, alpha, V)) #drop(rWishart(1, alpha, diag(1,p)))
  class_mu[[1]] = sapply(seq_len(n_class_trees), function(x) list(rMVN(mu = matrix(0, nrow=p), Q = kappa_class*diag(p))))
  class_phi = lapply(seq_len(n_class_trees), function(x) rep(class_mu[[1]][[x]], n))

  class_prior = lapply(seq_len(n_class_trees), function(x) class_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  class_accept = lapply(seq_len(n_class_trees), function(x) class_accept[x][[1]] = TRUE) #do the thing please

  new_R = diag(1, p)
  y[missing_index] = 0
  z = matrix(rep(1, n*p), nrow=n, ncol=p)
  z[missing_index] = -1

  Y = probit_predictors(x, y, include_x = include_x, include_y = include_y)

  new_omega = sim_omega(y = y, y_hat = Reduce("+", reg_phi), alpha = alpha, Vinv = Vinv)

  #####----- OUT-OF-SAMPLE PREDICTIONS -----#####
  if(predict){
    predict_change_id = vector(mode = "list", length = n_reg_trees)
    reg_phi_pred = vector(mode = "list", length = n_reg_trees)
    n_new = nrow(x_predict)
    new_y_post = vector(mode = "list", length = 0)
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

    #### -------------------- BART REGRESSION FOR DATA
    reg_mu[[i]] = vector(mode = "list", n_reg_trees)
    for(j in seq_len(n_reg_trees)) {
      ###----- Compute partial residuals -----###
      if (n_reg_trees==1) {
        partial_res_y = y
      } else {
        partial_res_y = y - Reduce("+", reg_phi[-j])
      }

      ###----- Set likelihood of first stump of each tree -----###
      if(i==1) reg_likely[[j]] = log_marginal_likelihood(node_partial_res = y, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      df = accepted_reg_trees[[j]]
      new_tree = propose_tree(df, x, min_node, max_attempt, i) # Propose new tree for tree j
      new_df = new_tree$new_df
      change_points = new_tree$change_points # Get change points for new tree
      decent_tree = new_tree$decent_tree
      # new_df = true_trees_data[[j]]
      # change_points = get_change_points(new_df, x)
      # decent_tree = TRUE
      # accept = TRUE

      accept = FALSE
      if(decent_tree){
        p2 = tree_priors(nodes = row.names(new_df), parents = unique(new_df$parent), depth = new_df$depth, prior_alpha, prior_beta)
        l2 = sum(sapply(split.data.frame(partial_res_y, change_points), log_marginal_likelihood, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        accept = MH(reg_prior[[j]], reg_likely[[j]], p2, l2)
      }

      if(accept){
        accepted_reg_trees[[j]] = new_df
        reg_prior[[j]] = p2
        reg_likely[[j]] = l2
        reg_change_id[[j]] = change_points
        partial_reg[[j]] = partial_res_y
      }

      ###----- Tree node updates -----###
      ##--Sample mu and compute phi, then store
      L_mu = sim_mu(reg_change_id[[j]], partial_res_y, kappa_reg, new_omega)
      L_phi = L_mu[reg_change_id[[j]],, drop=FALSE]
      reg_mu[[i]][[j]] = L_mu
      reg_phi[[j]] = L_phi # Used to calculate partial res.

      ###----- Out-of-Sample Predictions -----###
      if(predict){
        # all_points = get_change_points(accepted_reg_trees[[j]], rbind(x, x_predict))
        predict_change_id[[j]] = get_change_points(accepted_reg_trees[[j]], rbind(x, x_predict))[-c(1:n)] # assigns each x to a terminal node
        reg_phi_pred[[j]] = L_mu[predict_change_id[[j]],, drop=FALSE] # mu_j's for new x
      }

    } # End of j iterations

    ###----- BART updates -----###
    #--Get BART predictions
    y_hat = Reduce("+", reg_phi)

    #--Sample data precision
    new_omega = sim_omega(y, y_hat, alpha = alpha, Vinv = Vinv)
    kappa_reg = sim_kappa(mu = reg_mu[[i]], a = 16, b = 1/n_reg_trees)

    ###----- Probit BART -----###
    class_mu[[i]] = vector(mode = "list", n_class_trees)
    for(k in seq_len(n_class_trees)) {
      ###----- Compute partial residuals -----###
      if (n_class_trees==1) {
        partial_res_z = z
      } else {
        partial_res_z = z - Reduce("+", class_phi[-k])
      }

      ###----- Set likelihood of first stump of each tree -----###
      if(i==1) class_likely[[k]] = log_marginal_likelihood(node_partial_res = z, kappa = kappa_class, omega = new_R, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      df = accepted_class_trees[[k]]
      new_tree = propose_tree(df, Y, min_node, max_attempt, i, probit = TRUE, miss_row = miss_row) # Propose new tree for tree k
      new_df = new_tree$new_df
      change_points = new_tree$change_points # Get change points for new tree
      decent_tree = new_tree$decent_tree
      # new_df = true_trees_missing[[k]]
      # change_points = get_change_points(new_df, Y)
      # decent_tree = TRUE
      # accept = TRUE
      ###----- Metropolis-Hastings for accepting/rejecting proposed tree -----###
      accept = FALSE
      if(decent_tree){
        p2 = tree_priors(nodes = row.names(new_df), parents = unique(new_df$parent), depth = new_df$depth, prior_alpha, prior_beta)
        l2 = sum(sapply(split.data.frame(partial_res_z, change_points), log_marginal_likelihood, kappa = kappa_class, omega = new_R, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        l1 = sum(sapply(split.data.frame(partial_res_z, class_change_id[[k]]), log_marginal_likelihood, kappa = kappa_class, omega = new_R, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        accept = MH(class_prior[[k]], l1, p2, l2)
      }

      if(accept){
        accepted_class_trees[[k]] = new_df
        class_prior[[k]] = p2
        class_likely[[k]] = l2
        class_change_id[[k]] = change_points
      }

      ###----- Tree node updates -----###
      ##--Sample mu and compute phi, then store
      L_mu = sim_mu(class_change_id[[k]], partial_res_z, kappa_class, new_R)
      L_phi = L_mu[class_change_id[[k]],, drop=FALSE]
      class_mu[[i]][[k]] = L_mu
      class_phi[[k]] = L_phi # Used to calculate partial res.
    }

    #--Get classBART predictions
    z_hat = z = Reduce("+", class_phi)
    if(p==1){
      z = matrix(stats::rnorm(n, mean=z, sd=1), ncol=p, byrow=TRUE)
    } else {
      z = multi_rMVN(z, chol2inv(PD_chol(new_R)))
    }
    z[intersect(which(z<0), which(m==1))] = 0
    z[intersect(which(z>=0), which(m==0))] = 0
    kappa_class = sim_kappa(mu = class_mu[[i]], a = 16, b = 1/n_class_trees)

    if(include_y){ # If include_y==FALSE, then we are assuming MAR. Missing y's can be "imputed" from the model.
      #--Metropolis Hastings step for y_miss--#
      y_miss = update_y_miss_BART(x = x, y = y, z = z, z_hat = z_hat, y_hat = y_hat, n_trees = n_class_trees, R = new_R, Omega = new_omega, missing_index = missing_index, accepted_class_trees = accepted_class_trees, class_mu_i = class_mu[[i]], include_x = include_x, include_y = include_y, MH_sd = 0.1)
      y_miss_accept[i,] = y_miss$accept[missing_index]
      # y_miss = update_y_miss_BART(x = x, y = y, z = z, z_hat = z_hat, y_hat = y_hat, n_trees = n_class_trees, R = new_R, Omega = new_omega, missing_index = missing_index, accepted_class_trees = accepted_class_trees, class_mu_i = class_mu[[i]], include_x = include_x, include_y = include_y, MH_sd = 0.1)
      # y_miss_accept[i,] = rep(y_miss$accept, p)[missing_index]
      y[missing_index] = y_miss$y[missing_index]
    } else {
      y[missing_index] = multi_rMVN(y_hat, new_omega)[missing_index]
    }

    Y = probit_predictors(x, y, include_x = include_x, include_y = include_y, intercept = FALSE)

    ###----- Store posterior samples after burn-in, accounting for thinning -----###
    if(i > burn && i%%thin == 0){
      y_post = append(y_post, list(y_hat))
      pred = multi_rMVN(y_hat, new_omega)
      pred[missing_index] = y[missing_index]
      y_pred = append(y_pred, list(pred))
      omega_post = append(omega_post, list(new_omega))

      if(predict){
        new_pred_mean = Reduce("+", reg_phi_pred)
        new_predictions = multi_rMVN(new_pred_mean, new_omega)
        new_y_post = append(new_y_post, list(new_predictions))
      }
    }

  } # End of i iterations

  if(!predict) new_y_post = NA

  return(list(y_post = y_post, omega_post = omega_post,
              imputed = y, new_y_post = new_y_post, accepted_reg_trees = accepted_reg_trees, accepted_class_trees = accepted_class_trees,
              burn = burn, iters = iters, thin = thin,
              max_y = max_y, min_y = min_y,
              y_pred = y_pred))
}
