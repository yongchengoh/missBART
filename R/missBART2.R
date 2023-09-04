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
missBART2 = function(x, y, x_predict = c(), n_reg_trees = 150, n_class_trees = 50, burn = 1000, iters = 1000, thin = 3, predict = TRUE, MH_sd = 0.1,
                     tree_prior_params = tree_list(), hypers = hypers_list(),
                     scale = TRUE, include_x = TRUE, include_y = TRUE, show_progress = TRUE, progress_every = 10,
                     pdp_range = c(-0.5, 0.5), make_pdp = FALSE, mice_impute = FALSE, true_trees_data = NA, true_trees_missing = NA, z_true, true_change_points = NA, true_change_points_miss = NA, ...) {

  if(is.null(x_predict)) predict = FALSE
  y = as.matrix(y)
  x = as.matrix(x)
  if(!is.null(x_predict)) x_predict = as.matrix(x_predict)
  for(l in 1:ncol(x)){
    if(any(is.na(x[,l]))){
      x = cbind(x, 1-as.integer(is.na(x[,l])))
      if(predict) {
        x_predict = cbind(x_predict, 1-as.integer(is.na(x_predict[,l])))
        colnames(x_predict) = colnames(x)
      }
    }
  }

  missing_index = which(is.na(y))
  obs_index = which(!is.na(y))
  miss_row = apply(y, 1, function(x) any(is.na(x)))

  min_y = apply(y, 2, min, na.rm = TRUE) #min(y, na.rm = TRUE)
  max_y = apply(y, 2, max, na.rm = TRUE) #max(y, na.rm = TRUE)
  if(scale){
    y = t(apply(sweep(y, 2, min_y), 1, function(x) x/(max_y-min_y))) - 0.5
    if(nrow(y)==1) y = t(y)
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
  # alpha = p + 1
  # sample_t = 1/(apply(y, 2, sd, na.rm=TRUE))^2
  # V = diag(1/(summary(lm((y ~ x)))$sigma)^2, p)
  # V = -diag(sample_t, p)/(1.28*sqrt(2*max(4, alpha))-max(4, alpha)) #diag(1/(apply(y, 2, sd, na.rm=TRUE))^2/alpha, p)
  # Vinv = solve(V)

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
  z_post = vector(mode = "list", length = 0)
  y_impute = vector(mode = "list", length = 0)

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

  kappa_reg = 16*n_reg_trees  #2*sqrt(n_reg_trees)
  kappa_reg_list = c()

  reg_mu[[1]] = sapply(seq_len(n_reg_trees), function(x) list(rMVN(mu = matrix(0, nrow=p), Q = kappa_reg*diag(p))))
  reg_phi = lapply(seq_len(n_reg_trees), function(x) rep(reg_mu[[1]][[x]], n))

  reg_prior = lapply(seq_len(n_reg_trees), function(x) reg_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  reg_accept = lapply(seq_len(n_reg_trees), function(x) reg_accept[x][[1]] = TRUE)

  #####-------------------- SET INITIAL VALUES FOR PROBIT MODEL --------------------#####
  accepted_class_trees = lapply(seq_len(n_class_trees), function(x) accepted_class_trees[[x]] = df)
  class_change_id = lapply(seq_len(n_class_trees), function(x) class_change_id[[x]] = rep(1, n))

  kappa_class = (4/9)*(n_class_trees)

  class_mu[[1]] = sapply(seq_len(n_class_trees), function(x) list(rMVN(mu = matrix(0, nrow=p), Q = kappa_class*diag(p))))
  class_phi = lapply(seq_len(n_class_trees), function(x) rep(class_mu[[1]][[x]], n))

  class_prior = lapply(seq_len(n_class_trees), function(x) class_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  class_accept = lapply(seq_len(n_class_trees), function(x) class_accept[x][[1]] = TRUE) #do the thing please

  new_R = diag(1, p)

  if(mice_impute){
    imputed = as.matrix(mice::complete(mice::mice(cbind(y, x), print=FALSE)))[,1:p]
    y[missing_index] = imputed[missing_index]
  } else {
    y[missing_index] = 0
  }

  z = matrix(3, ncol = p, nrow = n)
  z[m==0] = -3

  Y = probit_predictors(x, y, include_x = include_x, include_y = include_y)

  nu = hypers$df
  # lambda = (sqrt(2/nu)*qnorm(1-hypers$q) + 1)/sample_t
  qchi = qchisq(1-hypers$q, nu)
  sigest = rep(0, p)
  for(i in 1:p){
    sigest[i] = summary(lm(y[,i]~x))$sigma
  }
  lambda = (sigest^2)*qchi/nu
  # curve(dgamma(x, shape = nu/2, rate = (nu*lambda/2), log=FALSE), from = 0, to = 4)
  # new_omega = sim_omega(y = y, y_hat = Reduce("+", reg_phi), alpha = alpha, Vinv = Vinv)
  new_omega = diag(rgamma(p, shape = nu/2, rate = nu*lambda/2), p) #sim_omega(y = y, y_hat = Reduce("+", reg_phi), nu = nu, lambda = lambda)
  curve(dgamma(x, shape = nu/2, rate = (nu*lambda/2), log=FALSE), from = 0, to = 2)


  #####----- OUT-OF-SAMPLE PREDICTIONS -----#####
  if(predict){
    predict_change_id = vector(mode = "list", length = n_reg_trees)
    reg_phi_pred = vector(mode = "list", length = n_reg_trees)
    n_new = nrow(x_predict)
    new_y_post = vector(mode = "list", length = 0)
  }

  if(p==1 & make_pdp){
    #####----- PDP PLOT -----#####
    if(include_x){
      pdp_list = pdp_param_mat_list(x = x, intercept = FALSE, y_range = pdp_range, include_x = include_x, n = n)
    } else {
      pdp_list = list(seq(pdp_range[1], pdp_range[2], length = 50))
    }
    pdp_phi = matrix(0, nrow = 1, ncol = 50)
    pdp_out = vector(mode = "list", length = total_iters)
    pdp_m = matrix(rep(m, 50), ncol=50)
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
      if(i==1) reg_likely[[j]] = log_marginal_likelihood(node_partial_res = partial_res_y, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      df = accepted_reg_trees[[j]]
      new_tree = propose_tree(df, x, min_node, max_attempt, i) # Propose new tree for tree j
      new_df = new_tree$new_df
      change_points = new_tree$change_points # Get change points for new tree
      decent_tree = new_tree$decent_tree
      # new_df = true_trees_data[[j]]
      # change_points = true_change_points[,j] #get_change_points(new_df, x)
      # decent_tree = TRUE
      # accept = TRUE

      accept = FALSE
      if(decent_tree){
        p2 = tree_priors(nodes = row.names(new_df), parents = unique(new_df$parent), depth = new_df$depth, prior_alpha, prior_beta)
        l2 = sum(sapply(split.data.frame(partial_res_y, change_points), log_marginal_likelihood, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        p1 = reg_prior[[j]]
        l1 = sum(sapply(split.data.frame(partial_res_y, reg_change_id[[j]]), log_marginal_likelihood, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)) #tree_likely[[j]]
        ratio = (p2 + l2) - (p1 + l1)
        accept = ratio >= 0 || - stats::rexp(1L) <= ratio
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
    # new_omega = sim_omega(y, y_hat, alpha = alpha, Vinv = Vinv)
    new_omega = sim_omega(y = y, y_hat = y_hat, nu = nu, lambda = lambda)
    # kappa_reg = sim_kappa(mu = reg_mu[[i]], a = 16, b = 1/n_reg_trees)

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
      if(i==1) class_likely[[k]] = log_marginal_likelihood(node_partial_res = partial_res_z, kappa = kappa_class, omega = new_R, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      df = accepted_class_trees[[k]]
      new_tree = propose_tree(df, Y, min_node, max_attempt, i, probit = TRUE, miss_row = miss_row) # Propose new tree for tree k
      new_df = new_tree$new_df
      change_points = new_tree$change_points # Get change points for new tree
      decent_tree = new_tree$decent_tree
      # new_df = true_trees_missing[[k]]
      # change_points = true_change_points_miss #get_change_points(new_df, Y)
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
        class_change_id[[k]] = change_points
        class_prior[[k]] = p2
        class_likely[[k]] = l2
        # print(new_df)
      }

      ###----- Tree node updates -----###
      ##--Sample mu and compute phi, then store
      L_mu = sim_mu(change_points = class_change_id[[k]], partial_res = partial_res_z, kappa = kappa_class, omega = new_R) #matrix(c(0.7927806, -1.5046567, 1.6511764), ncol=1)
      L_phi = L_mu[class_change_id[[k]],, drop=FALSE]
      class_mu[[i]][[k]] = L_mu
      class_phi[[k]] = L_phi # Used to calculate partial res.

      if(p==1 & make_pdp){
        ###----- PDP -----###
        if(!include_x){
          pdp_change_list = get_change_points(accepted_class_trees[[k]], rbind(Y, matrix(unlist(pdp_list), ncol=1)))[-c(1:n)]
          pdp_phi[1,] = pdp_phi[1,] + L_mu[pdp_change_list,,drop=FALSE]
          # predict_change_id[[j]] = get_change_points(accepted_reg_trees[[j]], rbind(x, x_predict))[-c(1:n)] # assigns each x to a terminal node
          # reg_phi_pred[[j]] = L_mu[predict_change_id[[j]],, drop=FALSE] # mu_j's for new x
        } else {
          pdp_change_list = lapply(pdp_list, function(pdp) get_change_points(accepted_class_trees[[k]], rbind(Y, pdp))[-c(1:n)])
          pdp_phi = pdp_phi + Reduce(rbind, lapply(pdp_change_list, function(x) L_mu[x]))
        }
      }
    }

    #--Get classBART predictions
    z_hat = z = Reduce("+", class_phi)
    if(p==1){
      if(make_pdp){
        pdp_z = pdp_phi #matrix(stats::rnorm(prod(dim(pdp_phi)), mean=pdp_phi, sd=1), ncol=ncol(pdp_phi))
        pdp_out[[i]] = pdp_z
      }
    } else {
      # z = multi_rMVN(z_hat, diag(1,p))
    }
    z[missing_index] = extraDistr::rtnorm(length(missing_index), mean = z_hat[missing_index], sd = 1, b = 0)
    z[obs_index] = extraDistr::rtnorm(length(obs_index), mean = z_hat[obs_index], sd = 1, a = 0)
    # kappa_class = 1 #sim_kappa(mu = class_mu[[i]], a = 16, b = 1/n_class_trees)

    if(include_y){ # If include_y==FALSE, then we are assuming MAR. Missing y's can be "imputed" from the model.
      #--Metropolis Hastings step for y_miss--#
      y_miss = update_y_miss_BART(x = x, y = y, z = z, z_hat = z_hat, y_hat = y_hat, n_trees = n_class_trees,
                                  R = new_R, Omega = new_omega, missing_index = missing_index,
                                  accepted_class_trees = accepted_class_trees, class_mu_i = class_mu[[i]],
                                  include_x = include_x, include_y = include_y, MH_sd = MH_sd, true_change_points = true_change_points_miss)
      y_miss_accept[i,] = y_miss$accept[missing_index]
      y[missing_index] = y_miss$y[missing_index]
    } else {
      y[missing_index] = multi_rMVN(y_hat, new_omega)[missing_index]
    }

    Y = probit_predictors(x, y, include_x = include_x, include_y = include_y, intercept = FALSE)

    for(k in 1:n_class_trees){
      class_change_id[[k]] = get_change_points(accepted_class_trees[[k]], Y)
      L_mu = class_mu[[i]][[k]]
      class_phi[[k]] = L_mu[class_change_id[[k]],, drop=FALSE]
    }

    # z_hat = z = Reduce("+", class_phi)
    # if(p==1){
    #   z = matrix(stats::rnorm(n, mean=z_hat, sd=1), ncol=p, byrow=TRUE)
    # } else {
    #   z = multi_rMVN(z_hat, diag(1,p))
    # }
    # z[intersect(which(z<0), which(m==1))] = 0
    # z[intersect(which(z>=0), which(m==0))] = 0

    ###----- Store posterior samples after burn-in, accounting for thinning -----###
    if(i > burn && i%%thin == 0){
      # y_post = append(y_post, list(unscale(y_hat, min_y, max_y)))
      # y_impute = append(y_impute, list(unscale(y[missing_index], min_y, max_y)))
      # omega_post = append(omega_post, list(1/sqrt(new_omega/(max_y - min_y)^2))) #list(1/sqrt(new_omega/(max_y - min_y)^2))
      if(scale){
        y_post = append(y_post, list(unscale(y_hat, min_y, max_y)))
        y_impute = append(y_impute, list(unscale(y[missing_index], min_y, max_y)))
        omega_post = append(omega_post, list(1/sqrt(new_omega/(max_y - min_y)^2)))
      } else {
        y_post = append(y_post, list((y_hat)))
        y_impute = append(y_impute, list(y[missing_index]))
        omega_post = append(omega_post, list(1/sqrt(new_omega)))
      }
      z_post = append(z_post, list(z))
      kappa_reg_list = c(kappa_reg_list, kappa_reg)

      if(predict){
        new_pred_mean = Reduce("+", reg_phi_pred)
        new_predictions = multi_rMVN(new_pred_mean, new_omega)
        new_y_post = append(new_y_post, list(unscale(new_predictions, min_y, max_y)))
      }
    }

    if(p>1 || !make_pdp) pdp_out = NA

  } # End of i iterations

  if(!predict) new_y_post = NA

  return(structure(list(y_post = y_post, omega_post = omega_post,
                        x = x, y_impute = y_impute, new_y_post = new_y_post,
                        accepted_reg_trees = accepted_reg_trees, accepted_class_trees = accepted_class_trees,
                        burn = burn, iters = iters, thin = thin,
                        max_y = max_y, min_y = min_y,
                        z_post = z_post,
                        y_pred = y_pred, pdp_out = pdp_out, kappa_reg = kappa_reg_list,
                        y_miss_accept = y_miss_accept, reg_mu = reg_mu, class_mu = class_mu), class = "bart"))
}

# print.bart <- function(bart_out, ...) {
#   if(!inherits(bart_out, "bart")) stop("x must be of class 'bart'")
#   y_post = bart_out$y_post
#   y = bart_out$imputed_y
#   x = bart_out$x
#   min_y = bart_out$min_y
#   min_x = bart_out$min_x
#
#   mean_y_post = Reduce("+", y_post)/length(y_post)
#   mean_y_post = unscale(mean_y_post, min = min_y, max = max_y, std = FALSE)
#   y = unscale(y, min = min_y, max = max_y, std = FALSE)
#
#   # for(i in 1:p){
#   #   for(j in c(1,6,7,16,22,25)){
#   #     print(ggplot(data = data.frame(x=x[,j], y=y[,i], m=m[,i]), aes(x, y, color=as.factor(m))) + geom_point(size = 0.9) + ylab(colnames(y)[i]) + xlab(colnames(x)[j]))
#   #   }
#   # }
#
#   bart_plot = list()
#   for(i in 1:p){
#     plot_data = data.frame(true = y[which(m[,i]==1),i], pred = mean_y_post[which(m[,i]==1),i])
#     min = min(plot_data)
#     max = max(plot_data)
#     cheat = data.frame(y_seq = seq(min, max, length=nrow(plot_data)))
#     bart_plot = ggplot(plot_data, aes(true, pred)) + geom_point() + geom_line(data=cheat, aes(y_seq, y_seq), colour="black", size=0.1)
#     print(bart_plot)
#   }
#
#   for(i in 1:p){
#     # for(j in c(1,6,7,16,22,25)){
#       # print(ggplot(data.frame(x=x[,j], y=mean_y_post[,i], m=m[,i]), aes(x, y, colour=factor(m))) + geom_point())
#       print(ggplot(data.frame(y = mean_y_post[,i], m = factor(m[,i])), aes(x = y, colour=m)) + geom_histogram(fill="white") + facet_grid(m ~ .))
#     # }
#   }
# }
