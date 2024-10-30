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
#' @param MH_sd standard deviation for MH proposal for missing Y
#' @param pdp_range range for partial dependence plots
#' @param make_pdp logical indicating whether to produce a partial dependence plot
#' @param mice_impute logical indicating whether to impute missing values via mice prior to prior calibration
#' @param true_trees_data true trees for BART component
#' @param true_trees_missing true trees for probit BART component
#' @param true_change_points true change points for BART trees
#' @param true_change_points_miss true change points for probit BART trees
#' @param ... Catches unused arguments
#'
#' @importFrom extraDistr "rtnorm"
#' @importFrom mice "complete" "mice"
#' @return a list containing BART predictions and imputed values
#' @export
#'
#' @examples
#' # x <- matrix(runif(6), ncol = 2)
#' # y <- matrix(runif(6), ncol = 2) %*% matrix(rnorm(4), ncol=2)
#' # bart_out <- missBART2(x, y, n_trees = 2, burn = 2,
#' #                       iters = 2, thin = 1, scale = FALSE)
missBART2 <- function(x, y, x_predict = c(), n_reg_trees = 100, n_class_trees = 100, burn = 1000, iters = 1000, thin = 2,
                      predict = TRUE, MH_sd = NULL, tree_prior_params = tree_list(...), hypers = hypers_list(...),
                      scale = TRUE, include_x = TRUE, include_y = TRUE, show_progress = TRUE, progress_every = 10,
                      pdp_range = c(-0.5, 0.5), make_pdp = FALSE, mice_impute = TRUE, true_trees_data = NA,
                      true_trees_missing = NA, true_change_points = NA, true_change_points_miss = NA, ...) {

  bart_img2()

  if(is.null(x_predict)) predict <- FALSE
  y <- as.matrix(y)
  x <- as.matrix(x)
  x_predict <- as.matrix(x_predict)

  missing_index <- which(is.na(y))
  obs_index <- which(!is.na(y))
  miss_row <- apply(is.na(y), 1, any)

  min_y <- apply(y, 2, min, na.rm = TRUE) #min(y, na.rm = TRUE)
  max_y <- apply(y, 2, max, na.rm = TRUE) #max(y, na.rm = TRUE)
  if(scale) {
    y <- t(apply(sweep(y, 2, min_y), 1, function(x) x/(max_y - min_y))) - 0.5
    if(nrow(y) == 1) y <- t(y)
  }

  #####-------------------- GET PARAMETERS --------------------#####
  p <- ncol(y) # No. of y variables
  q <- ncol(x) # No. of x variables
  r <- p * include_y + q * include_x
  Y_vars <- seq_len(r)
  n <- nrow(y)
  thinned <- iters
  total_iters <- burn + thin * iters

  m <- matrix(1, nrow=n, ncol=p)
  m[is.na(y)] <- 0

  #####-------------------- GET BART PRIOR PARAMETERS --------------------#####
  mu0 <- rep(hypers$mu0, p)
  kappa_reg <- ifelse(is.null(hypers$kappa), 4 * (stats::qnorm(0.975))^2 * n_reg_trees, hypers$kappa)
  nu <- hypers$df
  qchi <- stats::qchisq(1 - hypers$q, nu)
  sigest <- rep(0, p)
  for(i in seq_len(p)) {
    sigest[i] <- summary(stats::lm(y[,i]~x))$sigma
  }
  lambda <- (sigest^2) * qchi/nu

  if(p>1) {
    alpha <- ifelse(is.null(hypers$alpha), nu, hypers$alpha)
    if(is.null(hypers$V)) {
      V <- diag(1/(lambda * alpha), p)
      Vinv <- diag(lambda * alpha,  p)
    } else {
      V <- hypers$V
      Vinv <- solve(V)
    }
  }
  if(is.null(MH_sd)) MH_sd <- 0.5/p

  #####-------------------- FIRST IMPUTATION OF DATA --------------------#####
  if(mice_impute) {
    imputed <- as.matrix(mice::complete(mice::mice(cbind(y, x), print = FALSE))[,seq_len(p)])
    y[missing_index] <- imputed[missing_index]
  } else {
    y[missing_index] <- 0
  }

  #####-------------------- MISSING X --------------------#####
  if(any(is.na(x))){
    m_x <- data.frame(mx=matrix(1, nrow=nrow(x), ncol=ncol(x)))
    m_x[is.na(x)] <- 0
    x <- as.matrix(cbind(x, m_x))
    if(predict){
      m_x_pred <- data.frame(mx=matrix(1, nrow=nrow(x_predict), ncol=ncol(x_predict)))
      m_x_pred[is.na(x_predict)] <- 0
      x_predict <- cbind(x_predict, m_x_pred)
      x_predict <- as.matrix(x_predict)
    }
    q <- ncol(x)
    x_vars <- seq_len(q)
  }

  #####-------------------- GET TREE PRIOR PARAMETERS --------------------#####
  prior_alpha <- tree_prior_params$prior_alpha
  prior_beta <- tree_prior_params$prior_beta
  min_node <- max(tree_prior_params$min_node, p + 1)
  max_attempt <- tree_prior_params$max_attempt

  #####-------------------- CREATE STORAGE FOR REGRESSION BART --------------------#####
  accepted_reg_trees <- lapply(vector(mode = "list", length = n_reg_trees), as.list)
  reg_trees <- vector(mode = "list", length = 0)
  reg_change_id <- vector(mode = "list", length = n_reg_trees)

  reg_mu <- vector(mode = "list", length = thinned)
  reg_phi <- vector(mode = "list", length = n_reg_trees)
  omega_post <- vector(mode = "list", length = 0)

  reg_prior <- lapply(vector(mode = "list", length = n_reg_trees), as.list) # tree prior for accepted trees
  reg_likely <- lapply(vector(mode = "list", length = n_reg_trees), as.list) # likelihood for accepted trees
  reg_accept <- lapply(vector(mode = "list", length = n_reg_trees), as.list) # accept/reject status for all trees: reg_accept[[j]][[i]]

  #####-------------------- CREATE STORAGE FOR PROBIT BART --------------------#####
  R_post <- vector(mode = "list", length = 0)
  y_post <- vector(mode = "list", length = 0)
  y_pred <- vector(mode = "list", length = 0)
  z_post <- vector(mode = "list", length = 0)
  y_impute <- vector(mode = "list", length = 0)
  var_imp <- vector(mode = "list", length = 0)

  accepted_class_trees <- lapply(vector(mode = "list", length = n_class_trees), as.list)
  class_trees <- vector(mode = "list", length = 0)
  class_change_id <- vector(mode = "list", length = n_class_trees)

  class_mu <- vector(mode = "list", length = thinned)
  class_phi <- vector(mode = "list", length = n_class_trees)
  omega_post <- vector(mode = "list", length = 0)

  class_prior <- lapply(vector(mode = "list", length = n_class_trees), as.list) # tree prior for accepted trees
  class_likely <- lapply(vector(mode = "list", length = n_class_trees), as.list) # likelihood for accepted trees
  class_accept <- lapply(vector(mode = "list", length = n_class_trees), as.list) # accept/reject status for all trees: class_accept[[j]][[i]]
  class_moves <- NULL

  y_miss_accept <- matrix(nrow = total_iters, ncol = length(missing_index)) #vector(mode = "list", length = total_iters) # Metropolis-Hastings acceptance/rejection for missing y's

  partial_res_y <- matrix(nrow = n, ncol = p)
  partial_res_z <- matrix(nrow = n, ncol = p)

  partial_reg <- vector(mode = "list", length = n_reg_trees)

  #####-------------------- SET INITIAL VALUES FOR REGRESSION BART --------------------#####
  df <- data.frame(matrix(ncol = 8, nrow = 1))
  colnames(df) <- c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction", "NA_direction")
  df[1,] <- c(0,0,1,0,1,0,0,NA)

  accepted_reg_trees <- lapply(seq_len(n_reg_trees), function(x) accepted_reg_trees[[x]] = df)
  reg_change_id <- lapply(seq_len(n_reg_trees), function(x) reg_change_id[[x]] = rep(1, n))

  # kappa_reg_list <- c()
  reg_mu[[1]] <- sapply(seq_len(n_reg_trees), function(x) list(rMVN(mu = matrix(0, nrow=p), Q = diag(kappa_reg, p))))
  reg_phi <- lapply(seq_len(n_reg_trees), function(x) rep(reg_mu[[1]][[x]], n))

  reg_prior <- lapply(seq_len(n_reg_trees), function(x) reg_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  reg_accept <- lapply(seq_len(n_reg_trees), function(x) reg_accept[x][[1]] = TRUE)

  #####-------------------- SET INITIAL VALUES FOR PROBIT MODEL --------------------#####
  accepted_class_trees <- lapply(seq_len(n_class_trees), function(x) accepted_class_trees[[x]] = df)
  class_change_id <- lapply(seq_len(n_class_trees), function(x) class_change_id[[x]] = rep(1, n))

  kappa_class <- 4/9*(n_class_trees)
  class_mu[[1]] <- sapply(seq_len(n_class_trees), function(x) list(rMVN(mu = matrix(0, nrow=p), Q = diag(kappa_class, p))))
  class_phi <- lapply(seq_len(n_class_trees), function(x) rep(class_mu[[1]][[x]], n))

  class_prior <- lapply(seq_len(n_class_trees), function(x) class_prior[x][[1]] = log(node_priors(0, prior_alpha, prior_beta)))
  class_accept <- lapply(seq_len(n_class_trees), function(x) class_accept[x][[1]] = TRUE) #do the thing please

  z <- matrix(1, ncol = p, nrow = n)
  z[m == 0] <- -1

  new_R <- diag(p)
  new_omega <- diag(stats::rgamma(p, shape = nu/2, rate = nu * lambda/2), p)
  Y <- probit_predictors(x, y, include_x = include_x, include_y = include_y)

  #####----- OUT-OF-SAMPLE PREDICTIONS -----#####
  if(predict) {
    predict_change_id <- vector(mode = "list", length = n_reg_trees)
    reg_phi_pred <- vector(mode = "list", length = n_reg_trees)
    n_new <- nrow(x_predict)
    new_y_post <- vector(mode = "list", length = 0)
  }

  if(p == 1 &&  make_pdp) {
    #####----- PDP PLOT -----#####
    if(include_x) {
      pdp_list <- pdp_param_mat_list(x = x, intercept = FALSE, y_range = pdp_range, include_x = include_x, n = n)
    } else {
      pdp_list <- list(seq(pdp_range[1], pdp_range[2], length = 50))
    }
    pdp_phi <- matrix(0, nrow = 1, ncol = 50)
    pdp_out <- vector(mode = "list", length = total_iters)
    pdp_m <- matrix(rep(m, 50), ncol=50)
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

    #### -------------------- BART REGRESSION FOR DATA
    reg_mu[[i]] <- vector(mode = "list", n_reg_trees)
    for(j in seq_len(n_reg_trees)) {
      ###----- Compute partial residuals -----###
      if(n_reg_trees == 1) {
        partial_res_y <- y
      } else {
        partial_res_y <- y - Reduce("+", reg_phi[-j])
      }

      ###----- Set likelihood of first stump of each tree -----###
      if(i == 1) reg_likely[[j]] <- log_marginal_likelihood(node_partial_res = partial_res_y, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      df <- accepted_reg_trees[[j]]
      new_tree <- propose_tree(df, x, min_node, max_attempt, i) # Propose new tree for tree j
      new_df <- new_tree$new_df
      change_points <- new_tree$change_points # Get change points for new tree
      decent_tree <- new_tree$decent_tree
      # new_df <- true_trees_data[[j]]
      # change_points <- get_change_points(new_df, x) #true_change_points[,j]
      # decent_tree <- TRUE
      # accept <- TRUE

      accept <- FALSE
      if(decent_tree) {
        p2 <- tree_priors(new_df = new_df, prior_alpha, prior_beta)
        l2 <- sum(sapply(split.data.frame(partial_res_y, change_points), log_marginal_likelihood, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        p1 <- reg_prior[[j]]
        l1 <- sum(sapply(split.data.frame(partial_res_y, reg_change_id[[j]]), log_marginal_likelihood, kappa = kappa_reg, omega = new_omega, mu0 = mu0, Vinv = Vinv, alpha = alpha)) #tree_likely[[j]]
        ratio <- (p2 + l2) - (p1 + l1)
        accept <- ratio >= 0 || - stats::rexp(1L) <= ratio
      }

      if(accept) {
        accepted_reg_trees[[j]] <- new_df
        reg_prior[[j]] <- p2
        reg_likely[[j]] <- l2
        reg_change_id[[j]] <- change_points
        partial_reg[[j]] <- partial_res_y
      }

      ###----- Tree node updates -----###
      ##--Sample mu and compute phi, then store
      L_mu <- sim_mu(reg_change_id[[j]], partial_res_y, kappa_reg, new_omega)
      L_phi <- L_mu[reg_change_id[[j]],, drop=FALSE]
      reg_mu[[i]][[j]] <- L_mu
      reg_phi[[j]] <- L_phi # Used to calculate partial res.

      ###----- Out-of-Sample Predictions -----###
      if(predict) {
        # all_points <- get_change_points(accepted_reg_trees[[j]], rbind(x, x_predict))
        predict_change_id[[j]] <- get_change_points(accepted_reg_trees[[j]], rbind(x, x_predict))[-seq_len(n)] # assigns each x to a terminal node
        reg_phi_pred[[j]] <- L_mu[predict_change_id[[j]],, drop=FALSE] # mu_j's for new x
      }

    } # End of j iterations

    ###----- BART updates -----###
    #--Get BART predictions
    y_hat <- Reduce("+", reg_phi)

    #--Sample data precision
    if(p == 1) {
      new_omega <- sim_omega(y = y, y_hat = y_hat, nu = nu, lambda = lambda)
    } else {
      new_omega <- sim_omega(y = y, y_hat = y_hat, alpha = alpha, Vinv = Vinv)
    }
    # kappa_reg <- sim_kappa(mu = reg_mu[[i]], a = 16, b = 1/n_reg_trees)

    ###----- Probit BART -----###
    class_mu[[i]] <- vector(mode = "list", n_class_trees)

    for(k in seq_len(n_class_trees)) {

      ###----- Compute partial residuals -----###
      if(n_class_trees == 1) {
        partial_res_z <- z
      } else {
        partial_res_z <- z - Reduce("+", class_phi[-k])
      }

      ###----- Set likelihood of first stump of each tree -----###
      if(i == 1) class_likely[[k]] <- log_marginal_likelihood(node_partial_res = partial_res_z, kappa = kappa_class, omega = new_R, mu0 = mu0, Vinv = Vinv, alpha = alpha)

      ###----- Propose new tree -----###
      df <- accepted_class_trees[[k]]
      new_tree <- propose_tree(df, Y, min_node, max_attempt, i, probit = TRUE, miss_row = miss_row) # Propose new tree for tree k
      new_df <- new_tree$new_df
      change_points <- new_tree$change_points # Get change points for new tree
      decent_tree <- new_tree$decent_tree
      # new_df <- true_trees_missing[[k]]
      # change_points <- get_change_points(new_df, Y) #true_change_points_miss
      # decent_tree <- TRUE
      # accept <- TRUE

      ###----- Metropolis-Hastings for accepting/rejecting proposed tree -----###
      accept <- FALSE
      if(decent_tree) {
        p2 <- tree_priors(new_df = new_df, prior_alpha, prior_beta)
        l2 <- sum(sapply(split.data.frame(partial_res_z, change_points), log_marginal_likelihood, kappa = kappa_class, omega = new_R, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        l1 <- sum(sapply(split.data.frame(partial_res_z, class_change_id[[k]]), log_marginal_likelihood, kappa = kappa_class, omega = new_R, mu0 = mu0, Vinv = Vinv, alpha = alpha))
        accept <- MH(class_prior[[k]], l1, p2, l2)
      }

      if(accept) {
        accepted_class_trees[[k]] <- new_df
        class_change_id[[k]] <- change_points
        class_prior[[k]] <- p2
        class_likely[[k]] <- l2
      }

      ###----- Tree node updates -----###
      ##--Sample mu and compute phi, then store
      L_mu <- sim_mu(change_points = class_change_id[[k]], partial_res = partial_res_z, kappa = kappa_class, omega = new_R) #matrix(c(0.7927806, -1.5046567, 1.6511764), ncol=1)
      L_phi <- L_mu[class_change_id[[k]],, drop=FALSE]
      class_mu[[i]][[k]] <- L_mu
      class_phi[[k]] <- L_phi # Used to calculate partial res.

      if(p == 1 && make_pdp) {
        ###----- PDP -----###
        if(!include_x) {
          pdp_change_list <- get_change_points(accepted_class_trees[[k]], rbind(Y, matrix(unlist(pdp_list), ncol=1)))[-seq_len(n)]
          pdp_phi[1,] <- pdp_phi[1,] + L_mu[pdp_change_list,, drop=FALSE]
          # predict_change_id[[j]] <- get_change_points(accepted_reg_trees[[j]], rbind(x, x_predict))[-seq_len(n)] # assigns each x to a terminal node
          # reg_phi_pred[[j]] <- L_mu[predict_change_id[[j]],, drop=FALSE] # mu_j's for new x
        } else {
          pdp_change_list <- lapply(pdp_list, function(pdp) get_change_points(accepted_class_trees[[k]], rbind(Y, pdp))[-seq_len(n)])
          pdp_phi <- pdp_phi + Reduce(rbind, lapply(pdp_change_list, function(x) L_mu[x]))
        }
      }
    }

    #--Get classBART predictions
    z_hat <- z <- Reduce("+", class_phi)
    if(p == 1) {
      if(make_pdp) {
        pdp_z <- pdp_phi #matrix(stats::rnorm(prod(dim(pdp_phi)), mean=pdp_phi, sd=1), ncol=ncol(pdp_phi))
        pdp_out[[i]] <- pdp_z
      }
    } else {
      # z <- multi_rMVN(z_hat, diag(p))
    }
    z[missing_index] <- extraDistr::rtnorm(length(missing_index), mean = z_hat[missing_index], sd = 1, b = 0)
    z[obs_index] <- extraDistr::rtnorm(length(obs_index), mean = z_hat[obs_index], sd = 1, a = 0)
    # kappa_class <- 1 #sim_kappa(mu = class_mu[[i]], a = 16, b = 1/n_class_trees)

    if(include_y) { # If include_y==FALSE, then we are assuming MAR. Missing y's can be "imputed" from the model.
      #--Metropolis Hastings step for y_miss--#
      y_miss <- update_y_miss_BART(x = x, y = y, z = z, z_hat = z_hat, y_hat = y_hat, n_trees = n_class_trees,
                                   R = new_R, Omega = new_omega, missing_index = missing_index,
                                   accepted_class_trees = accepted_class_trees, class_mu_i = class_mu[[i]],
                                   true_change_points = true_change_points_miss, include_x = include_x, include_y = include_y, MH_sd = MH_sd)
      y_miss_accept[i,] <- matrix(rep(y_miss$accept, p), ncol = p)[missing_index]
      y[missing_index] <- y_miss$y[missing_index]
    } else {
      y[missing_index] <- multi_rMVN(y_hat, new_omega)[missing_index]
    }

    Y <- probit_predictors(x, y, include_x = include_x, include_y = include_y, intercept = FALSE)

    for(k in seq_len(n_class_trees)) {
      class_change_id[[k]] <- get_change_points(accepted_class_trees[[k]], Y)
      L_mu <- class_mu[[i]][[k]]
      class_phi[[k]] <- L_mu[class_change_id[[k]],, drop=FALSE]
    }

    ###----- Store posterior samples after burn-in, accounting for thinning -----###
    if(i > burn && i %% thin == 0) {
      if(scale) {
        y_post <- append(y_post, list(unscale(y_hat, min_y, max_y)))
        if(p == 1) {
          y_impute <- append(y_impute, list(unscale(y[missing_index], min_y, max_y)))
          omega_post <- append(omega_post, list(1/(new_omega/(max_y - min_y)^2))) # returns the residual variance on the original scale
        } else {
          y_impute <- append(y_impute, list(unscale(y, min_y, max_y)[missing_index]))
          # omega_post <- append(omega_post, list(chol2inv(chol(new_omega))))
          omega_post <- append(omega_post, list(diag(chol2inv(chol(new_omega))) * (max_y - min_y)^2)) # returns the residual covariance matrix on the original scale
        }
      } else {
        y_post <- append(y_post, list((y_hat)))
        y_impute <- append(y_impute, list(y[missing_index]))
        if(p == 1) {
          omega_post <- append(omega_post, list(1/(new_omega)))
        } else {
          omega_post <- append(omega_post, list(chol2inv(chol(new_omega))))
        }
      }
      z_post <- append(z_post, list(z))
      # kappa_reg_list <- c(kappa_reg_list, kappa_reg)

      var_imp <- append(var_imp, list(table(Reduce(append,lapply(accepted_class_trees, function(X) X[-1,4])))/2))

      reg_trees <- append(reg_trees, list(accepted_reg_trees))
      class_trees <- append(class_trees, list(accepted_class_trees))

      if(predict) {
        new_pred_mean <- Reduce("+", reg_phi_pred)
        new_predictions <- multi_rMVN(new_pred_mean, new_omega)
        new_y_post <- append(new_y_post, list(unscale(new_predictions, min_y, max_y)))
      }
    }

    if(p > 1 || !make_pdp) pdp_out <- NA

  } # End of i iterations

  if(!predict) new_y_post <- NA

  return(structure(list(y_post = y_post, omega_post = omega_post,
                        x = x, y_impute = y_impute, new_y_post = new_y_post,
                        reg_trees = reg_trees, class_trees = class_trees,
                        burn = burn, iters = iters, thin = thin,
                        max_y = max_y, min_y = min_y,
                        z_post = z_post,
                        y_pred = y_pred, pdp_out = pdp_out,
                        y_miss_accept = y_miss_accept,
                        reg_mu = reg_mu, class_mu = class_mu, MH_sd = MH_sd,
                        var_imp = var_imp), class = "bart"))
}
