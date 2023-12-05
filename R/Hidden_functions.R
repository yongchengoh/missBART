rMVNmu0 <- function(Q) {
  p  <- NCOL(Q)
  z  <- stats::rnorm(p)
  if(p == 1) z/sqrt(Q) else backsolve(PD_chol(Q), z)
}

PD_chol <- function(x, ...) tryCatch(chol(x, ...), error=function(e) {
  chol(Matrix::nearPD(x, base.matrix=TRUE)$mat, ...)
})

MH <- function(p1, l1, p2, l2, q1 = 0, q2 = 0){
  ratio = (p2 + l2 + q2) - (p1 + l1 + q1)
  accept = ratio >= 0 || - stats::rexp(1L) <= ratio
  return(accept)
}

sim_mu_node = function(node_partial_res, kappa, omega){
  n = nrow(node_partial_res)
  p = ncol(node_partial_res)
  Q = n*omega + diag(kappa, p)
  b = omega %*% colSums(node_partial_res)
  return(rMVNc(b, Q))
}

sim_mu <- function(change_points, partial_res, kappa, omega) {
  return(matrix((sapply(split.data.frame(partial_res, change_points), sim_mu_node, kappa = kappa, omega = omega)), nrow = length(unique(change_points)), byrow = TRUE))
}

# sim_omega = function(y, y_hat, alpha, Vinv){
#   n = nrow(y)
#   df = alpha + n
#   scale = crossprod(y-y_hat) + Vinv
#   scale = chol2inv(PD_chol(scale))
#   if(ncol(y)==1){
#     a = df/2
#     b = 0.5/scale
#     omega = stats::rgamma(1, a, b)
#   } else {
#     omega = stats::rWishart(1, df, scale)[,,1]
#   }
#   return(omega)
# }
sim_omega = function(y, y_hat, nu = NULL, lambda = NULL, alpha = NULL, Vinv = NULL){
  n = nrow(y)
  p = ncol(y)
  if(ncol(y)==1){
    if(is.null(nu) || is.null(lambda)) stop("Must specify values for nu and lambda")
    shape = (nu + n)/2
    rate = (colSums((y - y_hat)^2) + nu*lambda)/2
    omega = diag(stats::rgamma(p, shape, rate), p)
  } else {
  if(is.null(alpha) || is.null(Vinv)) stop("Must specify values for alpha and Vinv")
  df = alpha + n
  scale = chol2inv(chol(crossprod(y-y_hat) + Vinv))
  omega = stats::rWishart(1, df, scale)[,,1]
  }
  return(omega)
}

sim_kappa = function(mu, a, b){ #tree_mu[[i]]
  mu_mat = Reduce(rbind, mu)
  p = ncol(mu_mat)
  n_mu = nrow(mu_mat)
  shape = n_mu*p/2
  rate = b + sum(apply(mu_mat, 1, crossprod))/2
  return(stats::rgamma(1, shape, rate))
}

###--------------------------------- PROBIT UPDATES ---------------------------------###
probit_predictors = function(x, y, include_x, include_y, intercept = FALSE){
  n = nrow(x)
  if(intercept) Y = matrix(1,nrow=n, ncol=1) else Y = c()
  if(include_x & !include_y){
    Y = cbind(Y, x)
  } else if (!include_x & include_y){
    Y = cbind(Y, y)
  } else if (include_x & include_y){
    Y = cbind(Y, x, y)
  }
  return(Y)
}

#' update latent variable z
#'
#' @param Y
#' @param m
#' @param B
#' @param R
#'
#' @return
#' @export
#' @importFrom tmvtnorm "rtmvnorm"
#'
#' @examples
update_z = function(Y, m, B, R){
  n = nrow(Y)
  p = ncol(m)
  mu = Y %*% B
  # U = diag(1, n)
  z = matrix(nrow = n, ncol = p)
  if(p==1){
    z[m==0] = extraDistr::rtnorm(length(which(m==0)), mean = mu[m==0], sd = 1, a = -Inf, b = 0)
    z[m==1] = extraDistr::rtnorm(length(which(m==1)), mean = mu[m==1], sd = 1, a = 0, b = Inf)
  } else {
    m_lower = matrix(0, nrow = n, ncol = p)
    m_lower[m==0] = -Inf
    m_upper = matrix(0, nrow = n, ncol = p)
    m_upper[m==1] = Inf
    # for(i in 1:p){
    #   z[,i] = extraDistr::rtnorm(n, mean = mu[,i], sd = 1, a = m_lower[,i], b = m_upper[,i])
    # }
    for(i in 1:n){
      # z[i,] = tmvtnorm::rtmvnorm(n = 1, mean = mu[i,], sigma = R, lower = m_lower[i,], upper = m_upper[i,], algorithm = "gibbs")
      z[i,] = TruncatedNormal::rtmvnorm(n = 1, mu = mu[i,], sigma = R, lb = m_lower[i,], ub = m_upper[i,])
    }
  }

  return(z)
}

# MCMC sample of B, the matrix of probit parameters
update_B = function(y, z, Y, tau_b){
  r = ncol(Y)
  p = ncol(y)
  V = diag(p)
  U = chol2inv(chol(diag(tau_b, r) + t(Y)%*%Y))
  M = t(t(z) %*% Y %*% U)
  return(matrnorm(M, U, V))
  # b = colSums(sweep(Y, 1, z, "*"))
  # y_sum = matrix(0, ncol = r, nrow = r)
  # for(i in 1:nrow(y)){
  #   y_sum = y_sum + Y[i,] %*% t(Y[i,])
  # }
  # Q = diag(tau_b, r) + y_sum
  # mean = b %*% solve(Q)
  # return(matrix(rMVN(mu = mean, Q = Q), ncol = 1))
  # return(matrix(rMVNc(b = b, Q = Q), ncol = 1))
}

# MCMC sample of W, the latent variable in probit regression introduced by \cite{Talhouk, A., Doucet, A., & Murphy, K. (2012). Efficient Bayesian inference for multivariate probit models with sparse inverse correlation matrices. Journal of Computational and Graphical Statistics, 21(3), 739-757.}
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

# MCMC sample of \eqn{\Sigma} and \eqn{\gamma}, the latent variables in probit regression introduced by \cite{Talhouk, A., Doucet, A., & Murphy, K. (2012). Efficient Bayesian inference for multivariate probit models with sparse inverse correlation matrices. Journal of Computational and Graphical Statistics, 21(3), 739-757.}
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

# MCMC sample of R (correlation matrix) and B (probit parameters)
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

# MCMC sample of \eqn{\Psi} from \cite{Talhouk, A., Doucet, A., & Murphy, K. (2012). Efficient Bayesian inference for multivariate probit models with sparse inverse correlation matrices. Journal of Computational and Graphical Statistics, 21(3), 739-757.}
#' @importFrom CholWishart "rInvWishart"
update_Psi = function(nu, S, B, R){
  r = nrow(B)
  df = nu + r
  Sigma = S + B %*% chol2inv(PD_chol(R)) %*% t(B)
  return(rInvWishart(1, df, Sigma)[,,1])
}

# MCMC sample of missing data
# update_y_miss_reg <- function(x, y_hat, m, y, z, B, R, omega, include_x = TRUE, include_y = TRUE) {
#   M = apply(m, 2, as.logical)
#   n <- nrow(y_hat)
#   p <- ncol(y_hat)
#   q <- ncol(x)
#   if(include_x && !include_y) {
#     A   <- B[2:(q+1),, drop=FALSE]
#     B1  <- matrix(0, nrow=p, ncol=p)
#   } else if(!include_x && include_y) {
#     A   <- matrix(0, nrow=q, ncol=p)
#     B1  <- B[-1,, drop=FALSE]
#   } else if(include_x  && include_y) {
#     A   <- B[2:(q+1),, drop=FALSE]
#     B1  <- B[(q+2):(p+q+1),, drop=FALSE]
#   } else {
#     A   <- matrix(0, nrow=q, ncol=p)
#     B1  <- matrix(0, nrow=p, ncol=p)
#   }
#   b1    <- t(B[1,,drop=FALSE])
#   prec  <- omega + if(p == 1) tcrossprod(B1)/R else B1 %*% solve(R, t(B1))
#   Q     <- if(p == 1) sqrt(prec) else PD_chol(prec)
#   Oyhat <- tcrossprod(omega, y_hat)
#   zXA   <- if(include_x) z - x %*% A else z
#   BRinv <- B1 %*% chol2inv(PD_chol(R))
#   b     <- Oyhat + tcrossprod(BRinv, sweep(zXA, 2, b1, FUN="-", check.margin=FALSE))
#   Y     <- rMVNc(b, Q, is_chol=TRUE)
#   Y     <- if(p == 1) matrix(Y, ncol=p) else t(Y)
#   return(Y)
# }

update_y_miss_reg <- function(x, y_hat, m, z, B, R, omega, include_x = TRUE) {
  missing_index = which(m==0)
  # n_miss = length(missing_index)
  p = ncol(y_hat)
  # q = ncol(x)
  By = B[c((nrow(B)-p+1):nrow(B)),,drop=FALSE]
  A = B[-c((nrow(B)-p+1):nrow(B)),,drop=FALSE]
  X = if(include_x) cbind(rep(1,nrow(x)), x) else matrix(1, nrow = nrow(x), ncol = 1)
  R_inv = chol2inv(PD_chol(R))
  b = tcrossprod(omega, y_hat) + crossprod(t(By), tcrossprod(R_inv, (z - X %*% A))) #omega %*% t(y_hat) + By %*% solve(R) %*% t((z - X %*% A))
  Q = omega + By %*% R_inv %*% t(By) #crossprod(t(By), crossprod(t(R_inv), By)) #By %*% R_inv %*% By
  # print(Q)
  mu = t(b) %*% chol2inv(PD_chol(Q)) #if(p==1) t(b) %*% chol2inv(PD_chol(Q)) else (chol2inv(PD_chol(Q)) %*% t(b))
  # print((mu))
  y_miss = multi_rMVN(mean_mat = mu, precision = Q)[missing_index]

  # y_miss = matrix(ncol = p, nrow = n_miss)
  # for(i in 1:n_miss){
  #   miss_i = missing_index[i]
  #   X1 = if(include_x) matrix(c(1, x[miss_i,]), ncol=1) else matrix(1)
  #   b = y_hat[miss_i,] %*% omega + By %*% chol2inv(PD_chol(R)) %*% (z[miss_i,] - t(A) %*% X1)
  #   Q = omega + By %*% chol2inv(PD_chol(R)) %*% By
  #   y_miss[i,] = rMVNc(b = b, Q = Q)
  # }
  return(y_miss)
}

# update_y_miss_reg <- function(x, y_hat, m, y, z, B, R, omega, include_x = TRUE, include_y = TRUE) {
#   prec = omega + (B[2,])^2
#   mu = (y_hat%*%omega + (z - B[1,])%*%B[2,])/as.vector(prec) #((y_hat - B[1,]*B[2,])*as.vector(omega) + z*B[2,])/as.vector(prec)
#   Y = matrix(rnorm(nrow(mu), mu, 1/sqrt(prec)), ncol=1)
#   return(Y)
# }

# Get vector for sorting observations into terminal nodes
get_change_points <- function(df, x) {
  n = nrow(x)
  if(nrow(df) == 1) {
    change_points <- rep(1, n)
  } else {
    change_points <- c()
    terminal_nodes <- as.numeric(setdiff(row.names(df), df$parent))
    terminal_obs <- list()
    for(j in seq_along(terminal_nodes)) {
      y_index <- list()
      parent <- terminal_nodes[j]
      count <- 1
      while(parent != 1) {
        #Add in if/else for type of variable (numeric, integer, categorical)
        #Add in split if missing (not if m==0/1)
        type = class(x[,df$split_variable[parent]])
        # switch(EXPR = type,
        #        "numeric" = {
        #          y_index[[count]] = if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] <= df$split_value[parent]) else which(x[,df$split_variable[parent]] > df$split_value[parent])
        #          if(any(is.na(x[,df$split_variable[parent]]))){
        #            if(df$NA_direction[parent] == 1) y_index[[count]] = sort(c(y_index[[count]], which(is.na(x[,df$split_variable[parent]]))))
        #          }
        #        }, "factor" = {
        #          y_index[[count]] = if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] == df$split_value[parent]) else which(x[,df$split_variable[parent]] != df$split_value[parent])
        #        }) #End switch
        y_index[[count]] = if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] <= df$split_value[parent]) else which(x[,df$split_variable[parent]] > df$split_value[parent])
        if(any(is.na(x[,df$split_variable[parent]]))){
          if(df$NA_direction[parent] == 1) y_index[[count]] = sort(c(y_index[[count]], which(is.na(x[,df$split_variable[parent]]))))
        }
        parent <- df$parent[parent]
        count <- count + 1
      }
      terminal_obs[[j]] = Reduce(intersect, y_index)
      if(length(terminal_obs[[j]])==0) terminal_obs[[j]] = NA
    }
    terminal_obs = terminal_obs[!is.na(terminal_obs)]
    change_points = vector(length=n)

    for(k in 1:length(terminal_obs)){
      change_points[terminal_obs[[k]]] = k
    }
    # max_node = max(unique(change_points))
    # change_points[which(change_points==0)] = sample(c(max_node, max_node - 1), 1)
  }
  return(change_points)
}

update_y_miss_BART = function(x, y, z, z_hat, y_hat, n_trees, R, Omega, missing_index, accepted_class_trees, class_mu_i, include_x, include_y, MH_sd = 0.2, true_change_points){
  n = nrow(y)
  p = ncol(y)
  q = ncol(x)
  R_inv = if(p==1) 1 else chol2inv(PD_chol(R))
  n_miss = length(missing_index)

  p1 = -0.5*apply(y - y_hat, 1, function(x) crossprod(x, (Omega %*% x)))
  l1 = -0.5*apply(z - z_hat, 1, function(x) crossprod(x, (R_inv %*% x)))

  proposed_y = y
  proposed_y[missing_index] = stats::rnorm(n_miss, mean = y[missing_index], sd = MH_sd)
  Y = probit_predictors(x = x, y = proposed_y, include_x = include_x, include_y = include_y)
  z_hat = matrix(0, nrow=n, ncol=p)
  for(k in 1:n_trees){
    MH_change_id = get_change_points(accepted_class_trees[[k]], rbind(probit_predictors(x = x, y = y, include_x = include_x, include_y = include_y), Y))[-c(1:n)]
    # MH_change_id = true_change_points
    # MH_change_id = get_change_points(accepted_class_trees[[k]], Y)
    if(length(unique(MH_change_id)) != nrow(class_mu_i[[k]])) {
      print(MH_change_id)
      print(class_mu_i[[k]])
    }
    z_hat = z_hat + class_mu_i[[k]][MH_change_id,, drop=FALSE]
  }
  if(p==1){
    # z_hat = t(z_hat)
    Omega = matrix(Omega)
  }
  p2 = -0.5*apply(proposed_y - y_hat, 1, function(x) crossprod(x, (Omega %*% x))) #t(x) %*% Omega %*% x
  l2 = -0.5*apply(z - z_hat, 1, function(x) crossprod(x, (R_inv %*% x))) #t(x) %*% R_inv %*% x
  accept = mapply(MH, p1 = p1, l1 = l1, p2 = p2, l2 = l2)
  which_accepted = which(accept)
  y[which_accepted,] = proposed_y[which_accepted,]
  return(list(accept = accept, y = y, p2 = p2, l2 = l2))
}

# Computes the log marginal likelihood for accepting/rejecting BART trees
# log_marginal_likelihood <- function(node_partial_res, kappa, omega, mu0, Vinv, alpha) {
#   n = nrow(node_partial_res)
#   p = ncol(node_partial_res)
#   C = sum(apply(node_partial_res, 1, function(x) crossprod(crossprod(omega, x), x)))
#
#   vec = omega %*% colSums(node_partial_res) + kappa*mu0
#   mat = n*omega + diag(kappa, p)
#   vec2 = solve(mat, vec)
#   A = crossprod(crossprod(mat, vec2), vec2)
#
#   if(p==1){
#     det_omega = omega
#   } else {
#     det_omega = det(omega)
#   }
#
#   B = sum(diag(Vinv %*% omega))
#   return((n + alpha - p - 1)/2 * log(det_omega) - 0.5*(B + C - A))
# }
log_marginal_likelihood <- function(node_partial_res, kappa, omega, mu0, Vinv, alpha) {
  n = nrow(node_partial_res)
  p = ncol(node_partial_res)
  # C = sum(apply(node_partial_res, 1, function(x) crossprod(crossprod(omega, x), x)))
  C = sum(rowSums((node_partial_res %*% t(omega)) * node_partial_res))

  vec = omega %*% colSums(node_partial_res) + kappa*mu0
  mat = n*omega + diag(kappa, p)
  vec2 = solve(mat, vec)
  A = crossprod(vec2, vec)

  if(p==1){
    det_omega = omega
    det_mat = mat
  } else {
    det_omega = det(omega)
    det_mat = det(mat)
  }
  return(n*det_omega/2 + 0.5*kappa^p - 0.5 * det_mat - 0.5*(C - A))
}

# Compute tree priors at the node level
node_priors <- function(depth, prior_alpha, prior_beta) {
  return(prior_alpha * (1 + depth)^(-prior_beta))
}

# Compute tree priors
tree_priors <- function(nodes, parents, depth, prior_alpha, prior_beta) {
  depth <- as.numeric(depth)
  node_prob <- rep(NA, length=length(nodes))
  for(i in seq_along(nodes)) {
    node_prob[i] <- log(ifelse(nodes[i] %in% parents, node_priors(depth[i], prior_alpha, prior_beta), 1-node_priors(depth[i], prior_alpha, prior_beta)))
  }
  return(sum(node_prob))
}

# Unscale data (data was scaled to [-0.5,0.5])
#' Title
#'
#' @param scaled_val matrix of scaled values
#' @param min original min. value
#' @param max original max. value
#'
#' @return
#' @export
#'
#' @examples
unscale = function(scaled_val, min, max){
  if(is.matrix(scaled_val)){
    p = ncol(scaled_val)
    n = nrow(scaled_val)
    unscaled = matrix(nrow=n, ncol=p)
    for(j in 1:p){
      unscaled[,j] = (scaled_val[,j] + 0.5) * (max[j] - min[j]) + min[j]
    }
  } else {
    unscaled = (scaled_val + 0.5) * (max - min) + min
  }
  return(unscaled)
}

scale_bart = function(data, min = NULL, max = NULL){
  if(is.matrix(data)){
    if(is.null(min)) min = apply(data, 2, min, na.rm=TRUE)
    if(is.null(max)) max = apply(data, 2, max, na.rm=TRUE)
  } else {
    if(is.null(min)) min = min(data, na.rm=TRUE)
    if(is.null(max)) max = max(data, na.rm=TRUE)
  }
  scaled_val = (data - min)/(max - min) - 0.5
  return(scaled_val)
}

pdp_param_mat_list = function(x, y_range = c(-0.5, 0.5), grid_len = 20, intercept = FALSE, include_x = TRUE, n){ # Returns a list of size n (n pdp lines)
  y_grid = seq(y_range[1], y_range[2], length = grid_len)
  if(include_x){
    n = nrow(x)
    q = ncol(x)
    param_mat_list = vector(mode = "list", length = n)
    for(i in 1:n){
      param_mat_list[[i]] = cbind(matrix(rep(x[i,], grid_len), ncol = q, byrow = TRUE), y_grid)
      if(intercept) param_mat_list[[i]] = cbind(rep(1, grid_len), param_mat_list[[i]])
    }
  } else {
    if(missing(n)) {
      warning("Need to provide length of data")
      break
    }
    param_mat_list = matrix(y_grid, ncol=1)
    if(intercept) param_mat_list = cbind(rep(1, grid_len), param_mat_list[[i]])
  }
  return(param_mat_list)
}

plot_posterior = function(actual, post_list, q = c(0.025, 0.975), row_names = c()){
  if(missing(actual)){
    actual = Reduce("+", post_list)/length(post_list)
  }
  if(isSymmetric(actual)){
    post = Reduce(rbind, lapply(post_list, function(x) x[upper.tri(x, diag=TRUE)]))
    perc = apply(post, 2, quantile, q)
    predicted = colMeans(post)
    actual = matrix(actual[upper.tri(actual, diag = TRUE)], nrow=1)
  } else {
    perc = apply(Reduce(rbind, lapply(post_list, as.vector)), 2, quantile, q)
    predicted = as.vector(Reduce("+", post_list)/length(post_list))
  }
  data = data.frame("actual" = as.vector(actual), "predicted" = predicted, "lower" = perc[1,], "upper" = perc[2,])
  data$group = rep(seq(1:nrow(actual)), ncol(actual)) #rep(seq(1:(nrow(data)/2)), 2)
  if(is.null(row_names)){
    for(j in 1:ncol(actual)){
      for(i in 1:nrow(actual)){
        row_names = c(row_names, paste("b", i,j,sep = ""))
      }
    }
  }
  data$row.names = factor(row_names, levels = row_names)

  plot = ggplot(data, aes(x = factor(row.names))) +
    geom_point(aes(y = actual), col="black", size=1.3) +
    # geom_crossbar(aes(ymin = lower, ymax = upper, colour=as.factor(group)), width = 0.2) +
    geom_errorbar(aes(ymin = lower, ymax = upper, colour=as.factor(group)), width=0.5) +
    labs(x = "Parameters",
         y = "Values") +
    # scale_color_brewer(palette="RdYlGn") +
    theme_bw() +
    coord_flip()
  return(plot)
}

