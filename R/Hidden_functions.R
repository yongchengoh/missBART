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

sim_omega = function(y, y_hat, alpha, Vinv){
  n = nrow(y)
  df = alpha + n
  scale = crossprod(y-y_hat) + Vinv
  return(stats::rWishart(1, df, scale)[,,1])
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
update_z = function(Y, m, B, R){
  n = nrow(Y)
  mu = Y %*% B
  U = diag(1, n)
  z = matrnorm(mu, U, R)
  z[intersect(which(z<0), which(m==1))] = 0
  z[intersect(which(z>=0), which(m==0))] = 0
  return(z)
}

# MCMC sample of B, the matrix of probit parameters
update_B = function(x, y, z, tau_b = 0.01, include_x = TRUE, include_y = TRUE, center = TRUE){
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
  return(matrnorm(M, U, V))
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
        switch(EXPR = type,
               "numeric" = {
                 y_index[[count]] = if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] <= df$split_value[parent]) else which(x[,df$split_variable[parent]] > df$split_value[parent])
               }, "integer" = {
                 y_index[[count]] = if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] == df$split_value[parent]) else which(x[,df$split_variable[parent]] != df$split_value[parent])
               }, "factor" = {
                 y_index[[count]] = if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] == df$split_value[parent]) else which(x[,df$split_variable[parent]] != df$split_value[parent])
               }) #End switch
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
  }
  return(change_points)
}

# Computes the log marginal likelihood for accepting/rejecting BART trees
log_marginal_likelihood <- function(node_partial_res, kappa, omega) {
  n = nrow(node_partial_res)
  p = ncol(node_partial_res)
  rbar <- colMeans(node_partial_res)
  omega_mu = diag(kappa, p) + n*omega
  mu_mu = n * (chol2inv(PD_chol(omega_mu)) %*% omega %*% rbar)
  if(p==1){
    det_omega = omega
    det_omega_mu = omega_mu
  } else {
    det_omega = det(omega)
    det_omega_mu = det(omega_mu)
  }
  return(n/2 * log(det_omega) - 0.5 * log(det_omega_mu) - 0.5 * (t(mu_mu) %*% omega_mu %*% mu_mu + sum(apply(node_partial_res, 1, function(x) t(x) %*% omega %*% x))))
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
    node_prob[i] <- log(ifelse(nodes[i] %in% parents, node_priors(depth[i], prior_alpha, prior_beta), 1 - node_priors(depth[i], prior_alpha, prior_beta)))
  }
  return(sum(node_prob))
}
