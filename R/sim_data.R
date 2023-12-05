#---------- Simulating Friedman data ----------#
#' Simulate
#'
#' @param n number of observations
#' @param p number of responses
#' @param scale_par scale_par
#' @param omega_diag precision ~ rWishart(1, p+1, diag(omega_diag, p))
#'
#' @return
#' @export
#' @importFrom MASS mvrnorm
#'
#' @examples
#' data = sim_data_friedman(n = 1000, p = 3)
#' true_trees_data = data$true_trees
#' y = y_original = data$y
#' x = data$x
#' ome = data$ome
sim_data_friedman = function(n, p = 1, scale_par = 1, omega_diag = 1, Omega = NULL) {
  # Simulate some data using a multivariate version of Friedman
  # y = 10sin(πx1x2)+20(x3−0.5)2+10x4+5x5+ε
  X = matrix(NA, nrow = n, ncol = 5)
  for(i in 1:ncol(X)) X[,i] = stats::runif(n) #stats::rnorm(n, 0, 1)
  # pars = matrix(stats::rnorm(5 * p, sd = scale_par), ncol = p, nrow = 5) # 5 parameters on p dimensions
  y = mean = matrix(NA, ncol = p, nrow = n)
  if((is.null(Omega))){
    Omega = stats::rWishart(1, p+1, diag(omega_diag, p))[,,1]
  }
  if(p > 1) {
    err = mvrnorm(n, mu=rep(0, ncol(Omega)), Sigma = solve(Omega))
  } else {
    err = matrix(stats::rnorm(n, sd = 1/sqrt(Omega)), ncol = 1)
  }
  for(j in 1:p) {
    # mean[,j] = pars[1,j]*sin(X[,1]*X[,2]) + pars[2,j] * (X[,3]-0.5)^2 + pars[3,j] * X[,4] + pars[5,j] * X[,5]
    mean[,j] = 10*sin(X[,1]*X[,2]) + 20 * (X[,3]-0.5)^2 + 10 * X[,4] + 5 * X[,5]
    y[,j] = mean[,j] + err[,j]
  }

  return(list(y = y, x = X, ome = Omega, mean = mean, p = p, q = 5))
}

#' Simulating data from trees
#'
#' @param n number of observations
#' @param p number of responses
#' @param q number of covariates
#' @param min_x min_x
#' @param max_x max_x
#' @param trees number of trees
#' @param ... Catches unused arguments
#'
#' @return
#' @export
#' @importFrom
#'
#' @examples
#' data = sim_data_trees(n = 100, p = 3, q = 4, trees = 6)
#' y = y_original = data$y
#' x = data$x
#' ome = data$ome
#' true_trees_data = data$true_trees
sim_data_trees = function(n, p, q, min_x = 0, max_x = 1, trees = 1, ome = NULL, kappa = NULL, splits = NULL, ...){
  x = matrix(stats::runif(n*q, min = min_x, max = max_x), ncol = q)
  sum_mu = matrix(0, ncol = p, nrow = n)
  if(is.null(ome)) ome = stats::rWishart(1, p+1, diag(1,p))[,,1]
  # ome = ifelse(is.null(ome), stats::rWishart(1, p+1, diag(1,p))[,,1], ome)  #stats::rWishart(1, p+1, diag(1,p))[,,1]
  if(!all(dim(ome) == p) || !is.matrix(ome)) stop("ome must be a pxp matrix")
  if(is.null(kappa)) kappa = 16*trees
  true_trees = vector(mode = "list", length = trees)
  n_splits = splits
  for(i in 1:trees){
    df1 <- data.frame(matrix(ncol = 8, nrow = 1))
    colnames(df1) = c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction", "NA_direction")
    df1[1,] <- c(0,0,1,0,1,0,0,NA)
    if(is.null(splits)) n_splits = sample(seq(1,5), 1)
    mu = multi_rMVN(matrix(0, ncol=p, nrow = n_splits+1), kappa*diag(1,p))
    for(j in 1:n_splits){
      new_tree = propose_tree(df1, x, min_node = 20, max_attempt = 10, i = 2)
      df1 = new_tree$new_df
    }
    true_trees[[i]] = df1
    sum_mu = sum_mu + mu[new_tree$change_points,,drop=FALSE]
  }
  y = multi_rMVN(mean_mat = sum_mu, precision = ome)
  return(list(y = y, x = x, ome = ome, true_trees = true_trees))
}

#' Simulating missing values from a probit regression model
#'
#' @param x regression covariates
#' @param y regression response
#' @param include_x logical. Include x in probit model?
#' @param include_y logical. Include y in probit model?
#' @param min_missing_prop minimum proportion of missingness
#' @param max_missing_prop maximum proportion of missingness
#' @param ... Catches unused arguments
#'
#' @return
#' @export
#'
#' @examples
sim_missing = function(x, y, include_x = FALSE, include_y = FALSE, min_missing_prop = 0.6, max_missing_prop = 0.9, ...){
  p = ncol(y)
  q = ncol(x)
  n = nrow(y)
  r = 1
  if(include_x & !include_y){
    r = 1 + q
  } else if (!include_x & include_y){
    r = 1 + p
  } else if (include_x & include_y){
    r = 1 + p + q
  }
  corR = diag(1, p)
  psd = FALSE
  while(!psd){
    corR[upper.tri(corR)] = sample(seq(-1,1,length=100), sum(seq_len(p-1)))
    corR[lower.tri(corR)] = t(corR)[lower.tri(corR)]
    psd = all(eigen(corR)$values >= 0)
  }

  Psi_1 = rInvWishart(1, r+1, diag(1,r))[,,1]

  for(seed in sample(seq(1000,100000), size = 10000)){
    set.seed(seed)

    B = matrnorm(matrix(0, nrow=r, ncol=p), Psi_1, corR)
    if(include_x & !include_y){
      phi = cbind(rep(1, n), x) %*% B
    } else if (!include_x & include_y){
      phi = cbind(rep(1,n), y) %*% B
    } else if (include_x & include_y){
      phi = cbind(rep(1,n), x, y) %*% B
    } else {
      phi = matrix(B, nrow = n, ncol = p, byrow = TRUE)
    }
    z_mod = multi_rMVN(phi, solve(corR))
    m = matrix(1, nrow=n, ncol=p)
    m[z_mod<=0] = 0

    mis_prop = colSums(m)/n
    if(all(mis_prop > min_missing_prop) & all(mis_prop < max_missing_prop)) break
  }

  y[m==0] = NA
  missing_prop = colSums(m)/n
  missing_index = which(is.na(y))
  obs_index = which(!is.na(y))

  return(list(B = B, m = m, missing_y = y, missing_prop = missing_prop, missing_index = missing_index, obs_index = obs_index, corR = corR, z_mod = z_mod))
}

#' Simulate missing data from a BART model
#'
#' @param x regression covariates
#' @param y regression response
#' @param trees number of trees
#' @param include_x logical. Include x in probit model?
#' @param include_y logical. Include y in probit model?
#' @param min_missing_prop minimum proportion of missingness
#' @param max_missing_prop maximum proportion of missingness
#' @param ... Catches unused arguments
#'
#' @return
#' @export
#'
#' @examples
sim_missing_trees = function(x, y, trees = 1, include_x = FALSE, include_y = TRUE, min_missing_prop = 0.6, max_missing_prop = 0.9, ...){
  p = ncol(y)
  q = ncol(x)
  n = nrow(x)
  if(include_x & !include_y){
    r = q
    Y = x
  } else if (!include_x & include_y){
    r = p
    Y = y
  } else if (include_x & include_y){
    r = p + q
    Y = cbind(x, y)
  }

  if(!exists('min_node')) min_node = round(n/100)

  kappa = 4/9 * trees
  sum_mu = matrix(0, ncol = p, nrow = n)
  true_trees = vector(mode = "list", length = trees)

  # decent = FALSE
  # decent_tree = c()
  # while(!decent){
  for(seed in sample(seq(1000,100000),size = 10000)){
    set.seed(seed)
    for(i in 1:trees){
      df2 = data.frame(matrix(ncol = 8, nrow = 1))
      colnames(df2) = c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction", "NA_direction")
      df2[1,] = c(0,0,1,0,1,0,0,NA)
      n_splits = 2
      # n_splits = sample(seq(2, 3), 1)
      for(j in 1:n_splits){
        new_tree = propose_tree(df2, Y, min_node = min_node, max_attempt = 10, i = 2)
        df2 = new_tree$new_df
      }
      # mu = multi_rMVN(matrix(0, ncol=p, nrow = n_splits+1), kappa*diag(1,p))
      mu = matrix(c(1, -0.7, 1), ncol=1)
      true_trees[[i]] = df2
      sum_mu = sum_mu + mu[new_tree$change_points,,drop=FALSE]
    }
    z_mod = multi_rMVN(mean_mat = sum_mu, precision = diag(1,p))
    m = matrix(1, nrow=n, ncol=p)
    m[z_mod<=0] = 0
    mis_prop = colSums(m)/n

    if(all(mis_prop > min_missing_prop) & all(mis_prop < max_missing_prop)) break
  }

  y[m==0] = NA
  missing_prop = colSums(m)/n
  missing_index = which(is.na(y))
  obs_index = which(!is.na(y))
  miss_row = apply(y, 1, function(x) any(is.na(x)))

  return(list(m = m, missing_y = y, missing_prop = missing_prop, missing_index = missing_index, obs_index = obs_index, true_trees = true_trees, z = z_mod))
}

