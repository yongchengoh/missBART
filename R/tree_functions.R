#' Title
#'
#' @param df dataframe
#' @param x covariates
#' @param x_vars vector
#' @param min_node minimum number of observations that should fall into a single terminal node.
#' @param max_attempt maximum number of attempts to find a suitable tree
#' @param i iteration number
#'
#' @return
#' @export
#'
#' @examples
propose_tree = function(df, x, x_vars, min_node, max_attempt = 100, i) {
  n = nrow(x)
  if(nrow(df) == 1) {
    #-If tree is only a stump, only grow it.
    MOVE = "GROW"
    grow_variable = sample(x_vars, 1)
    grow_node = 1
    if(inherits(class(x[, grow_variable]), "factor")) {
      lower = upper = NA
    } else {
      lower = min(x[,grow_variable])
      upper = max(x[,grow_variable])
    }
    df$split_variable[1] = grow_variable
  } else if(i<=5){
    MOVE = "GROW"
  } else if(max(df$depth) == 1) {
    #-If tree is not a stump, we can randomly choose to grow or prune the tree
    MOVE = sample(c("GROW", "PRUNE"), 1)
  } else {
    MOVE = sample(c("GROW", "PRUNE", "CHANGE", "SWAP"), 1)
  } # End grow/prune/change/swap

  decent_tree = FALSE
  attempt = 1
  while(isFALSE(decent_tree) && attempt <= max_attempt) {
    ### Updating the proposed tree, "new_df"
    switch(EXPR=MOVE,
           GROW= {
             if(nrow(df) != 1) {
               #-If GROW, we grow the tree from a randomly picked terminal node.
               terminal_nodes = as.numeric(setdiff(row.names(df), df$parent))
               grow_node = sample(terminal_nodes, 1)
               bad_grow_variable = TRUE
               while(bad_grow_variable){
                 grow_variable = sample(x_vars, 1)
                 type = class(x[,grow_variable])
                 parent = grow_node # This is NOT the parent node. Used to find the splitting boundary for current splitting variable (depends on previous splits on the same variable)
                 while(df$split_variable[parent] != grow_variable) {
                   parent = df$parent[parent]
                   if(parent == 0) {
                     break
                   }
                 } #End while(df$split_variable[parent] != grow_variable)

                 switch(EXPR = type,
                        "numeric" = {
                          if(parent != 0) {
                            lower = ifelse(df$direction[parent] == 0, df$lower[parent], df$split_value[parent])
                            upper = ifelse(df$direction[parent] == 0, df$split_value[parent], df$upper[parent])
                          } else {
                            lower = min(x[,grow_variable], na.rm = TRUE)
                            upper = max(x[,grow_variable], na.rm = TRUE)
                          }
                          new_point = stats::runif(1, lower, upper)
                        }, "integer" = {
                          lower = 0
                          upper = 1
                          new_point = ifelse(grow_variable %in% df$split_variable, NA, sample(x[,grow_variable], 1))
                        }, "factor" = {
                          lower = NA
                          upper = NA
                          split_set = setdiff(unique(x[,grow_variable]), unique(df$split_value[df$split_variable==grow_variable]))
                          new_point = ifelse(length(split_set == 1), NA, sample(split_set, 1)) #setequal(unique(x[,grow_variable]), unique(df$split_value[df$split_variable==grow_variable]))
                        }) # End switch
                 if(!is.na(new_point)) {bad_grow_variable = FALSE}
               } #End while(bad_grow_variable)
             } else { #This is if nrow(df) == 1
               new_point <- ifelse(class(x[,grow_variable])=="numeric", stats::runif(1, lower, upper), sample(unique(x[,grow_variable]), 1))
             }
             # new_point = sample(x[which(x[,grow_variable]<upper & x[,grow_variable]>=lower),grow_variable], 1)
             new_depth = df$depth[grow_node] + 1
             new_df = rbind(df, c(grow_node, lower, upper, grow_variable, new_point, new_depth, 0), c(grow_node, lower, upper, grow_variable, new_point, new_depth, 1))
           }, PRUNE= {
             #-If tree only has one level (one split point), we prune back to a "stump" and the full dataset is used.
             #-Otherwise, we prune off a randomly selected pair of terminal nodes belonging to the same parent
             parents = unique(df$parent[df$parent != 0])
             children = c()
             prune_parent = c()
             #-Iterate through list of parent nodes, find the 2 children of that parent
             # If either child is a parent, we cannot prune that parent.
             for(k in seq_along(parents)) {
               children = which(df$parent == parents[k])
               if(isFALSE(children[1] %in% parents || children[2] %in% parents)) {
                 prune_parent = c(prune_parent, df$parent[children[1]])
               }
             }
             if(length(prune_parent) == 1) {
               prune = prune_parent
             } else {
               prune = sample(prune_parent, size=1)
             }
             new_df = df[!df$parent==prune,] # remove the children of the selected parent
             # Rename the parents
             if(nrow(new_df) > 1) {
               for(k in 2:nrow(new_df)) {
                 new_df$parent[k] = which(row.names(new_df) == new_df$parent[k])
               }
               row.names(new_df) = seq_len(nrow(new_df))
             }
           }, CHANGE= {
             # If CHANGE, we choose an internal node (i.e. a parent node) and change the split variable & value by resampling uniformly
             change_node = sample(unique(df$parent[df$parent != 0]), 1)
             change_variable = sample(x_vars, 1)
             parent = change_node
             while(df$split_variable[parent] != change_variable) {
               parent = df$parent[parent]
               if(parent == 0) {
                 break
               }
             }
             if(parent != 0) {
               lower = ifelse(df$direction[parent] == 0, df$lower[parent], df$split_value[parent])
               upper = ifelse(df$direction[parent] == 0, df$split_value[parent], df$upper[parent])
             } else {
               lower = min(x[,change_variable])
               upper = max(x[,change_variable])
             }
             change_value = stats::runif(1, lower, upper)
             # change_value = sample(x[which(x[,change_variable]<upper & x[,change_variable]>=lower),change_variable],1)
             new_df = df
             new_df$lower[new_df$parent == change_node] = lower
             new_df$upper[new_df$parent == change_node] = upper
             new_df$split_variable[new_df$parent == change_node] = change_variable
             new_df$split_value[new_df$parent == change_node] = change_value
             if(change_node == 1) {
               new_df$split_variable[change_node] = change_variable
             }
           }, SWAP= {
             sampling = TRUE
             # attempts = 1
             while(isTRUE(sampling)) {
               swap_node_parent = sample(unique(df$parent[df$parent != 0]), 1)
               swap_node_child = sample(which(df$parent == swap_node_parent), 1)
               sampling = !(swap_node_child %in% unique(df$parent[df$parent != 0]))
             }
             new_df = df
             new_df[which(df$parent == swap_node_parent), 2:5] = df[df$parent == swap_node_child,  2:5]
             new_df[which(df$parent == swap_node_child),  2:5] = df[df$parent == swap_node_parent, 2:5]
           } # End SWAP
    ) # End switch

    change_points = get_change_points(new_df, x)
    empty_terminal = length(unique(change_points)) != length(as.numeric(setdiff(row.names(new_df), new_df$parent)))
    decent_tree = isFALSE(empty_terminal) && all(table(change_points) >= min_node)
    # if(probit & decent_tree) decent_tree = !Reduce(any, lapply(split(miss_row, change_points), function(x) all(isTRUE(x))))
    # decent_tree = all(tabulate(change_points) >= min_node)
    attempt = attempt + 1
  } # End while loop
  return(list("new_df"=new_df, "MOVE"=MOVE, "change_points"=change_points, "decent_tree"=decent_tree))
}

#' Get vector for sorting observations into terminal nodes
#'
#' @param df dataframe
#' @param x covariates
#'
#' @return a vector
#' @export
#'
#' @examples
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

#' Computes the log marginal likelihood for accepting/rejecting BART trees
#'
#' @param node_partial_res matrix
#' @param kappa scalar
#' @param omega matrix
#'
#' @return
#' @export
#'
#' @examples
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

#' Compute tree priors at the node level
#'
#' @param depth scalar
#' @param prior_alpha scalar
#' @param prior_beta scalar
#'
#' @return
#' @export
#'
#' @examples
node_priors <- function(depth, prior_alpha, prior_beta) {
  return(prior_alpha * (1 + depth)^(-prior_beta))
}

#' Compute tree priors
#'
#' @param nodes nodes of trees
#' @param parents parents of nodes
#' @param depth current depth of tree
#' @param prior_alpha hyperparameter
#' @param prior_beta hyperparameter
#'
#' @return
#' @export
#'
#' @examples
tree_priors <- function(nodes, parents, depth, prior_alpha, prior_beta) {
  depth <- as.numeric(depth)
  node_prob <- rep(NA, length=length(nodes))
  for(i in seq_along(nodes)) {
    node_prob[i] <- log(ifelse(nodes[i] %in% parents, node_priors(depth[i], prior_alpha, prior_beta), 1 - node_priors(depth[i], prior_alpha, prior_beta)))
  }
  return(sum(node_prob))
}
