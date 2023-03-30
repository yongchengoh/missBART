#' A function for proposing an updated tree structure
#'
#' @param df dataframe
#' @param x covariates
#' @param min_node minimum number of observations that should fall into a single terminal node.
#' @param max_attempt maximum number of attempts to find a suitable tree
#' @param i iteration number
#' @param probit logical. Indicates whether tree is for continuous or binary outcomes.
#' @param miss_row index of missing rows
#' @param cat_list list of used categorical subset
#' @importFrom purrr "rbernoulli"
#'
#' @return A dataframe containing details of a new tree structure
#' @export
#'
#' @examples
#' # Create a dataframe representing the root node of a tree
#' df <- data.frame(matrix(ncol = 7, nrow = 1))
#' colnames(df) <- c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction")
#' df[1,] <- c(0,0,1,0,1,0,0)
#' x <- matrix(stats::runif(9), ncol=3)
#' propose_tree(df, x, min_node = 1, max_attempt = 1, i = 1)
propose_tree = function(df, x, min_node, max_attempt = 50, i, probit = FALSE, miss_row = NA, cat_list = NA) {
  n = nrow(x)
  x_vars = seq_len(ncol(x))
  if(nrow(df) == 1) {
    #-If tree is only a stump, only grow it.
    MOVE = "GROW"
    # grow_variable = sample(x_vars, 1)
    # grow_node = 1
    # if(class(x[, grow_variable])=="factor") {
    #   lower = upper = NA
    # } else {
    #   lower = min(x[,grow_variable], na.rm = TRUE)
    #   upper = max(x[,grow_variable], na.rm = TRUE)
    # }
    # df$split_variable[1] = grow_variable
  } else if(i<=5){
    MOVE = "GROW"
  } else if(max(df$depth) == 3) {
    #-If tree is not a stump, we can randomly choose to grow or prune the tree
    MOVE = sample(c("GROW", "PRUNE"), 1)
  } else {
    # MOVE = sample(c("GROW", "PRUNE", "CHANGE", "SWAP"), 1)
    MOVE = sample(c("GROW", "PRUNE"), 1)
  } # End grow/prune/change/swap

  decent_tree = FALSE
  attempt = 1
  while(isFALSE(decent_tree) && attempt <= max_attempt) {
    ### Updating the proposed tree, "new_df"
    switch(EXPR=MOVE,
           GROW= {
             # if(nrow(df) != 1) {
               #-If GROW, we grow the tree from a randomly picked terminal node.
               terminal_nodes = as.numeric(setdiff(row.names(df), df$parent))
               grow_node = sample(terminal_nodes, 1)
               bad_grow_variable = TRUE
               while(bad_grow_variable){
                 grow_variable = sample(x_vars, 1)
                 # type = class(x[,grow_variable])
                 parent = grow_node # This is NOT the parent node. Used to find the splitting boundary for current splitting variable (depends on previous splits on the same variable)
                 while(df$split_variable[parent] != grow_variable) {
                   parent = df$parent[parent]
                   if(parent == 0) {
                     break
                   }
                 } #End while(df$split_variable[parent] != grow_variable)

                 if(is.numeric(x[,grow_variable])){
                    if(parent != 0) {
                      lower = ifelse(df$direction[parent] == 0, df$lower[parent], df$split_value[parent])
                      upper = ifelse(df$direction[parent] == 0, df$split_value[parent], df$upper[parent])
                    } else {
                      lower = min(x[,grow_variable], na.rm = TRUE)
                      upper = max(x[,grow_variable], na.rm = TRUE)
                    }
                    # x_set = x[which(x[,grow_variable] < upper & x[,grow_variable] >= lower), grow_variable]
                    # new_point = ifelse(length(x_set)==0, NA, sample(x_set, 1))
                    new_point = stats::runif(1, lower, upper)
                 } else {
                    #unique(x[,grow_variable])[as.logical(rbern(length(unique(x[,grow_variable]))))]
                    categories = unique(x[,grow_variable]) #We still need to check which node we're splitting on, whether the variable has already been split on
                    subset = categories[rbernoulli(length(categories), p = 0.5)]
                    lower = NA
                    upper = NA
                    # unique(x[,grow_variable])[as.logical(rbern(length(unique(x[,grow_variable]))))]
                    # split_set = setdiff(unique(x[,grow_variable]), unique(df$split_value[df$split_variable==grow_variable]))
                    # new_point = sample(split_set,sample(1:length(split_set), 1)) #ifelse(length(split_set) == 1, NA, sample(split_set, 1)) #setequal(unique(x[,grow_variable]), unique(df$split_value[df$split_variable==grow_variable]))
                 }
                 if(!is.na(new_point)) {bad_grow_variable = FALSE}
               } #End while(bad_grow_variable)
             # } else { #This is if nrow(df) == 1
             #   new_point <- ifelse(class(x[,grow_variable])=="numeric", sample(unique(x[,grow_variable]), 1), sample(unique(x[,grow_variable]), 1))
             # }
             new_depth = df$depth[grow_node] + 1
             if(any(is.na(x[,grow_variable]))) {
               NA_split_direc = sample(c(0,1), 2)
             } else {
               NA_split_direc = rep(NA, 2)
             }
             new_df = rbind(df, c(grow_node, lower, upper, grow_variable, new_point, new_depth, 0, NA_split_direc[1]), c(grow_node, lower, upper, grow_variable, new_point, new_depth, 1, NA_split_direc[2]))
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
               lower = min(x[,change_variable], na.rm = TRUE)
               upper = max(x[,change_variable], na.rm = TRUE)
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
    if(all(probit, decent_tree)) decent_tree = !Reduce(any, lapply(split(miss_row, change_points), function(x) all(isTRUE(x))))
    # decent_tree = all(tabulate(change_points) >= min_node)
    attempt = attempt + 1
  } # End while loop
  return(list("new_df"=new_df, "MOVE"=MOVE, "change_points"=change_points, "decent_tree"=decent_tree))
}


