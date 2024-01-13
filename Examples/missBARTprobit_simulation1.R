rm(list=ls())
devtools::load_all()

## SIMULATE DATA
# set.seed(5555)
set.seed(10000)
n <- 500
data <- sim_data_trees(n, 1, 1, ome = diag(0.5, 1)) #mlbench::mlbench.friedman1(n, sd = 5)
y <- y_original = matrix(data$y, ncol = 1)
x <- data$x
true_sd <- 1/sqrt(data$ome)
true_sd
range(y)

# y <- y_original <- scale_bart(y)

# SIMULATE MISSINGNESS
missing <- sim_missing_probit(x, y, include_x = TRUE, include_y = FALSE, min_missing_prop = 0.7, max_missing_prop = 0.9)
m <- missing$m
y <- missing$missing_y
missing$B

range(y_original)
range(y, na.rm=TRUE)

missing_index <- which(is.na(y))
obs_index <- which(!is.na(y))

plot(x[obs_index], y_original[obs_index], ylim = range(y_original))
points(x[missing_index], y_original[missing_index], col = "red")
# points(x, y_hat)

# TRAIN/TEST SPLIT
split_prop <- 0.6
n_train <- round(n * split_prop)
n_test <- n - n_train
train_id <- sample(seq_len(n), n_train)
test_id <- seq_len(n)[-train_id]
x_train <- x[train_id,,drop=FALSE]
y_train <- y[train_id,, drop=FALSE]
x_test <- x[test_id,,drop=FALSE]
y_test <- y[test_id,,drop=FALSE]

## MODEL 1
model1 <- missBARTprobit(x = x_train, y = y_train, predict = TRUE, n_trees = 90,
                         burn = 500, iters = 1000, thin = 2, x_predict = x_test,
                         show_progress = TRUE, scale = TRUE, make_pdp = FALSE,
                         pdp_range = c(pdp_range1, pdp_range2),
                         include_x = TRUE, include_y = FALSE,
                         mice_impute = TRUE,
                         hypers = hypers_list(df = 10, q = 0.75),
                         true_trees_data = list(fix_tree))

plot(unlist(model1$omega_post), type = "l", ylab = "sd - BART + probit regression", main = "Linear Data, Tree Missingness")
abline(c(mean(unlist(model1$omega_post)), 0))
# abline(c(1/(range(data$y)[2] - range(data$y)[1]), 0), col = "red")
abline(c(true_sd, 0), col = "red")
mean(unlist(model1$omega_post))

# GET OOS RMSE & PLOT POSTERIOR PREDICTIONS
new_y_post <- Reduce(cbind, model1$new_y_post)
mean_new_y_post <- rowMeans(new_y_post)
pred_int <- apply(new_y_post, 1, quantile, c(0.25, 0.75))
model1_data <- data.frame(pred = mean_new_y_post, lower = pred_int[1,], upper = pred_int[2,], true = y_original[test_id])

rmse <- sum(sqrt((model1_data$true - model1_data$pred)^2)/nrow(model1_data))

sum(apply(model1_data, 1, function(x) x[4] <= x[3] & x[4] >= x[2]))/length(test_id)

ggplot(data = model1_data, aes(true, pred)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  xlim(min(model1_data), max(model1_data)) +
  ylim(min(model1_data), max(model1_data))

quantile(unlist(model1$omega_post), c(0.25, 0.75))
true_sd <= quantile(unlist(model1$omega_post), 0.75) & true_sd >= quantile(unlist(model1$omega_post), 0.25)

# plot(y_original[missing_index], Reduce("+", model1$y_impute)/length(model1$y_impute), xlim = range(y_original), ylim = range(y_original), ylab = "y_missBART[missing_index]")
# abline(coef = c(0,1))
# sqrt(sum((y_original[missing_index] - Reduce("+", model1$y_impute)/length(model1$y_impute))^2)/n)
#
# ## MODEL 2
# devtools::load_all()
# model2 <- missBART2(x, y, predict = FALSE, n_reg_trees = 100, n_class_trees = 100, scale = TRUE,
#                     burn = 1000, iters = 1000, thin = 2, x_predict = x,
#                     show_progress = TRUE, progress_every = 10,
#                     pdp_range = c(pdp_range1, pdp_range2), make_pdp = FALSE,
#                     include_x = TRUE, include_y = TRUE, MH_sd = 0.1,
#                     mice_impute = TRUE,
#                     hypers = hypers_list(df = 10, q = 0.75),
#                     true_trees_data = list(fix_tree),
#                     true_trees_missing = list(fix_tree_miss),
#                     z_true = matrix(z_true, ncol=1),
#                     true_change_points_miss = get_change_points(df_miss, y_original), true_change_points = matrix(get_change_points(df, x), ncol=1))
# plot(unlist(model2$omega_post), type = "l", ylab = "sd - BART + probit BART", main = "Linear Data, Tree Missingness")
# abline(c(mean(unlist(model2$omega_post)), 0))
# abline(c(5/(range(data$y)[2] - range(data$y)[1]), 0), col = "red")
# mean(unlist(model2$omega_post))
#
# colMeans(Reduce(rbind,lapply(model2$class_mu, unlist)))
# colMeans(posterior$mcmc[[1]][,c(5:7)])
# fix_mu_miss
# plot(Reduce("+", model2$z_post)/length(model2$z_post), z_true)
#
# plot(Reduce("+", model2$y_impute)/length(model2$y_impute), y_jags)
# abline(coef = c(0,1))
#
# Reduce("+", lapply(model2$reg_mu, function(x) Reduce(cbind, x)))/length(model2$reg_mu)
# colMeans(posterior$mcmc[[1]][,c(2:4)])
#
# plot(y_original[missing_index], Reduce("+", model2$y_impute)/length(model2$y_impute), xlim = range(y_original), ylim = range(y_original), ylab = "y_missBART2[missing_index]")
# abline(coef = c(0,1))
# sqrt(sum((y_original[missing_index] - Reduce("+", model2$y_impute)/length(model2$y_impute))^2)/n)
