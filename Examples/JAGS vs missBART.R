rm(list=ls())
devtools::load_all()

## SIMULATE DATA
# set.seed(5555)
set.seed(12345)
n = 500
a = 1
b = 2
prec = 1
x = x_original = matrix(runif(n = n, min = 0, max = 3), nrow = n)
y = y_original = a + b*x + rnorm(n, 0, sqrt(1/prec))
x_onehot = cbind(x<=1, x>1 & x<=2, x>2)*1

## SCALE DATA & SET SD/PRECISION HYPERPARAMETERS
min_y = min(y,na.rm=TRUE)
max_y = max(y,na.rm=TRUE)
y = y_original = (y - min_y)/(max_y - min_y) - 0.5

## FIX TREE
df = data.frame(matrix(ncol = 8, nrow = 5))
colnames(df) = c("parent", "lower", "upper", "split_variable", "split_value", "depth", "direction", "NA_direction")
df[1,] = c(0,0,1,0,1,0,0,NA)
df[2,] = c(1, 0, 3, 1, 1, 1, 0, NA)
df[3,] = c(1, 0, 3, 1, 1, 1, 1, NA)
df[4,] = c(3, 1, 3, 1, 2, 2, 0, NA)
df[5,] = c(3, 1, 3, 1, 2, 2, 1, NA)
fix_tree = df
fix_mu = matrix(c(mean(y[x<=1], na.rm = TRUE), mean(y[x>1 & x<=2], na.rm = TRUE), mean(y[x>2], na.rm = TRUE)), ncol = 1)

## SIMULATE MISSINGNESS
gam = 0.5
del = 0
ome = -1
p_m = pnorm(gam + del*(x) + ome*(y))
m = matrix(rbinom(n, 1, p_m), ncol = 1)
y[m==0] = NA
plot(y_original, m)
points(y_original, p_m)
range(y, na.rm=TRUE)
range(y_original)

missing_index = which(is.na(y))
obs_index = which(!is.na(y))

plot(x[obs_index], y[obs_index], ylim = range(y_original))
points(x[missing_index], y_original[missing_index], col = "red")

nu = 10
sigest = summary(lm(y~x))$sigma
qchi = qchisq(1-0.75, nu)
lambda = (sigest^2)*qchi/nu

## RUN JAGS
library(runjags)
modelString = "model{
  for (i in 1:n){
    y[i] ~ dnorm(mu[i], 1/sd^2)
    m[i] ~ dbin(p[i],1)
    probit(p[i]) = gamma + omega*(y[i])
  }

  mu = x_onehot %*% mu_vec
  sd = 1/sqrt(tau)
  y_miss = y[missing_index]

  # Priors
  for(j in 1:3){
    mu_vec[j] ~ dnorm(0, 16)
  }
  gamma ~ dnorm(0, 0.01)
  omega ~ dnorm(0, 0.01)
  tau ~ dgamma(nu/2, lambda*nu/2)
}"
my_data = list("y" = as.vector(y),
               "n" = n, "m" = as.vector(m), "x_onehot" = x_onehot,
               "nu" = nu, "lambda" = lambda, "missing_index" = missing_index)
posterior = run.jags(model = modelString,
                     data = my_data,
                     monitor = c("gamma", "omega", "sd", "mu_vec", "y_miss"),
                     n.chains = 1,
                     burnin = 2000,
                     sample = 1000,
                     thin = 1)

# modelString = "model{
#   for (i in 1:n){
#     y[i] ~ dnorm(alpha + beta*x[i], 1/sd^2)
#     m[i] ~ dbin(p[i],1)
#     probit(p[i]) = gamma + omega*(y[i])
#   }
#
#   sd = 1/sqrt(tau)
#   y_miss = y[missing_index]
#
#   # Priors
#   alpha ~ dnorm(0, 0.01)
#   beta ~ dnorm(0, 0.01)
#   gamma ~ dnorm(0, 0.01)
#   omega ~ dnorm(0, 0.01)
#   tau ~ dgamma(nu/2, lambda*nu/2)
# }"
#
# my_data = list("y" = as.vector(y),
#                "n" = n, "m" = as.vector(m), "x" = as.vector(x),
#                "nu" = nu, "lambda" = lambda, "missing_index" = missing_index)
#
# posterior = run.jags(model = modelString,
#                      data = my_data,
#                      monitor = c("gamma", "omega", "sd", "alpha", "beta", "y_miss"),
#                      n.chains = 1,
#                      burnin = 500,
#                      sample = 1000,
#                      thin = 2)


## Plot sd on original scale
post_sd = posterior$mcmc[[1]][,3]*(max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
plot(post_sd)
mean(post_sd)

## Plot missing y's on original scale
y_jags = (colMeans(posterior$mcmc[[1]][,-c(1:6)]) + 0.5) * (max(y, na.rm=TRUE) - min(y, na.rm=TRUE)) + min(y, na.rm=TRUE)
plot(y_original[missing_index], y_jags)
abline(coef = c(0,1))

## MODEL 1
devtools::load_all()
model1 = missBARTprobit(x, y, predict = FALSE, n_trees = 50,
                        burn = 2000, iters = 1000, thin = 1, x_predict = x,
                        show_progress = TRUE, scale = FALSE, make_pdp = FALSE,
                        pdp_range = c(pdp_range1, pdp_range2),
                        include_x = FALSE, include_y = TRUE,
                        mice_impute = TRUE,
                        hypers = hypers_list(df = 10, q = 0.75),
                        true_trees_data = list(fix_tree))

plot(unlist(model1$omega_post)*(max_y - min_y), type = "l", ylab = "sd - BART + probit regression", main = "Linear Data, Probit Missingness")
abline(c(mean(unlist(model1$omega_post)*(max_y - min_y)), 0))
abline(c(1, 0), col = "red")
mean(unlist(model1$omega_post)*(max_y - min_y))

plot(y_original[missing_index], Reduce("+", model1$y_impute)/length(model1$y_impute), ylab = "y_missBART[missing_index]")
abline(coef = c(0,1))
sqrt(sum((y_original[missing_index] - Reduce("+", model1$y_impute)/length(model1$y_impute))^2)/n)

Reduce("+", model1$B_post)/length(model1$B_post)
colMeans(posterior$mcmc[[1]][,c(1,2)])

## MODEL 2
model2 = missBART2(x, y, predict = FALSE, n_reg_trees = 50, n_class_trees = 50,
                    burn = 2000, iters = 1000, thin = 1, x_predict = x,
                    scale = FALSE, show_progress = TRUE, progress_every = 10,
                    pdp_range = c(pdp_range1, pdp_range2), make_pdp = FALSE,
                    include_x = FALSE, include_y = TRUE, MH_sd = 0.2,
                    mice_impute = TRUE,
                    hypers = hypers_list(df = 10, q = 0.75),
                    true_trees_data = list(fix_tree))
plot(unlist(model2$omega_post)*(max_y - min_y), type = "l", ylab = "sd - BART + probit BART", main = "Linear Data, Probit Missingness")
abline(c(mean(unlist(model2$omega_post)*(max_y - min_y)), 0))
abline(c(1, 0), col = "red")
mean(unlist(model2$omega_post)*(max_y - min_y))

plot(y_original[missing_index], Reduce("+", model2$y_impute)/length(model2$y_impute), xlim = c(-1, 10), ylim = c(-1, 10), ylab = "y_missBART2[missing_index]")
abline(coef = c(0,1))
sqrt(sum((y_original[missing_index] - Reduce("+", model2$y_impute)/length(model2$y_impute))^2)/n)
