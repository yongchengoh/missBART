# devtools::install_github("yongchengoh/missBART")
# library(missBART)

data = sim_data_friedman(n = 1000, p = 5)
y = data$y
x = data$x
n = nrow(y)
p = ncol(y)
bart_out = mvBART(x, y, n_trees = 90, burn = 500, iters = 1000, thin = 2, predict = FALSE)
y_post = bart_out$y_pred
mean_y_post = Reduce("+", y_post)/length(y_post)
mean_y_post = unscale(mean_y_post, min = bart_out$min_y, max = bart_out$max_y)

library(ggplot2)
for(j in 1:ncol(y)){
  min = min(y[,j], mean_y_post[,j])
  max = max(y[,j], mean_y_post[,j])
  cheat = data.frame(y_seq = seq(min, max, length=n))
  plot_data = data.frame(y_original = y[,j], y = mean_y_post[,j])
  plot = ggplot(plot_data, aes(y_original, y)) +
    geom_point(size=1, alpha = 0.8) +
    theme_bw() +
    ylab(paste("Predicted Y",j,sep="")) +
    xlab(paste("Data Y",j,sep="")) +
    ggtitle("Multivariate BART") +
    geom_line(data=cheat, aes(y_seq, y_seq), colour="black", size=0.5)
  print(plot)
}


