# devtools::install_github("yongchengoh/missBART")
library(missBART)

set.seed(12345678)
#----- SIMULATE COMPLETE DATA
n <- 1000
p <- 3
data <- sim_data_friedman(n = n, p = p)
true_trees_data <- data$true_trees
y <- y_original <- data$y
x <- data$x
ome <- data$ome

#----- SIMULATE MISSINGNESS IN DATA
missing <- sim_missing_probit(x = x, y = y, include_x = TRUE, include_y = TRUE, max_missing_prop = 0.99, min_node = 10, min_missing_prop = 0.7, intercept = TRUE)
missing$missing_prop
y <- missing$missing_y
missing_index <- missing$missing_index
obs_index <- missing$obs_index
m <- missing$m
B <- missing$B
corR <- missing$corR

#----- PLOT DATA
for(j in 1:ncol(y)) {
  plot <- ggplot(data.frame(x = x[,1], y = y_original[,j]), aes(x, y, color=factor(m[,j]))) + geom_point() + ggtitle(paste("In-sample y", j, sep=""))
  print(plot)
}

model1 <- missBARTprobit(x, y, predict = FALSE)
y_pred <- model1$y_pred

mean_y_post <- Reduce("+", y_pred)/length(y_pred)
mean_y_post <- unscale(mean_y_post, min = model1$min_y, max = model1$max_y)

BART_probit_plot <- vector(mode = "list", length = p)
BART_probit_facet <- vector(mode = "list", length = p)
cheat <- vector("list", length = p)
colours <- c("#FD7446FF", "#709AE1FF", "#1A9993FF", "#FD8CC1FF", "#FED439FF")
for(j in seq_len(p)) {
  min <- min(y_original[,j], mean_y_post[,j])
  max <- max(y_original[,j], mean_y_post[,j])
  cheat[[j]] <- data.frame(y_seq = seq(min, max, length=n))
  plot_data <- data.frame(y_original = y_original[,j], y_post = mean_y_post[,j], Observed = factor(m[,j]))
  plot <- ggplot(plot_data, aes(y_original, y_post)) +
    geom_point(aes(color=Observed, shape = Observed), size=0.8, alpha = 0.8) +
    scale_shape_manual(name = "", labels = c("Missing", "Observed"),values=c(17, 19)) +
    scale_colour_manual(name = "", labels = c("Missing", "Observed"), values = c(colours[1], colours[3])) +
    theme_bw() +
    ylab(paste("Predicted Y",j,sep="")) +
    xlab(paste("Data Y",j,sep="")) +
    ggtitle("Multivariate BART and Probit Regression") +
    geom_line(data=cheat[[j]], aes(y_seq, y_seq), colour="black", size=0.5)

  BART_probit_plot[[j]] <- plot
  print(plot)

  facet_labels <- c("Missing", "Observed")
  names(facet_labels) <- c('0', '1')
  facet_plot <- plot + facet_grid(as.factor(m[,j]), labeller = as_labeller(facet_labels)) #vars(test$m)
  BART_probit_facet[[j]] <- facet_plot
  print(facet_plot)
}
