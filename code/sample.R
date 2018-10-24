rm(list = ls(all = TRUE))
source("model.R")
data(iris)
x = as.matrix(iris[, 1:4])
fit = gaussian_mixture_model(x = x, K = 3)
z_true = iris[, 5]
table(fit$z_map, z_true)
