setwd("./R")

source("metrics.R")
source("help_functions.R")
source('crf_km.R')

list.of.packages <- c("ggplot2", "grf", "quantregForest", "randomForestSRC", "survival")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library(ggplot2)
library(quantregForest)
library(randomForestSRC)
library(survival)
library(grf)

# ------------------------------------------------------- #

attach(mtcars)
op <- par(mfrow=c(2,2), mar=c(2,2,1,1)+0.1, oma = c(0,0,0,0) + 0.1, pty="s")

# Load in the data
n <- 300
n_test <- 300

# training data
Xtrain <- Xtrain <- sort(runif(n = n, min = 0, max = 2*pi))
sigma <- 0.3
Ttrain <- sin(Xtrain) + rnorm(n, mean = 0, sd = sigma) + 2.5
ctrain <- sin(Xtrain) + rexp(n = n, rate = 0.2) + 1
Ytrain <- pmin(Ttrain, ctrain)
censorInd <- 1*(Ttrain <= ctrain)
data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)

# plot training data
plot(Xtrain, Ytrain, cex = 0.2)
points(Xtrain[!censorInd], Ytrain[!censorInd], type = 'p', col = 'red', cex = 0.3)
points(Xtrain[!censorInd], Ttrain[!censorInd], type = 'p', col = 'green', cex = 0.3)

# test data
Xtest <- sort(runif(n = n_test, min = 0, max = 2*pi))
Ytest <- sin(Xtest) + rnorm(n_test, mean = 0, sd = sigma) + 2.5
data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))
colnames(data_train) <- c('x', 'y', 'ind')
colnames(data_test) <- c('x', 'y', 'ind')

# parameters
ntree = 500
taus <- c(0.1, 0.5, 0.9)
nodesize.crf <- 10
nodesize.grf <- 10

one_run = function(ntree, tau) {
  # build generalizedForest model
  # get quantiles
  Yc <- crf.km(y ~ x, ntree = ntree, nodesize = nodesize.crf, data_train = data_train, data_test = data_test,
               yname = 'y', iname = 'ind', tau = tau, method = "grf", calibrate_taus = taus)
  
  # RSF
  v.rsf <- rfsrc(Surv(y, ind) ~ ., data = data_train, ntree = ntree, nodesize = nodesize.grf)
  surv.rsf <- predict(v.rsf, newdata = data_test)
  Yrsf <- find_quantile(surv = surv.rsf, max_value = max(data_train$y), tau = tau)
  
  # plot.survival(surv.rsf)
  
  # # generalized random forest (Stefan's)
  # # latent T ~ x
  grf_qf_latent <- quantile_forest(data_train[,1,drop=FALSE], Ttrain, quantiles = tau, num.trees = ntree, min.node.size = nodesize.grf)
  Ygrf_latent <- predict(grf_qf_latent, data_test[,1,drop=FALSE], quantiles = tau)
  # # censored Y ~ x
  grf_qf <- quantile_forest(data_train[,1,drop=FALSE], Ytrain, quantiles = tau, num.trees = ntree, min.node.size = nodesize.grf)
  Ygrf <- predict(grf_qf, data_test[,1,drop=FALSE], quantiles = tau)
  
  # Meinshasen
  #qrf_latent <- quantregForest(x=data_train[,1,drop=FALSE], y=Ttrain, nodesize=nodesize.crf, ntree=ntree)
  #Yqrf_latent <- predict(qrf_latent, data_test[,1,drop=FALSE], what = tau)
  
  # comparison
  plot(Xtest, Ytest, cex = 0.04, xlab = 'x', ylab = 'y', xlim=c(0, 6), ylim=c(0,4))
  quantiles <- sin(Xtest) + qnorm(tau, 0, sigma) + 2.5
  lines(Xtest, quantiles, col = 'black', cex = 2)
  lines(Xtest, Ygrf, col = 'black', lty = 2, cex = 1)
  lines(Xtest, Ygrf_latent, col = 'green', type = 'b', pch = 18, lty = 1, cex = .5)
  lines(Xtest, Yrsf, col = 'blue', lty = 5, cex = 1)
  lines(Xtest, Yc$predicted, col='red', lty = 5, cex = 1)
  
  # Add a legend
  legend(0.1, 1.8, legend=c("true quantile", "grf", "grf-oracle", "rsf", "crf"),
         lty=c(1, 2, 1, 5, 5), cex=c(1,1,1,1,1), pch = c(-1,-1,18, -1, -1), col = c('black', 'black', 'green', 'blue', 'red'))
  title(main = paste("tau =", tau))
}

for (tau in taus){
  one_run(ntree, tau)
}