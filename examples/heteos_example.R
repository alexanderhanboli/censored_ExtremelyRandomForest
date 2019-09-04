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
op <- par(mfrow=c(2,2), mar=c(4,4,1,1)+0.1, oma = c(0,0,0,0) + 0.1, pty="s")

# Create the data
n <- 2000
n_test <- 500
p <- 40

# training data
Xtrain <- matrix(runif(n = n*p, min = -1, max = 1), nrow = n, ncol = p)
Ttrain <- rnorm(n = n, mean = 0, sd = 1 + 1*(Xtrain[,1]>0))
#ctrain <- rep(1000000000, n)
ctrain <- rexp(n = n, rate = 0.05) - 2
Ytrain <- pmin(Ttrain, ctrain)
censorInd <- 1*(Ttrain <= ctrain)
data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
# plot training data
plot(Xtrain[,1], Ytrain, cex = 0.2)
points(Xtrain[!censorInd,1], Ytrain[!censorInd], type = 'p', col = 'red', cex = 0.3)
points(Xtrain[!censorInd,1], Ttrain[!censorInd], type = 'p', col = 'green', cex = 0.3)

# test data
Xtest <- matrix(runif(n = n_test*p, min = -1, max = 1), nrow = n_test, ncol = p)
Ytest <- rnorm(n = n_test, mean = 0, sd = 1 + 1*(Xtest[,1]>0))
data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))

# column names
xnam <- paste0('x', 1:p)
colnames(data_train) <- c(xnam, 'y', 'status')
colnames(data_test) <- c(xnam, 'y', 'status')

# parameters
ntree = 1500
taus <- c(0.1, 0.5, 0.9)
nodesize.crf <- 60
nodesize.qrf <- 60
nodesize.grf <- 60

one_run = function(ntree, tau, nodesize) {
  # build censored Extreme Forest model
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  Yc <- crf.km(fmla, ntree = ntree, nodesize = nodesize.crf, data_train = data_train, data_test = data_test, 
               yname = 'y', iname = 'status', tau = tau, method = "grf", splitrule = "extratrees")
  Yc <- Yc$predicted
  
  # generalized random forest (Stefan's)
  # grf_qf_latent <- quantile_forest(data_train[,1:p,drop=FALSE], Ttrain, quantiles = tau, 
  # num.trees = ntree, min.node.size = nodesize)
  # Ygrf_latent <- predict(grf_qf_latent, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  grf_qf <- quantile_forest(data_train[,1:p,drop=FALSE], Ytrain, quantiles = tau, 
                            num.trees = ntree, min.node.size = nodesize.grf)
  Ygrf <- predict(grf_qf, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  # quantile random forest (Meinshasen)
  # qrf_latent <- quantregForest(x=Xtrain, y=Ttrain, nodesize=nodesize, ntree=ntree)
  # Yqrf_latent <- predict(qrf_latent, Xtest, what = tau)
  
  qrf <- quantile_forest(data_train[,1:p,drop=FALSE], Ytrain, quantiles = tau, 
                         num.trees = ntree, min.node.size = nodesize.qrf, regression.splitting = TRUE)
  Yqrf <- predict(qrf, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  # comparison
  plot(Xtest[,1], Ytest, cex = 0.04, xlab = 'x', ylab = 'y', ylim = c(-3, 3))
  Xsort <- sort(Xtest[,1], index.return=TRUE)
  X1 <- Xsort$x
  Xindex <- Xsort$ix
  quantiles <- qnorm(tau, 0, 1 + 1*(X1>0))
  lines(X1, quantiles, col = 'black', cex = 2)
  lines(X1, Ygrf[Xindex], col = 'green', lty = 5, cex = 1)
  # lines(Xtest[,1], Ygrf_latent, col = 'black', type = 'b', pch = 18, lty = 1, cex = .5)
  lines(X1, Yqrf[Xindex], col = 'blue', lty = 5, cex = 1)
  lines(X1, Yc[Xindex], col='red', lty = 5, cex = 1)
  
  # Add a legend
  legend(-1, 3, legend=c("true quantile", "grf", "qrf", "crf"),
         lty=c(1, 5, 5, 5), cex=0.8, pch = c(-1,-1, -1, -1), col = c('black', 'green', 'blue', 'red'))
  title(main = paste("tau =", tau))
}

for (tau in taus){
  one_run(ntree, tau, nodesize)
}