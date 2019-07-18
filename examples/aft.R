setwd("../R")

source("metrics.R")
source("help_functions.R")
source('crf_km.R')
source("crf.R")

list.of.packages <- c("ggplot2", "quantregForest", "randomForestSRC", "survival")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library(ggplot2)
library(quantregForest)
library(randomForestSRC)
library(survival)

# help functions
get.event.info <- function(obj, subset = NULL) {
  ## survival case
  if (grepl("surv", obj$family)) {
    if (!is.null(obj$yvar)) {
      if (is.null(subset)) {
        subset <- (1:nrow(cbind(obj$yvar)))
      }
      r.dim <- 2
      time <- obj$yvar[subset, 1]
      cens <- obj$yvar[subset, 2]
      ## censoring must be coded coherently
      if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
        stop("for survival families censoring variable must be coded as a non-negative integer")
      }
      ## Extract the unique event types.
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- sort(unique(event))
    }
    ##everything else
    else {
      r.dim <- 0
      event <- event.type <- cens <- cens <- time <- NULL
    }
    ## Set grid of time points.
    time.interest <- obj$time.interest
  }
  else {
    ## NULL for other families
    if ((obj$family == "regr+") | (obj$family == "class+")) {
      r.dim <- dim(obj$yvar)[2]
    }
    else {
      r.dim <- 1
    }
    event <- event.type <- cens <- time.interest <- cens <- time <- NULL
  }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest, time = time, r.dim = r.dim))
}

find_quantile <- function(surv, max_value, tau){
  b_size <- dim(surv$survival)[1]
  event.info <- get.event.info(surv)
  bin <- event.info$time.interest
  quantiles <- rep(0, b_size)
  for (i in 1:b_size){
    check <- surv$survival[i,] >= 1-tau
    j <- max(sum(check), 1)
    quantile <- bin[j]
    quantiles[i] <- quantile
  }

  return(quantiles)
}

# ------------------------------------------------------- #

attach(mtcars)
op <- par(mar=c(4,4,1,1)+0.1, oma = c(0,0,0,0) + 0.1, pty="s")

# Load in the data
n <- 300
n_test <- 300

# training data
Xtrain <- sort(runif(n = n, min = 0, max = 2))
sigma <- 0.3
Ttrain <- exp(Xtrain + rnorm(n, mean = 0, sd = sigma))
ctrain <- rexp(n = n, rate = 0.08)
Ytrain <- pmin(Ttrain, ctrain)
censorInd <- 1*(Ttrain <= ctrain)
data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
# plot training data
plot(Xtrain, Ytrain, cex = 0.2)
points(Xtrain[!censorInd], Ytrain[!censorInd], type = 'p', col = 'red', cex = 0.3)
points(Xtrain[!censorInd], Ttrain[!censorInd], type = 'p', col = 'green', cex = 0.3)

# test data
Xtest <- sort(runif(n = n_test, min = 0, max = 2))
Ytest <- exp(Xtest + rnorm(n_test, mean = 0, sd = sigma))
data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))

colnames(data_train) <- c('x', 'y', 'ind')
colnames(data_test) <- c('x', 'y', 'ind')

# build an AFT model
# aft <- survreg(Surv(y, ind) ~ ., data_train, dist='weibull', scale=1)
# Yaft <- predict(aft, data_test)

# build generalizedForest model
# get quantiles
tau <- 0.8
nodesize <- 5
Yc <- crf.km(y ~ x, ntree = 1000, nodesize = 6*nodesize, data_train = data_train, data_test = data_test,
          yname = 'y', iname = 'ind',
          tau = tau)

# RSF
v.rsf <- rfsrc(Surv(y, ind) ~ ., data = data_train, ntree = 1000)
surv.rsf <- predict(v.rsf, newdata = data_test)
Yrsf <- find_quantile(surv = surv.rsf, max_value = max(data_train$y), tau = tau)

plot.survival(surv.rsf)

# without debiasing
#Yc_nodebias <- crf(y ~ x, ntree = 2000, nodesize = nodesize, data_train = data_train, data_test = data_test,
#          yname = 'y', iname = 'ind',
#          tau = tau, fixed_censoring = FALSE, fixed_c = NULL, debias = FALSE)

# # generalized random forest (Stefan's)
# # latent T ~ x
grf_qf_latent <- quantile_forest(data_train[,1,drop=FALSE], Ttrain, quantiles = tau, num.trees = 1000, min.node.size = nodesize)
Ygrf_latent <- predict(grf_qf_latent, Xtest, quantiles = tau)
# # censored Y ~ x
grf_qf <- quantile_forest(data_train[,1,drop=FALSE], Ytrain, quantiles = tau, num.trees = 1000, min.node.size = nodesize)
Ygrf <- predict(grf_qf, Xtest, quantiles = tau)

# # survival forest
# surv_rf <- rfsrc(Surv(y, ind) ~ x, data = data_train, ntree = 1000, nodesize = 100)
# Ysurv <- predict(surv_rf, newdata = data_test)$predicted

# Meinshasen
#qrf_latent <- quantregForest(x=data_train[,1,drop=FALSE], y=Ttrain, nodesize=3*nodesize, ntree=1000)
#Yqrf_latent <- predict(qrf_latent, data_test[,1,drop=FALSE], what = tau)

# comparison
plot(Xtest, Ytest, cex = 0.04, xlab = 'x', ylab = 'y')
quantiles <- exp(Xtest + qnorm(tau, 0, sigma))
lines(Xtest, quantiles, col = 'black', cex = 2)
lines(Xtest, Yc$predicted, col='red', lty = 5, cex = 1)
lines(Xtest, Ygrf, col = 'black', lty = 2, cex = 1)
lines(Xtest, Ygrf_latent, col = 'black', type = 'b', pch = 18, lty = 1, cex = .5)
lines(Xtest, Yrsf, col = 'blue', lty = 5, cex = 1)

# Add a legend
legend(0.1, 10, legend=c("true quantile", "cRF", "gRF", "gRF-oracle", "rsf"),
       lty=c(1, 5, 2, 1, 5), cex=0.8, pch = c(-1,-1,-1, 18, -1), col = c('black', 'red', 'black', 'black', 'blue'))
title(main = paste("tau =", tau))
