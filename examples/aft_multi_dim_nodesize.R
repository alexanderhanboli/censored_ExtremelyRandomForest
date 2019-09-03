setwd("./R")

source("metrics.R")
source("help_functions.R")
source("crf_km.R")

library(quantregForest)
library(ggplot2)
library(grf)
library(randomForestSRC)


# Load in the data
n <- 500
n_test <- 300
p <- 5
ntree <- 1000

one_run <- function(n, n_test, p, tau, nodesize, ntree) {
  
  # training data
  Xtrain <- matrix(runif(n = n*p, min = 0, max = 2), nrow = n, ncol = p)
  beta <- c(0.1,0.2,0.3,0.4,0.5)
  sigma <- 0.3
  Ttrain <- exp(Xtrain%*%beta + rnorm(n, mean = 0, sd = sigma))
  ctrain <- rexp(n = n, rate = 0.05)
  Ytrain <- pmin(Ttrain, ctrain)
  censorInd <- (Ttrain <= ctrain)
  print(paste("censoring level is", 1-mean(censorInd)))
  data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
  
  # test data
  Xtest <- matrix(runif(n = n_test*p, min = 0, max = 2), nrow = n_test, ncol = p)
  quantile_test <- exp(Xtest%*%beta + qnorm(tau, 0, sigma))
  Ytest <- exp(Xtest%*%beta + rnorm(n_test, mean = 0, sd = sigma))
  data_test <- cbind.data.frame(Xtest, Ytest, rep(TRUE, n_test))
  # column names
  xnam <- paste0('x', 1:p)
  colnames(data_train) <- c(xnam, 'y', 'status')
  colnames(data_test) <- c(xnam, 'y', 'status')
  
  # build censored Extreme Forest model
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  Yc <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
               yname = 'y', iname = 'status', tau = tau, method = "ranger", splitrule = "extratrees")$predicted
  
  # generalized random forest (Stefan's)
  grf_qf_latent <- quantile_forest(data_train[,1:p,drop=FALSE], Ttrain, quantiles = tau, 
                                   num.trees = ntree, min.node.size = nodesize)
  Ygrf_latent <- predict(grf_qf_latent, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  grf_qf <- quantile_forest(data_train[,1:p,drop=FALSE], Ytrain, quantiles = tau, 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf <- predict(grf_qf, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  # quantile random forest (Meinshasen)
  qrf_latent <- quantregForest(x=Xtrain, y=Ttrain, nodesize=nodesize, ntree=ntree)
  Yqrf_latent <- predict(qrf_latent, Xtest, what = tau)
  
  qrf <- quantregForest(x=Xtrain, y=Ytrain, nodesize=nodesize, ntree=ntree)
  Yqrf <- predict(qrf, Xtest, what = tau)
  
  # RSF
  v.rsf <- rfsrc(Surv(y, status) ~ ., data = data_train, ntree = ntree, nodesize = nodesize)
  surv.rsf <- predict(v.rsf, newdata = data_test)
  Ysurv <- find_quantile(surv = surv.rsf, max_value = max(data_train$y), tau = tau)
  
  # results
  return(
    list(
      'crf'=metrics(data_test$y, Yc, quantile_test, tau),
      'qrf'=metrics(data_test$y, Yqrf, quantile_test, tau),
      'qrf_latent'=metrics(data_test$y, Yqrf_latent, quantile_test, tau),
      'grf'=metrics(data_test$y, Ygrf, quantile_test, tau),
      'surv'=metrics(data_test$y, Ysurv, quantile_test, tau),
      'grf_latent'=metrics(data_test$y, Ygrf_latent, quantile_test, tau),
      
      'crf_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Yc),
      'qrf_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Yqrf),
      'qrf_latent_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Yqrf_latent),
      'grf_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Ygrf),
      'surv_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Ysurv),
      'grf_latent_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Ygrf_latent)
    )
  )
}

B = 120
tau <- 0.2
nodesize <- 0
mse_result <- list('nodesize' = rep(NA,B), 'crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
mad_result <- list('nodesize' = rep(NA,B), 'crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
quantile_result <- list('nodesize' = rep(NA,B), 'crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
cindex_result <- list('nodesize' = rep(NA,B), 'crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'surv'=rep(NA,B))

for (t in 1:B) {
  print(t)
  
  if (t%%10 == 1) {
    nodesize <- nodesize + 5
  }
  print(paste0("nodesize is ", nodesize))
  tmp <- one_run(n, n_test, p, tau, nodesize, ntree)
  
  print("saving results...")
  mse_result$nodesize[t] <- nodesize
  mse_result$crf[t] <- tmp$crf['mse']
  mse_result$qrf[t] <- tmp$qrf['mse']
  mse_result$qrf_oracle[t] <- tmp$qrf_latent['mse']
  mse_result$grf[t] <- tmp$grf['mse']
  mse_result$rsf[t] <- tmp$surv['mse']
  mse_result$grf_oracle[t] <- tmp$grf_latent['mse']
  
  mad_result$nodesize[t] <- nodesize
  mad_result$crf[t] <- tmp$crf['mad']
  mad_result$qrf[t] <- tmp$qrf['mad']
  mad_result$qrf_oracle[t] <- tmp$qrf_latent['mad']
  mad_result$grf[t] <- tmp$grf['mad']
  mad_result$rsf[t] <- tmp$surv['mad']
  mad_result$grf_oracle[t] <- tmp$grf_latent['mad']
  
  quantile_result$nodesize[t] <- nodesize
  quantile_result$crf[t] <- tmp$crf['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$qrf_oracle[t] <- tmp$qrf_latent['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  quantile_result$rsf[t] <- tmp$surv['quantile_loss']
  quantile_result$grf_oracle[t] <- tmp$grf_latent['quantile_loss']
  
  cindex_result$nodesize[t] <- nodesize
  cindex_result$crf[t] <- tmp$crf_c
  cindex_result$rsf[t] <- tmp$surv_c
  cindex_result$qrf[t] <- tmp$qrf_c
  cindex_result$qrf_oracle[t] <- tmp$qrf_latent_c
  cindex_result$grf[t] <- tmp$grf_c
  cindex_result$grf_oracle[t] <- tmp$grf_latent_c
}

# plot
require(reshape2)
dd <- melt(as.data.frame(quantile_result), id.vars = 'nodesize')
dd.agg <- aggregate(value ~ nodesize + variable, dd, function(x) c(mean = mean(x), sd = sd(x)))
dd.agg$mean <- dd.agg[-1][[2]][,1]
dd.agg$sd <- dd.agg[-1][[2]][,2]
dd.agg$value <- NULL
ggplot(data = dd.agg, aes(x=nodesize, y=mean, colour=variable)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.05)) +
  geom_line() +
  geom_point() +
  labs(x = "Nodesize", y = paste("Quantile loss, tau =", tau), fill = "Nodesize")
ggsave(paste0("aft_multiD_quantile_loss_nodesize_tau_0", 10*tau, ".pdf"), width = 5, height = 5, path = "../examples/figs")