#setwd("./R")

source("metrics.R")
source("help_functions.R")
source("crf_km.R")

library(quantregForest)
library(ggplot2)
library(grf)
library(randomForestSRC)


# Load in the data
n <- 1000
n_test <- 200
p <- 20

one_run <- function(n, n_test, p, tau, nodesize, ntree, rate) {
  
  # training data
  Xtrain <- matrix(runif(n = n*p, min = 0, max = 2), nrow = n, ncol = p)
  sigma <- 0.3
  Ttrain <- exp(Xtrain[,1] + rnorm(n, mean = 0, sd = sigma))
  ctrain <- rexp(n = n, rate = rate)
  Ytrain <- pmin(Ttrain, ctrain)
  censorInd <- (Ttrain <= ctrain)
  print(paste("censoring level is", 1-mean(censorInd)))
  data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
  
  # test data
  Xtest <- matrix(runif(n = n_test*p, min = 0, max = 2), nrow = n_test, ncol = p)
  quantile_test <- exp(Xtest[,1] + qnorm(tau, 0, sigma))
  Ytest <- exp(Xtest[,1] + rnorm(n_test, mean = 0, sd = sigma))
  data_test <- cbind.data.frame(Xtest, Ytest, rep(TRUE, n_test))
  
  # assign column names
  xnam <- paste0('x', 1:p)
  colnames(data_train) <- c(xnam, 'y', 'status')
  colnames(data_test) <- c(xnam, 'y', 'status')
  
  # build crf models
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  Yc.qrf <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                   yname = 'y', iname = 'status', tau = tau, method = "grf", calibrate_taus = tau, 
                   reg.split = TRUE)$predicted
  
  Yc.grf <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                   yname = 'y', iname = 'status', tau = tau, method = "grf", calibrate_taus = tau)$predicted
  
  # generalized random forest (Stefan's)
  grf_qf_latent <- quantile_forest(data_train[,1:p,drop=FALSE], Ttrain, quantiles = tau, 
                                   num.trees = ntree, min.node.size = nodesize)
  Ygrf_latent <- predict(grf_qf_latent, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  grf_qf <- quantile_forest(data_train[,1:p,drop=FALSE], Ytrain, quantiles = tau, 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf <- predict(grf_qf, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  # quantile random forest (Meinshasen)
  qrf_latent <- quantile_forest(data_train[,1:p,drop=FALSE], Ttrain, quantiles = tau, 
                                num.trees = ntree, min.node.size = nodesize, regression.splitting = TRUE)
  Yqrf_latent <- predict(qrf_latent, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  qrf <- quantile_forest(data_train[,1:p,drop=FALSE], Ytrain, quantiles = tau, 
                         num.trees = ntree, min.node.size = nodesize, regression.splitting = TRUE)
  Yqrf <- predict(qrf, data_test[,1:p,drop=FALSE], quantiles = tau)
  
  # RSF
  v.rsf <- rfsrc(Surv(y, status) ~ ., data = data_train, ntree = ntree, nodesize = nodesize)
  surv.rsf <- predict(v.rsf, newdata = data_test)
  Ysurv <- find_quantile(surv = surv.rsf, max_value = max(data_train$y), tau = tau)
  
  # results
  return(
    list(
      'censoring_level'=1-mean(censorInd),
      'crf_quantile'=metrics(data_test$y, Yc.qrf, quantile_test, tau),
      'crf_generalized'=metrics(data_test$y, Yc.grf, quantile_test, tau),
      'qrf'=metrics(data_test$y, Yqrf, quantile_test, tau),
      'qrf_latent'=metrics(data_test$y, Yqrf_latent, quantile_test, tau),
      'grf'=metrics(data_test$y, Ygrf, quantile_test, tau),
      'surv'=metrics(data_test$y, Ysurv, quantile_test, tau),
      'grf_latent'=metrics(data_test$y, Ygrf_latent, quantile_test, tau),
      
      'crf_quantile_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Yc.qrf),
      'crf_generalized_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Yc.grf),
      'qrf_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Yqrf),
      'qrf_latent_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Yqrf_latent),
      'grf_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Ygrf),
      'surv_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Ysurv),
      'grf_latent_c'=randomForestSRC:::get.cindex(data_test$y, data_test$status, Ygrf_latent)
    )
  )
}

B = 80
tau <- 0.3
nodesize <- 30
rate <- 0.08
ntree <- 1000
mse_result <- list('ntree'=rep(NA,B), 'rate' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
mad_result <- list('ntree'=rep(NA,B), 'rate' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
quantile_result <- list('ntree'=rep(NA,B), 'rate' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
cindex_result <- list('ntree'=rep(NA,B), 'rate' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
censoring_level <- rep(NA, B)

for (t in 1:B) {
  print(t)
  
  if (t%%10 == 1) {
    rate <- rate + 0.02
  }
  tmp <- one_run(n, n_test, p, tau, nodesize, ntree, rate)
  
  censoring_level[t] <- tmp$censoring_level
  
  print("saving results...")
  #mse_result$ntree[t] <- ntree
  mse_result$rate[t] <- rate
  mse_result$crf_quantile[t] <- tmp$crf_quantile['mse']
  mse_result$crf_generalized[t] <- tmp$crf_generalized['mse']
  mse_result$qrf[t] <- tmp$qrf['mse']
  mse_result$grf[t] <- tmp$grf['mse']
  mse_result$rsf[t] <- tmp$surv['mse']
  mse_result$grf_oracle[t] <- tmp$grf_latent['mse']
  mse_result$qrf_oracle[t] <- tmp$qrf_latent['mse']
  
  #mad_result$ntree[t] <- ntree
  mad_result$rate[t] <- rate
  mad_result$crf_quantile[t] <- tmp$crf_quantile['mad']
  mad_result$crf_generalized[t] <- tmp$crf_generalized['mad']
  mad_result$qrf[t] <- tmp$qrf['mad']
  mad_result$grf[t] <- tmp$grf['mad']
  mad_result$rsf[t] <- tmp$surv['mad']
  mad_result$grf_oracle[t] <- tmp$grf_latent['mad']
  mad_result$qrf_oracle[t] <- tmp$qrf_latent['mad']
  
  #quantile_result$ntree[t] <- ntree
  quantile_result$rate[t] <- rate
  quantile_result$crf_quantile[t] <- tmp$crf_quantile['quantile_loss']
  quantile_result$crf_generalized[t] <- tmp$crf_generalized['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  quantile_result$rsf[t] <- tmp$surv['quantile_loss']
  quantile_result$grf_oracle[t] <- tmp$grf_latent['quantile_loss']
  quantile_result$qrf_oracle[t] <- tmp$qrf_latent['quantile_loss']
  
  #cindex_result$ntree[t] <- ntree
  cindex_result$rate[t] <- rate
  cindex_result$crf_quantile[t] <- tmp$crf_quantile_c
  cindex_result$crf_generalized[t] <- tmp$crf_generalized_c
  cindex_result$rsf[t] <- tmp$surv_c
  cindex_result$qrf[t] <- tmp$qrf_c
  cindex_result$grf[t] <- tmp$grf_c
  cindex_result$grf_oracle[t] <- tmp$grf_latent_c
  cindex_result$qrf_oracle[t] <- tmp$qrf_latent_c
}

# plot
require(reshape2)
dd <- melt(as.data.frame(quantile_result), id.vars = 'rate')
dd.agg <- aggregate(value ~ rate + variable, dd, function(x) c(mean = mean(x), sd = sd(x)/sqrt(10)))
dd.agg$mean <- dd.agg[-1][[2]][,1]
dd.agg$sd <- dd.agg[-1][[2]][,2]
dd.agg$value <- NULL
ggplot(data = dd.agg, aes(x=rate, y=mean, colour=variable)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.01, position=position_dodge(0.0)) +
  geom_line() +
  geom_point() +
  labs(x = expression(lambda), y = paste("Quantile loss, tau =", tau), fill = "rate")
ggsave(paste0("run3_aft_multiD_quantile_loss_rate_tau_0", 10*tau, ".pdf"), width = 5, height = 5, path = "../examples/figs")