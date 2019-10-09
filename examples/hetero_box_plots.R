setwd("./R")

source("metrics.R")
source("help_functions.R")
source("crf_km.R")

library(quantregForest)
library(ggplot2)
library(grf)
library(randomForestSRC)


# Load in the data
n <- 2000
n_test <- 200
p <- 40
ntree <- 2000

one_run <- function(n, n_test, p, tau, nodesize, ntree) {
  
  # training data
  # training data
  Xtrain <- matrix(runif(n = n*p, min = -1, max = 1), nrow = n, ncol = p)
  Ttrain <- 10 + rnorm(n = n, mean = 0, sd = 1 + 1*(Xtrain[,1]>0))
  ctrain <- rexp(n = n, rate = 0.05) + 8
  Ytrain <- pmin(Ttrain, ctrain)
  censorInd <- 1*(Ttrain <= ctrain)
  print(paste("censoring level is", 1-mean(censorInd)))
  data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
  
  # test data
  # test data
  Xtest <- matrix(runif(n = n_test*p, min = -1, max = 1), nrow = n_test, ncol = p)
  quantile_test <- 10 + qnorm(tau, 0, 1 + 1*(Xtest[,1]>0))
  Ytest <- 10 + rnorm(n = n_test, mean = 0, sd = 1 + 1*(Xtest[,1]>0))
  data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))
  
  # column names
  xnam <- paste0('x', 1:p)
  colnames(data_train) <- c(xnam, 'y', 'status')
  colnames(data_test) <- c(xnam, 'y', 'status')
  
  # build censored Extreme Forest model
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  Yc.qrf <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                   yname = 'y', iname = 'status', tau = tau, method = "ranger", splitrule = "extratrees")$predicted
  
  Yc.grf <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                   yname = 'y', iname = 'status', tau = tau, method = "grf", calibrate_taus = c(0.1, 0.5, 0.9))$predicted
  
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
      'crf_quantile'=metrics(data_test$y, Yc.qrf, quantile_test, tau),
      'crf_generalized'=metrics(data_test$y, Yc.grf, quantile_test, tau),
      'qrf'=metrics(data_test$y, Yqrf, quantile_test, tau),
      'qrf_latent'=metrics(data_test$y, Yqrf_latent, quantile_test, tau),
      'grf'=metrics(data_test$y, Ygrf, quantile_test, tau),
      'surv'=metrics(data_test$y, Ysurv, quantile_test, tau),
      'grf_latent'=metrics(data_test$y, Ygrf_latent, quantile_test, tau),
      
      'crf_quantile_c'=1-randomForestSRC:::get.cindex(data_test$y, data_test$status, Yc.qrf),
      'crf_generalized_c'=1-randomForestSRC:::get.cindex(data_test$y, data_test$status, Yc.grf),
      'qrf_c'=1-randomForestSRC:::get.cindex(data_test$y, data_test$status, Yqrf),
      'qrf_latent_c'=1-randomForestSRC:::get.cindex(data_test$y, data_test$status, Yqrf_latent),
      'grf_c'=1-randomForestSRC:::get.cindex(data_test$y, data_test$status, Ygrf),
      'surv_c'=1-randomForestSRC:::get.cindex(data_test$y, data_test$status, Ysurv),
      'grf_latent_c'=1-randomForestSRC:::get.cindex(data_test$y, data_test$status, Ygrf_latent)
    )
  )
}

B = 20
tau <- 0.1
nodesize <- 80
mse_result <- list('nodesize' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
mad_result <- list('nodesize' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
quantile_result <- list('nodesize' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))
cindex_result <- list('nodesize' = rep(NA,B), 'crf_quantile'=rep(NA,B), 'crf_generalized'=rep(NA,B), 'qrf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf'=rep(NA,B), 'grf_oracle'=rep(NA,B), 'rsf'=rep(NA,B))

for (t in 1:B) {
  print(t)
  
  tmp <- one_run(n, n_test, p, tau, nodesize, ntree)
  
  print("saving results...")
  
  quantile_result$crf_quantile[t] <- tmp$crf_quantile['quantile_loss']
  quantile_result$crf_generalized[t] <- tmp$crf_generalized['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  quantile_result$rsf[t] <- tmp$surv['quantile_loss']
  quantile_result$grf_oracle[t] <- tmp$grf_latent['quantile_loss']
  quantile_result$qrf_oracle[t] <- tmp$qrf_latent['quantile_loss']
  
  cindex_result$crf_quantile[t] <- tmp$crf_quantile_c
  cindex_result$crf_generalized[t] <- tmp$crf_generalized_c
  cindex_result$rsf[t] <- tmp$surv_c
  cindex_result$qrf[t] <- tmp$qrf_c
  cindex_result$grf[t] <- tmp$grf_c
  cindex_result$grf_oracle[t] <- tmp$grf_latent_c
  cindex_result$qrf_oracle[t] <- tmp$qrf_latent_c
}

# plot
# boxplot
require(reshape2)
# cindex results
dd <- as.data.frame(cindex_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("C-index, tau =", tau), fill = "Model")
ggsave(paste0("hetero_cindex_boxplots_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "../examples/figs")

# quantile loss results
dd <- as.data.frame(quantile_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("Quantile loss, tau =", tau), fill = "Model")
ggsave(paste0("hetero_quantile_loss_boxplots_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "../examples/figs")
