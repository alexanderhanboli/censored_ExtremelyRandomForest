source("metrics.R")
source("help_functions.R")
source("crf_km.R")

library(quantregForest)
library(ggplot2)
library(grf)
library(randomForestSRC)
library(timereg)
library(mlbench)

## data, formula specifications
dataset <- 'BostonHousing'

if (dataset == 'BostonHousing') {
  data("BostonHousing", package = "mlbench")
  data.all <- na.omit(BostonHousing)
  mean.data <- mean(data.all$medv)
  
  yid <- which(colnames(data.all) == 'medv')
  colnames(data.all)[yid] <- 'y'
  
  censoring <- rexp(n = nrow(data.all), rate = 0.25/mean.data)
  data.all$truth <- data.all$y
  data.all$status <- 1*(data.all$y <= censoring)
  data.all$y <- pmin(data.all$y, censoring)
  xnam <- colnames(data.all)[-which(colnames(data.all) %in% c('status', 'y', 'truth'))]
  ynam <- 'y'
  inam <- 'status'
  for (col in colnames(data.all)) {
    data.all[[col]] <- as.numeric(data.all[[col]])
  }
  print(mean(data.all[[inam]]))
}



one_run <- function(b, data, xnam, ynam, inam, tau, nodesize, ntree) {
  cat("\n bootstrap step:", b, "\n")
  
  train <- sample(1:nrow(data), floor(nrow(data)*0.8), replace = FALSE)
  data_train <- data[train, ]
  data_test <- data[-train, ]
  #data_test <- data[-train, ][data[-train, inam] == 1, ]
  print(paste0("num of training data: ", nrow(data_train)))
  print(paste0("num of test data: ", nrow(data_test)))
  
  # build generalizedForest model
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  print("doing crf...")
  Yc.qrf <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                   yname = ynam, iname = inam, xname = xnam, tau = tau, method = "ranger", splitrule = "variance")$predicted
  print("doing crf...")
  Yc.grf <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                   yname = ynam, iname = inam, xname = xnam, tau = tau, method = "grf", calibrate_taus = tau)$predicted
  
  
  # generalized random forest (Stefan's)
  print("doing grf...")
  grf_qf <- quantile_forest(data_train[,xnam], data_train[,ynam], quantiles = tau, 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf <- predict(grf_qf, data_test[,xnam], quantiles = tau)
  
  grf_qf <- quantile_forest(data_train[,xnam], data_train[,'truth'], quantiles = tau, 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf.oracle <- predict(grf_qf, data_test[,xnam], quantiles = tau)
  
  # quantile random forest (Meinshasen)
  print("doing qrf...")
  qrf <- quantregForest(x=data_train[,xnam], y=data_train[,ynam], nodesize=nodesize, ntree=ntree)
  Yqrf <- predict(qrf, data_test[,xnam], what = tau)
  
  qrf <- quantregForest(x=data_train[,xnam], y=data_train[,'truth'], nodesize=nodesize, ntree=ntree)
  Yqrf.oracle <- predict(qrf, data_test[,xnam], what = tau)
  
  # survival forest
  print("doing rsf...")
  v.rsf <- rfsrc(as.formula(paste("Surv(y, status) ~ ", paste(xnam, collapse= "+"))), data = data_train, ntree = ntree, nodesize = nodesize)
  surv.rsf <- predict(v.rsf, newdata = data_test)
  Ysurv <- find_quantile(surv = surv.rsf, max_value = max(data_train$y), tau = tau)
  
  # results
  return(
    list(
      'crf.qrf'=metrics(data_test[['truth']], Yc.qrf, data_test[['truth']], tau),
      'crf.grf'=metrics(data_test[['truth']], Yc.grf, data_test[['truth']], tau),
      'qrf'=metrics(data_test[['truth']], Yqrf, data_test[['truth']], tau),
      'qrf.oracle'=metrics(data_test[['truth']], Yqrf.oracle, data_test[['truth']], tau),
      'grf'=metrics(data_test[['truth']], Ygrf, data_test[['truth']], tau),
      'grf.oracle'=metrics(data_test[['truth']], Ygrf.oracle, data_test[['truth']], tau),
      'rsf'=metrics(data_test[['truth']], Ysurv, data_test[['truth']], tau),
      
      'crf.qrf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Yc.qrf),
      'crf.grf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Yc.grf),
      'qrf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Yqrf),
      'qrf.c.oracle'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Yqrf.oracle),
      'rsf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Ysurv),
      'grf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Ygrf),
      'grf.c.oracle'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Ygrf.oracle)
    )
  )
}

B <- 80
nodesize <- 0
ntree <- 1000
tau <- 0.1
quantile_result <- list('nodesize' = rep(NA,B), 'crf.grf'=rep(NA,B), 'crf.qrf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'rsf'=rep(NA,B), 'grf.oracle'=rep(NA,B), 'qrf.oracle'=rep(NA,B))
cindex_result <- list('nodesize' = rep(NA,B), 'crf.grf'=rep(NA,B), 'crf.qrf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'rsf'=rep(NA,B), 'grf.oracle'=rep(NA,B), 'qrf.oracle'=rep(NA,B))

print(mean(data.all[[inam]]))

for (t in 1:B) {
  print(t)
  
  if (t%%10 == 1) {
    nodesize <- nodesize + 5
  }
  print(paste0("nodesize is ", nodesize))
  
  tmp <- one_run(t, data.all, xnam, ynam, inam, tau, nodesize, ntree)
  
  quantile_result$nodesize[t] <- nodesize
  quantile_result$crf.grf[t] <- tmp$crf.grf['quantile_loss']
  quantile_result$crf.qrf[t] <- tmp$crf.qrf['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$qrf.oracle[t] <- tmp$qrf.oracle['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  quantile_result$grf.oracle[t] <- tmp$grf.oracle['quantile_loss']
  quantile_result$rsf[t] <- tmp$rsf['quantile_loss']
  
  cindex_result$nodesize[t] <- nodesize
  cindex_result$crf.grf[t] <- tmp$crf.grf.c
  cindex_result$crf.qrf[t] <- tmp$crf.qrf.c
  cindex_result$qrf[t] <- tmp$qrf.c
  cindex_result$qrf.oracle[t] <- tmp$qrf.c.oracle
  cindex_result$rsf[t] <- tmp$rsf.c
  cindex_result$grf[t] <- tmp$grf.c
  cindex_result$grf.oracle[t] <- tmp$grf.c.oracle
} 


# plot
require(reshape2)
dd <- melt(as.data.frame(quantile_result), id.vars = 'nodesize')
dd.agg <- aggregate(value ~ nodesize + variable, dd, function(x) c(mean = mean(x), sd = sd(x)/sqrt(10)))
dd.agg$mean <- dd.agg[-1][[2]][,1]
dd.agg$sd <- dd.agg[-1][[2]][,2]
dd.agg$value <- NULL
ggplot(data = dd.agg, aes(x=nodesize, y=mean, colour=variable)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.05)) +
  geom_line() +
  geom_point() +
  labs(x = "Nodesize", y = paste("Quantile loss, tau =", tau), fill = "Nodesize")
ggsave(paste0(dataset, "_quantile_loss_nodesize_tau_0", 10*tau, ".pdf"), width = 5, height = 5, path = "../examples/figs")