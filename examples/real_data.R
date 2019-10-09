source("metrics.R")
source("help_functions.R")
source("crf_km.R")

library(quantregForest)
library(ggplot2)
library(grf)
library(randomForestSRC)
library(timereg)
library(mlbench)
library(survival)

## data, formula specifications
dataset <- 'pbc'

if (dataset == 'pbc') {
  data(pbc, package = "randomForestSRC")
  data.all <- na.omit(pbc) ##remove NA's
  xnam <- colnames(data.all)[-which(colnames(data.all) %in% c('status', 'days'))]
  yid <- which(colnames(data.all) == 'days')
  colnames(data.all)[yid] <- 'y'
  ynam <- 'y'
  data.all$y <- as.numeric(data.all$y)
  inam <- 'status'
  print(mean(data.all[[inam]]))
  #surv.f <- as.formula(Surv(y, status) ~ .)
} else if (dataset == 'sTRACE') {
  data(sTRACE, package = "timereg")
  data.all <- na.omit(sTRACE)
  data.all$status <- 1*(data.all$status == 9)
  xnam <- c("age", "sex", "diabetes", "chf", "vf")
  yid <- which(colnames(data.all) == 'time')
  colnames(data.all)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.all[[inam]]))
} else if (dataset == 'csl') {
  data(csl, package = "timereg")
  data.all <- na.omit(csl)
  data.all$status <- data.all$dc
  data.all$dc <- NULL
  xnam <- c("prot", "age", "sex", "diabetes", "chf", "vf")
  yid <- which(colnames(data.all) == 'time')
  colnames(data.all)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.all[[inam]]))
} else if (dataset == 'diabetes') {
  data(diabetes, package = "timereg")
  data.all <- na.omit(diabetes)
  data.all$id <- NULL
  xnam <- colnames(data.all)[-which(colnames(data.all) %in% c('status', 'time'))]
  yid <- which(colnames(data.all) == 'time')
  colnames(data.all)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.all[[inam]]))
} else if (dataset == 'bmt') {
  data(bmt, package = "timereg")
  data.all <- na.omit(bmt)
  data.all$status <- 1*(data.all$cause > 0)
  data.all$cause <- NULL
  xnam <- colnames(data.all)[-which(colnames(data.all) %in% c('status', 'time'))]
  yid <- which(colnames(data.all) == 'time')
  colnames(data.all)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.all[[inam]]))
} else if (dataset == 'follic') {
  data(follic, package = "randomForestSRC")
  data.all <- na.omit(follic)
  data.all$status <- 1*(data.all$status > 1)
  data.all$chn <- 1*(data.all$ch == 'Y')
  data.all$rtn <- 1*(data.all$rt == 'Y')
  data.all$ch <- NULL
  data.all$rt <- NULL
  xnam <- colnames(data.all)[-which(colnames(data.all) %in% c('status', 'time'))]
  yid <- which(colnames(data.all) == 'time')
  colnames(data.all)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.all[[inam]]))
} else if (dataset == 'lung') {
  data(lung, package = "survival")
  data.all <- na.omit(lung)
  data.all$inst <- NULL
  data.all$status <- 1*(data.all$status > 1)
  xnam <- colnames(data.all)[-which(colnames(data.all) %in% c('status', 'time'))]
  yid <- which(colnames(data.all) == 'time')
  colnames(data.all)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
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
  
  # quantile random forest (Meinshasen)
  print("doing qrf...")
  qrf <- quantregForest(x=data_train[,xnam], y=data_train[,ynam], nodesize=nodesize, ntree=ntree)
  Yqrf <- predict(qrf, data_test[,xnam], what = tau)
  
  # survival forest
  print("doing rsf...")
  v.rsf <- rfsrc(as.formula(paste("Surv(y, status) ~ ", paste(xnam, collapse= "+"))), data = data_train, ntree = ntree, nodesize = nodesize)
  surv.rsf <- predict(v.rsf, newdata = data_test)
  Ysurv <- find_quantile(surv = surv.rsf, max_value = max(data_train$y), tau = tau)
  
  # results
  return(
    list(
      'crf.qrf'=metrics(data_test[[ynam]], Yc.qrf, data_test[[ynam]], tau),
      'crf.grf'=metrics(data_test[[ynam]], Yc.grf, data_test[[ynam]], tau),
      'qrf'=metrics(data_test[[ynam]], Yqrf, data_test[[ynam]], tau),
      'grf'=metrics(data_test[[ynam]], Ygrf, data_test[[ynam]], tau),
      'rsf'=metrics(data_test[[ynam]], Ysurv, data_test[[ynam]], tau),
      
      'crf.qrf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Yc.qrf),
      'crf.grf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Yc.grf),
      'qrf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Yqrf),
      'rsf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Ysurv),
      'grf.c'=1-randomForestSRC::get.cindex(data_test[[ynam]], data_test[[inam]], Ygrf)
    )
  )
}

B <- 160
nodesize <- 0
ntree <- 1000
tau <- 0.9
quantile_result <- list('nodesize' = rep(NA,B), 'crf.grf'=rep(NA,B), 'crf.qrf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'rsf'=rep(NA,B))
cindex_result <- list('nodesize' = rep(NA,B), 'crf.grf'=rep(NA,B), 'crf.qrf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'rsf'=rep(NA,B))

for (t in 1:B) {
  print(t)
  
  if (t%%20 == 1) {
    nodesize <- nodesize + 5
  }
  print(paste0("nodesize is ", nodesize))
  
  tmp <- one_run(t, data.all, xnam, ynam, inam, tau, nodesize, ntree)
  
  quantile_result$nodesize[t] <- nodesize
  quantile_result$crf.grf[t] <- tmp$crf.grf['quantile_loss']
  quantile_result$crf.qrf[t] <- tmp$crf.qrf['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  quantile_result$rsf[t] <- tmp$rsf['quantile_loss']
  
  cindex_result$nodesize[t] <- nodesize
  cindex_result$crf.grf[t] <- tmp$crf.grf.c
  cindex_result$crf.qrf[t] <- tmp$crf.qrf.c
  cindex_result$qrf[t] <- tmp$qrf.c
  cindex_result$rsf[t] <- tmp$rsf.c
  cindex_result$grf[t] <- tmp$grf.c
} 


# plot
require(reshape2)
dd <- melt(as.data.frame(cindex_result), id.vars = 'nodesize')
dd.agg <- aggregate(value ~ nodesize + variable, dd, function(x) c(mean = mean(x), sd = sd(x)/sqrt(50)))
dd.agg$mean <- dd.agg[-1][[2]][,1]
dd.agg$sd <- dd.agg[-1][[2]][,2]
dd.agg$value <- NULL
ggplot(data = dd.agg, aes(x=nodesize, y=mean, colour=variable)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.05)) +
  geom_line() +
  geom_point() +
  labs(x = "Nodesize", y = paste("C-index, tau =", tau), fill = "Nodesize")
ggsave(paste0(dataset, "_c_index_nodesize_tau_0", 10*tau, ".pdf"), width = 5, height = 5, path = "../examples/figs")

# dd <- melt(as.data.frame(quantile_result), id.vars = 'nodesize')
# dd.agg <- aggregate(value ~ nodesize + variable, dd, function(x) c(mean = mean(x), sd = sd(x)/sqrt(10)))
# dd.agg$mean <- dd.agg[-1][[2]][,1]
# dd.agg$sd <- dd.agg[-1][[2]][,2]
# dd.agg$value <- NULL
# ggplot(data = dd.agg, aes(x=nodesize, y=mean, colour=variable)) + 
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.05)) +
#   geom_line() +
#   geom_point() +
#   labs(x = "Nodesize", y = paste("Quantile loss, tau =", tau), fill = "Nodesize")
# ggsave(paste0(dataset, "_quantile_loss_nodesize_tau_0", 10*tau, ".pdf"), width = 5, height = 5, path = "../examples/figs")