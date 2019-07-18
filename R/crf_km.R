source("metrics.R")
source("help_functions.R")
source("Csurv.R")
source("crf.R")
source("cranger.R")

list.of.packages <- c("survival", "randomForest")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library(survival)
library(randomForest)

crf.km <- function(fmla, ntree, nodesize, data_train, data_test, yname, iname,
                tau, method = "randomForest", splitrule = "extratrees", nnb = FALSE) {
  # build Forest model
  if (method == "randomForest") {
    rf <- randomForest(fmla, data = data_train, ntree = ntree, nodesize=nodesize)
    # get proximity matrix
    proxMtx <- rf.getWeights(rf, data_train, data_test, yname, iname)
  } else if (method == "ranger") {
    rf <- ranger(fmla, data = data_train, num.trees = ntree, min.node.size=nodesize, splitrule=splitrule)
    # get proximity matrix
    proxMtx <- ranger.getWeights(rf, data_train, data_test, yname, iname)
  }
  # censor forest
  n <- nrow(data_test)
  ntrain <- nrow(data_train)
  Yc <- rep(NA, n) # to store new predictions
  Ytrain <- data_train[[yname]]
  censorInd <- data_train[[iname]]
  right_censoring <- TRUE
  # find minimum
  for (r in 1:NROW(data_test)) {
    # C survival estimate
    if (nnb) {
      boot.idx <- proxMtx[r, ] > 1/ntrain
      C.boot <- Ytrain[boot.idx]
      i.boot <- 1 - censorInd[boot.idx]
      C.km <- survfit(Surv(C.boot, i.boot) ~ 1, type = 'kaplan-meier')
      C.surv <- stepfun(C.km$time, c(1, C.km$surv))
    } else {
      boot.idx <- proxMtx[r, ] > 0
    }

    max.uncensored <- max(Ytrain[boot.idx & censorInd==1])
    #min.all <- min(Ytrain[boot.idx])
    #candidateY <- seq(min.all, max.uncensored, length.out = 1000)
    candidateY <- Ytrain[boot.idx & Ytrain<=max.uncensored]

    Yc[r] <- candidateY[1]
    min_loss <- 10

    denom <- sapply(Ytrain, function(x) {1*(Ytrain >= x)%*%proxMtx[r, ]})
    denom[denom == 0] <- 1
    base <- (1 - proxMtx[r, ]/denom)^(1 - censorInd)
    for (lambda in candidateY) {
      if (right_censoring) {
        if (nnb) {
          kappa <- C.surv(lambda)
        } else {
          kappa <- C.surv(lambda, Ytrain, base)
          #print(kappa)
        }
        loss1 <- proxMtx[r, ]%*%(1*(Ytrain > lambda))
        loss <- abs((1-tau)*kappa - loss1)
      } else {
        #
      }
      if (loss < min_loss) {
        Yc[r] <- lambda
        min_loss <- loss
      }
    }
  }
  return(
    list(
      'predicted'=Yc
    )
  )
}
