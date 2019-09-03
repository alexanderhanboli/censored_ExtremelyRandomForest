# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

generate_data <- function(n = 50, p = 60, mu = NULL, beta0 = NULL, 
                              whichCov = "identity",
                              toeplitzCoef = 0.3){
  library(mvtnorm)
  n <- n
  p <- p
  if(is.null(mu)) mu <- rep(0, p)
  
  cov.kind.pool = c("random","toeplitz","fixed", "identity") # types of different cov matrices
  if(whichCov %in% cov.kind.pool) {cov.kind <- whichCov}
  else if(whichCov %in% c(1,2,3,4)) {cov.kind <- cov.kind.pool[whichCov]}
  else stop("Currently do not support such covariance matrix type. Please enter 1,2,3,4 or one of 'toeplitz', 'random', 'fixed', 'identity'")
  
  if (cov.kind == "random"){
    Sigma.factor <- matrix(rnorm(p*p)/sqrt(p), nrow = p)
    Sigma <- t(Sigma.factor)%*%Sigma.factor  
  }else if (cov.kind == "toeplitz"){
    Sigma <- toeplitz(toeplitzCoef^(0:(p-1))) # set 0: identity matrix
  }else if (cov.kind == "identity"){
    Sigma <- diag(rep(1,p))
  }else{
    Sigma <- matrix(rep(fixedCoef, p*p), nrow = p)
    diag(Sigma) <- 1
  }
  
  # generate x
  x <- rmvnorm(n = n, mean = mu, sigma = Sigma)
  #x <- matrix(rnorm(n*p), nrow = n)
  if(is.null(beta0)) {
    beta0 <- rep(0, p)
    beta0[1:5] = c(1, 1, 1, 1, 1)
  }
  if(length(beta0) != p){
    stop(paste("Please enter beta0 with correct length: ", p))
  }
  
  y <- x%*%beta0 + rnorm(n)
  
  return(list(
    'x' = x,
    'y' = y,
    'n' = n,
    'p' = p,
    'beta0' = beta0
  ))
}


quantile_loss <- function(u, tau) {
  return (u*(tau - (u < 0)))
}

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