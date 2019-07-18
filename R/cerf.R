cerf = function(data, y_name, delta_name, ntrees=100, maxnodes=31, nodesize=10, as.value = TRUE,
                keplan.meier=FALSE, k.neigh=10, Kth.tau.cond.quantile=FALSE,
                plot=TRUE){



  y = data[ , y_name]
  delta_data = data[ , delta_name]

  formula_rf = reformulate(names(data[ ,!(names(data) %in% c(y_name, delta_name)), drop=F]),
                           y_name)

  rf = randomForest(formula=formula_rf, data=data,
                    ntree=ntrees, maxnodes = maxnodes, nodesize=nodesize)

  pred = predict(rf,data[ ,!(names(data) %in% c(y_name, delta_name)), drop=F], nodes = T)
  nodes = attributes(pred)$nodes #
  nodesize = nodes
  weight = nodes
  quant_dist=c()
  quants = seq(0.001,0.999,length.out = 100)
  weight_x=c()


  for (k in 1:ntrees){
    # get nodesize based on number of unique terminal nodes for each tree
    mapping=plyr::count(nodes[,k])
    nodesize[,k] = plyr::mapvalues(nodes[,k], from=mapping$x, to=mapping$freq)
  }

  main_calc = function(u, sequ){
    u=u
    # Calculate weights
    for (k in 1:ncol(nodes)){
      # check if x0 in the same node as other Xs
      o=outer(nodes[u,k],nodes[,k],'==')
      # final weight
      weight[,k] = o/nodesize[,k] #sumCol = 1
    }
    weight_x = matrixStats::rowSums2(weight)/ncol(nodes) # normalize wrt to ntree

    g=sapply(sequ,fun,weight_x=weight_x) #distribution of quantiles estimates based on tau

    if(as.value == TRUE){
      quant_dist=quantile(y,g)
      quant_dist
    }else{
      g
    }

  }

  surv_min = function(q,tau,delta_data,y,weight_x,k.neigh){


    cond = quantile(y,q)

    #survival function - Beran estimator
    if(keplan.meier==FALSE){
      ys = y[y<=cond & weight_x!=0]
      ws = weight_x[y<=cond & weight_x!=0]
      delta = delta_data[y<=cond & weight_x!=0]

      product = sapply(1:length(ys),function(i) (1- ws[i]/sum(weight_x[y>=ys[i]]))^(1-delta[i]))
      if(is.list(product)) product=1 #if failed
      product = prod(product)
    }

    # survival function - Keplan-Meier estimator usin knn
    # k.neigh - k-nearest neigbours
    # sort vectors according to weight_x and get first k

    if(keplan.meier==TRUE){
      ys = y[order(weight_x, decreasing = 1) & y<=cond & weight_x!=0]
      delta = delta_data[order(weight_x, decreasing = 1) & y<=cond & weight_x!=0]
      if(length(ys)>k.neigh){
        ys = ys[1:k.neigh]
        delta = delta[1:k.neigh]
      }

      product = sapply(1:length(ys), function(i) (1 - 1/ (sum(y>=ys[i])))^(1-delta[i]))
      if(is.list(product)) product=1 #if failed
      product = prod(product)
    }

    s_curr = abs((1-tau)*product - sum(weight_x[y>cond]))
    if(length(s_curr)==0) s_curr = 1
    s_curr

  }

  fun = function(tau,weight_x){
    t=tau
    s = sapply(quants,surv_min,tau=t,delta_data=delta_data,y=y,weight_x=weight_x,k.neigh=k.neigh)
    quants[which.min(s)]
  }

  if(Kth.tau.cond.quantile!=FALSE & Kth.tau.cond.quantile<=1000){
    sequ = seq(0.01,0.99,length.out = 30) #tau
    output = main_calc(Kth.tau.cond.quantile,sequ)

  } else{
    sequ = c(0.2,0.5,0.8)
    #sequ=0.5
    output = sapply(1:nrow(data),main_calc, sequ = sequ)

  }

  if(plot==TRUE){
    plot(seq(0,1,length.out = length(output)),
         sort(output),
         type = 'l',
         ylab='quantile',
         xlab='response', lwd=1.5, col='#4C00FFFF')
    xaxis = seq(0,1,length.out = length(y))
    if(Kth.tau.cond.quantile==FALSE){
      if(as.value==TRUE){
        lines(xaxis,sort(y), col='#00FF4DFF')
      } else {
        lines(xaxis,quantile(y,xaxis), col='#00FF4DFF')
      }
    }
  }
  return(output)

}
