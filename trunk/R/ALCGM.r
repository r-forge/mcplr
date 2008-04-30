setClass("AlGcmNominal",
  contains="GcmNominal",
  lambda="vector",
  w="matrix"
)
setClass("AlGcmInterval",
  contains="GcmInterval",
  lambda="vector",
  w="matrix"
)
setMethod("fit",signature(object="AlGcmInterval"),
  function(object,ntimes=NULL,...) {
    if(is.null(ntimes)) ntimes <- nrow(object@y)
    lt <- length(ntimes)
    et <- cumsum(ntimes)
    bt <- c(1,et[-lt]+1)
    for(case in 1:lt) {
      repl <- getReplication(object,ntimes=ntimes,case=case,...)
      x <- repl$x
      y <- repl$y
      pars <- repl$parameters
      fit <- algcminterval.fit(x=x,y=y,parameters=pars,distance=object@distance,similarity=object@similarity,sampling=object@sampling)
      object@weight[[case]] <- fit$weight
      object@lambda[bt[case]:et[case]] <- fit$lambda
      object@w[bt[case]:et[case],] <- fit$weight
    }
    return(object)
  }
)

algcminterval.fit <- function(x,y,parameters,distance,similarity,sampling,...) {
  grad <- function(pred,x,y,p,d,s,pq,pr,v) {
    nt <- nrow(x)
    nx <- ncol(x)
    out <- vector()
    N <- sum(p*s)
    for(i in 1:nx) {
      out[i] <- -(pq/pr)*((pred-y[nt,])/N)*sum((y[-nt,]-pred)*p[-nt]*s[-nt]*d[-nt]^((pq-pr)/pr)*exp(v[i])*(x[nt,i]-x[-nt,i])^pr)
    }
    out
  }
  nt <- nrow(y)
  nx <- ncol(x)
  v <- matrix(0,nrow=nt,ncol=nx)
  weight <- matrix(0,nrow=nt,ncol=nt)
  pr <- parameters$r
  pq <- parameters$q
  if(is.null(pr)) if(attr(distance,"name") == "euclidian") pr <- 2 else pr <- 1
  if(is.null(pq)) if(attr(similarity,"name") == "gaussian") pq <- 2 else pq <- 1
  v[1,] <- v[2,] <- log(parameters$w*parameters$lambda^(pq/pr))
  for(i in 2:nt) {
    dis <- distance(x=x,parameters=parameters,...)[i,]
    sim <- similarity(distance=dis,parameters=parameters,...)
    sam <- sampling(nt=i,parameters=parameters,...)[,i]
    w <- sim*sam
    if(sum(w)==0) {
      w <- rep(1,length(w))
      sim <- rep(1e-10,length(sim))
    }
    w <- w/sum(w)
    w[is.nan(w)] <- 0 # FIX ME
    weight[1:i,i] <- w
    pred <- sum(y[1:i]*w)
    if(i < nt) {
      v[i+1,] <- v[i,] - parameters$eta*grad(pred,as.matrix(x[1:i,]),as.matrix(y[1:i,]),sam,dis,sim,pq,pr,v[i,])
      #v[i+1,][v[i+1,]<0] <- 1e-10
      parameters$lambda <- sum(exp(v[i+1,]))
      parameters$w <- exp(v[i+1,])/parameters$lambda
      parameters$lambda <- parameters$lambda^(pr/pq)
    }
    lambda <- rowSums(exp(v))^(pr/pq)
    w <- exp(v)/rowSums(exp(v))
    #w <- t(t(w)/colSums(w))
  }
  return(list(weight=weight,w=w,lambda=lambda))
}

algcminterval.fit.old <- function(x,y,parameters,distance,similarity,sampling,...) {
  grad <- function(pred,x,y,p,d,s,pq,pr) {
    nt <- nrow(x)
    nx <- ncol(x)
    out <- vector()
#    for(i in 1:nx) {
#      out[i] <- -sum((pred-y[nt,])*p[-nt]*s[-nt]*y[-nt,]*d[-nt]^(pq/pr)*(pq/pr)*(x[nt,i]-x[-nt,i])^pr)
#    }
    N <- sum(p*s)
    for(i in 1:nx) {
      out[i] <- -(pq/pr)*((pred-y[nt,])/N)*sum((y[-nt,]-pred)*p[-nt]*s[-nt]*d[-nt]^((pq-pr)/pr)*(x[nt,i]-x[-nt,i])^pr)
    }
    out
  }
  nt <- nrow(y)
  nx <- ncol(x)
  v <- matrix(0,nrow=nt,ncol=nx)
  weight <- matrix(0,nrow=nt,ncol=nt)
  pr <- parameters$r
  pq <- parameters$q
  if(is.null(pr)) if(attr(distance,"name") == "euclidian") pr <- 2 else pr <- 1
  if(is.null(pq)) if(attr(similarity,"name") == "gaussian") pq <- 2 else pq <- 1
  v[1,] <- v[2,] <- parameters$w*parameters$lambda^(pq/pr)
  for(i in 2:nt) {
    #dis <- distance(x=x,parameters=parameters,...)[i,]
    #sim <- similarity(distance=dis,parameters=parameters,...)
    #sam <- sampling(nt=i,parameters=parameters,...)
    dis <- distance(x=x,parameters=parameters)[1:i,i]
    sim <- similarity(distance=dis,parameters=parameters)
    sam <- sampling(nt=i,parameters=parameters)[,i]
    w <- sim*sam
    if(sum(w)==0) {
      w <- rep(1,length(w))
      sim <- rep(1e-10,length(sim))
    }
    w <- w/sum(w)
    w[is.nan(w)] <- 0 # FIX ME
    weight[1:i,i] <- w
    pred <- sum(y[1:i]*w)
    if(i < nt) {
      v[i+1,] <- v[i,] - parameters$eta*grad(pred,as.matrix(x[1:i,]),as.matrix(y[1:i,]),sam,dis,sim,pq,pr)
      v[i+1,][v[i+1,]<0] <- 1e-10
      parameters$lambda <- sum(v[i+1,])
      parameters$w <- v[i+1,]/parameters$lambda
      parameters$lambda <- parameters$lambda^(pr/pq)
    }
    #w <- t(t(w)/colSums(w))
  }
  return(list(weight=weight,v=v))
}


# with exponential parametrization
