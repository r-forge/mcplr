################################################################################
### R_W.r
### Rescorla Wagner model
###
### (c) 2007, M. Speekenbrink
################################################################################

setClass("RescorlaWagner",
  contains="LearningModel",
  representation(
    weights="matrix"
  )
)
setMethod("setPars",signature(object="RescorlaWagner"),
  function(object,pars,unconstrained=FALSE,...,rval=c("object","parameters")) {
    sP.u <- function(pars,fixl) {
      if(!is.null(pars$alpha) && !fixl$alpha) pars$alpha <- exp(pars$alpha)
      if(!is.null(pars$beta) && !fixl$beta) pars$beta <- exp(pars$beta)
      pars
    }
    #object <- callNextMethod()
    parl <- callNextMethod(object=object,pars=pars,rval="parameters",unconstrained=unconstrained,...)
    if(length(object@parStruct@fix)>0) fixl <- relist(object@parStruct@fix,skeleton=object@parameters) else fixl <- relist(rep(FALSE,length(unlist(object@parameters))),skeleton=object@parameters)
    if(unconstrained) {
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) parl[[case]] <- sP.u(parl[[case]],fixl[[case]])
      } else {
        parl <- sP.u(parl,fixl)
      }
    }
    switch(rval,
      object = {
        object@parameters <- parl
        object
      },
      parameters = parl)
  }
)

setMethod("getPars",signature(object="RescorlaWagner"),
  function(object,which="all",unconstrained=FALSE,...) {
    gP.u <- function(pars) {
      pars$alpha <- log(pars$alpha)
      pars$beta <- log(pars$beta)
      pars
    }
    if(unconstrained) {
      pars <- object@parameters
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) pars[[case]] <- gP.u(pars[[case]])
      } else {
        pars <- gP.u(pars)
      }
      pars <- unlist(pars)
      if(which=="free") {
        if(length(object@parStruct@fix)>0) {
          pars <- pars[!object@parStruct@fix]
        }
      }
    } else {
      #pars <- callNextMethod()
      pars <- callNextMethod(object=object,which=which,unconstrained=unconstrained,...)
    }
    return(pars)
  }
)
setMethod("fit",signature(object="RescorlaWagner"),
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        fit <- R_W.fit(x=x,y=y,alpha=pars$alpha,beta=pars$beta,lambda=pars$lambda,ws=pars$ws)
        object@weights[object@nTimes@bt[case]:object@nTimes@et[case],] <- fit$weights
      }
    } else {
      pars <- object@parameters
      fit <- R_W.fit(x=object@x,y=object@y,alpha=pars$alpha,beta=pars$beta,lambda=pars$lambda,ws=pars$ws)
      object@weights <- fit$weights
    }
    return(object)
  }
)

R_W.fit <- function(y,x,alpha,beta,lambda,ws) {
  nx <- ncol(x)
  nt <- nrow(x)
  ny <- ncol(y)
  #if(length(alpha)==1) alpha <- rep(alpha,nx)
  if(length(beta)==1) beta <- rep(beta,ny)
  cid <- matrix(1:(nx*ny),nrow=nx)
  weight <- matrix(,nrow=nt,ncol=length(cid))
  weight[1,] <- w <- ws
  for(i in 1:(nt-1)) {
    for(k in 1:ny) {
      weight[i+1,cid[,k]] <- weight[i,cid[,k]] + alpha*sum(c(y[i,k],1-y[i,k])*beta[k,])*(sum(c(y[i,k],1-y[i,k])*lambda[k,])-t(weight[i,cid[,k]])%*%x[i,])*x[i,]
    }
  }
  return(list(weights=weight))
}

setMethod("predict",signature(object="RescorlaWagner"),
  function(object,...) {
    ny <- ncol(object@y)
    nt <- nrow(object@y)
    nx <- ncol(object@x)
    cid <- matrix(1:(nx*ny),nrow=nx)
    pred <- matrix(nrow=nt,ncol=ny)
    for(i in 1:ny) pred[,i] <- rowSums(object@x*object@weights[,cid[,i]])
    return(pred)
  }
)

setMethod("plot",signature(x="RescorlaWagner",y="missing"),
  function(x,y,...) {
    nx <- ncol(x@x)
    ny <- ncol(x@y)
    nw <- ncol(x@weights)
    nt <- nrow(x@weights)
    dat <- data.frame(
      id=factor(rep(rep(1:x@nTimes@cases,x@nTimes@n),nx*ny),labels="series",1:x@nTimes@cases),
      trial=rep(unlist(sapply(x@nTimes@n,seq,from=1,by=1)),nx*ny),
      w=as.vector(x@weights),
      xid=factor(rep(rep(1:nx,each=nt),ny),labels=paste("x",1:nx,sep="")),
      yid=factor(rep(1:ny,each=nt*nx),labels=paste("y",1:ny,sep="")))
    if(x@nTimes@cases>1) {
      plot <- xyplot(w~trial|yid*id,groups=xid,data=dat,type="l",as.table=TRUE,auto.key=TRUE,...)
    } else {
      plot <- xyplot(w~trial|yid,groups=xid,data=dat,type="l",as.table=TRUE,auto.key=TRUE,...)
    }
    plot
  }
)

RescorlaWagner <- function(formula,parameters=list(alpha=.1,beta=c(1,1),lambda=c(1,0),ws=0),data,subset,fixed=list(alpha=FALSE,beta=TRUE,lambda=TRUE,ws=TRUE),parStruct,remove.intercept=FALSE,base=NULL,ntimes=NULL,replicate=TRUE) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=remove.intercept) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=remove.intercept)
  x <- dat$x
  y <- dat$y
  nx <- ncol(x)
  nt <- nrow(x)
  ny <- ncol(y)
  parfill <- function(parameters) {
    if(is.null(parameters$alpha)) parameters$alpha <- .1
    if(is.null(parameters$beta)) parameters$beta <- c(1,1)
    if(is.null(parameters$lambda)) parameters$lambda <- c(1,0)
    if(is.null(parameters$ws)) parameters$ws <- 0
    if(!is.matrix(parameters$lambda)) {
      if(length(parameters$lambda)!=2 && length(parameters$lambda)!=2*ny) stop("lambda should have length 2 or 2*ny")
      parameters$lambda <- matrix(parameters$lambda,ncol=2,nrow=ny,byrow=T)
    } else {
      if(ncol(parameters$lambda)!=2) stop("lambda should be an ny*2 matrix")
      if(nrow(parameters$lambda)==1) parameters$lambda <- matrix(parameters$lambda,ncol=2,nrow=ny,byrow=T)
      if(nrow(parameters$lamda)!=ny) stop("lambda should be an ny*2 matrix")
    }
    # initialize alpha
    if(length(parameters$alpha)!=1 && length(parameters$alpha)!=nx) stop("alpha must have length 1 or nx")
    # initialize beta
    if(!is.matrix(parameters$beta)) {
      if(length(parameters$beta)!=2 && length(parameters$beta)!=2*ny) stop("beta should have length 2 or 2*ny")
      parameters$beta <- matrix(parameters$beta,ncol=2,nrow=ny,byrow=T)
    } else {
      if(ncol(parameters$beta)!=2) stop("beta should be an ny*2 matrix")
      if(nrow(parameters$beta)==1) parameters$beta <- matrix(parameters$beta,ncol=2,nrow=ny,byrow=T)
      if(nrow(parameters$beta)!=ny) stop("beta should be an ny*2 matrix")
    }
    #if(length(parameters$beta)!=1 && length(parameters$beta)!=ny) stop("beta must have length 1 or ny")
    # intialize ws
    if(length(parameters$ws)!=1 && length(parameters$ws)!=nx*ny) stop("ws must have length 1 or nx*ny")
    parameters
  }
  if(is.null(ntimes) | length(ntimes)==1 | replicate) {
    # intialize lambda
    parameters <- parfill(parameters)
  } else {
    # setup a parlist
    nrep <- length(ntimes)
    # check structure of supplied list
    if(all(unlist(lapply(parameters,is.list))) && length(parameters)==nrep) {
      for(i in 1:nrep) {
        parameters[[i]] <- parfill(parameters[[i]])
      }
    } else {
      parameters <- parfill(parameters)
      parameters <- rep(list(parameters),nrep)
    }
  }
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
  if(missing(parStruct)) {
    parStruct <- ParStruct(parameters=parameters,replicate=replicate,
      fixed = fixed,
      ntimes = {if(missing(ntimes)) NULL else ntimes})
  }
  
  cid <- matrix(1:(nx*ny),nrow=nx)
  
  mod <- new("RescorlaWagner",
    x = x,
    y = y,
    weights = matrix(nrow=nt,ncol=length(cid)),
    parameters = parameters,
    parStruct=parStruct,
    nTimes=nTimes)
  mod <- fit(mod)
  mod
}