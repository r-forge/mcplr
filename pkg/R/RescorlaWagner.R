################################################################################
### R_W.r
### Rescorla Wagner model
###
### (c) 2007, M. Speekenbrink
################################################################################

setClass("RescorlaWagner",
  contains="LearningModel",
  representation(
    weights="array"
  )
)

setClass("ContinuousRescorlaWagner",
  contains="RescorlaWagner"
)

setMethod("canRepar",signature(object="RescorlaWagner"),
  function(object,...) {
    repar <- TRUE
    if(!is(object@parStruct@constraints,"S4")) repar <- FALSE # returns "S4" when it is a virtual class
    #if(is(object@parStruct@constraints,"LinConstraintsList") || is(object@parStruct@constraints,"BoxConstraintsList")) {
    #  res <- FALSE
    #}
    # TODO: check whether \lambda is fixed
    return(repar)
  }
)

setMethod("setTransPars",signature(object="RescorlaWagner"),
  function(object,pars,...,rval=c("object","parameters")) {
    sP.u <- function(pars,fixl) {
      if(!is.null(pars$alpha) && !fixl$alpha) pars$alpha <- exp(pars$alpha)
      if(!is.null(pars$beta) && !fixl$beta) pars$beta <- exp(pars$beta)
      pars
    }
    #object <- callNextMethod()
    parl <- getPars(object=object,pars=pars,rval="parameters",...)
    if(length(object@parStruct@fix)>0) fixl <- relist(object@parStruct@fix,skeleton=object@parStruct@parameters) else fixl <- relist(rep(FALSE,length(unlist(object@parStruct@parameters))),skeleton=object@parStruct@parameters)
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) parl[[case]] <- sP.u(parl[[case]],fixl[[case]])
      } else {
        parl <- sP.u(parl,fixl)
      }
    switch(rval,
      object = {
        object@parStruct@parameters <- parl
        object
      },
      parameters = parl)
  }
)


setMethod("getTransPars",signature(object="RescorlaWagner"),
  function(object,which="all",...) {
    gP.u <- function(pars) {
      pars$alpha <- log(pars$alpha)
      pars$beta <- log(pars$beta)
      pars
    }
      pars <- object@parStruct@parameters
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
    return(pars)
  }
)
setMethod("runm",signature(object="RescorlaWagner"),
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        runm <- R_W.runm(x=x,y=y,alpha=pars$alpha,beta=pars$beta,lambda=pars$lambda,ws=pars$ws)
        object@weights[,,object@nTimes@bt[case]:object@nTimes@et[case]] <- runm$weights
      }
    } else {
      pars <- object@parStruct@parameters
      runm <- R_W.runm(x=object@x,y=object@y,alpha=pars$alpha,beta=pars$beta,lambda=pars$lambda,ws=pars$ws)
      object@weights <- runm$weights
    }
    return(object)
  }
)

setMethod("runm",signature(object="ContinuousRescorlaWagner"),
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        runm <- C_R_W.runm(x=x,y=y,alpha=pars$alpha,beta=pars$beta,lambda=pars$lambda,ws=pars$ws)
        object@weights[,,object@nTimes@bt[case]:object@nTimes@et[case]] <- runm$weights
      }
    } else {
      pars <- object@parStruct@parameters
      runm <- C_R_W.runm(x=object@x,y=object@y,alpha=pars$alpha,beta=pars$beta,lambda=pars$lambda,ws=pars$ws)
      object@weights <- runm$weights
    }
    return(object)
  }
)

R_W.runm <- function(y,x,alpha,beta,lambda,ws) {
  nx <- ncol(x)
  nt <- nrow(x)
  ny <- ncol(y)
  #if(length(alpha)==1) alpha <- rep(alpha,nx)
  if(length(lambda) == 1) lambda <- rep(lambda,length=ny)
  if(length(beta)==1) beta <- matrix(beta,nrow=ny,ncol=2)
  if(is.vector(beta)) beta <- matrix(beta,nrow=ny,ncol=2)
  
  #cid <- matrix(1:(nx*ny),nrow=nx)
  weights <- array(0.0,dim=c(nx,ny,nt))
  weights[,,1] <- w <- ws
  for(i in 1:(nt-1)) {
    for(k in 1:ny) {
      if(y[i,k] == 0) bet <- beta[k,1] else bet <- beta[k,2]
      weights[,k,i+1] <- weights[,k,i] + alpha*bet*(lambda[k]*y[i,k]-t(weights[,k,i])%*%x[i,])*x[i,]
    }
  }
  return(list(weights=weights))
}

C_R_W.runm <- function(y,x,alpha,beta,lambda,ws) {
  nx <- ncol(x)
  nt <- nrow(x)
  ny <- ncol(y)
  if(length(alpha)==1) alpha <- rep(alpha,nx)
  if(length(alpha)!=nx) stop("alpha needs to have length 1 or nx")
  if(length(beta)==1) beta <- rep(beta,ny)
  if(length(alpha)!=nx) stop("beta needs to have length 1 or ny")
  if(length(lambda)==1) y <- y*lambda else {
    if(length(lambda)!=ny) stop("lambda needs to have length 1 or ny") 
    y <- t(t(y)*lambda)
  }
  eta <- outer(alpha,beta)
  weights <- array(0.0,dim=c(nx,ny,nt))
  weights[,,1] <- ws
  ypred <- array(0.0,dim=c(ny,nt))
  weights <- array(.C("slfn",
	  y=as.double(y),
	  ny=as.integer(ny), 
	  x=as.double(x), 
	  nx=as.integer(nx), 
	  bt=as.integer(1),
	  et=as.integer(nt),
	  lt=as.integer(1), 
	  eta=as.double(eta), 
	  actfun = as.integer(1),
	  w = as.double(weights),
	  ypred = as.double(ypred),
	  PACKAGE="mcplR"
  )$w,dim=c(nx,ny,nt))
  
  return(list(weights=weights))
}
  #if(is.vector(beta)) beta <- matrix(beta,nrow=ny,ncol=2)
#  cid <- matrix(1:(nx*ny),nrow=nx)
#  weight <- matrix(,nrow=nt,ncol=length(cid))
#  weight[1,] <- w <- ws
#  for(i in 1:(nt-1)) {
#    for(k in 1:ny) {
#      #if(y[i,k] == 0) bet <- beta[k,1] else bet <- beta[k,2]
#      weight[i+1,cid[,k]] <- weight[i,cid[,k]] + alpha*beta[k]*(lambda[k]*y[i,k]-t(weight[i,cid[,k]])%*%x[i,])*x[i,]
#    }
#  }
#  return(list(weights=weight))
#}

setMethod("predict",signature(object="RescorlaWagner"),
  function(object,...) {
    ny <- ncol(object@y)
    nt <- nrow(object@y)
    nx <- ncol(object@x)
    pred <- matrix(0.0,ncol=ny,nrow=nt)
    for(i in 1:ny) {
      pred[,i] <- rowSums(t(object@weights[,i,])*object@x)
    }
    #cid <- matrix(1:(nx*ny),nrow=nx)
    #pred <- matrix(nrow=nt,ncol=ny)
    #if(!is.matrix(object@weights)) for(i in 1:ny) pred[,i] <- colSums(t(object@x)*object@weights[,i,]) else pred <- colSums(t(object@x)*object@weights)
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
      plot <- xyplot(w~trial|yid*id,groups=dat$xid,data=dat,type="l",as.table=TRUE,auto.key=TRUE,...)
    } else {
      plot <- xyplot(w~trial|yid,groups=dat$xid,data=dat,type="l",as.table=TRUE,auto.key=TRUE,...)
    }
    plot
  }
)

RescorlaWagner <- function(formula,parameters=list(alpha=.1,beta=1,lambda=1,ws=0),data,subset,fixed=list(alpha=FALSE,beta=TRUE,lambda=TRUE,ws=TRUE),parStruct,remove.intercept=FALSE,base=NULL,ntimes=NULL,replicate=TRUE) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=remove.intercept) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=remove.intercept)
  x <- dat$x
  y <- dat$y
  if(!all(y %in% c(0,1))) {
  	warning("The criterion variable is not dichotomous. Will create a Rescorla Wagner model for continuous criterion.")
  	cont <- TRUE
  } else cont <- FALSE
  nx <- ncol(x)
  nt <- nrow(x)
  ny <- ncol(y)
  parfill <- function(parameters) {
    if(is.null(parameters$alpha)) parameters$alpha <- .1
    if(is.null(parameters$beta)) {
    	#if(!cont) parameters$beta <- c(1,1) else parameters$beta <- 1
    	parameters$beta <- 1
    }
    if(is.null(parameters$lambda)) parameters$lambda <- 1
    if(is.null(parameters$ws)) parameters$ws <- 0
    # initialize alpha
    if(length(parameters$alpha)!=1 && length(parameters$alpha)!=nx) stop("alpha must have length 1 or nx")
    # initialize beta
    if(!cont) {
		  if(!is.matrix(parameters$beta)) {
		    if(!(length(parameters$beta) %in%c(1,2)) && length(parameters$beta)!=2*ny) stop("beta should have length 1, 2 or 2*ny")
		    if(length(parameters$beta) > 1) parameters$beta <- matrix(parameters$beta,ncol=2,nrow=ny,byrow=T)
		  } else {
		    if(ncol(parameters$beta)!=2) stop("beta should be an ny*2 matrix")
		    if(nrow(parameters$beta)==1) parameters$beta <- matrix(parameters$beta,ncol=2,nrow=ny,byrow=T)
		    if(nrow(parameters$beta)!=ny) stop("beta should be an ny*2 matrix")
		  }
    } else {
    	if(length(parameters$beta)!=ny && length(parameters$beta)!=2) stop("beta should have length 1 or ny")
    	parameters$beta <- as.vector(parameters$beta)
    	#}
    }
    #if(length(parameters$beta)!=1 && length(parameters$beta)!=ny) stop("beta must have length 1 or ny")
    # intialize ws
    if(length(parameters$lambda) != ny && length(parameters$lambda) != 1) stop("lambda should have length 1 or ny")
    if(length(parameters$ws)!=1 && length(parameters$ws)!=nx*ny) stop("ws must have length 1 or nx*ny")
    parameters <- parameters[c("alpha","beta","lambda","ws")] # put in correct order
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
    parStruct <- ParStruct(
      parameters=parameters,
      replicate=replicate,
      fixed = fixed,
      ntimes = {if(missing(ntimes)) NULL else ntimes})
  }
  
  cid <- matrix(1:(nx*ny),nrow=nx)
  
  if(!cont) {
		mod <- new("RescorlaWagner",
		  x = x,
		  y = y,
		  weights = array(0.0,dim=c(nx,ny,nt)),
		  parStruct=parStruct,
		  nTimes=nTimes)
  } else {
		mod <- new("ContinuousRescorlaWagner",
		  x = x,
		  y = y,
		  weights = array(0.0,dim=c(nx,ny,nt)),
		  parStruct=parStruct,
		  nTimes=nTimes)
  }  
  mod <- runm(mod)
  mod
}
