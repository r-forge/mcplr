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
        object@weight[object@nTimes@bt[case]:object@nTimes@et[case],] <- fit$weights
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

rescorlaWagner <- function(formula,parameters=list(alpha=.1,beta=c(1,1),lambda=c(1,0),ws=0),data,alpha,beta,ws=0,intercept=TRUE,base=NULL,ntimes=NULL,replicate=TRUE,fixed,parStruct,subset) {
  if(!intercept) remi <- TRUE else remi <- FALSE
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=remi) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=remi)
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
  if(is.null(ntimes) | replicate) {
    # intialize lambda
    parameters <- parfill(parameters)
  } else {
    # setup a parlist
    nrep <- length(ntimes)
    # check structure of supplied list
    if(all(lapply(parameters,is.list)) && length(parameters)==nrep) {
      for(i in 1:nrep) {
        parameters[[i]] <- parfill(parameters[[i]])
      }
    } else {
      parameters <- parfill(parameters)
      parameters <- rep(list(parameters),nrep)
    }
  }

  if(missing(parStruct)) {
    parStruct <- new("ParStruct")
    if(replicate) parStruct@replicate <- TRUE else parStruct@replicate <- FALSE
    if(is.null(ntimes) | replicate) {
      fix <- rep(FALSE,length(unlist(parameters)))
      fix <- relist(fix,skeleton=parameters)
    } else {
      fix <- rep(FALSE,length(unlist(parameters[[1]])))
      fix <- relist(fix,skeleton=parameters[[1]])
    }
    if(!missing(fixed)) {
      for(i in 1:length(fix)) {
        if(!is.null(fixed[[names(fix)[i]]]) && fixed[[names(fix)[i]]]) fix[[i]] <- rep(TRUE,length(fix[[i]]))
      }
    }
    if(is.null(ntimes) | replicate) {
      parStruct@fix <- unlist(fix)
    } else {
      fix <- rep(list(fix,length(ntimes)))
      parStruct@fix <- unlist(fix)
    }
  }

  if(is.null(ntimes)) ntimes <- nrow(y)
  nTimes <- nTimes(ntimes)
    
  mod <- new("RescorlaWagner",
    x = x,
    y = y,
    weights = matrix(nrow=nt,ncol=ny),
    parameters = parameters,
    parStruct=parStruct,
    nTimes=nTimes)
  mod <- fit(mod)
  mod
}