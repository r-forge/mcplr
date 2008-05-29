# Single Layer Feedforward Network
setClass("SLFN",
  contains="LearningModel",
  representation(
    weight="array",
    activation="character",
    actfun="function",
    gradient="function",
    window.size="numeric")
)

setMethod("setPars",signature(object="SLFN"),
  function(object,pars,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    sP.u <- function(pars,fixl) {
      if(!is.null(pars$eta) && !fixl$eta) pars$eta <- exp(pars$eta)
      if(!is.null(pars$alpha) && !fixl$alpha) pars$alpha <- exp(pars$alpha)
      if(!is.null(pars$beta) && !fixl$beta) pars$beta <- ifelse(pars$beta>600,1,exp(pars$beta)/(1+exp(pars$beta)))
      pars
    }
    #object <- callNextMethod()
    parl <- callNextMethod(object=object,pars=pars,internal=internal,...,rval="parameters")
    if(length(object@parStruct@fix)>0) fixl <- relist(object@parStruct@fix,skeleton=object@parameters) else fixl <- relist(rep(FALSE,length(unlist(object@parameters))),skeleton=object@parameters)
    if(internal && is.null(object@parStruct@constraints)) {
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
    #object@parameters
    #return(object)
  }
)

setMethod("getPars",signature(object="SLFN"),
  function(object,which="all",internal=FALSE,...) {
    gP.u <- function(pars) {
      pars$alpha <- log(pars$alpha)
      pars$eta <- log(pars$eta)
      pars$beta <- log(pars$beta/(1-pars$beta))
      pars
    }
    if(internal && is.null(object@parStruct@constraints)) {
      pars <- object@parameters
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) pars[[case]] <- gP.u(pars[[case]])
      } else {
        pars <- gP.u(pars)
      }
      pars <- unlist(as.relistable(pars))
      if(which=="free") {
        if(length(object@parStruct@fix)>0) {
          pars <- pars[!object@parStruct@fix]
        }
      }
    } else {
      pars <- callNextMethod(object=object,which=which,internal=internal,...)
      #pars <- callNextMethod()
    }
    return(pars)
  }
)
setMethod("fit","SLFN",
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        #x <- object@x[bt[case]:et[case],]
        #y <- object@y[bt[case]:et[case],]
        #eta <- object@parameters$eta[case,]
        #alpha <- object@parameters$alpha[case,]
        #beta <- object@parameters$beta[case,]
        #ws <- object@parameters$ws[case,]
        fit <- slfn.fit(x=x,y=y,eta=pars$eta,alpha=pars$alpha,beta=pars$beta,ws=pars$ws,grad=object@gradient,window.size=object@window.size)
        object@weight[object@nTimes@bt[case]:object@nTimes@et[case],] <- fit$weight
      }
    } else {
      fit <- slfn.fit(x=object@x,y=object@y,eta=object@parameters$eta,alpha=object@parameters$alpha,beta=object@parameters$beta,ws=object@parameters$ws,grad=object@gradient,window.size=object@window.size)
      object@weight <- fit$weight
    }
    return(object)
  }
)
slfn.fit <- function(y,x,ws,eta,alpha,beta,grad,window.size=0) {
  if(!all(alpha >= 0)) stop("negative alpha not allowed")
  if(!all(beta >= 0)) stop("negative beta not allowed")
  if(!all(beta <= 1)) stop("beta must be maximum 1")
  if(!all(eta >= 0)) stop("negative eta not allowed")
  if(window.size < 0) stop("negative window.size not allowed")
  #x <- as.matrix(x)
  nt <- nrow(x)
  nx <- ncol(x)
  ny <- ncol(y)
  if(length(ws)!=nx & length(ws) != nx*ny & length(ws)!=1) warning("ws should have length 1, ncol(x), or ncol(x)*ncol(y)")
  if(length(eta)!=nx & length(eta) != nx*ny & length(eta)!=1) warning("eta should have length 1, ncol(x), or ncol(x)*ncol(y)")
  if(length(alpha)!=nx & length(alpha) != nx*ny & length(alpha)!=1) warning("alpha should have length 1, ncol(x), or ncol(x)*ncol(y)")
  if(length(beta)!=nx & length(beta)!= nx*ny & length(beta)!=1) warning("beta should have length 1, ncol(x), or ncol(x)*ncol(y), or number of cues")
  xnames <- dimnames(x)[[2]]
  
  #if(length(ws)==1) ws <- rep(ws,length=nx)
  ws <- matrix(ws,ncol=ny,nrow=nx)
  #cid <- matrix(1:(nx*ny),nrow=nx)
  weight <- array(,dim=c(nt,nx,ny)) #nrow=nt,ncol=length(cid))
  weight[1,,] <- w <- ws
  for(i in 1:(nt-1)) {
    j <- ifelse(i-window.size > 0, i-window.size,1)
    #for(k in 1:ny) weight[i+1,cid[,k]] <- weight[i,cid[,k]] - (eta/(i)^alpha)*(1/abs(i-j+1))*grad(y=y[j:i,k],x=x[j:i,],w=weight[i,cid[,k]]) + ifelse(i>1,beta*(weight[i,cid[,k]]-as.numeric(weight[i-1,cid[,k]])),0)
    weight[i+1,,] <- weight[i,,] - (eta/(i)^alpha)*(1/abs(i-j+1))*grad(y=y[j:i,],x=x[j:i,],w=weight[i,,]) + ifelse(i>1,beta*(weight[i,,]-weight[i-1,,]),0)
  }
  #if(ny==1) dimnames(weight)[[2]] <- xnames
  return(list(weight=weight))
}

slfn <- function(formula,parameters=list(eta=.01,alpha=0,beta=0,ws=0),type=c("linear","logistic"),data,window.size=0,intercept=TRUE,base=NULL,ntimes=NULL,replicate=T,subset) {
  lin <- function(x) {
    x
  }
  logis <- function(x) {
    1/(1+exp(-x))
  }
  grad.logis <- function(y,x,w) {
    if(!is.matrix(x)) x <- matrix(x,ncol=length(x))
    if(!is.matrix(y)) y <- matrix(y,ncol=length(y))
    t(x)%*%(logis(x%*%w) - y)
#    
#    if(!is.matrix(x)) {
#      return(as.numeric(x*(logis(x%*%w) - y)))
#    } else return(as.numeric(t(x)%*%(logis(x%*%w) - y)))
  }
  grad.lin <- function(y,x,w) {
    if(!is.matrix(x)) x <- matrix(x,ncol=length(x))
    if(!is.matrix(y)) y <- matrix(y,ncol=length(y))
    t(x)%*%(x%*%w - y)
#    {
#      return(as.numeric(x*(x%*%w - y)))
#    } else {
#      tmp <- t(matrix(x,nrow=ncol(y),ncol=300,byrow=T)*as.vector((x%*%w - y)))
#      return(as.numeric(t(x)%*%(x%*%w - y)))
#    }
  }
  if(!intercept) remi <- TRUE else remi <- FALSE
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=remi) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=remi)
  x <- dat$x
  y <- dat$y
  fun <- switch(type[1],
    logistic = logis,
    linear = lin,
    lin
  )
  grad <- switch(type[1],
    logistic = grad.logis,
    linear = grad.lin,
    grad.lin)
  #if(missing(ws)) ws <- rep(0,ncol(x))
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$eta)) parameters$eta <- .01
    if(is.null(parameters$alpha)) parameters$alpha <- 0
    if(is.null(parameters$beta)) parameters$beta <- 0
    if(is.null(parameters$ws)) parameters$ws <- 0
  } else {
    if(length(parameters)==0) {
      nrep <- length(ntimes)
      parameters$eta <- .01
      parameters$alpha <- 0
      parameters$beta <- 0
      parameters$ws <- 0
      parameters <- rep(list(parameters),nrep)
    } else warning("there is no validity check for the given parameters when combined with ntimes and replicate=FALSE \n Please make sure the supplied list is valid")
  }
  if(is.null(ntimes)) ntimes <- nrow(y)
  nTimes <- nTimes(ntimes)
  new("SLFN",x=x,y=y,parameters=parameters,nTimes=nTimes,activation=type[1],actfun=fun,gradient=grad,window.size=window.size)
}

setMethod("predict",signature(object="SLFN"),
  function(object,type="link",...) {
    ny <- ncol(object@y)
    nt <- nrow(object@y)
    nx <- ncol(object@x)
    #cid <- matrix(1:(nx*ny),nrow=nx)  
    pred <- apply(array(as.numeric(object@x)*as.numeric(object@weight),dim=c(nt,nx,ny)),c(1,3),sum)
    #matrix(nrow=nt,ncol=ny)      
    #for(i in 1:ny) pred[,i] <- rowSums(object@x*object@weight[,cid[,i]])
    if(type=="response") pred <- object@actfun(pred)
    return(pred)
  }
)

