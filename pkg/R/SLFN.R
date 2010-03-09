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
    if(length(object@parStruct@fix)>0) fixl <- relist(object@parStruct@fix,skeleton=object@parStruct@parameters) else fixl <- relist(rep(FALSE,length(unlist(object@parStruct@parameters))),skeleton=object@parStruct@parameters)
    if(internal && is.null(object@parStruct@constraints)) {
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) parl[[case]] <- sP.u(parl[[case]],fixl[[case]])
      } else {
        parl <- sP.u(parl,fixl)
      }
    }
    switch(rval,
      object = {
        object@parStruct@parameters <- parl
        object
      },
      parameters = parl)
    #object@parStruct@parameters
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
      pars <- object@parStruct@parameters
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
setMethod("runm","SLFN",
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        #x <- object@x[bt[case]:et[case],]
        #y <- object@y[bt[case]:et[case],]
        #eta <- object@parStruct@parameters$eta[case,]
        #alpha <- object@parStruct@parameters$alpha[case,]
        #beta <- object@parStruct@parameters$beta[case,]
        #ws <- object@parStruct@parameters$ws[case,]
        runm <- slfn.runm(x=x,y=y,eta=pars$eta,alpha=pars$alpha,beta=pars$beta,ws=pars$ws,grad=object@gradient,window.size=object@window.size)
        object@weight[object@nTimes@bt[case]:object@nTimes@et[case],,] <- runm$weight
      }
    } else {
      runm <- slfn.runm(x=object@x,y=object@y,eta=object@parStruct@parameters$eta,alpha=object@parStruct@parameters$alpha,beta=object@parStruct@parameters$beta,ws=object@parStruct@parameters$ws,grad=object@gradient,window.size=object@window.size)
      object@weight <- runm$weight
    }
    return(object)
  }
)
slfn.runm <- function(y,x,ws,eta,alpha,beta,grad,window.size=0) {
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
    if(i>1 & any(beta > 0)) {
    	weight[i+1,,] <- weight[i,,] - (eta/(i)^alpha)*(1/abs(i-j+1))*grad(y=y[j:i,],x=x[j:i,],w=weight[i,,]) + beta*(weight[i,,]-weight[i-1,,])
    } else {
    	weight[i+1,,] <- weight[i,,] - (eta/(i)^alpha)*(1/abs(i-j+1))*grad(y=y[j:i,],x=x[j:i,],w=weight[i,,])
    }
  }
  #if(ny==1) dimnames(weight)[[2]] <- xnames
  return(list(weight=weight))
}

SLFN <- function(formula,parameters=list(eta=.01,alpha=0,beta=0,ws=0),type=c("linear","logistic"),data,subset,fixed,parStruct,window.size=0,remove.intercept=FALSE,base=NULL,ntimes=NULL,replicate=T) {
  type <- match.arg(type)
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
  }
  grad.lin <- function(y,x,w) {
    if(!is.matrix(x)) x <- matrix(x,ncol=length(x))
    if(!is.matrix(y)) y <- matrix(y,ncol=length(y))
    t(x)%*%(x%*%w - y)
  }

  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=remove.intercept) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=remove.intercept)
  x <- dat$x
  y <- dat$y
  fun <- switch(type,
    logistic = logis,
    linear = lin,
    lin
  )
  grad <- switch(type,
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
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
  if(missing(parStruct)) {
    parStruct <- ParStruct(parameters=parameters,replicate=replicate,
      fixed = if(missing(fixed)) NULL else fixed,
      ntimes = if(missing(ntimes)) NULL else ntimes)
  }
  
  mod <- new("SLFN",
    x=x,
    y=y,
    weight=array(dim=c(sum(nTimes@n),ncol(x),ncol(y))),
    parStruct=parStruct,
    nTimes=nTimes,
    activation=type[1],
    actfun=fun,
    gradient=grad,
    window.size=window.size
  )
  mod <- runm(mod)
  mod
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

setMethod("plot",signature(x="SLFN",y="missing"),
  function(x,y,...) {
    nx <- ncol(x@x)
    ny <- ncol(x@y)
    nw <- prod(dim(x@weight)[c(2,3)])
    nt <- dim(x@weight)[1]
    dat <- data.frame(
      id=factor(rep(rep(1:x@nTimes@cases,x@nTimes@n),nx*ny),labels="series",1:x@nTimes@cases),
      trial=rep(unlist(sapply(x@nTimes@n,seq,from=1,by=1)),nx*ny),
      w=as.vector(x@weight),
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
