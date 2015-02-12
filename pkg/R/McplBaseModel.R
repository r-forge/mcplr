setClass("McplBaseModel",
  representation(
    x="matrix",
    y="matrix",
    parStruct="ParStruct",
    nTimes="NTimes"
  )
)

setMethod("getPars",signature(object="McplBaseModel"),
  function(object,...) {
    getPars(object=object@parStruct,...)
  }
)

setMethod("getFreePars",signature(object="McplBaseModel"),
 function(object,...) {
   getFreePars(object=object@parStruct,...)
 }
)

setMethod("setPars",signature(object="McplBaseModel"),
  function(object,pars,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    if(rval == "object") {
      object@parStruct <- setPars(object=object@parStruct,pars=pars,...,rval=rval)
      return(object)
    } else {
      return(setPars(object=object@parStruct,pars=pars,...,rval=rval))
    }
  }
)

setMethod("setFreePars",signature(object="McplBaseModel",pars="numeric"),
  function(object,pars,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    if(rval == "object") {
      object@parStruct <- setFreePars(object=object@parStruct,pars=pars,...,rval=rval)
      return(object)
    } else {
      return(setFreePars(object=object@parStruct,pars=pars,...,rval=rval))
    }
  }
)

setMethod("getConstraints",signature(object="McplBaseModel"),
  function(object,...) {
    getConstraints(object@parStruct,...) 
  }
)

setMethod("logLik",signature=(object="McplBaseModel"),
  function(object,discount=0,...) {
    if(discount > 0) {
      discount <- as.numeric(mapply(seq,from=object@nTimes@bt,to=object@nTimes@bt-1 + discount))
      out <- sum(log(dens(object)[-discount]))
      nobs <- sum(object@nTimes@n) - length(discount)
    } else {
      out <- sum(log(dens(object)))
      nobs <- sum(object@nTimes@n)
    }
    attr(out,"nobs") <- nobs
    attr(out,"df") <- length(getFreePars(object,...))
    class(out) <- "logLik"
    out   
  }
)

setMethod("AIC",signature(object="McplBaseModel"),
  function(object,npar,...,k=2) {
    logL <- logLik(object,...)
    if(missing(npar)) npar <- length(getFreePars(object,...))#return(AIC(logL,...,k=k)) else {
    return(-2*logL + k*npar)
  }
)

setMethod("AICc",signature(object="McplBaseModel"),
  function(object,npar,...) {
    logL <- logLik(object,...)
    if(missing(npar)) npar <- length(getFreePars(object,...))
    nobs <- attr(logL,"nobs")
    if(is.null(nobs)) nobs <- nrow(object@y)
    -2*logL + (2*nobs*npar)/(nobs-npar-1)
  }
)
#setMethod("BIC",signature(object="McplBaseModel"),
#  function(object,npar,...) {
#    logL <- logLik(object)
#    if(missing(npar)) npar <- length(getPars(object,which="free",...))
#    nobs <- attr(logL,"nobs")
#    if(is.null(nobs)) nobs <- nrow(object@y)
#    -2*logL + npar*log(nobs)
#  }
#)
setMethod("Rsq",signature(object="McplBaseModel"),
  function(object,...) {
    p <- predict(object,type="response",...)
    SSt <- sum((object@y - mean(object@y))^2)
    SSe <- sum((object@y - p)^2)
    1-SSe/SSt
  }
)

setMethod("getReplication",signature(object="McplBaseModel"),
  function(object,case,ntimes=NULL,...) {
    if(is.null(ntimes)) {
      lt <- object@nTimes@cases
      bt <- object@nTimes@bt
      et <- object@nTimes@et
    } else {
      lt <- length(ntimes)
      et <- cumsum(ntimes)
      bt <- c(1,et[-lt]+1)
    }
    x <- object@x[bt[case]:et[case],]
    y <- object@y[bt[case]:et[case],]
    if(!is.matrix(x)) x <- as.matrix(x)
    if(!is.matrix(y)) y <- as.matrix(y)
    #if(!object@parStruct@replicate) parameters <- object@parStruct@parameters[[case]] else parameters <- object@parStruct@parameters
    parameters <- getPars(object@parStruct,replication=case)
    return(list(x=x,y=y,parameters=parameters))
  }
)

setMethod("show",signature(object="McplBaseModel"),
  function(object) {
    print(object)
  }       
)

setMethod("print",signature(x="McplBaseModel"),
  function(x,...) {
    cat("Parameters:\n")
    print(x@parStruct)
  }
)

setMethod("summary",signature(object="McplBaseModel"),
  function(object,fits=TRUE,...) {
    #cat("Parameters:\n")
    print(object@parStruct)
    #cat("\n")
#     x <- object
#     #cat("ResponseModel, class:",is(x),"\n\n")
#     pars <- getPars(x,which="all",...)
#     attr(pars,"skeleton") <- NULL
#     fx <- x@parStruct@fix
#     if(sum(fx)>0) {
#     #if(length(fx)>0) {
#       frpars <- pars[!fx]
#       fxpars <- pars[fx]
#     } else {
#       frpars <- pars
#       fxpars <- NULL
#     }
#     if(length(frpars)>0) {
#       cat("Estimated parameters:\n")
#       print(frpars)
#       cat("\n")
#     }
#     if(length(fxpars>0)) {
#       cat("Fixed parameters:\n")
#       fpars <- getPars(x,which="all",...)
#       print(fpars[x@parStruct@fix])
#       cat("\n")
#     }
    if(fits) {
      mf <- list()
      try({
        mf$logLik=logLik(x,...)
        mf$AIC=AIC(x,...)
        mf$BIC=BIC(x,...)
        #mf$AICc=AICc(x,...)
      },silent=TRUE)
      if(length(mf>0)) {
        cat("Model fit:\n")
        print(unlist(mf))
        cat("\n")
      }
    }
  }
)

setMethod("fit",signature(object="McplBaseModel"),
  function(object,method="Nelder-Mead",...) {
    MLoptfun <- function(pars,object,...) {
      object <- setFreePars(object,pars,...,rval="object")
      object <- runm(object,...)
      -logLik(object)
    }
    LSoptfun <- function(pars,object,...) {
      object <- setFreePars(object,pars,...,rval="object")
      object <- runm(object,...)
      sum((predict(object,type="response",...)-object@y)^2)
    }
    pars <- getFreePars(object,...)
    if(hasMethod("logLik",is(object))) optfun <- MLoptfun else optfun <- LSoptfun
    constraints <- getConstraints(object,...)
    if(is(constraints,"Unconstrained")) {
      if(length(pars) == 1 & method=="Nelder-Mead") method <- "BFGS"
      opt <- optim(par=pars,fn=optfun,object=object,method=method,...)  
    } else if(is(constraints,"BoxContraints")) {
      opt <- optim(par=pars,fn=optfun,object=object,method="L-BFGS-B",min=constraints@min,max=constraints@max,...)  
    } else if(is(constraints,"LinearConstraints")) {
      opt <- constrOptim(theta=pars,f=optfun,grad=NULL,ui=constraints@Amat,ci=constraints@bvec,...)
    } else {
      stop("Cannot determine constraints of this MCPL model")
    }
    object <- setFreePars(object,opt$par,...,rval="object")
    object <- runm(object,...)
    object
  }
)

setMethod("has.runm",signature(object="McplBaseModel"),
  function(object,...) {
    FALSE
  }
)

setMethod("has.lFr",signature(object="McplBaseModel"),
  function(object,...) {
    FALSE
  }
)

setMethod("has.rFl",signature(object="McplBaseModel"),
  function(object,...) {
    FALSE
  }
)
