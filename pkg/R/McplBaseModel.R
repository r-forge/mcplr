setClass("McplBaseModel",
  representation(
    x="matrix",
    y="matrix",
    parStruct="ParStruct",
    nTimes="NTimes"
  )
)
setMethod("getPars",signature(object="McplBaseModel"),
  # Note:
  # if one element in parameters corresponding to id[i] is free, so
  #   is this parameter!
  function(object,...) {
    getPars(object=object@parStruct,...)
  }
)

setMethod("setPars",signature(object="McplBaseModel"),
  function(object,pars,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    if(rval == "object") {
      object@parStruct@parameters <- setPars(object=object@parStruct,pars=pars,...,rval="parameters")
      return(object)
    } else {
      return(setPars(object=object@parStruct,pars=pars,...,rval="parameters"))
    }
  }
)
setMethod("AIC",signature(object="McplBaseModel"),
  function(object,npar,...,k=2) {
    logL <- logLik(object,...)
    if(missing(npar)) return(AIC(logL,...,k=k)) else {
      nobs <- attr(logL,"nobs")
      if(is.null(nobs)) nobs <- nrow(object@y)
      return(-2*logL + k*npar)
    }
  }
)

setMethod("AICc",signature(object="McplBaseModel"),
  function(object,npar,...) {
    logL <- logLik(object,...)
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
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
    if(!object@parStruct@replicate) parameters <- object@parStruct@parameters[[case]] else parameters <- object@parStruct@parameters
    return(list(x=x,y=y,parameters=parameters))
  }
)
setMethod("show",signature(object="McplBaseModel"),
  function(object) {

  }
)

setMethod("summary",signature(object="McplBaseModel"),
  function(object,fits=TRUE,...) {
    x <- object
    #cat("ResponseModel, class:",is(x),"\n\n")
    pars <- getPars(x,which="all",...)
    attr(pars,"skeleton") <- NULL
    fx <- x@parStruct@fix
    if(sum(fx)>0) {
    #if(length(fx)>0) {
      frpars <- pars[!fx]
      fxpars <- pars[fx]
    } else {
      frpars <- pars
      fxpars <- NULL
    }
    if(length(frpars)>0) {
      cat("Estimated parameters:\n")
      print(frpars)
      cat("\n")
    }
    if(length(fxpars>0)) {
      cat("Fixed parameters:\n")
      fpars <- getPars(x,which="all",...)
      print(fpars[x@parStruct@fix])
      cat("\n")
    }
    if(fits) {
      mf <- list()
      try({
        mf$logLik=logLik(x,...)
        mf$BIC=BIC(x,...)
        mf$AICc=AICc(x,...)
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
      object@parStruct@parameters <- setPars(object,pars,...,rval="parameters",internal=TRUE)
      object <- runm(object,...)
      -logLik(object)
    }
    LSoptfun <- function(pars,object,...) {
      object@parStruct@parameters <- setPars(object,pars,...,rval="parameters",internal=TRUE)
      object <- runm(object,...)
      sum((predict(object,type="response",...)-object@y)^2)
    }
    pars <- getPars(object,which="free",...,internal=TRUE)
    if(hasMethod("logLik",is(object))) optfun <- MLoptfun else optfun <- LSoptfun
    if(!is.null(object@parStruct@constraints)) {
      switch(is(object@parStruct@constraints),
        "LinConstraintsList" = {
          A <- object@parStruct@constraints@Amat
          b <- object@parStruct@constraints@bvec
          opt <- constrOptim(theta=pars,f=optfun,grad=NULL,ui=A,ci=b,object=object,...)
          object@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
        },
        "BoxConstraintsList" = {
          opt <- optim(pars,fn=optfun,method="L-BFGS-B",object=object,min=object@parStruct@constraints@min,max=object@parStruct@constraints@min,...)
          object@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
        },
        {
          warning("This extension of ConstraintsList is not implemented; using default optimisation.")
          opt <- optim(pars,fn=optfun,method=method,object=object,...)
          object@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
        }
       )
    } else {
      opt <- optim(pars,fn=optfun,method=method,object=object,...)
      object@parStruct@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
    }
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

setMethod("canRepar",signature(object="McplBaseModel"),
  function(object,...) {
    FALSE
  }
)
