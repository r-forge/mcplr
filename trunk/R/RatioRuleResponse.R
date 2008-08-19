setClass("RatioRuleResponse",
  contains="ResponseModel",
  representation(
    transformation = "function"
  )
)

setMethod("estimate",signature(object="RatioRuleResponse"),
  function(object,...) {
    optfun <- function(par,object,...) {
      object@parameters <- setPars(object,par,rval="parameters",...)
      -sum(logLik(object,...))
    }
    mf <- match.call()
    pstart <- unlist(object@parameters)
    #if(length(pstart)!=1) stop("Ratio Rule response must have a single parameter")
    if(!is.null(mf$CML.method) || !is.null(mf$method)) {
      if(!is.null(mf$CML.method)) mf$method <- mf$CML.method
      mf$par <- pstart
      mf$fn <- optfun
      mf$object <- object
      mf[[1]] <- as.name("optim")
      opt <- eval(mf,parent.frame())
      #opt <- optim(log(pstart),fn=optfun,object=object,...) else
    } else {
      opt <- optim(pstart,fn=optfun,object=object,...)
    }
    #object <- setPars(object,exp(opt$par))
    object@parameters <- setPars(object,opt$par,rval="parameters",...)
    object <- fit(object,...)
    object
  }
)

setMethod("predict",signature(object="RatioRuleResponse"),
  function(object,...) {
    #beta <- object@parameters$beta
    out <- object@transformation(object,...)
    out <- out/rowSums(out)
    #out <- apply(object@x,1,function(x) exp(x)/sum(exp(x)))
    if(!is.matrix(out)) out <- matrix(out,ncol=1)
    out
  }
)

setMethod("logLik",signature(object="RatioRuleResponse"),
  function(object,eps=.Machine$double.eps,...) {
    pred <- predict(object,type="response",...)
    nt <- NROW(pred)
    if(ncol(as.matrix(pred))==1) {
      pred <- rowSums(cbind(pred*object@y,(1-pred)*(1-object@y)))
    } else {
      pred <- rowSums(object@y*pred)
    }
    pred[pred > 1-eps] <- 1-eps
    pred[pred < eps] <- eps
    #if(!is.null(discount)) {
    #  discount <-  as.numeric(sapply(object@nTImes@bt,"+",discount))
    #  LL <- sum(log(pred[-discount]))
    #  attr(LL,"nobs") <- length(pred)-length(discount)
    #} else {
      LL <- sum(log(pred))
    #}
    attr(LL,"nobs") <- nt
    attr(LL,"df") <- length(getPars(object,which="free"))
    class(LL) <- "logLik"
    LL
  }
)

RatioRuleResponse.trans.exp <- function(object,...) {
  exp(object@parameters$beta*object@x)
}

RatioRuleResponse.trans.none <- function(object,...) {
  object@x
}

RatioRuleResponse <- function(formula,parameters=list(beta=1),transformation=c("exponential","none"),
                        data,ntimes=NULL,replicate=TRUE,fixed,
                        parStruct,subset) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset) else 
    dat <- mcpl.prepare(formula,data)

  y <- dat$y
  x <- dat$x
  
  if(!is.function(transformation)) {
    transformation <- match.arg(transformation)
    trans <- switch(transformation,
      exponential = RatioRuleResponse.trans.exp,
      none = RatioRuleResponse.trans.none)
  } else {
    trans <- transformation
  }
  
  parfill <- function(parameters) {
    #pars <- list()
    if(!is.list(parameters)) parameters <- as.list(parameters)
    if(is.null(parameters$beta)) parameters$beta <- 1
    parameters
  }

  if(is.null(ntimes) | replicate) {
    # intialize
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
    tfix <- NULL
    if(!missing(fixed)) tfix <- fixed
    parStruct <- ParStruct(parameters,replicate=replicate,
                    fixed=tfix,ntimes=ntimes)
  }

  if(is.null(ntimes)) ntimes <- nrow(y)
  nTimes <- nTimes(ntimes)
    
  mod <- new("RatioRuleResponse",
    x = x,
    y = y,
    parameters = parameters,
    parStruct=parStruct,
    nTimes=nTimes,
    transformation=trans)
  mod <- fit(mod)
  mod                     
                        
}

