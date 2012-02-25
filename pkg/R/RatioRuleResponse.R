setClass("RatioRuleResponse",
  contains="ResponseModel",
  representation(
    transformation = "function"
  )
)

setMethod("fit",signature(object="RatioRuleResponse"),
  function(object,...) {
    optfun <- function(par,object,...) {
      object@parStruct@parameters <- setPars(object,par,rval="parameters",...)
      -sum(logLik(object,...))
    }
    mf <- match.call()
    pstart <- unlist(object@parStruct@parameters)
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
    object@parStruct@parameters <- setPars(object,opt$par,rval="parameters",...)
    object <- runm(object,...)
    object
  }
)

setMethod("predict",signature(object="RatioRuleResponse"),
  function(object,...) {
    #beta <- object@parStruct@parameters$beta
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

setMethod("simulate",signature(object="RatioRuleResponse"),
	function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) set.seed(seed)
    if(missing(times)) {
      pr <- predict(object,type=response)
    } else {
      pr <- predict(object,type=response)[times,]
    }
    nt <- nrow(pr)
    response <- t(apply(pr,1,rmultinom,n=1,size=1))
    
		#if(nsim > 1) response <- matrix(response,ncol=nsim)
		#response <- as.matrix(response)
		
		object@y <- as.matrix(response)
		if(!missing(times)) {
      object@x <- object@x[rep(times,nsim),]
      ntim <- rep(0,length=object@nTimes@cases)
  		for(i in 1:length(ntim)) {
  		  ntim[i] <- sum(seq(object@nTimes@bt[i],object@nTimes@et[i]) %in% times)
      }
      warning("simulation with a times argument may result in wrong parStruct argument; please check parameters.")
      object@parStruct <- rep(object@parStruct,times=nsim)
    } else {
      object@x <- object@x[rep(1:nrow(object@x),nsim),]
      ntim <- object@nTimes@n
      object@parStruct <- rep(object@parStruct,times=nsim)
    }
    ntim <- rep(ntim,nsim)
    object@nTimes <- nTimes(ntim)
    if(!is.null(seed)) {
      set.seed(old.seed)
    }
		return(object)
	}
)

RatioRuleResponse.trans.exp <- function(object,...) {
  if(object@parStruct@replicate) {
    return(exp(object@parStruct@parameters$beta*object@x))
  } else {
    beta <- getPars(object,"beta")
    beta <- rep(beta,each=object@nTimes@n)
    return(exp(beta*object@x))
  }
}

RatioRuleResponse.trans.pow <- function(object,...) {
    if(object@parStruct@replicate) {
      return(object@x^object@parStruct@parameters$beta)
    } else {
      beta <- getPars(object,"beta")
      beta <- rep(beta,each=object@nTimes@n)
      return(object@x^beta)
    }
}

RatioRuleResponse.trans.none <- function(object,...) {
  object@x
}

RatioRuleResponse <- function(formula,parameters=list(beta=1),transformation=c("exponential","power","none"),
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
      power = RatioRuleResonse.trans.pow,
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
    parStruct=parStruct,
    nTimes=nTimes,
    transformation=trans)
  mod <- runm(mod)
  mod                     
                        
}

