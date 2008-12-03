setClass("RatioRuleResponse",
  contains="ResponseModel",
  representation(
    family="ANY"
  )
)
setMethod("estimate",signature(object="RatioRuleResponse"),
  function(object,...) {
    mf <- match.call()
    pstart <- unlist(object@parameters)
    if(length(pstart)!=1) stop("Ratio Rule response must have a single parameter")
    if(ncol(object@y) == 1) {
      # use glm.fit
      fit <- glm.fit(x=object@x,y=object@y,family=object@family,...)
      pars <- as.relistable(object@parameters)
      pars$beta <- fit$coefficients
      #object <- setPars(object,unlist(pars))
      object@parameters <- setPars(object,unlist(pars),rval="parameters",...)
      object <- fit(object,...)
    } else {
      optfun <- function(par,object,...) {
        beta <- exp(par[1])
        p <- predict(object,type="response")
        p <- p^beta
        p <- p/rowSums(p)
        -sum(log(rowSums(p*object@y)))
      }
      if(!is.null(mf$CML.method) || !is.null(mf$method)) {
        if(!is.null(mf$CML.method)) mf$method <- mf$CML.method
        mf$par <- log(pstart)
        mf$fn <- optfun
        mf$object <- object
        mf[[1]] <- as.name("optim")
        opt <- eval(mf,parent.frame())
        #opt <- optim(log(pstart),fn=optfun,object=object,...) else
      } else {
        opt <- optim(log(pstart),fn=optfun,method="BFGS",object=object,...)
      }
      #object <- setPars(object,exp(opt$par))
      object@parameters <- setPars(object,exp(opt$par),rval="parameters",...)
      object <- fit(object,...)
    }
    object
  }
)

setMethod("predict",signature(object="RatioRuleResponse"),
  function(object,...) {
    beta <- object@parameters$beta
    out <- object@family$linkinv(beta*object@x)
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

RatioRuleResponse <- function(formula,parameters=list(beta=1),
                        data,base=NULL,ntimes=NULL,replicate=TRUE,fixed,
                        parStruct,family,subset) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base) else 
    dat <- mcpl.prepare(formula,data,base=base)

  y <- dat$y
  
  if(!is.null(base)) {
    #y <- y[,base]
    if(missing(family)) if(ncol(y)==1) family <- binomial() else family <- multinomial()
  } else {
    if(missing(family)) if(ncol(y)==2) family <- binomial() else family <- multinomial()
  }
  x <- dat$x
  
  parfill <- function(parameters) {
    pars <- list()
    if(is.null(parameters$beta)) pars$beta <- 1 else pars$beta <- parameters$beta
    pars
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
    family=family)
  mod <- fit(mod)
  mod
}

