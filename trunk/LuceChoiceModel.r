setClass("LuceChoiceModel",
  contains="ResponseModel",
  representation(
    family="ANY"
  )
)
setMethod("estimate",signature(object="LuceChoiceModel"),
  function(object,...) {
    mf <- match.call()
    pstart <- unlist(object@parameters)
    if(length(pstart)!=1) stop("Luce Choice Rule must have a single parameter")
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
        beta <- exp(par)
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
setMethod("predict",signature(object="LuceChoiceModel"),
  function(object,...) {
    beta <- object@parameters$beta
    #out <- object@family$linkinv(beta*object@x)
    out <- t(apply(object@x,1,function(x) exp(x)/sum(exp(x))))
    if(!is.matrix(out)) out <- matrix(out,ncol=1)
    out
  }
)

setMethod("logLik",signature(object="LuceChoiceModel"),
  function(object,eps=.Machine$double.eps,...) {
    pred <- predict(object,type="response",...)
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
    LL
  }
)
