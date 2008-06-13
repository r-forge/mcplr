setClass("McplModelList",
  representation(
    mcplModels="list"
  )
)
mcplModelList <- function(list) {
  new("McplModelList",
    mcplModels=list)
}
setMethod("logLik",signature(object="McplModelList"),
  function(object,...) {
    LL <- vector("numeric")
    for(i in 1:length(object@mcplModels)) {
      LL[i] <- logLik(object@mcplModels[[i]],...)
    }
    sum(LL)
  }
)
setMethod("getPars",signature(object="McplModelList"),
  function(object,which="all",...) {
    pars <- list()
    for(i in 1:length(object@mcplModels)) {
      pars[[i]] <- getPars(object@mcplModels[[i]],which=which,...)
    }
    pars <- as.relistable(pars)
    unlist(pars)
  }
)
setMethod("predict",signature(object="McplModelList"),
  function(object,type="link",...) {
    p <- predict(object@mcplModels[[1]],type=type,...)
    for(i in 2:length(object@mcplModels)) {
      if(is.matrix(p)) p <- rbind(p,predict(object@mcplModels[[i]],type=type,...)) else p <- c(p,predict(object@mcplModels[[i]],type=type,...))
    }
    p
  }
)


setMethod("AIC",signature(object="McplModelList"),
  function(object,npar,...,k=2) {
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    nobs <- 0
    for(i in 1:length(object@mcplModels)) {
      nobs <- nobs + nrow(object@mcplModels[[i]]@responseModel@y)
    }
    logL <- logLik(object,...)
    AIC(logL=logL,npar=npar,nobs=nobs,...)
  }
)

setMethod("AICc",signature(object="McplModelList"),
  function(object,npar,...) {
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    nobs <- 0
    for(i in 1:length(object@mcplModels)) {
      nobs <- nobs + nrow(object@mcplModels[[i]]@responseModel@y)
    }
    logL <- logLik(object,...)
    AICc(logL=logL,npar=npar,nobs=nobs,...)
  }
)
setMethod("BIC",signature(object="McplModelList"),
  function(object,npar,...) {
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    nobs <- 0
    for(i in 1:length(object)) {
      nobs <- nobs + nrow(object@mcplModels[[i]]@responseModel@y)
    }
    logL <- logLik(object,...)
    BIC(logL=logL,npar=npar,nobs=nobs,...)
  }
)
setMethod("RSquare",signature(object="McplModelList"),
  function(object,...) {
    warning("RSquare is based on responseModel")
    p <- as.matrix(predict(object@mcplModels[[1]]@responseModel,type="response"))
    y <- object@mcplModels[[1]]@responseModel@y
    for(i in 2:length(object@mcplModels)) {
      p <- rbind(p,predict(object@mcplModels[[i]]@responseModel,type="response"))
      y <- rbind(y,object@mcplModels[[i]]@responseModel@y)
    }
    SSt <- sum((y - mean(y))^2)
    SSe <- sum((y - p)^2)
    1-SSe/SSt
  }
)