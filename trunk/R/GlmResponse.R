# GLM Response Model
### A response model which can be fitted by glm.fit
setClass("GlmResponse",
  contains="ResponseModel",
  representation(
    family="ANY"
    #sigma="matrix"
  )
)


GlmResponse <- function(formula,data,family=gaussian(),parameters=list(),ntimes=NULL,replicate=TRUE,fixed,base=NULL,
                        parStruct,subset) {
  if(!missing(subset)) dat <- mcpl.prepare(formula=formula,data=data,subset=subset,base=base) else
    dat <- mcpl.prepare(formula=formula,data=data,base=base)
  y <- dat$y
  x <- dat$x
  
  parfill <- function(parameters) {
    if(family$family=="gaussian") if(is.null(parameters$sd)) parameters$sd <- 1 else parameters$sd <- parameters$sd
    if(is.null(parameters$coefficients)) parameters$coefficients <- rep(0,NCOL(x))
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

  mod <- new("GlmResponse",
    x = x,
    y = y,
    parameters = parameters,
    parStruct=parStruct,
    nTimes=nTimes,
    family=family)
  mod <- fit(mod)
  mod
}

setMethod("estimate",signature(object="GlmResponse"),
	function(object,...) {
    # y = response
    pars <- object@parameters
    fit <- glm.fit(x=object@x,y=object@y,family=object@family,start=pars$coefficients)
    pars$coefficients <- fit$coefficients
    if(object@family$family=="gaussian") pars$sd <- sqrt(sum((object@y - predict(object))^2)/(length(object@y)-1))
    #object <- setpars(object,unlist(pars))
    object@parameters <- setPars(object,unlist(pars),rval="parameters",...)
    object <- fit(object,...)
    return(object)
	}
)
setMethod("predict","GlmResponse",
	function(object,type="link") {
    # y = response
    if(NCOL(object@x) > 1) mu <- object@x%*%as.matrix(object@parameters$coefficients) else mu <- object@parameters$coefficients*object@x
    if(type=="link") return(mu) else {
      if(type=="response") {
        return(object@family$linkinv(mu))
      }
    }
	}
)
setMethod("logLik","GlmResponse",
  function(object) {
    switch(object@family$family,
      gaussian = {
        mu <- predict(object,type=response)
        sum(dnorm(x=object@y,mean=mu,sd=object@parameters$sd,log=TRUE))
      },
      binomial = {
        p <- predict(object,type="response")
        sum(dbinom(x=object@y,size=1,prob=p,log=TRUE))
      },
      poisson = {
        lambda <- predict(object,type="response")
        sum(dpois(x=object@y,lambda=lambda,log=TRUE))
      },
      Gamma = {
        shape <- predict(object,type="response")
        sum(dgamma(x=object@y,shape=shape,log=TRUE))
      },
      stop("family",object@family$family,"not implemented (yet)")
    )
  }
)