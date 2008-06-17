# GLM Response Model
### A response model which can be fitted by glm.fit
setClass("GlmResponse",
  contains="ResponseModel",
  representation(
    family="ANY"
    #sigma="matrix"
  )
)
setMethod("estimate",signature(object="GlmResponse"),
	function(object) {
    # y = response
    pars <- object@parameters
    fit <- glm.fit(x=object@x,y=object@y,family=object@family,start=pars$coefficients)
    pars$coefficients <- fit$coefficients
    if(object@family$family=="gaussian") pars$sd <- sqrt(sum((object@y - predict(object))^2)/(length(object@y)-1))
    #object <- setpars(object,unlist(pars))
    object@parameters <- setpars(object,unlist(pars),rval="parameters",...)
    object <- fit(object,...)
    return(object)
	}
)
setMethod("predict","GlmResponse",
	function(object,type="link") {
    # y = response
    mu <- crossprod(object@x,object@parameters$coefficients)
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
        dnorm(x=object@y,mean=mu,sd=object@parameters$sd,log=TRUE)
      },
      binomial = {
        p <- predict(object,type="response")
        dbinom(x=object@y,size=1,prob=p,log=TRUE)
      },
      poisson = {
        lambda <- predict(object,type="response")
        dpois(x=object@y,lambda=lambda,log=TRUE)
      },
      Gamma = {
        shape <- predict(object,type="response")
        dgamma(x=object@y,shape=shape,log=TRUE)
      },
      stop("family",object@family$family,"not implemented (yet)")
    )
  }
)