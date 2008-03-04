# GLM Response Model
### A response model which can be fitted by glm.fit
setClass("GlmResponseModel",
  contains="ResponseModel",
  representation(
    family="ANY",
    sigma="matrix"
  )
)
setMethod("estimate",signature(object="GlmResponseModel"),
	function(object) {
    # y = response
    pars <- object@parameters
    fit <- glm.fit(x=object@x,y=object@y,family=object@family)
    pars$coefficients <- fit$coefficients
    #object <- setpars(object,unlist(pars))
    object@parameters <- setpars(object,unlist(pars),rval="parameters",...)
    object <- fit(object,...)
    return(object)
	}
)
setMethod("predict","GlmResponseModel",
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
setMethod("logLik","GlmResponseModel",
  function(object) {
    switch(object@family$family,
      gaussian = {
        mu <- predict(object,type=response)
        dnorm(x=object@y,mean=mu,sd=object@sigma,log=TRUE)
      },
      binomial = {
        p <- predict(object,type="response")
        dbinom(x=object@y,size=1,prob=p,log=TRUE)
      },
      stop("family",object@family$family,"not implemented (yet)")
    )
  }
)