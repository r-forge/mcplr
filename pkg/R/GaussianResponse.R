setClass("GaussianResponse",
  contains="ResponseModel"
)
setMethod("estimate",signature(object="GaussianResponse"),
	function(object) {
    pars <- object@parameters
    pars$sd <- sd(object@x-object@y)
    object@parameters <- setPars(object,unlist(pars),rval="parameters",...)
    object <- fit(object,...)
    return(object)
	}
)
setMethod("predict","GaussianResponse",
	function(object,type="link") {
    object@x
  }
)
setMethod("logLik","GaussianResponse",
  function(object) {
    mu <- predict(object,type=response)
    sum(dnorm(x=object@y,mean=mu,sd=object@parameters$sd,log=TRUE))
    #sum(pnorm(x=object@y+.5,mean=mu,sd=object@parameters$sd) - pnorm(x=object@y+.5,mean=mu,sd=object@parameters$sd))
  }
)

gaussianResponse <- function(formula,parameters=list(sd=1),data,subset,ntimes=NULL) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,remove.intercept=TRUE) else dat <- mcpl.prepare(formula,data,remove.intercept=TRUE)
  y <- dat$y
  x <- dat$x
  if(is.null(ntimes)) ntimes <- nrow(y)
  nTimes <- nTimes(ntimes)
  new("GaussianResponse",
    y = y,
    x = x,
    parameters=parameters,
    nTimes=nTimes)
}