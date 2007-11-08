setClass("GaussianModel",
  contains="ResponseModel"
)
setMethod("estimate",signature(object="GaussianModel"),
	function(object) {
    pars <- object@parameters
    pars$sd <- sd(object@x-object@y)
    object@parameters <- setPars(object,unlist(pars),rval="parameters",...)
    object <- fit(object,...)
    return(object)
	}
)
setMethod("predict","GaussianModel",
	function(object,type="link") {
    object@x
  }
)
setMethod("logLik","GaussianModel",
  function(object) {
    mu <- predict(object,type=response)
    sum(dnorm(x=object@y,mean=mu,sd=object@parameters$sd,log=TRUE))
  }
)

gaussianModel <- function(formula,parameters=list(sd=1),data,subset,ntimes=NULL) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,remove.intercept=TRUE) else dat <- mcpl.prepare(formula,data,remove.intercept=TRUE)
  y <- dat$y
  x <- dat$x
  if(is.null(ntimes)) ntimes <- nrow(y)
  nTimes <- nTimes(ntimes)
  new("GaussianModel",
    y = y,
    x = x,
    parameters=parameters,
    nTimes=nTimes)
}