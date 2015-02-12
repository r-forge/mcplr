setClass("GaussianResponse",
  contains="ResponseModel"
)

setMethod("fit",signature(object="GaussianResponse"),
  function(object,...) {
    pars <- object@parStruct@parameters
    pars$sd <- sd(object@x-object@y)
    object@parStruct@parameters <- setPars(object,unlist(pars),rval="parameters",...)
    object <- runm(object,...)
    return(object)
  }
)
setMethod("predict","GaussianResponse",
	function(object,type="link") {
    object@x
  }
)

setMethod("dens","GaussianResponse",
  function(object,...) {
    sd <- getPars(object,...)$sd
    mu <- predict(object,type="response",...)
    dnorm(x=object@y,mean=mu,sd=sd)
    #sum(pnorm(x=object@y+.5,mean=mu,sd=object@parStruct@parameters$sd) - pnorm(x=object@y+.5,mean=mu,sd=object@parStruct@parameters$sd))
  }
)

# setMethod("logLik","GaussianResponse",
#   function(object) {
#     mu <- predict(object,type=response)
#     sum(dnorm(x=object@y,mean=mu,sd=object@parStruct@parameters$sd,log=TRUE))
#     #sum(pnorm(x=object@y+.5,mean=mu,sd=object@parStruct@parameters$sd) - pnorm(x=object@y+.5,mean=mu,sd=object@parStruct@parameters$sd))
#   }
# )

setMethod("simulate",signature(object="GaussianResponse"),
	function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) set.seed(seed)
    if(missing(times)) {
      pr <- predict(object,type=response)
    } else {
      pr <- predict(object,type=response)[times,]
    }
    nt <- nrow(pr)
    response <- rnorm(nt*nsim,mean=pr,sd=object@parStruct@parameters$sd)
    
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

GaussianResponse <- function(formula,parameters=list(sd=1),fixed,parStruct,data,subset,ntimes=NULL,replicate=TRUE) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,remove.intercept=TRUE) else dat <- mcpl.prepare(formula,data,remove.intercept=TRUE)
  y <- dat$y
  x <- dat$x
  
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$sd)) {
      parameters$sd <- 1
    } else {
      if(length(parameters$sd) > 1) {
        warning("sd should have length 1. Only first value will be used")
        parameters$sd <- parameters$sd[1]
      }
    }
  } else {
    # setup a parlist
    if(length(parameters)==0) {
      nrep <- length(ntimes)
      parameters$sd <- 1
      parameters <- rep(list(parameters),nrep)
    } else warning("there is no validity check for the given parameters when combined with ntimes and replicate=FALSE \n Please make sure the supplied list is valid")
  }
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
  if(missing(parStruct)) {
    parStruct <- ParStruct(parameters=parameters,replicate=replicate,
      fixed = if(missing(fixed)) NULL else fixed,
      ntimes = if(missing(ntimes)) NULL else ntimes)
  }
  
  new("GaussianResponse",
    y = y,
    x = x,
    parStruct=parStruct,
    nTimes=nTimes)
}
