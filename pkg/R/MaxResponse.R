setClass("MaxResponse",
  contains="ResponseModel"
)

setMethod("predict",signature(object="MaxResponse"),
  function(object,...) {
    #beta <- object@parStruct@parameters$beta
    #out <- object@transformation(object,...)
    p <- 1/(1+exp(-object@parStruct@parameters$beta))
    if(length(beta) > 1) {
      N <- object@parStruct@n
      if(length(p) == length(n)) p <- rep(p,each=n) else stop("parameter beta should have length 1 or number of replications") 
    }
    out <- (object@x==apply(object@x,1,max))*(1-p) + (1-(object@x==apply(object@x,1,max)))*p
    #out <- out/rowSums(out)
    #out <- apply(object@x,1,function(x) exp(x)/sum(exp(x)))
    if(!is.matrix(out)) out <- matrix(out,ncol=1)
    out
  }
)

setMethod("logLik",signature(object="MaxResponse"),
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
    LL <- sum(log(pred))
    attr(LL,"nobs") <- nt
    attr(LL,"df") <- length(getPars(object,which="free"))
    class(LL) <- "logLik"
    LL
  }
)

MaxResponse <- function(formula,data,parameters=list(beta=2.944439),
                        ntimes=NULL,replicate=TRUE,fixed,
                        parStruct,subset) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset) else dat <- mcpl.prepare(formula,data)

  y <- dat$y
  x <- dat$x
  
  
  parfill <- function(parameters) {
    #pars <- list()
    if(!is.list(parameters)) parameters <- as.list(parameters)
    if(is.null(parameters$beta)) parameters$beta <- 2.944439
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
    
  mod <- new("MaxResponse",
    x = x,
    y = y,
    parStruct=parStruct,
    nTimes=nTimes)
  #mod <- runm(mod)
  mod                   
}

