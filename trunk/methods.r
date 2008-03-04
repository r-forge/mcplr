# Generics
setGeneric("nPar",function(object) standardGeneric("nPar"))
setGeneric("setPar",function(object,par,...) standardGeneric("setPar"))
setGeneric("getPar",function(object,...) standardGeneric("getPar"))
setGeneric("setX",function(object,psi,rho,...) standardGeneric("setX"))
setGeneric("fit",function(object,...) standardGeneric("fit"))

# npar

setMethod("nPar","lrModel",
	function(object) {
		return(length(unlist(object@parameters)))
	}
)

setMethod("nPar","mcplModel",
	function(object) {
    parid <- unique(unlist(object@parid))
		return(sum(parid!=0))
	}
)

# setpar

setMethod("setPar","lrModel",
  function(object,par) {
    object@parameters <- par
    object <- fit(object)
  }
)

relist <- function(vec,original,j=1) {
    # turns vector 'vec' into a list of identical structure to list 'original'
    # usage e.g.
    #   tmp <- unlist(list)
    #   tmp <- tmp + 3
    #   list <- relist(tmp,list)
    if(length(vec) != length(original)) stop("vec and original must have same length")
    x <- original
    if (length(x) == 0) return(x)
    if (is.list(x)) for(i in seq_along(x)) {
        rl <- relist(vec,x[[i]],j=j)
        x[[i]] <- rl
        attr(x[[i]],"j") <- NULL
        j <- attr(rl,"j")
    } else {
        dm <- dim(x)
        x <- vec[j:(j+(length(x)-1))]
        dim(x) <- dm
        j <- j+length(x)
        if(j < length(vec)) attr(x,"j") <- j
    }
    x
}

setMethod("setPar","mcplModel",
    # sets parameters of slots learnModel and respModel using the parid list
    # then re-fits both using new parameters
    function(object,par) {
        varpar <- unlist(rapply(object@parid, function(x) par[x],how="replace"))
        fixpar <- unlist(list(object@learnModel@parameters,object@respModel@parameters))
        allpar <- (varpar==0)*fixpar + varpar
        allpar <- relist(allpar,object@parid)
        object@learnModel <- setPar(object@learnModel,allpar[[1]])
        object@respModel <- setPsi(object@respModel,object@learnModel@psi)
        object@respModel <- setPar(object@respModel,allpar[[2]])
        object
    }
)

# getPar

setMethod("getPar","mcplModel",
  function(object) {
#    parid <- unlist(object@parid)
#    pars <- vector()
#    npar <- nPar(object)
#    for(i in 1:npar) {
#      pars[i] <-
#    unlist(object@parameters)
    object@parameters
  }
)

# setX
setMethod("setX","respModel",
  function(object,psi,rho) {
    object@x <- cbind(psi,rho)
    object
  }
)

# fit

setMethod("fit","respModel",
  function(object,...) {
    fit <- glm.fit(x=object@x,y=object@y,family=object@family,intercept=FALSE,...)
    object@par$coefficients <- fit$coefficients
    return(object)
  }
)

setMethod("fit","mcplModel",
  function(object,method="BFGS",criterion="logLik",...) {
    optfun <- function(par,obj,pidx,...) {
      obj <- setPar(obj,par) # this also re-fits the models
      switch(criterion,
        logLik = -logLik(obj),
        LS = SS(obj),
        -logLik(obj)
      )
    }
    pstart <- getPar(obj)
    fit <- optim(pstart,fn=optfun,method=method,...)
    object <- setpar(object,fit$par)
    object
  }
)