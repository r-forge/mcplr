setClass("ParStruct",
  representation(
    parameters="list",
    id="integer",
    fix="logical",
    constraints="ConstraintsList",
    replicate="logical"
  )
)

ParStruct <- function(parameters,replicate=TRUE,fixed=NULL,ntimes=NULL,
                constraints=NULL) {         
  parStruct <- new("ParStruct")
  parStruct@parameters <- parameters
  if(replicate) parStruct@replicate <- TRUE else parStruct@replicate <- FALSE

  if(!is.null(fixed)) {
    if(is.null(ntimes) | length(ntimes)==1 | replicate) {
      if(is.list(fixed)) fix <- fixedListToVec(fixed,parameters) else fix <- fixed
    } else {
      cases <- length(ntimes)
      fix <- vector("list",length=cases)
      if(is.list(fixed)) {
        if(is.null(names(fixed)) && length(fixed) == cases) {
          for(i in 1:cases) {
            if(is.list(fixed[[i]])) {
              fix[[i]] <- fixedListToVec(fixed[[i]],parameters)
            } else {
              fix[[i]] <- fixed[[i]]
            }
          }
        } else {
           for(i in 1:cases) {
            fix[[i]] <- fixedListToVec(fixed,parameters[[i]])
          }
        }
      } else {
        if(length(fixed) == length(unlist(parameters))) {
          fix <- fixed
        } else {
          fix <- rep(fixed,length=length(unlist(parameters)))
        }
      }
    }
  } else {
    fix <- rep(FALSE,length=length(unlist(parameters)))
  }
  parStruct@fix <- unlist(fix)
  if(!is.null(constraints)) parStruct@constraints <- constraints
  parStruct
}

setMethod("getPars",signature(object="ParStruct"),
  # Note:
  # if one element in parameters corresponding to id[i] is free, so
  #   is this parameter!
  function(object,which="all",internal=FALSE,...) {
    pars <- unlist(object@parameters)
    switch(which,
      free = {
        parid <- object@id
        fix <- object@fix
        if(length(parid)>0) {
          # get unique parameters
          if(length(fix)>0) parid.v <- unique(parid[!fix]) else parid.v <- unique(parid)
          parid.n <- length(parid.v)
          newpars <- vector()
          for(i in 1:parid.n) newpars[i] <- pars[which(parid==parid.v[i])[1]]
          pars <- newpars
        } else {
          if(length(fix)>0) {
            if(length(fix) == length(pars)) pars <- pars[!fix] else stop("length of fix in parameterList does not match number of parameters")
          }
        }
        if(length(pars)==0) pars <- NULL
        pars
      },
      all = pars,
      pars
    )
    return(pars)
  }
)
setMethod("setPars",signature(object="ParStruct"),
  function(object,pars,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    oldpars <- unlist(object@parameters)
    if(length(pars) > 0) {
      if(length(pars)!=length(oldpars)) {
        parid <- object@id
        fix <- object@fix
        if(length(parid)>0) {
          if(length(fix)>0) parid.v <- unique(parid[!fix]) else parid.v <- unique(parid)
          parid.n <- length(parid.v)
          if(length(pars)!=parid.n) stop("length of parid does not match length of pars")
          for(i in 1:parid.n) oldpars[which(parid==parid.v[i])] <- pars[i]
        } else {
          if(length(fix)>0) {
            if(sum(!fix)!=length(pars)) stop("parid not given and length of pars does not equal number of nonfixed parameters")
            oldpars[!fix] <- pars
          } else {
            stop("cannot work with par of this length in setPar")
          }
        }
      } else {
        oldpars <- pars
      }
      object@parameters <- relist(oldpars,skeleton=object@parameters)
    }
    switch(rval,
      object = object,
      parameters = object@parameters)
    #object <- fit(object)
    #return(object)
    #newpars <- relist(oldpars,skeleton=object@parameters)
    #newpars
  }
)

setMethod("rep",signature(x="ParStruct"),
  function(x,times=1) {
    if(x@replicate) return(x)
    x@parameters <- rep(x@parameters,times)
    x@id <- rep(x@id,times)
    x@fix <- rep(x@fix,times)
    x@constraintsList <- rep(x@constraintsList,times=times)
    return(x)
  }
)
