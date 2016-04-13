setClass("ParStruct",
  representation(
    parameters="numeric",
    skeleton="list",
    repeated="logical" # is this a list 
  )
)

setMethod("print",signature(x="ParStruct"),
  function(x,...) {
    cat("Parameters\n")
    print(unlist(getPars(x,...)))
    cat("\n")
  }
)

setMethod("show",signature(object="ParStruct"),
  function(object) {
    print(object)
  }
)
                  
setMethod("getPars",signature(object="ParStruct",name="ANY",replication="ANY"),
  function(object,name,replication,...) {
    tmp <- relist(object@parameters,skeleton=object@skeleton)
    if(object@repeated) {
      if(missing(name)) {
        if(missing(replication)) {
          return(tmp)
        } else {
          return(tmp[[replication]])
        }
      } else {
        if(missing(replication)) {
          # this is annoying; have to create a list with replications of name?
          out <- vector("list",length=length(tmp))
          for(i in 1:length(tmp)) {
            out[[i]] <- tmp[[i]][name]
          }
          return(out)
        } else {
          return(tmp[[replication]][name])
        }
      }
    } else {
      if(missing(name)) {
        return(tmp)
      } else if(is(name,"character")) {
        return(tmp[name])
      } else {
        stop("name argument should be character")
      }
    }
  }
)

setMethod("setPars",signature(object="ParStruct",pars="list"),
  # setPars is the recommended user function...
  function(object,pars,...,rval=c("object","parameters")) {
    # should check structure of list ...
    rval <- match.arg(rval)
    tmp <- relist(object@parameters,skeleton=object@skeleton)
    for(name in names(pars)) {
      tmp[[name]] <- pars[[name]]
    }
    object@parameters <- unlist(tmp) 
    switch(rval,
     object = object,
     parameters = object@parameters
    )
  }
)

setMethod("getFreePars",signature(object="ParStruct"),
  function(object,replication,...) {  
    return(object@parameters)
  }
)

## This function is called by fit; should be hidden from user
## always takes a numeric vector as pars
setMethod("setFreePars",signature(object="ParStruct",pars="numeric"),
  function(object,pars,...,rval=c("object","parameters")) { 
    rval <- match.arg(rval)
    oldpars <- object@parameters
    if(length(pars) > 0) {
      # all parameters are free..
      if(length(pars)!=length(oldpars)) stop("length of parid does not match length of pars")
      #oldpars <- pars
      object@parameters <- pars#relist(oldpars,skeleton=object@parameters)
    }
    switch(rval,
     object = object,
     parameters = object@parameters
    )
  }
)

setMethod("getConstraints",signature(object="ParStruct"),
  function(object,...) {
  return(new("Unconstrained",
    dim=length(object@parameters)))
  }
)  

setClass("SimpleConstraintsParStruct",
  contains="ParStruct",
  representation(
    id="integer", # This can be longer than parameters...
    fix="logical",
    min="numeric",
    max="numeric"
  )
)

setMethod("getPars",signature(object="SimpleConstraintsParStruct",name="ANY",replication="ANY"),
  function(object,name,replication,...) {
    object@parameters <- object@parameters[object@id]
    callNextMethod()
  }
)

setMethod("getConstraints",signature(object="SimpleConstraintsParStruct"),
  function(object,...) {
    if(length(unique(object@id)) == length(object@parameters)) {
      ident <- FALSE
    } else ident <- TRUE
    if(all(object@min == -Inf) & all(object@max == Inf)) {
      box <- FALSE
    } else box <- TRUE
    if(ident & box) {
      ## need to construct a linear constraints?
      return(new("BoxConstraints",
                 min=object@min[!object@fix],
                 max=object@max[!object@fix]))
      #stop("Not implemented yet")
    } else if(box) {
      # need to construct an identity constraints
      return(new("BoxConstraints",
                 min=object@min[!object@fix],
                 max=object@max[!object@fix]))
    } else {
      return(new("Unconstrained",
                 dim=length(object@parameters)))
    }
    

  }
)

setMethod("getFreePars",signature(object="SimpleConstraintsParStruct"),
  function(object,...) {
    pars <- object@parameters
    if(length(pars)==0) return(NULL)
    #parid <- object@id
    fix <- object@fix
    pars[!fix]
#     if(length(parid)>0) {
#      # get unique parameters
#      if(length(fix)>0) parid.v <- unique(parid[!fix]) else parid.v <- unique(parid)
#      parid.n <- length(parid.v)
#      newpars <- vector()
#      for(i in 1:parid.n) newpars[i] <- pars[which(parid==parid.v[i])[1]]
#      pars <- newpars
#     } else {
#      if(length(fix)>0) {
#        if(length(fix) == length(pars)) pars <- pars[!fix] else stop("length of fix in parameterList does not match number of parameters")
#      }
#     }
#     if(length(pars)==0) pars <- NULL
#     #pars
#     #)
#     return(pars)
   }     
)

setMethod("setFreePars",signature(object="SimpleConstraintsParStruct",pars="numeric"),
  function(object,pars,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    oldpars <- getFreePars(object,...)
    if(length(pars) == length(oldpars)) {
     object@parameters[!object@fix] <- pars
    } else {
      stop("number of parameters provided does not match the number of parameters in the object")
    }
#     
#     oldpars <- unlist(getPars(object,...))
#     if(length(pars) > 0) {
#       if(length(pars)!=length(oldpars)) {
#         #warning("This is currently buggy!")
#         parid <- object@id
#         fix <- object@fix
#         if(length(parid)>0) {
#           if(length(fix)>0) parid.v <- unique(parid[!fix]) else parid.v <- unique(parid)
#           parid.n <- length(parid.v)
#           if(length(pars)!=parid.n) stop("length of parid does not match length of pars")
#           for(i in 1:parid.n) oldpars[which(parid==parid.v[i])] <- pars[i]
#         } else {
#           if(length(fix)>0) {
#             if(sum(!fix)!=length(pars)) stop("parid not given and length of pars does not equal number of nonfixed parameters")
#             oldpars[!fix] <- pars
#           } else {
#             stop("cannot work with par of this length in setPar")
#           }
#         }
#       } else {
#         oldpars <- pars
#       }
#       object@parameters <- oldpars#relist(oldpars,skeleton=object@parameters)
#     }
    switch(rval,
       object = object,
       parameters = object@parameters
    )
  }
)

setMethod("print",signature(x="SimpleConstraintsParStruct"),
  function(x,...) {
    fix <- x@fix[x@id]
    if(sum(!fix) > 0) {
      cat("Free parameters\n")
      print(unlist(getPars(x,...)[!fix]))
      cat("\n")
    }
    if(sum(fix) > 0) {
    #cat("\n")
      cat("Fixed parameters\n")
      print(unlist(getPars(x,...)[fix]))
    }
    cat("\n")
  }
)

# setClass("RepSimpleConstraintsParStruct",
#   contains="SimpleConstraintsParStruct",
#   representation(
#     cases="integer"
#   )
# )

setClass("LinearConstraintsParStruct",
  contains="ParStruct",
  representation(
    Amat="matrix",
    bvec="numeric"
  )
)

ParStruct <- function(parameters,replicate=TRUE,fixed=NULL,min=NULL,max=NULL,id=NULL,ntimes=NULL) {    
  if(is.null(ntimes) | length(ntimes)==1 | replicate) {
    parStruct <- new("ParStruct",repeated=FALSE)
  } else {
    #parStruct <- new("RepParStruct")
    parStruct <- new("ParStruct",repeated=TRUE)
    #parStruct@cases <- as.integer(length(ntimes))
  }
  parStruct@parameters <- unlist(parameters)
  parStruct@skeleton <- parameters
  #if(replicate) parStruct@replicate <- TRUE else parStruct@replicate <- FALSE

  if(any(c(!is.null(fixed),!is.null(min),!is.null(max),!is.null(id)))) {
    if(is.null(ntimes) | length(ntimes)==1 | replicate) {
      parStruct@repeated <- FALSE
      parStruct <- as(parStruct,"SimpleConstraintsParStruct")
    } else {
      #parStruct <- as(parStruct,"RepSimpleConstraintsParStruct")
      parStruct@repeated <- TRUE
      parStruct <- as(parStruct,"SimpleConstraintsParStruct")
      #parStruct@cases <- as.integer(length(ntimes))
    }
    if(!is.null(fixed)) {
      if(is.null(ntimes) | length(ntimes)==1 | replicate) {
        if(is.list(fixed)) fix <- fixedListToVec(fixed,parameters) else fix <- fixed
      } else {
        cases <- length(ntimes)
        #fix <- vector("list",length=cases)
        if(is.list(fixed)) {
          fix <- fixedRepListToVec(fixed,parameters,cases)
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
    if(!is.null(id)) {
      if(is.null(ntimes) | length(ntimes)==1 | replicate) {
        if(is.list(id)) id <- idListToVec(id,parameters)
      } else {
        cases <- length(ntimes)
        if(is.list(id)) {
          id <- idRepListToVec(id,parameters,cases)
        } else {
          if(length(id) != length(unlist(parameters))) {
            id <- rep(id,length=length(unlist(parameters)))
          }
        }
      }
    } else {
      id <- seq(from=1,to=length(unlist(parameters)))
    }
    parStruct@id <- unlist(id)
    
    if(!is.null(min)) {
      if(is.null(ntimes) | length(ntimes)==1 | replicate) {
        if(is.list(min)) min <- minListToVec(min,parameters)
      } else {
        cases <- length(ntimes)
        if(is.list(min)) {
          min <- minRepListToVec(min,parameters,cases)
        } else {
          if(length(min) != length(unlist(parameters))) {
            min <- rep(min,length=length(unlist(parameters)))
          }
        }
      }
    } else {
      min <- rep(-Inf,length=length(unlist(parameters)))
    }
    parStruct@min <- unlist(min)
    
    if(!is.null(max)) {
      if(is.null(ntimes) | length(ntimes)==1 | replicate) {
        if(is.list(max)) max <- maxListToVec(max,parameters)
      } else {
        cases <- length(ntimes)
        if(is.list(max)) {
          max <- maxRepListToVec(max,parameters,cases)
        } else {
          if(length(max) != length(unlist(parameters))) {
            max <- rep(max,length=length(unlist(parameters)))
          }
        }
      }
    } else {
      max <- rep(Inf,length=length(unlist(parameters)))
    }
    parStruct@max <- unlist(max)
  }
  #if(!is.null(constraints)) parStruct@constraints <- constraints else parStruct@constraints <- new("NoConstraints")
  parStruct
}


setMethod("rep",signature(x="ParStruct"),
  function(x,times=1) {
    if(x@replicate) return(x)
    x@parameters <- rep(x@parameters,times)
    x@skeleton <- rep(x@skeleton,times)
    #x@id <- rep(x@id,times)
    #x@fix <- rep(x@fix,times)
    #x@constraints <- rep(x@constraints,times=times)
    return(x)
  }
)

setMethod("rep",signature(x="SimpleConstraintsParStruct"),
  function(x,times=1) {
    if(x@replicate) return(x)
    x@parameters <- rep(x@parameters,times)
    x@skeleton <- rep(x@skeleton,times)
    x@id <- rep(x@id,times)
    x@fix <- rep(x@fix,times)
    #x@constraints <- rep(x@constraints,times=times)
    return(x)
  }
)

#setMethod("getConstraints",signature(object="ParStruct"),
#  function(object,...) {
#    getConstraints(object@parStruct,...)
#    #return(object@parStruct@constraints)  
#  }
#)


#setMethod("getFreePars",signature(object="ANY"),
#  function(object,...) {
#    getPars(object,which="free",...)
#  }         
#)

#setMethod("setFreePars",signature(object="ANY",pars="numeric"),
#  function(object,pars,calledBy,...) {
#    setPars(object=object,pars=pars,calledBy=calledBy,....,rval="object")
#  }
#)