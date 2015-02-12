# TransParStruct: Transformable parameters

setClass("TransParStruct",
  contains="ParStruct",
  representation(
    transform="function",
    #inv.transform="function",
    fix="logical"#,
    #tfix="logical"
  )
)

TransParStruct <- function(parameters,transform,replicate=TRUE,fixed,ntimes) {
  tmp <- ParStruct(parameters=parameters,replicate=replicate,fixed=fixed,
                   ntimes=ntimes)
  parStruct <- new("TransParStruct",
    parameters=tmp@parameters,
    skeleton=tmp@skeleton,
    repeated=tmp@repeated,
    fix=tmp@fix,
  #parStruct@fix <- fixed
    transform = transform
  )
  parStruct
}

setMethod("getPars",signature(object="TransParStruct",name="ANY",replication="ANY"),
  function(object,name,which="all",...) {
    #cb <- match.arg(calledBy)
    #if(cb == "fit") {
    #  object <- object@transform(object,...)
    #} #else tmp <- object
    
    # transform the list of parameters
    tmp <- object@transform(relist(object@parameters,skeleton=object@skeleton),...)
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

setMethod("getFreePars",signature(object="TransParStruct"),
  function(object,...) {
    #cb <- match.arg(calledBy)
    #object <- object@transform(object,...)
    #unlist(object@parameters)[!object@tfix]
    return(object@parameters[!object@fix])
  }
)

#setMethod("getPars",signature(object="TransParStruct",name="ANY"),
#  function(object,which=c("all","free"),calledBy=c("user","fit"),...) {
#    callNextMethod()
    #which <- match.arg(which)
    #cb <- match.arg(calledBy)
    #if(cb == "fit") {
    #  object <- object@transform(object,...)
    #}
    #object <- as(object,"ParStruct")
    #getFreePars(object=object,which=which,calledBy=cb,...)
#  }
#)

setMethod("setPars",signature(object="TransParStruct",pars="list"),
    function(object,pars,calledBy=c("user","fit"),...,rval=c("object","parameters")) {
      stop("setPars for TransParStruct not yet supported")
      #cb <- match.arg(calledBy)
      #rval <- match.arg(rval)
      #if(rval == "object") object <- callNextMethod() else object@parameters <- callNextMethod()
      #if(cb == "fit") object <- object@inv.transform(object,...)
      #if(rval=="object") return(object) else return(object@parameters)
    }
          
)

setMethod("setFreePars",signature(object="TransParStruct",pars="numeric"),
    function(object,pars,...,rval=c("object","parameters")) {      
      rv <- match.arg(rval)
      if(length(pars)!=sum(!object@fix)) stop("length of pars provided does not match free parameters")
      object@parameters[!object@fix] <- pars
      if(rv=="object") return(object) else return(object@parameters)
    }
)


#       
#       if(cb == "fit") {
#         # transform
#         #object <- object@inv.transform(object,...)
#         if(is.list(pars)) {
#           # should not really happen...
#           #object@pars <- object@inv.transform(pars,...)
#           pars <- unlist(pars)
#         } #else {
#           # try to relist it
#         object@parameters <- relist(pars,skeleton=object@parameters)
#         object <- object@inv.transform(object,...)
#         pars <- unlist(pars)
#         #}
#       }
#       callNextMethod(object,pars,...,rval=rval)
#    }
#)

# setMethod("setPars",signature(object="RepTransParStruct"),
#   function(object,pars,calledBy=c("user","fit"),...,rval=c("object","parameters")) {
#     cb <- match.arg(calledBy)
#     rval <- match.arg(rval)
#     if(cb == "fit") {
#       # transform
#       #object <- object@inv.transform(object,...)
#       if(is.list(pars)) {
#         # should not really happen...
#         pars <- object@inv.transform(pars,...)
#         pars <- unlist(pars)
#       } else {
#         # try to relist it
#         pars <- relist(pars,skeleton=object@parameters)
#         for(i in 1:object@cases) {
#           pars[[i]] <- object@inv.transform(pars[[i]],...)
#         }
#         pars <- unlist(pars)
#       }
#     }
#     callNextMethod(object,pars,...,rval=rval)
#   }
# )

# setMethod("getConstraints",signature(object="TransParStruct"),
#   function(object,...) {
#     return(new("Unconstrained",
#                dim=length(getFreePars(object,...))))
#   }       
# )