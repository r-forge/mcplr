# learningModel
### x = cues
### y = criterion
### parameters
### contraints = list

setClass("ConstraintsList",
#  contains="list"
)

setClass("BoxConstraintsList",
  contains="ConstraintsList",
  representation(
    min="numeric",
    max="numeric"
  )
)

setClass("LinConstraintsList",
  contains="ConstraintsList",
  representation(
    Amat = "matrix",
    bvec = "numeric"
  )
)

setClass("ParStruct",
  representation(
    id="integer",
    fix="logical",
    constraints="ConstraintsList",
    replicate="logical" # switch to indicate identical parameters for each nrep
  )
)

ParStruct <- function(parameters,replicate=TRUE,fixed=NULL,ntimes=NULL,
                constraints=NULL) {
  parStruct <- new("ParStruct")
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


setClass("NTimes",
  representation(
    n="integer",
    cases="integer",
    bt="integer",
    et="integer"
  )
)

nTimes <- function(ntimes) {
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)
  out <- new("NTimes",
    n=as.integer(ntimes),
    cases=as.integer(lt),
    bt=as.integer(bt),
    et=as.integer(et)
  )
  out
}


setMethod("BIC",signature(object="missing"),
  function(object,logL,npar,nobs,...) {
    -2*logL + npar*log(nobs)
  }
)
setMethod("AIC",signature(object="missing"),
  function(object,logL,npar,nobs,...,k=2) {
    -2*logL + k*npar
  }
)
setMethod("AICc",signature(object="missing"),
  function(object,logL,npar,nobs,...) {
    -2*logL + (2*nobs*npar)/(nobs-npar-1)
  }
)





