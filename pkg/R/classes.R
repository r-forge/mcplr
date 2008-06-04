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
    if(is.null(ntimes) | replicate) {
      fix <- rep(FALSE,length(unlist(parameters)))
      fix <- relist(fix,skeleton=parameters)
    } else {
      fix <- rep(FALSE,length(unlist(parameters[[1]])))
      fix <- relist(fix,skeleton=parameters[[1]])
    }
    if(!is.null(fixed)) {
      for(i in 1:length(fix)) {
        if(!is.null(fixed[[names(fix)[i]]]) && fixed[[names(fix)[i]]]) fix[[i]] <- rep(TRUE,length(fix[[i]]))
      }
    }
    if(is.null(ntimes) | replicate) {
      parStruct@fix <- unlist(fix)
    } else {
      fix <- rep(list(fix,length(ntimes)))
      parStruct@fix <- unlist(fix)
    }
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





