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





