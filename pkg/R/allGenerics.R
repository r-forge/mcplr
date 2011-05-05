.First.lib <- function(lib, pkg) {
  require(stats4)
}

.Last.lib <- function(libpath) {}

# Guess what: all generics

setGeneric("fit",function(object,...) standardGeneric("fit"))
setGeneric("runm",function(object,...) standardGeneric("runm"))
#setGeneric("estimate",function(object,...) standardGeneric("estimate"))
setGeneric("getPars",function(object,name,...) standardGeneric("getPars"))
setGeneric("setPars",function(object,pars,...,rval=c("object","parameters")) standardGeneric("setPars"))
setGeneric("getReplication",function(object,case,...) standardGeneric("getReplication"))
setGeneric("lFr",function(x,y,...) standardGeneric("lFr"))
setGeneric("rFl",function(x,y,...) standardGeneric("rFl"))
#setGeneric("BIC",function(object,...) standardGeneric("BIC"))
setGeneric("AICc",function(object,...) standardGeneric("AICc"))
setGeneric("Rsq",function(object,...) standardGeneric("Rsq"))
setGeneric("is.unconstrained",function(object,...) standardGeneric("is.unconstrained"))
setGeneric("simulate", function(object,nsim=1,seed=NULL, ...) standardGeneric("simulate"))
setGeneric("has.runm", function(object,...) standardGeneric("has.runm"))
setGeneric("has.lFr", function(object,...) standardGeneric("has.lFr"))
setGeneric("has.rFl", function(object,...) standardGeneric("has.rFl"))
setGeneric("canRepar", function(object,...) standardGeneric("canRepar"))
setGeneric("getTransPars",function(object,...) standardGeneric("getTransPars"))
setGeneric("setTransPars",function(object,pars,...,rval=c("object","parameters")) standardGeneric("setTransPars"))
