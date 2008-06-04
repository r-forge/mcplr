.First.lib <- function(lib, pkg) {
  require(stats4)
}

.Last.lib <- function(libpath) {}

# Guess what: all generics

setGeneric("fit",function(object,...) standardGeneric("fit"))
setGeneric("estimate",function(object,...) standardGeneric("estimate"))
setGeneric("getPars",function(object,...) standardGeneric("getPars"))
setGeneric("setPars",function(object,pars,...,rval=c("object","parameters")) standardGeneric("setPars"))
setGeneric("getReplication",function(object,case,...) standardGeneric("getReplication"))
setGeneric("lFr",function(x,y,...) standardGeneric("lFr"))
setGeneric("rFl",function(x,y,...) standardGeneric("rFl"))
#setGeneric("BIC",function(object,...) standardGeneric("BIC"))
setGeneric("AICc",function(object,...) standardGeneric("AICc"))
setGeneric("RSquare",function(object,...) standardGeneric("RSquare"))
setGeneric("is.unconstrained",function(object,...) standardGeneric("is.unconstrained"))
