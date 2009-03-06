setGeneric("showPars",function(object,...) standardGeneric("showPars"))
setMethod("showPars",signature(object="ParStruct"),
  function(object,...) {
    print(object@parameters,...)
  }
)

setMethod("showPars",signature(object="McplBaseModel"),
  function(object,...) {
    showPars(object@parStruct,...)
  }
)