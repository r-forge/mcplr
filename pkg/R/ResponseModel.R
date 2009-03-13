setClass("ResponseModel",
  contains="McplBaseModel"
)

setMethod("runm",signature(object="ResponseModel"),
  function(object,...) {
    object
  }
)

setMethod("summary",signature(object="ResponseModel"),
  function(object,...) {
    cat("Response Model, class:",is(object)[1],"\n\n")
    callNextMethod(object=object,...)
  }
)

setMethod("show",signature(object="ResponseModel"),
  function(object) {
    cat("Response Model, class:",is(object)[1],"\n\n")
    callNextMethod()
  }
)


