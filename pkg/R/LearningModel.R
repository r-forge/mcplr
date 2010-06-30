setClass("LearningModel",
  contains="McplBaseModel"
)

setMethod("summary",signature(object="LearningModel"),
  function(object,...) {
    cat("Learning Model, class:",is(object)[1],"\n\n")
    callNextMethod(object=object,...)
  }
)

setMethod("show",signature(object="LearningModel"),
  function(object) {
    cat("Learning Model, class:",is(object)[1],"\n\n")
    callNextMethod()
  }
)

setMethod("has.lFr",signature(object="LearningModel"),
  function(object,...) {
    TRUE
  }
)
