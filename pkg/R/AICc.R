setMethod("AICc", signature(object="logLik"),
  function(object, ...) {
    -2 * c(object) + 2*attr(object, "df") * (attr(object, "nobs")/(attr(object, "nobs") - attr(object,"df") - 2) )
  }
)

setMethod("AICc", signature(object="ANY"),
          function(object, ...) AICc(object=logLik(object, ...)) )