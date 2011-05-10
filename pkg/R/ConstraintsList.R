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
