mcpl.prepare <- function(formula,data,subset,base=NULL,remove.intercept=FALSE) {
  # Returns a matrix of numeric input variables and output variable
  if (missing(data)) {
      data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf,"any")
#  if(is.factor(eval(attr(mt,"variables"),envir=mf)[[1]])) {
#    Y <- model.response(mf)
#    Y <- model.matrix(~Y-1)
#    if(!is.null(base)) Y <- Y[,-base]
#  } else {
#    Y <- model.response(mf)
#  }
  if(is.factor(Y)) {
    Y <- model.matrix(~Y-1)
    if(!is.null(base)) Y <- Y[,-base]
  }
  if(is.vector(Y)) Y <- matrix(Y,ncol=1)
  X <- if(!is.empty.model(mt)) {
    model.matrix(mt, mf, contrasts)
  } else matrix(, NROW(Y), 0)
  Xn <- as.matrix(X)
  if(remove.intercept & "(Intercept)"%in%dimnames(Xn)[[2]]) Xn <- Xn[,-1]
  x.names <- dimnames(Xn)[[2]]
  #y.names <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
  return(list(x=Xn,y=Y,x.names=x.names))
}

