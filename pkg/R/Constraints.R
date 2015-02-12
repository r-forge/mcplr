setClass("Constraints",
  contains="VIRTUAL"
)

setClass("Unconstrained",
  contains="Constraints",
  representation(
    dim="integer"
  )
)

setClass("BoxConstraints",
  contains="Constraints",
  representation(
    min="numeric",
    max="numeric"
  )
)

setClass("LinearConstraints",
  contains="Constraints",
  representation(
    Amat = "matrix",
    bvec = "numeric"
  )
)

setMethod("combineConstraints",signature(x="Unconstrained",y="Unconstrained"),
  function(x,y) {
    return(new("Unconstrained",
               dim = x@dim+y@dim))
  }
)

setMethod("combineConstraints",signature(x="Unconstrained",y="BoxConstraints"),
  function(x,y) {
    min <- c(rep(-Inf,x@dim),y@min)
    max <- c(rep(Inf,x@dim),y@max)
    return(new("BoxConstraints",
     min=min,
     max=max)
    )
  }
)

setMethod("combineConstraints",signature(x="BoxConstraints",y="Unconstrained"),
  function(x,y) {
    min <- c(x@min,rep(-Inf,y@dim))
    max <- c(x@max,rep(Inf,y@dim))
    return(new("BoxConstraints",
       min=min,
       max=max)
    )
  }
)

setMethod("combineConstraints",signature(x="Unconstrained",y="LinearConstraints"),
  function(x,y) {
    A <- cbind(matrix(0,ncol=x@dim,nrow=nrow(y@Amat)))
    return(new("LinearConstraints",
      Amat=A,
      bvec=y@bvec)
    )
  }
)

setMethod("combineConstraints",signature(x="LinearConstraints",y="Unconstrained"),
  function(x,y) {
    A <- cbind(matrix(0,ncol=y@dim,nrow=nrow(x@Amat)))
    return(new("LinearConstraints",
         Amat=A,
         bvec=x@bvec)
    )
  }
)

setMethod("combineConstraints",signature(x="BoxConstraints",y="BoxConstraints"),
  function(x,y) {
    min <- c(x@min,y@min)
    max <- c(x@max,y@max)
    return(new("BoxConstraints",
               min=min,
               max=max)
    )
  }
)

setMethod("combineConstraints",signature(x="BoxConstraints",y="LinearConstraints"),
  function(x,y) {
    min <- x@min
    max <- x@max
    npar <- length(min)
    minMat <- matrix(0.0,ncol=npar,nrow=npar)
    diag(minMat) <- as.numeric(min > -Inf)
    maxMat <- matrix(0.0,ncol=npar,nrow=npar)
    diag(maxMat) <- -as.numeric(max < Inf)
    A1 <- rbind(minMat[rowSums(minMat) > 0,],maxMat[rowSums(maxMat) > 0,])
    b1 <- c(min[min > -Inf],-max[max < Inf])
    A2 <- y@Amat
    b2 <- y@bvec
    Amat <- bdiag(list(A1,A2))
    bvec <- c(b1,b2)
    return(new("LinearConstraints",
               Amat=Amat,
               bvec=bvec)
    )
  }
)

setMethod("combineConstraints",signature(x="LinearConstraints",y="BoxConstraints"),
  function(x,y) {
    min <- y@min
    max <- y@max
    npar <- length(min)
    minMat <- matrix(0.0,ncol=npar,nrow=npar)
    diag(minMat) <- as.numeric(min > -Inf)
    maxMat <- matrix(0.0,ncol=npar,nrow=npar)
    diag(maxMat) <- -as.numeric(max < Inf)
    A2 <- rbind(minMat[rowSums(minMat) > 0,],maxMat[rowSums(maxMat) > 0,])
    b2 <- c(min[min > -Inf],-max[max < Inf])
    A1 <- x@Amat
    b1 <- x@bvec
    Amat <- bdiag(list(A1,A2))
    bvec <- c(b1,b2)
    return(new("LinearConstraints",
               Amat=Amat,
               bvec=bvec)
    )
  }
)

setMethod("combineConstraints",signature(x="LinearConstraints",y="LinearConstraints"),
  function(x,y) {
    Amat <- bdiag(list(x@Amat,y@Amat))
    bvec <- c(x@bvec,y@bvec)
    new("LinearConstraints",
        Amat = Amat,
        bvec=bvec)
  }          
)

# setMethod("getConstraints",signature(object="McplModel"),
#           function(object,...) {
#             lCon <- getConstraints(object@learningModel)
#             rCon <- getConstraints(object@responseModel)
#             if(is(lCon,"NoConstraints") & is(rCon,"NoConstraints")) return(new("NoConstraints"))
#             if(is(lCon,"LinearConstraints") | is(rCon,"LinearConstraints")) {
#               npar <- length(unlist(getPars(object@learningModel)))
#               if(is(lCon,"LinearConstraints")) {
#                 A1 <- lCon@Amat
#                 b1 <- lCon@bvec
#               } else if(is(lCon,"BoxConstraints")) {
#                 # minimum first
#                 minMat <- matrix(0.0,ncol=npar,nrow=npar)
#                 diag(minMat) <- as.numeric(lCon@min > -Inf)
#                 maxMat <- matrix(0.0,ncol=npar,nrow=npar)
#                 diag(maxMat) <- -as.numeric(lCon@max < Inf)
#                 A1 <- rbind(minMat[rowSums(minMat) > 0,],maxMat[rowSums(maxMat) > 0,])
#                 b1 <- c(lCon@min[lCon@min > -Inf],-lCon@max[lCon@max < Inf])
#               } else if(is(lCon,"NoConstraints")) {
#                 A1 <- matrix(0,ncol=npar,nrow=0)
#                 b1 <- vector(,length=0)
#               }  else {
#                 stop("Cannot determine constraints")
#               }
#               npar <- length(unlist(getPars(object@responseModel)))
#               if(is(rCon,"LinearConstraints")) {
#                 A2 <- rCon@Amat
#                 b2 <- rCon@bvec
#               } else if(is(rCon,"BoxConstraints")) {
#                 # minimum first
#                 minMat <- matrix(0.0,ncol=npar,nrow=npar)
#                 diag(minMat) <- as.numeric(rCon@min > -Inf)
#                 maxMat <- matrix(0.0,ncol=npar,nrow=npar)
#                 diag(maxMat) <- -as.numeric(rCon@max < Inf)
#                 A2 <- rbind(minMat[rowSums(minMat) > 0,],maxMat[rowSums(maxMat) > 0,])
#                 b2 <- c(rCon@min[rCon@min > -Inf],-rCon@max[rCon@max < Inf])
#               } else if(is(rCon,"NoConstraints")) {
#                 A2 <- matrix(0,ncol=npar,nrow=0)
#                 b2 <- vector(,length=0)
#               }  else {
#                 stop("Cannot determine constraints")
#               }
#               Amat <- bdiag(list(A1,A2))
#               bvec <- c(b1,b2)
#               return(new("LinearConstraints",
#                          Amat = Amat,
#                          bvec = bvec)
#               )
#             } else if(is(lCon,"BoxConstraints") | is(rCon,"BoxConstraints")) {
#               # no linear constraints, so can just use box constraints if necessary
#               npar <- length(unlist(getPars(object@learningModel)))
#               if(is(lCon,"BoxConstraints")) {
#                 # minimum first
#                 min1 <- lCon@min
#                 max1 <- lCon@max
#               } else if(is(lCon,"NoConstraints")){
#                 min1 <- rep(-Inf,npar)
#                 max1 <- rep(Inf,npar)
#               } else {
#                 stop("Cannot determine constraints")
#               }
#               npar <- length(unlist(getPars(object@responseModel)))
#               if(is(rCon,"BoxConstraints")) {
#                 # minimum first
#                 min2 <- rCon@min
#                 max2 <- rCon@max
#               } else if(is(rCon,"NoConstraints")){
#                 min2 <- rep(-Inf,npar)
#                 max2 <- rep(Inf,npar)
#               } else {
#                 stop("Cannot determine constraints")
#               }  
#             }
#           }
# )
