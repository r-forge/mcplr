# Radial Basis Function Network
# a simple RBFN with fixed kernel locations, but estimable kernel width.
# RBFN uses radial basis for x only
# RBFN2 uses radial basis for x and y

setClass("RBFN",
  contains="SLFN",
  representation(
    phi = "matrix",
    basis.function="function",
    basis.dimension="numeric",
    basis.location="list"
  )
)
setClass("RBFN2",
  contains="RBFN",
  representation(
    y.phi = "matrix",
    y.basis.function="function",
    y.basis.dimension="numeric",
    y.basis.location="list"
  )
)
setMethod("setPars",signature(object="RBFN"),
  function(object,pars,unconstrained=FALSE,...,rval=c("object","parameters")) {
    sP.u <- function(pars,fixl) {
#      if(!is.null(pars$eta) && !fixl$eta) pars$eta <- exp(pars$eta)
#      if(!is.null(pars$alpha) && !fixl$alpha) pars$alpha <- exp(pars$alpha)
#      if(!is.null(pars$beta) && !fixl$beta) pars$beta <- exp(pars$beta)/(1+exp(pars$beta))
      if(!is.null(pars$basis.scale) && !fixl$basis.scale) {
        if(!inherits(pars$basis.scale,"pdMat")) pars$basis.scale <- exp(pars$basis.scale)
      }
      if(!is.null(pars$y.basis.scale) && !fixl$y.basis.scale) {
        if(!inherits(pars$y.basis.scale,"pdMat")) pars$y.basis.scale <- exp(pars$y.basis.scale)
      }
      pars
    }
    parl <- callNextMethod(object=object,pars=pars,rval="parameters",unconstrained=unconstrained,...)
    if(length(object@parStruct@fix)>0) fixl <- relist(object@parStruct@fix,skeleton=object@parameters) else fixl <- relist(rep(FALSE,length(unlist(object@parameters))),skeleton=object@parameters)
    if(unconstrained) {
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) parl[[case]] <- sP.u(parl[[case]],fixl[[case]])
      } else {
        parl <- sP.u(parl,fixl)
      }
    }
    switch(rval,
      object = {
        object@parameters <- parl
        object
      },
      parameters = parl)
    #object@parameters
    #return(object)
  }
)

setMethod("getPars",signature(object="RBFN"),
  function(object,which="all",unconstrained=FALSE,...) {
    gP.u <- function(pars) {
#      if(!is.null(pars$alpha)) pars$alpha <- log(pars$alpha)
#      if(!is.null(pars$eta)) pars$eta <- log(pars$eta)
#      if(!is.null(pars$beta)) pars$beta <- log(pars$beta/(1-pars$beta))
      if(!is.null(pars$basis.scale)) {
        if(!inherits(pars$basis.scale,"pdMat")) pars$basis.scale <- log(pars$basis.scale)
      }
      if(!is.null(pars$y.basis.scale)) {
        if(!inherits(pars$y.basis.scale,"pdMat")) pars$y.basis.scale <- log(pars$y.basis.scale)
      }
      pars
    }
    if(unconstrained) {
      parl <- object@parameters
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) parl[[case]] <- gP.u(parl[[case]])
      } else {
        parl <- gP.u(parl)
      }
      object@parameters <- parl
    }
    pars <- callNextMethod(object=object,which=which,unconstrained=unconstrained,...)
    return(pars)
  }
)

setMethod("fit","RBFN",
  function(object,...) {
    x2phi.additive <- function(x,fun,loc,sca) {
      out <- vector()
      nx <- ncol(x)
      for(i in 1:nx) {
        out <- cbind(out,apply(as.matrix(loc[[i]]),1,fun,x=x[,i],scale=sca[[i]]))
      }
      out
    }
    x2phi.configural <- function(x,fun,loc,sca) {
      apply(loc,1,fun,x=x,scale=sca)
    }
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        intercept <- FALSE
        # transform x to basis phi
        if(all(x[,1]==1)) {
          intercept <- TRUE
          x <- x[,-1]
        }
        if(object@basis.dimension > 1) {
          object@phi[object@nTimes@bt[case]:object@nTimes@et[case],] <- bx <- x2phi.additive(x=x,fun=object@basis.function,loc=object@basis.location,sca=object@parameters$basis.scale)
        } else {
          object@phi[object@nTimes@bt[case]:object@nTimes@et[case],] <- bx <- x2phi.configural(x=x,fun=object@basis.function,loc=object@basis.location[[1]],sca=object@parameters$basis.scale)
        }
#        bx <- matrix(nrow=nrow(x),ncol=object@basis.n)
#        for(i in 1:object@basis.n) {
#          bx[,i] <- apply(x,1,object@basis.function,location=object@basis.location[[i]],scale=object@basis.scale[[i]])
#        }
        if(intercept) bx <- cbind(1,bx) # NOTE: phi does not contain intercept!!!
        pars <- repl$parameters
        fit <- slfn.fit(x=bx,y=y,eta=pars$eta,alpha=pars$alpha,beta=pars$beta,ws=pars$ws,grad=object@gradient,window.size=object@window.size)
        object@weight[object@nTimes@bt[case]:object@nTimes@et[case],,] <- fit$weight
      }
    } else {
      x <- object@x
      intercept <- FALSE
      if(all(x[,1]==1)) {
        intercept <- TRUE
        x <- x[,-1]
      }
      if(object@basis.dimension > 1) {
        object@phi <- bx <- x2phi.additive(x=x,fun=object@basis.function,loc=object@basis.location,sca=object@parameters$basis.scale)
      } else {
        object@phi <- bx <- x2phi.configural(x=x,fun=object@basis.function,loc=object@basis.location[[1]],sca=object@parameters$basis.scale)
      }
      if(intercept) bx <- cbind(1,bx)
      fit <- slfn.fit(x=bx,y=object@y,eta=object@parameters$eta,alpha=object@parameters$alpha,beta=object@parameters$beta,ws=object@parameters$ws,grad=object@gradient,window.size=object@window.size)
      object@weight <- fit$weight
    }
    return(object)
  }
)

setMethod("fit","RBFN2",
  function(object,...) {
    x2phi.additive <- function(x,fun,loc,sca) {
      out <- vector()
      nx <- ncol(x)
      for(i in 1:nx) {
        out <- cbind(out,apply(as.matrix(loc[[i]]),1,fun,x=x[,i],scale=sca[[i]]))
      }
      out
    }
    x2phi.configural <- function(x,fun,loc,sca) {
      apply(loc,1,fun,x=x,scale=sca)
    }
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        intercept <- FALSE
        # transform x to basis phi
        if(all(x[,1]==1)) {
          intercept <- TRUE
          x <- x[,-1]
        }
        if(object@basis.dimension > 1) {
          object@phi[object@nTimes@bt[case]:object@nTimes@et[case],] <- bx <- x2phi.additive(x=x,fun=object@basis.function,loc=object@basis.location,sca=object@parameters$basis.scale)
        } else {
          object@phi[object@nTimes@bt[case]:object@nTimes@et[case],] <- bx <- x2phi.configural(x=x,fun=object@basis.function,loc=object@basis.location[[1]],sca=object@parameters$basis.scale)
        }
        if(intercept) bx <- cbind(1,bx) # NOTE: phi does not contain intercept!!!
        # transform y to basis phi
        if(object@y.basis.dimension > 1) {
          object@y.phi[object@nTimes@bt[case]:object@nTimes@et[case],] <- phiy <- x2phi.additive(x=y,fun=object@y.basis.function,loc=object@y.basis.location,sca=object@parameters$y.basis.scale)
        } else {
          object@y.phi[object@nTimes@bt[case]:object@nTimes@et[case],] <- phiy <- x2phi.configural(x=y,fun=object@y.basis.function,loc=object@y.basis.location[[1]],sca=object@parameters$y.basis.scale)
        }
        pars <- repl$parameters
        fit <- slfn.fit(x=bx,y=phiy,eta=pars$eta,alpha=pars$alpha,beta=pars$beta,ws=pars$ws,grad=object@gradient,window.size=object@window.size)
        object@weight[object@nTimes@bt[case]:object@nTimes@et[case],,] <- fit$weight
      }
    } else {
      x <- object@x
      y <- object@y
      # transform x to basis phi
      intercept <- FALSE
      if(all(x[,1]==1)) {
        intercept <- TRUE
        x <- x[,-1]
      }
      if(object@basis.dimension > 1) {
        object@phi <- bx <- x2phi.additive(x=x,fun=object@basis.function,loc=object@basis.location,sca=object@parameters$basis.scale)
      } else {
        object@phi <- bx <- x2phi.configural(x=x,fun=object@basis.function,loc=object@basis.location[[1]],sca=object@parameters$basis.scale)
      }
      if(intercept) bx <- cbind(1,bx)
      # transform y to basis phi
      if(object@y.basis.dimension > 1) {
        object@y.phi <- phiy <- x2phi.additive(x=y,fun=object@y.basis.function,loc=object@y.basis.location,sca=object@parameters$y.basis.scale)
      } else {
        object@y.phi <- phiy <- x2phi.configural(x=y,fun=object@y.basis.function,loc=object@y.basis.location[[1]],sca=object@parameters$y.basis.scale)
      }
      fit <- slfn.fit(x=bx,y=phiy,eta=object@parameters$eta,alpha=object@parameters$alpha,beta=object@parameters$beta,ws=object@parameters$ws,grad=object@gradient,window.size=object@window.size)
      object@weight <- fit$weight
    }
    return(object)
  }
)
rbfn.basis <- function(name=c("gaussian","exponential")) {
  name <- match.arg(name)
  switch(name,
    exponential = function(x,location,scale) {
      if(is.matrix(x) && ncol(x) > 1) {
        x <- abs(t(x)-location)
        out <- t(exp(-rowSums(scale*x)))
      } else {
        out <- exp(-scale*abs(x-location))
      }
      out
    },
    gaussian = function(x,location,scale) {
      if(is.matrix(x) && ncol(x) > 1) {
        if(!is.matrix(scale)) scale <- scale*diag(ncol(x))
        out <- exp(-.5*mahalanobis(x,location,scale))
        #iscale <- as.matrix(solve(scale))
        #dx <- sweep(x, 2, location)
        #out <- .5*rowSums((dx %*% iscale) * dx)
        #logdet <- sum(log(eigen(scale, symmetric = TRUE, only.values = TRUE)$values))
        #x <- t(x)-location
        
        #out <- apply(x,2,function(x) exp(-.5*(t(x)%*%iscale%*%x)))
      } else {
        out <- exp(-.5*((x-location)/scale)^2)
      }
      out
    }
  )
}





rbfn <- function(formula,parameters=list(eta=.01,alpha=0,beta=0,ws=0),basis=rbfn.basis(name="gaussian"),basis.location,configural=TRUE,basis.grid.n=FALSE,type=c("linear","logistic"),data,window.size=0,intercept=TRUE,base=NULL,ntimes=NULL,replicate=TRUE,subset) {
  # parameters <- list with at least "eta" and "basis.scale"
  # possible other values include alpha, beta, ws
  cl <- match.call()
  sf <- match.call(expand.dots = FALSE)
  m <- match(names(formals(slfn)), names(sf), 0)
  sf <- sf[c(1, m)]
  sf[[1]] <- as.name("slfn")
  sf <- eval(sf, parent.frame())
  rbf <- as(sf,"RBFN")
  rbf@basis.function <- basis
  if(intercept) x <- rbf@x[,-1] else x <- rbf@x
  if(missing(basis.location)) {
    if(configural) {
      if(basis.grid.n) {
        #use an equally spaced grid
        basis.location <- list()
        for(i in 1:ncol(x)) {
          basis.location[[i]] <- seq(min(x[,i]),max(x[,i]),length=round(basis.grid.n^(1/ncol(x))))
        }
        basis.location <- expand.grid(basis.location)
      } else {
        # use x values as nodes
        basis.location = list(unique(x))
      }
    } else {
      if(basis.grid.n) {
        #use an equally spaced grid
        basis.location <- list()
        for(i in 1:ncol(x)) {
          basis.location[[i]] <- seq(min(x[,i]),max(x[,i]),length=round(basis.grid.n/ncol(x)))
        }
      } else {
        # use x values as nodes
        for(i in 1:ncol(x)) {
          basis.location[[i]] <- unique(x[order(x[,i]),i])
        }
      }
    }
  }
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$basis.scale)) {
      basis.scale <- list()
      for(i in 1:ncol(x)) {
        if(configural) tmp <- basis.location[[1]][,i] else tmp <- basis.location[[i]]
        basis.scale[[i]] <- .6*mean(diff(tmp[order(tmp)]))
      }
      if(configural) basis.scale <- nlme::pdDiag(diag(unlist(basis.scale)))
      rbf@parameters$basis.scale <- basis.scale
    }
  } else {
    for(i in 1:length(rbf@parameters)) {
      basis.scale <- list()
      for(i in 1:ncol(x)) {
        if(configural) tmp <- basis.location[[1]][,i] else tmp <- basis.location[[i]]
        basis.scale[[i]] <- .6*mean(diff(tmp[order(tmp)]))
      }
      if(configural) basis.scale <- nlme::pdDiag(diag(unlist(basis.scale)))
      rbf@parameters[[i]]$basis.scale <- basis.scale
    }
  }
  rbf@basis.location  <- basis.location
  if(configural) rbf@basis.dimension <- 1 else {
    rbf@basis.dimension <- ncol(x)
  }
  if(configural) nphi <- nrow(rbf@basis.location[[1]]) else nphi <- sum(unlist(lapply(rbf@basis.location),length))
  rbf@phi <- matrix(,nrow=nrow(x),ncol=nphi)
  rbf
}

rbfn2 <- function(formula,parameters=list(eta=.01,alpha=0,beta=0,ws=0),basis=rbfn.basis(name="gaussian"),basis.location,configural=TRUE,basis.grid.n=FALSE,y.basis=rbfn.basis(name="gaussian"),y.basis.location,y.configural=TRUE,y.basis.grid.n=FALSE,type=c("linear","logistic"),data,window.size=0,intercept=TRUE,base=NULL,ntimes=NULL,replicate=TRUE,subset) {
  cl <- match.call()
  sf <- match.call(expand.dots = FALSE)
  m <- match(names(formals(rbfn)), names(sf), 0)
  sf <- sf[c(1, m)]
  sf[[1]] <- as.name("rbfn")
  sf <- eval(sf, parent.frame())
  rbf <- as(sf,"RBFN2")
  rbf@y.basis.function <- y.basis
  y <- rbf@y
  if(missing(y.basis.location)) {
    if(y.configural) {
      if(y.basis.grid.n) {
        #use an equally spaced grid
        y.basis.location <- list()
        for(i in 1:ncol(y)) {
          y.basis.location[[i]] <- seq(min(y[,i]),max(y[,i]),length=round(y.basis.grid.n^(1/ncol(y))))
        }
        y.basis.location <- expand.grid(y.basis.location)
      } else {
        # use x values as nodes
        y.basis.location = list(as.matrix(unique(y[order(y[,1]),])))
      }
    } else {
      if(y.basis.grid.n) {
        #use an equally spaced grid
        y.basis.location <- list()
        for(i in 1:ncol(x)) {
          y.basis.location[[i]] <- seq(min(y[,i]),max(y[,i]),length=round(y.basis.grid.n/ncol(y)))
        }
      } else {
        # use x values as nodes
        for(i in 1:ncol(y)) {
          y.basis.location[[i]] <- unique(y[order(y[,i]),i])
        }
      }
    }
  }
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$y.basis.scale)) {
      y.basis.scale <- list()
      for(i in 1:ncol(y)) {
        if(y.configural) tmp <- y.basis.location[[1]][,i] else tmp <- y.basis.location[[i]]
        y.basis.scale[[i]] <- .6*mean(diff(tmp[order(tmp)]))
      }
      if(y.configural) {
        if(ncol(y)>1) y.basis.scale <- nlme::pdDiag(diag(unlist(y.basis.scale))) else y.basis.scale <- unlist(y.basis.scale)
      }
      rbf@parameters$y.basis.scale <- y.basis.scale
    }
  } else {
    if(is.null(parameters$y.basis.scale)) {
      for(i in 1:length(rbf@parameters)) {
        y.basis.scale <- list()
        for(i in 1:ncol(x)) {
          if(configural) tmp <- y.basis.location[[1]][,i] else tmp <- y.basis.location[[i]]
          y.basis.scale[[i]] <- .6*mean(diff(tmp[order(tmp)]))
        }
        if(y.configural) {
          if(ncol(y)>1) y.basis.scale <- nlme::pdDiag(diag(unlist(y.basis.scale))) else y.basis.scale <- unlist(y.basis.scale)
        }
        rbf@parameters[[i]]$y.basis.scale <- y.basis.scale
      }
    }
  }
  rbf@y.basis.location  <- y.basis.location
  if(y.configural) rbf@y.basis.dimension <- 1 else {
    rbf@y.basis.dimension <- ncol(y)
  }
  if(y.configural) nphi <- nrow(rbf@y.basis.location[[1]]) else nphi <- sum(unlist(lapply(rbf@y.basis.location),length))
  rbf@y.phi <- matrix(,nrow=nrow(y),ncol=nphi)
  rbf
}

setMethod("predict",signature(object="RBFN"),
  function(object,type="link",...) {
    ny <- ncol(object@y)
    nt <- nrow(object@y)
    nphi <- ncol(object@phi)
    if(ncol(object@weight)==ncol(object@phi)+1) {
      phi <- cbind(1,object@phi)
      nphi <- nphi+1
    } else phi <- object@phi
    #cid <- matrix(1:(nphi*ny),nrow=nphi)
#    pred <- matrix(nrow=nt,ncol=ny)
#    for(i in 1:ny) pred[,i] <- rowSums(object@phi*object@weight[,cid[,i]])
#    if(type=="response") pred <- object@actfun(pred)
    #cid <- matrix(1:(nx*ny),nrow=nx)  
    pred <- apply(array(as.numeric(object@phi)*as.numeric(object@weight),dim=c(nt,nphi,ny)),c(1,3),sum)
    #matrix(nrow=nt,ncol=ny)      
    #for(i in 1:ny) pred[,i] <- rowSums(object@x*object@weight[,cid[,i]])
    if(type=="response") pred <- object@actfun(pred)
    #return(pred)
    return(pred)
  }
)
  
setMethod("predict",signature(object="RBFN2"),
  function(object,type="link",...) {
    ny <- ncol(object@y)
    nt <- nrow(object@y)
    nphi <- ncol(object@phi)
    nphiy <- ncol(object@y.phi)
    if(ncol(object@weight)==ncol(object@phi)+1) {
      # add intercept
      phi <- cbind(1,object@phi)
      nphi <- nphi+1
    } else phi <- object@phi
    #cid <- matrix(1:(nphi*nphiy),nrow=nphi)
    #pred <- matrix(nrow=nt,ncol=nphiy)
    pred <- apply(array(as.numeric(object@phi)*as.numeric(object@weight),dim=c(nt,nphi,nphiy)),c(1,3),sum)
#    for(i in 1:nphiy) pred[,i] <- rowSums(object@phi*object@weight[,cid[,i]])
    if(type=="response") pred <- object@actfun(pred)
#    if(type=="response") {
#      out <- matrix(nrow=nt,ncol=ny)
#      tpred <- object@actfun(pred)
#      if(object@y.basis.dimension==1) {
#        tpred <- apply(tpred,1,function(x) {
#          if(!all(x >=0)) warning("activations should be >= 0")
#          if(sum(x)>0) out <- x/sum(x) else out <- rep(1/length(x),length(x))
#        })
#        for(i in 1:ny) out[,i] <- colSums(tpred*object@y.basis.location[[1]][,i])
#      } else {
        # TODO
#        stop("THIS IS TO DO!")
#      }
#      pred <- out
    #}
    return(pred)
  }
)

setMethod("lFr",signature(x="RBFN",y="ResponseModel"),
  function(x,y,...) {
    y@x <- predict(x,type="response",...)
    y
  }
)