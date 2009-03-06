################################################################################
### GCM.r
### Generalized Context Model (GCM, after Nosofsky, 1986)
###
### This implementation allows for a continuous criterion and includes
###     additional `sampling' functions, such that more recent exemplars
###     obtain more weight (with exponential or power law forgetting)
###
### (c) 2007, M. Speekenbrink
################################################################################

setClass("GCM",
  contains="LearningModel",
  representation(
    weights="list", # each element is a lower triangular matrix with (normalized) weight of each y[t]
    distance="function", # returns an nt*nt matrix with distances
    similarity="function", # computes similarities from distance matrix
    sampling="function" # computes nt*nt matrix with sampling weights (column wise!)
  )
)
setClass("GCMnominal",
  contains="GCM"
)
setClass("GCMinterval",
  contains="GCM"
)

setMethod("is.unconstrained",signature(object="GCM"),
  function(object,...) {
    if(is(object@parStruct@constraints,"LinConstraintsList") || is(object@parStruct@constraints,"BoxConstraintsList")) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
)

setMethod("fit",signature(object="GCM"),
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        fit <- gcm.fit(x=x,y=y,parameters=pars,distance=object@distance,similarity=object@similarity,sampling=object@sampling)
        object@weights[[case]] <- fit$weights
      }
    } else {
      fit <- gcm.fit(x=object@x,y=object@y,parameters=object@parStruct@parameters,distance=object@distance,similarity=object@similarity,sampling=object@sampling)
      object@weights[[1]] <- fit$weights
    }
    return(object)    
  }
)

gcm.fit <- gcm.cont.fit <- gcm.discr.fit <- function(x,y,parameters,distance,similarity,sampling,...) {
  dis <- distance(x=x,parameters=parameters,...)
  sim <- similarity(distance=dis,parameters=parameters,...)
  sam <- sampling(nt=nrow(x),parameters=parameters,...)
  w <- sim*sam
  w <- apply(w,2,function(x) x/sum(x))
  w[is.nan(w)] <- 0 # FIX ME!
#  if(!all(colSums(w[,-1])==1)) {
#    warning("problem in mixture weights: do not all sum to 1")
#    ids <- which(colSums(w[,-1])==0)+1
#    w[,ids] <- 1
#    w[lower.tri(w,diag=TRUE)] <- 0
#    w[,ids] <- apply(as.matrix(w[,ids]),2,function(x) x/sum(x))
#    w[is.nan(w)] <- 0 # FIX ME!!
#    if(!all(colSums(w[,-1])==1)) warning("problem remains")
#  }
#  w[,colSums(w)==0]
  return(list(weights=w))
}

gcm.fit.unconstrained <- function(x,y,parameters,distance,similarity,sampling,...) {
  parameters$lambda <- sum(parameters$w)
  parameters$w <- parameters$w/parameters$lambda
  if(is.null(parameters$r)) {
    if(attr(distance,"name") == "euclidian") parameters$r <- 2 else parameters$r <- 1
  }
  if(is.null(parameters$q)) {
    if(attr(similarity,"name") == "gaussian") parameters$q <- 2 else parameters$q <- 1
  }
  parameters$lambda <- parameters$lambda^(parameters$r/parameters$q)
    dis <- distance(x=x,parameters=parameters,...)
  sim <- similarity(distance=dis,parameters=parameters,...)
  sam <- sampling(nt=nrow(x),parameters=parameters,...)
  w <- sim*sam
  w <- apply(w,2,function(x) x/sum(x))
  return(list(weights=w))
}

setMethod("setPars",signature(object="GCM"),
  function(object,pars,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    if(!internal || !is.null(object@parStruct@constraints)) {
      parl <- callNextMethod(object=object,pars=pars,internal=internal,rval="parameters",...)
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases){
          if(length(parl[[case]]$w) == ncol(object@x) - 1) {
            parl[[case]]$w <- c(parl[[case]]$w,1-sum(parl[[case]]$w))
          }
        }
      } else {
        if(length(parl$w) == ncol(object@x) - 1) {
          parl$w <- c(parl$w,1-sum(parl$w))
        }
      }
    } else {
      # use reparametrization
      sP.u <- function(pars,fix) {
        pr <- pars$r
        pq <- pars$q
        if(is.null(pr)) if(attr(object@distance,"name") == "euclidian") pr <- 2 else pr <- 1
        if(is.null(pq)) if(attr(object@distance,"name") == "gaussian") pq <- 2 else pq <- 1
        pars$w <- exp(pars$w)
        pars$w[pars$w == Inf] <- 1e+90
        pars$lambda <- sum(pars$w)
        pars$w <- pars$w/pars$lambda
        pars$lambda <- pars$lambda^(pr/pq)
        if(!is.null(pars$gamma)) pars$gamma <- exp(pars$gamma)
        if(!is.null(fix$sdy) && !fix$sdy) pars$sdy <- exp(pars$sdy)
        pars
      }
      fix <- fixold <- object@parStruct@fix
      if(length(fix)==0) fix <- rep(FALSE,length(unlist(object@parStruct@parameters)))
      fix <- relist(fix,skeleton=object@parStruct@parameters)
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) fix[[case]]$lambda <- TRUE
      } else {
        fix$lambda <- TRUE
      }
      object@parStruct@fix <- unlist(fix)
      #object <- callNextMethod(object=object,pars=pars,...)
      #pars <- object@parStruct@parameters
      parl <- callNextMethod(object=object,pars=pars,internal=internal,rval="parameters",...)
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases){
          parl[[case]] <- sP.u(parl[[case]],fix[[case]])
        }
      } else {
        parl <- sP.u(parl,fix)
      }
      object@parStruct@fix <- fixold
    }
    switch(rval,
      object = {
        object@parStruct@parameters <- parl
        object},
      parameters = parl)
  }
)

setMethod("getPars",signature(object="GCM"),
  function(object,which="all",internal=FALSE,...) {
    gP.u <- function(pars) {
      pr <- pars$r
      pq <- pars$q
      if(is.null(pr)) if(attr(object@distance,"name") == "euclidian") pr <- 2 else pr <- 1
      if(is.null(pq)) if(attr(object@distance,"name") == "gaussian") pq <- 2 else pq <- 1
      pars$lambda <- pars$lambda^{pq/pr}
      pars$w <- log(pars$lambda*pars$w)
      pars$lambda <- NULL
      if(!is.null(pars$gamma)) pars$gamma <- log(pars$gamma)
      if(!is.null(pars$sdy)) pars$sdy <- log(pars$sdy)
      pars
    }
    if(internal && is.null(object@parStruct@constraints)) {
      pars <- object@parStruct@parameters
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases){
          pars[[case]] <- gP.u(pars[[case]])
        }
      } else {
        pars <- gP.u(pars)
      }
      pars <- unlist(pars)
      if(which=="free") {
        if(length(object@parStruct@fix)>0) {
          fixl <- relist(object@parStruct@fix,skeleton=object@parStruct@parameters)
          if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
            for(case in 1:object@nTimes@cases){
              fixl[[case]]$lambda <- NULL
            }
          } else {
            fixl$lambda <- NULL
          }
          fix <- unlist(fixl)
          pars <- pars[!fix]
        }
      }
    } else {
      pars <- callNextMethod(object=object,which=which,...)
    }
    return(pars)
  }
)

setMethod("predict",signature(object="GCM"),
  function(object,...) {
    pred <- vector()
    for(case in 1:object@nTimes@cases) {
      pred <- rbind(pred,apply(as.matrix(object@y[object@nTimes@bt[case]:object@nTimes@et[case],]),2,function(x,weights) {
        colSums(x*weights)
      },weights=object@weights[[case]]))
    }
    pred
  }
)

setMethod("logLik",signature(object="GCMnominal"),
  function(object,discount=1,eps=.Machine$double.eps,...) {
    discount <- unlist(lapply(object@nTimes@bt-1,"+",discount))
    discount <- discount[discount>0]
    pred <- predict(object,type="response",...)
    pred <- rowSums(object@y*pred)[-discount]
    pred[pred > 1-eps] <- 1-eps
    pred[pred < eps] <- eps
    LL <- sum(log(pred))
    out <- LL
    nobs <- length(pred)
    attr(out,"nobs") <- nobs
    attr(out,"df") <- length(getPars(object,which="free"))
    class(out) <- "logLik"
    out
  }
)

setMethod("logLik",signature(object="GCMinterval"),
  function(object,discount=1,...) {
    LL <- vector("double")
    nobs <- 0
    for(case in 1:object@nTimes@cases) {
      L <- outer(as.vector(object@y[object@nTimes@bt[case]:object@nTimes@et[case],]),as.vector(object@y[object@nTimes@bt[case]:object@nTimes@et[case],]),dnorm,sd=object@parStruct@parameters$sdy)
      weights <- object@weights[[case]]
      zw <- which(colSums(weights)==0)
      LL[case] <- sum(log(colSums(weights*L)[-discount]))
      nobs <- nobs + ncol(weights) - length(zw[,-discount]) - length(discount)
      #zwt <- zwt + length(zw)
    }
    out <- sum(LL)
    attr(out,"nobs") <- nobs
    attr(out,"df") <- length(getPars(object,which="free"))
    class(out) <- "logLik"
    out
  }
)

setMethod("estimate",signature(object="GCM"),
  function(object,method="Nelder-Mead",...) {
    optfun <- function(pars,object,repar,...) {
      #object <- setPars(object,pars,unconstrained=unconstrained)
      object@parStruct@parameters <- setPars(object,pars,repar=repar,...,rval="parameters")
      object <- fit(object,...)
      -logLik(object)
    }
    if(!is.unconstrained(object,...)) {
      pars <- getPars(object,which="free",internal=TRUE,,...)
      object@parStruct@parameters <- switch(is(object@parStruct@constraints)[1],
        "LinConstraintsList" = {
          A <- object@parStruct@constraints@Amat
          b <- object@parStruct@constraints@bvec
          opt <- constrOptim(theta=pars,f=optfun,grad=NULL,ui=A,ci=b,object=object,repar=FALSE,...)
          setPars(object,opt$par,internal=TRUE,...,rval="parameters")
        },
        "BoxConstraintsList" = {
          opt <- optim(pars,fn=optfun,method="L-BFGS-B",object=object,min=object@parStruct@constraints@min,max=object@parStruct@constraints@min,repar=FALSE,...)
          setPars(object,opt$par,internal=TRUE,...,rval="parameters")
        },
        {
          warning("This extension of ConstraintsList is not implemented; using default optimisation.")
          opt <- optim(pars,fn=optfun,method=method,object=object,repar=FALSE,...)
          setPars(object,opt$par,internal=TRUE,...,rval="parameters")
        }
      )
    } else {
      pars <- getPars(object,which="free",repar=TRUE,...)
      opt <- optim(pars,fn=optfun,method=method,object=object,repar=TRUE,...)
      object@parStruct@parameters <- setPars(object,opt$par,internal=TRUE,...,rval="parameters")
    }
    object <- fit(object,...)
    object
  }
)

gcm.distance <- function(type="cityblock") {
  dis.city <- function(x,parameters,...) {
    nc <- ncol(x)
    nr <- nrow(x)
    w <- parameters$w
    dis <- matrix(0,ncol=nr,nrow=nr)
    for(i in 1:nc) {
      dis <- dis + w[i]*abs(outer(x[,i],x[,i],"-"))
    }
    dis
  }
  dis.eucl <- function(x,parameters,...) {
    nc <- ncol(x)
    nr <- nrow(x)
    w <- parameters$w
    dis <- matrix(0,ncol=nr,nrow=nr)
    for(i in 1:nc) {
      dis <- dis + w[i]*abs(outer(x[,i],x[,i],"-"))^2
    }
    dis^(1/2)
  }
  dis.mink <- function(x,parameters,...) {
    nc <- ncol(x)
    nr <- nrow(x)
    r <- parameters$r
    w <- parameters$w
    dis <- matrix(0,ncol=nr,nrow=nr)
    for(i in 1:nc) {
      dis <- dis + w[i]*abs(outer(x[,i],x[,i],"-"))^r
    }
    dis^(1/r)
  }
  fun <- switch(type,
    cityblock = dis.city,
    euclidian = dis.eucl,
    minkowski = dis.mink,
    dis.city
  )
  attr(fun,"name") <- type
  fun
}

gcm.similarity <- function(type="exponential") {
  sim.exp <- function(distance,parameters,...) {
    lambda <- parameters$lambda
    sim <- exp(-lambda*distance)
    sim
  }
  sim.gauss <- function(distance,parameters,...) {
    lambda <- parameters$lambda
    sim <- exp(-lambda*distance^2)
    sim
  }
  sim.gen <- function(distance,parameters,...) {
    lambda <- parameters$lambda
    q <- parameters$q
    sim <- exp(-lambda*distance^q)
    sim
  }
  fun <- switch(type,
    exponential = sim.exp,
    gaussian = sim.gauss,
    general = sim.gen,
    sim.exp
  )
  attr(fun,"name") <- type
  fun
}

gcm.sampling <- function(type="uniform") {
  sam.exp <- function(nt,parameters,...) {
    k <- 1:nt
    gamma <- parameters$gamma
    sam <- exp(outer(k,k-1,"-")*gamma)
    sam[lower.tri(sam,diag=TRUE)] <- 0
#    sam <- apply(sam,2,function(x) x/sum(x)) # normalize
#    sam[is.nan(sam)] <- 0
    sam
  }
  sam.pow <- function(nt,parameters,...) {
    k <- 1:nt
    gamma <- parameters$gamma
    sam <- (-outer(k,k,"-"))^(-gamma)
    sam[lower.tri(sam,diag=TRUE)] <- 0
#    sam <- apply(sam,2,function(x) x/sum(x)) # normalize
#    sam[is.nan(sam)] <- 0
    sam
  }
  sam.unif <- function(nt,parameters,...) {
    sam <- matrix(1,ncol=nt,nrow=nt)
    sam[lower.tri(sam,diag=TRUE)] <- 0
#    sam <- apply(sam,2,function(x) x/sum(x)) # normalize
#    sam[is.nan(sam)] <- 0
    sam
  }
  fun <- switch(type,
    exponential = sam.exp,
    power = sam.pow,
    uniform = sam.unif,
    sam.unif
  )
  attr(fun,"name") <- type
  fun
}

GCM <- function(formula,level=c("nominal","interval"),distance=c("cityblock","euclidian","minkowski"),similarity=c("exponential","gaussian","general"),sampling=c("uniform","power","exponential"),parameters=list(w=NULL,lambda=1),fixed,parStruct,data,subset,ntimes=NULL,replicate=TRUE,base=NULL,remove.intercept=FALSE) {
  level <- match.arg(level)
  if(!is.function(distance)) distance <- gcm.distance(match.arg(distance))
  if(!is.function(similarity)) similarity <- gcm.similarity(match.arg(similarity))
  if(!is.function(sampling)) sampling <- gcm.sampling(match.arg(sampling))
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=remove.intercept) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=remove.intercept)
  x <- dat$x
  y <- dat$y
  nw <- ncol(x)
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$w)) {
      parameters$w <- rep(1,nw)
      parameters$w <- parameters$w/sum(parameters$w)
    } else {
      if(length(parameters$w) < nw) stop("w must have length of at least",ncol(x))
      if(length(parameters$w) > nw) warning("supplied w has more elements than required",ncol(x))
      parameters$w <- parameters$w[1:nw]
    }
    if(is.null(parameters$lambda)) parameters$lambda <- 1
    if(is.null(parameters$r)) {
      if(attr(distance,"name") == "minkowski") {
        parameters$r <- 1
      }
    }
    if(is.null(parameters$q)) {
      if(attr(distance,"name")=="minkowski") {
          parameters$q <- 1
      }
    }
    if(is.null(parameters$gamma)) {
      if(attr(sampling,"name") != "uniform") {
        if(is.null(parameters$gamma)) parameters$gamma <- 1
      }
    }
    if(level=="interval") {
      if(is.null(parameters$sdy)) parameters$sdy <- 1
    }
  } else {
    # setup a parlist
    if(length(parameters)==0) {
      nrep <- length(ntimes)
      parameters$w <- rep(1,nw)
      parameters$w <- parameters$w/sum(parameters$w)
      parameters$lambda <- 1
      if(attr(distance,"name") == "minkowski") {
        parameters$r <- 1
      }
      if(attr(distance,"name")=="general") {
          parameters$q <- 1
      }
      if(attr(sampling,"name") != "uniform") {
        if(is.null(parameters$gamma)) parameters$gamma <- 1
      }
      if(level=="interval") {
        if(is.null(parameters$sdy)) parameters$sdy <- 1
      }
      parameters <- rep(list(parameters),nrep)
    } else warning("there is no validity check for the given parameters when combined with ntimes and replicate=FALSE \n Please make sure the supplied list is valid")
  }
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
  if(missing(parStruct)) {
    parStruct <- ParStruct(parameters=parameters,replicate=replicate,
      fixed = if(missing(fixed)) NULL else fixed,
      ntimes = if(missing(ntimes)) NULL else ntimes)
  }
  
  if(level=="nominal") {
    mod <- new("GCMnominal",
      x=x,
      y=y,
      parStruct=parStruct,
      nTimes=nTimes,
      distance=distance,
      similarity=similarity,
      sampling=sampling
    )
  }
  if(level=="interval") {
    mod <- new("GCMinterval",
      x=x,
      y=y,
      parStruct=parStruct,
      nTimes=nTimes,
      distance=distance,
      similarity=similarity,
      sampling=sampling
    )
  }
  mod <- fit(mod)
  mod
}

setMethod("lFr",signature(x="GCMinterval",y="GaussianMixtureResponse"), 
  function(x,y,...) {
    for(case in 1:x@nTimes@cases) {
      y@weights[[case]] <- x@weights[[case]]
      y@means[[case]] <- t(matrix(x@y[x@nTimes@bt[case]:x@nTimes@et[case],],nrow=nrow(as.matrix(x@y[x@nTimes@bt[case]:x@nTimes@et[case],])),ncol=ncol(as.matrix(x@y[x@nTimes@bt[case]:x@nTimes@et[case],]))))
    }
    y
  }
)