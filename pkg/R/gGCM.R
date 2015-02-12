setClass("gGCM",
  contains="McplModel")

setClass("gGCMnominal",
  contains="gGCM")
  
setClass("gGCMinterval",
  contains="gGCM")
  
setClass("gGcmLearning",
  contains="GCMlearning",
  representation(
    weights="list", # each element is a lower triangular matrix with (normalized) weight of each y[t]
    distance="function", # returns an nt*nt matrix with distances
    similarity="function", # computes similarities from distance matrix
    sampling="function" # computes nt*nt matrix with sampling weights (column wise!)
  )
)

setClass("gGcmLearningNominal",
  contains="gGcmLearning"
)
setClass("gGcmLearningInterval",
  contains="gGcmLearning"
)

setMethod("runm",signature(object="gGcmLearning"),
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        runm <- ggcm.runm(x=x,y=y,parameters=pars,distance=object@distance,similarity=object@similarity,sampling=object@sampling)
        object@weights[[case]] <- runm$weights
      }
    } else {
      runm <- ggcm.runm(x=object@x,y=object@y,parameters=object@parStruct@parameters,distance=object@distance,similarity=object@similarity,sampling=object@sampling)
      object@weights[[1]] <- runm$weights
    }
    return(object)    
  }
)

ggcm.runm <- ggcm.cont.runm <- ggcm.discr.runm <- function(x,y,parameters,distance,similarity,sampling,...) {
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

ggcm.runm.unconstrained <- function(x,y,parameters,distance,similarity,sampling,...) {
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

setMethod("predict",signature(object="gGcmLearning"),
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


setMethod("logLik",signature(object="gGcmLearningNominal"),
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

setMethod("logLik",signature(object="gGcmLearningInterval"),
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


setMethod("fit",signature(object="gGcmLearning"),
  function(object,method="Nelder-Mead",...) {
    optfun <- function(pars,object,repar,...) {
      #object <- setPars(object,pars,unconstrained=unconstrained)
      object@parStruct@parameters <- setPars(object,pars,repar=repar,...,rval="parameters")
      object <- runm(object,...)
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
    object <- runm(object,...)
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

gGCM <- function(learning,response,level=c("nominal","interval"),distance=c("cityblock","euclidian","minkowski"),similarity=c("exponential","gaussian","general"),sampling=c("uniform","power","exponential"),parameters=list(w=NULL, lambda=1,r=1,q=1,gamma=NULL,theta=NULL,sdy=NULL,sdr=NULL),fixed,data,subset,ntimes=NULL,replicate=TRUE,base=NULL,remove.intercept=FALSE) {
  level <- match.arg(level)
  if(!is.function(distance)) {
  	distance <- match.arg(distance)
  	if(!is.null(parameters$r)) {
  		dist <- switch(parameters$r,
  		  "cityblock",
  		  "euclidian")
  		if(is.null(dist)) dist <- "minkowski"
  		if(dist != distance) warning("mismatch between distance argument and parameter r; will use",dist,"distance function")
  		distance <- dist
  	}
  	distance <- gcm.distance(distance)
  }
  if(!is.function(similarity)) {
    similarity <- match.arg(similarity)
  	if(!is.null(parameters$q)) {
  		sim <- switch(parameters$q,
  		  "exponential",
  		  "gaussian")
  		if(is.null(sim)) sim <- "general"
  		if(sim != similarity) warning("mismatch between similarity argument and parameter q; will use",sim,"similarity function")
  		similarity <- sim
  	}
  	similarity <- gcm.similarity(similarity)
  }
  if(!is.function(sampling)) sampling <- gcm.sampling(match.arg(sampling))
  if(!missing(subset)) {
    dat <- mcpl.prepare(learning,data,subset,base=base,remove.intercept=remove.intercept)
    rdat <- mcpl.prepare(response,data,subset,base=base,remove.intercept=remove.intercept)
  } else {
    dat <- mcpl.prepare(learning,data,base=base,remove.intercept=remove.intercept)
    rdat <- mcpl.prepare(response,data,base=base,remove.intercept=remove.intercept)
  }
  x <- dat$x
  y <- dat$y
  resp <- rdat$y
  nw <- ncol(x)
  lpars <- list()
  rpars <- list()
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$w)) {
      lpars$w <- rep(1,nw)
      lpars$w <- lpars$w/sum(lpars$w)
    } else {
      if(length(parameters$w) < nw) stop("w must have length of at least",ncol(x))
      if(length(parameters$w) > nw) warning("supplied w has more elements than required",ncol(x))
      lpars$w <- parameters$w[1:nw]
      if(sum(lpars$w)!=1) lpars$w <- lpars$w/sum(lpars$w)
    }
    if(is.null(parameters$lambda)) lpars$lambda <- 1 else lpars$lambda <- parameters$lambda
    if(is.null(parameters$r)) {
      if(attr(distance,"name") == "minkowski") {
        lpars$r <- 1
      }
    } else lpars$r <- parameters$r
    if(is.null(parameters$q)) {
      if(attr(distance,"name")=="minkowski") {
          lpars$q <- 1
      }
    } else lpars$q <- parameters$q
    if(is.null(parameters$gamma)) {
      if(attr(sampling,"name") != "uniform") {
        if(is.null(parameters$gamma)) lpars$gamma <- 1 else lpars$gamma <- parameters$gamma
      }
    }
    if(level=="interval") {
      if(is.null(parameters$sdy)) lpars$sdy <- 1 else lpars$sdy <- parameters$sdy
      if(is.null(parameters$sdr)) rpars$sd <- 1 else rpars$sd <- parameters$sdr
    }
    if(level=="nominal") {
      if(is.null(parameters$beta)) rpars$beta <- 1 else rpars$beta <- parameters$beta
    }
  } else {
    # setup a parlist
    if(length(parameters)==0) {
      nrep <- length(ntimes)
      lpars$w <- rep(1,nw)
      lpars$w <- lpars$w/sum(lpars$w)
      lpars$lambda <- 1
      if(attr(distance,"name") == "minkowski") {
        lpars$r <- 1
      }
      if(attr(distance,"name")=="general") {
          lpars$q <- 1
      }
      if(attr(sampling,"name") != "uniform") {
        if(is.null(parameters$gamma)) lpars$gamma <- 1
      }
      if(level=="interval") {
        if(is.null(parameters$sdy)) lpars$sdy <- 1
        if(is.null(parameters$sdr)) rpars$sd <- 1
      }
      lpars <- rep(list(lpars),nrep)
      rpars <- rep(list(rpars),nrep)
    } else {
      warning("there is no validity check for the given parameters when combined with ntimes and replicate=FALSE \n Please make sure the supplied list is valid")
      lpars <- parameters[["learning"]]
      rpars <- parameters[["response"]]
    }
  }
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
#  if(missing(parStruct)) {
#    parStruct <- ParStruct(parameters=parameters,replicate=replicate,
#        fixed = if(missing(fixed)) NULL else fixed,
#        ntimes = if(missing(ntimes)) NULL else ntimes
#      )
#  } else {
#    if(!is.list(parStruct
#  }

  if(!missing(fixed)) {
    if(is.list(fixed)) {
      lfixed <- fixed[which(names(fixed) %in% names(lpars))]
      rfixed <- fixed[which(names(fixed) %in% names(rpars))]
    } else {
      if(length(fixed) != unlist(c(rpars,lpars))) stop("argument fixed does not have the correct length")
      lfixed <- fixed[1:length(unlist(lfixed))]
      rfixed <- fixed[(length(lfixed) + 1):(length(lfixed) + length(rfixed))]
    }
  } else {
    lfixed <- rep(FALSE,length(unlist(lpars)))
    rfixed <- rep(FALSE,length(unlist(rpars)))
  }

  lParStruct <- ParStruct(parameters=lpars,replicate=replicate,
        fixed = if(length(lfixed) == 0) NULL else lfixed,
        ntimes = if(missing(ntimes)) NULL else ntimes
  )
  rParStruct <- ParStruct(parameters=rpars,replicate=replicate,
        fixed = if(missing(rfixed)) NULL else rfixed,
        ntimes = if(missing(ntimes)) NULL else ntimes
  )
  if(level=="nominal") {
    lmod <- new("gGcmLearningNominal",
      x=x,
      y=y,
      parStruct=lParStruct,
      nTimes=nTimes,
      distance=distance,
      similarity=similarity,
      sampling=sampling
    )
    lmod <- runm(lmod)
    rmod <- new("ExpRatioRuleResponse",
        x = predict(lmod),
        y = resp,
#         transformation = function(object,...) {
#           beta <- getPars(object,name="beta")$beta
#           if(object@parStruct@repeated) {
#             beta <- rep(beta,each=object@nTimes@n)
#           } else {
#             if(length(beta) > 1) warning("beta parameters has length > 1; using first value only")
#             beta <- beta[1]
#           }
#           return(exp(beta*object@x))
#         },
        parStruct=rParStruct,
        nTimes=nTimes
      )
      tmod <- new("gGCMnominal",
        learningModel = lmod,
        responseModel = rmod)
  }
  if(level=="interval") {
    lmod <- new("gGcmLearningInterval",
      x=x,
      y=y,
      parStruct=lParStruct,
      nTimes=nTimes,
      distance=distance,
      similarity=similarity,
      sampling=sampling
    )
    lmod <- runm(lmod)
    rmod <- new("GaussianResponse",
        x = predict(lmod),
        y = resp,
        parStruct=rParStruct,
        nTimes=nTimes
    )
    tmod <- new("gGCMinterval",
      learningModel = lmod,
      responseModel = rmod)
  }
  tmod
}

setMethod("lFr",signature(x="gGCMinterval",y="GaussianMixtureResponse"), 
  function(x,y,...) {
    for(case in 1:x@nTimes@cases) {
      y@weights[[case]] <- x@weights[[case]]
      y@means[[case]] <- t(matrix(x@y[x@nTimes@bt[case]:x@nTimes@et[case],],nrow=nrow(as.matrix(x@y[x@nTimes@bt[case]:x@nTimes@et[case],])),ncol=ncol(as.matrix(x@y[x@nTimes@bt[case]:x@nTimes@et[case],]))))
    }
    y
  }
)

setMethod("logLik",signature(object="gGCMnominal"),
  function(object,from=c("responseModel","learningModel"),discount=1,...) {
    from <- match.arg(from)
    LL <- switch(from,
      responseModel = {
        rmod <- object@responseModel
        discount <- unlist(lapply(rmod@nTimes@bt-1,"+",discount))
        discount <- discount[discount>0]
        if(is.matrix(rmod@x)) rmod@x <- rmod@x[-discount,] else {
          if(is.vector(rmod@x)) rmod@x <- rmod@x[-discount] else warning("discount in gGCM only works when response is a vector or matrix")
        }
        if(is.matrix(rmod@y)) rmod@y <- rmod@y[-discount,] else {
          if(is.vector(rmod@y)) rmod@y <- rmod@y[-discount] else warning("discount in gGCM only works when response is a vector or matrix")
        }
        logLik(rmod,...)
      },
      learningModel = logLik(object@learningModel,discount=discount,...))
    attr(LL,"df") <- switch(from,
      responseModel = length(getPars(object,which="free",...)),
      learningModel = length(getPars(object@learningModel,which="free",...)))
    return(LL)
  })
