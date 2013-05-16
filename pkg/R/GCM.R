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

# TODO:
# set discount in logLik

setClass("GCM",
  contains="McplModel")

setClass("GCMlearning",
  contains="LearningModel")
  
setClass("GCMresponse",
  contains="ResponseModel")
  

setMethod("is.unconstrained",signature(object="GCM"),
  function(object,...) {
    chck <- TRUE
    if(is(object@learningModel@parStruct@constraints,"LinConstraintsList") || is(object@learningModel@parStruct@constraints,"BoxConstraintsList")) {
      chck <- FALSE
    }
    if(is(object@responseModel@parStruct@constraints,"LinConstraintsList") || is(object@responseModel@parStruct@constraints,"BoxConstraintsList")) {
      chck <- FALSE
    }
    return(chck)
  }
)

setMethod("canRepar",signature(object="GCMlearning"),
  function(object,...) {
    repar <- TRUE
    if(!is(object@parStruct@constraints,"S4")) repar <- FALSE
    #if(is(object@parStruct@constraints,"LinConstraintsList") || is(object@parStruct@constraints,"BoxConstraintsList")) {
    #  repar <- FALSE
    #}
    # TODO: check whether \lambda is fixed
    return(repar)
  }
)

setMethod("canRepar",signature(object="GCMresponse"),
          function(object,...) {
            repar <- TRUE
            if(!is(object@parStruct@constraints,"S4")) repar <- FALSE
            return(repar)
          }
)

setMethod("canRepar",signature(object="GCM"),
          function(object,...){
            return(all(canRepar(object@learningModel),canRepar(object@responseModel)))
          }
)

setMethod("runm",signature(object="GCM"),
  function(object,...) {
    object@responseModel@x <- predict(object,...)
    object
  }
)

setMethod("setTransPars",signature(object="GCMlearning"),
  function(object,pars,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
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
    switch(rval,
      object = {
        object@parStruct@parameters <- parl
        object},
      parameters = parl)
  }
)

setMethod("setTransPars",signature(object="GCMresponse"),
          function(object,pars,...,rval=c("object","parameters")) {
            sP.u <- function(pars,fixl) {
              if(!is.null(pars$gamma) && !fixl$gamma) pars$gamma <- exp(pars$gamma)
              pars
            }
            parl <- getPars(object=object,pars=pars,rval="parameters",...)
            if(length(object@parStruct@fix)>0) fixl <- relist(object@parStruct@fix,skeleton=object@parStruct@parameters) else fixl <- relist(rep(FALSE,length(unlist(object@parStruct@parameters))),skeleton=object@parStruct@parameters)
            if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
              for(case in 1:object@nTimes@cases) parl[[case]] <- sP.u(parl[[case]],fixl[[case]])
            } else {
              parl <- sP.u(parl,fixl)
            }
            switch(rval,
                   object = {
                     object@parStruct@parameters <- parl
                     object
                   },
                   parameters = parl)
          }
)


setMethod("getTransPars",signature(object="GCMlearning"),
  function(object,which="all",...) {
    gP.u <- function(pars) {
          pr <- pars$r
          pq <- pars$q
          #if(is.null(pr)) if(attr(object@distance,"name") == "euclidian") pr <- 2 else pr <- 1
          #if(is.null(pq)) if(attr(object@distance,"name") == "gaussian") pq <- 2 else pq <- 1
          pars$lambda <- pars$lambda^{pq/pr}
          pars$w <- log(pars$lambda*pars$w)
          pars$lambda <- NULL
          if(!is.null(pars$gamma)) pars$gamma <- log(pars$gamma)
          if(!is.null(pars$sdy)) pars$sdy <- log(pars$sdy)
          pars
        }
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
    return(pars)
  }
)


setMethod("getTransPars",signature(object="GCMresponse"),
          function(object,which="all",...) {
            gP.u <- function(pars) {
              pars$gamma <- log(pars$gamma)
              pars
            }
            pars <- object@parStruct@parameters
            if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
              for(case in 1:object@nTimes@cases) pars[[case]] <- gP.u(pars[[case]])
            } else {
              pars <- gP.u(pars)
            }
            pars <- unlist(pars)
            if(which=="free") {
              if(length(object@parStruct@fix)>0) {
                pars <- pars[!object@parStruct@fix]
              }
            }
            return(pars)
          }
)

setMethod("setTransPars",signature(object="GCM"),
  function(object,pars,...) {
    object@learningModel <- setTransPars(object@learningModel,pars,...)
    object@responseModel <- setTransPars(object@responseModel,pars,...)
    object
  }
)

setMethod("getTransPars",signature(object="GCM"),
          function(object,which="all",...) {
            pars <- list()
            pars[[1]] <- getTransPars(object@learningModel,which=which,...)
            pars[[2]] <- getTransPars(object@responseModel,which=which,...)
            pars <- as.relistable(pars)
            unlist(pars)
          }
)

setMethod("predict",signature(object="GCM"),
  function(object,...) {
    x <- t(object@learningModel@x)
    y <- t(object@learningModel@y)
    ny <- nrow(y)
    nx <- nrow(x)
    bt <- object@learningModel@nTimes@bt
    lt <- object@learningModel@nTimes@cases
    et <- object@learningModel@nTimes@et
    nt <- sum(object@learningModel@nTimes@n)
    parameters <- object@learningModel@parStruct@parameters
    w <- parameters$w
    r <- parameters$r
    q <- parameters$q
    lambda <- parameters$lambda
    gamma <- object@responseModel@parStruct@parameters$gamma
    dist <- vector("double",length=nt)
    sim <- vector("double",length=ny)
    ypred <- vector("double",length=ny*nt)
    
    out <- .C("gcm_nominal",
      y = as.integer(y),
      ny = as.integer(ny),
      x = as.double(x),
      nx = as.integer(nx),
      bt = as.integer(bt),
      et = as.integer(et),
      lt = as.integer(lt),
      w = as.double(w),
      r = as.double(r),
      q = as.double(q),
      lambda = as.double(lambda),
      gamma = as.double(gamma),
      dist = as.double(dist),
      sim = as.double(sim),
      ypred = as.double(ypred),
      PACKAGE="mcplR"
    )
    return(matrix(out$ypred,nrow=nt,ncol=ny,byrow=TRUE))
  }
)

setMethod("logLik",signature=(object="GCM"),
  function(object,discount=1,...) {
    discount <- unlist(lapply(object@learningModel@nTimes@bt-1,"+",discount))
    discount <- discount[discount>0]
    out <- sum(log(rowSums(object@responseModel@x[-discount,]*object@responseModel@y[-discount,])))
    nobs <- sum(object@learningModel@nTimes@n) - length(discount)
    attr(out,"nobs") <- nobs
    attr(out,"df") <- length(getPars(object,which="free"))
    class(out) <- "logLik"
    out   
  }
)

setMethod("fit",signature(object="GCM"),
  function(object,method="Nelder-Mead",...) {
    optfun <- function(pars,object,...) {
      #object <- setPars(object,pars,unconstrained=unconstrained)
      if(canRepar(object,...)) object <- setTransPars(object,pars,...) else object <- setPars(object,pars,...)
      object <- runm(object,...)
      -logLik(object)
    }
    if(!is.unconstrained(object,...)) {
      stop("Constraints are not (yet) implemented for GCM; use gGCM for (partial) support of parameter constraints")  
    } else {
      if(canRepar(object,...)) pars <- getTransPars(object,which="free",...) else pars <- getPars(object,which="free",...) 
      opt <- optim(pars,fn=optfun,method=method,object=object,...)
      if(canRepar(object,...)) object <- setTransPars(object,opt$par,...) else object <- setPars(object,opt$par,...) 
    }
    object <- runm(object,...)
    object
  }
)

GCM <- function(learning,response,parameters=list(w=NULL, lambda=1,r=1,q=1,gamma=NULL),fixed,data,subset,ntimes=NULL,replicate=TRUE,remove.intercept=FALSE) {

  if(!missing(subset)) {
    dat <- mcpl.prepare(learning,data,subset,base=NULL,remove.intercept=remove.intercept)
    rdat <- mcpl.prepare(response,data,subset,base=NULL,remove.intercept=remove.intercept)
  } else {
    dat <- mcpl.prepare(learning,data,base=NULL,remove.intercept=remove.intercept)
    rdat <- mcpl.prepare(response,data,base=NULL,remove.intercept=remove.intercept)
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
    if(is.null(parameters$r)) lpars$r <- 1 else lpars$r <- parameters$r
    if(is.null(parameters$q)) lpars$q <- 1 else lpars$q <- parameters$q
    if(is.null(parameters$gamma)) rpars$gamma <- 1 else rpars$gamma <- parameters$gamma
  } else {
    # setup a parlist
    if(length(parameters)==0) {
      nrep <- length(ntimes)
      lpars$w <- rep(1,nw)
      lpars$w <- lpars$w/sum(lpars$w)
      lpars$lambda <- 1
      lpars$r <- 1
      lpars$q <- 1
      rpars$gamma <- 1 
      lpars <- rep(list(lpars),nrep)
      rpars <- rep(list(rpars),nrep)
    } else {
      warning("there is no validity check for the given parameters when combined with ntimes and replicate=FALSE \n Please make sure the supplied list is valid")
      lpars <- parameters[["learning"]]
      rpars <- parameters[["response"]]
    }
  }
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
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
        fixed = if(missing(lfixed)) NULL else lfixed,
        ntimes = if(missing(ntimes)) NULL else ntimes
  )
  rParStruct <- ParStruct(parameters=rpars,replicate=replicate,
        fixed = if(missing(rfixed)) NULL else rfixed,
        ntimes = if(missing(ntimes)) NULL else ntimes
  )
  lmod <- new("GCMlearning",
    x=x,
    y=y,
    parStruct=lParStruct,
    nTimes=nTimes
  )
  #lmod <- runm(lmod)
  rmod <- new("GCMresponse",
    #x = predict(lmod),
    y = resp,
    parStruct=rParStruct,
    nTimes=nTimes
  )
  tmod <- new("GCM",
  learningModel = lmod,
  responseModel = rmod)
  tmod <- runm(tmod)
  tmod
}

