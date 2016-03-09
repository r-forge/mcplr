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

setMethod("runm",signature(object="GCM"),
  function(object,...) {
    object@responseModel@x <- predict(object,...)
    object
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
    parameters <- getPars(object@learningModel)#@parStruct@parameters
    w <- parameters$w
    r <- parameters$r
    q <- parameters$q
    lambda <- parameters$lambda
    gamma <- getPars(object@responseModel)$gamma
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
    out <- matrix(out$ypred,nrow=nt,ncol=ny,byrow=TRUE)
    out[rowSums(out) == 0,] <- 1/ncol(out)
    # should check for 0-s
    return(out)
  }
)

setMethod("dens",signature=(object="GCM"),
  function(object,...) {
    dens(object@responseModel,...)
  }
)

setMethod("dens",signature=(object="GCMresponse"),
  function(object,...) {
    rowSums(object@x*object@y,...)
  }
)

setMethod("logLik",signature=(object="GCM"),
  function(object,discount=0,...) {
    if(discount > 0) {
      discount <- as.numeric(mapply(seq,from=object@learningModel@nTimes@bt,to=object@learningModel@nTimes@bt-1 + discount))
    }
    out <- sum(log(dens(object)[-discount]))
    nobs <- sum(object@learningModel@nTimes@n) - length(discount)
    attr(out,"nobs") <- nobs
    attr(out,"df") <- length(getFreePars(object,...))
    class(out) <- "logLik"
    out   
  }
)



GCM <- function(learning,response,parameters=list(w=NULL,lambda=1,r=1,q=1,gamma=NULL),fixed,data,subset,ntimes=NULL,replicate=TRUE,remove.intercept=FALSE) {

#   if(!missing(subset)) {
#     dat <- mcpl.prepare(learning,data,subset,base=NULL,remove.intercept=remove.intercept)
#     rdat <- mcpl.prepare(response,data,subset,base=NULL,remove.intercept=remove.intercept)
#   } else {
#     dat <- mcpl.prepare(learning,data,base=NULL,remove.intercept=remove.intercept)
#     rdat <- mcpl.prepare(response,data,base=NULL,remove.intercept=remove.intercept)
#   }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("learning", "data", "subset","remove.intercept"), names(mf), 0L)
  mf <- as.list(mf[m])
  names(mf) <- c("formula",names(mf)[-1])
  dat <- do.call("mcpl.prepare",mf)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("response", "data", "subset","remove.intercept"), names(mf), 0L)
  mf <- as.list(mf[m])
  names(mf) <- c("formula",names(mf)[-1])
  rdat <- do.call("mcpl.prepare",mf)
  
  x <- dat$x
  y <- dat$y
  resp <- rdat$y
  nw <- ncol(x)
  lparfill <- function(parameters) {
    lpars <- list()
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
    lpars <- lpars[c("w","lambda","r","q")]
    lpars
  }
  rparfill <- function(parameters) {
    rpars <- list()
    if(is.null(parameters$gamma)) rpars$gamma <- 1 else rpars$gamma <- parameters$gamma
    rpars <- rpars["gamma"]
    rpars
  }
  if(is.null(ntimes) | replicate) {
    lpars <- lparfill(parameters)
    rpars <- rparfill(parameters)
  } else {
    nrep <- length(ntimes)
    # setup a parlist
    if(all(unlist(lapply(parameters,is.list))) && length(parameters)==nrep) {
      for(i in 1:nrep) {
        lpars[[i]] <- lparfill(parameters[[i]])
        rpars[[i]] <- rparfill(parameters[[i]])
      }
    } else {
      warning("there is no validity check for the given parameters when combined with ntimes and replicate=FALSE \n Please make sure the supplied list is valid")
      lpars <- lparfill(parameters)
      rpars <- rparfill(parameters)
      #parameters <- rep(list(parameters),nrep)
    }
  }
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  if(!missing(fixed)) {
    if(is.list(fixed)) {
      lfixed <- fixed[which(names(fixed) %in% names(lpars))]
      lfixed <- unlist(fixedListToVec(lfixed,lpars))
      rfixed <- fixed[which(names(fixed) %in% names(rpars))]
      rfixed <- unlist(fixedListToVec(rfixed,rpars))
    } else {
      if(length(fixed) != unlist(c(rpars,lpars))) stop("argument fixed does not have the correct length")
      lfixed <- fixed[1:length(unlist(lfixed))]
      rfixed <- fixed[(length(lfixed) + 1):(length(lfixed) + length(rfixed))]
    }
  } else {
    lfixed <- rep(FALSE,length(unlist(lpars)))
    rfixed <- rep(FALSE,length(unlist(rpars)))
  }
  lFixList <- relist(lfixed,skeleton=lpars)
  if(all(lFixList$w == FALSE) & lFixList$lambda == FALSE) {
    # should be transformable
    pars <- lpars
    pr <- pars$r
    pq <- pars$q
    pars$lambda <- pars$lambda^{pq/pr}
    pars$w <- log(pars$lambda*pars$w)
    pars$lambda <- NULL    
    
    lParStruct <- TransParStruct(
      parameters=pars,replicate=replicate,
      fixed = if(missing(lfixed)) NULL else lfixed,
      ntimes = if(missing(ntimes)) NULL else ntimes,
#       transform = function(pars,...) {
#         #pars <- object@parameters
#         #pars <- parStruct@parameters
#         pr <- pars$r
#         pq <- pars$q
#         #if(is.null(pr)) if(attr(object@distance,"name") == "euclidian") pr <- 2 else pr <- 1
#         #if(is.null(pq)) if(attr(object@distance,"name") == "gaussian") pq <- 2 else pq <- 1
#         pars$lambda <- pars$lambda^{pq/pr}
#         pars$w <- log(pars$lambda*pars$w)
#         pars$lambda <- NULL
#         #if(!is.null(pars$gamma)) pars$gamma <- log(pars$gamma)
#         #if(!is.null(pars$sdy)) pars$sdy <- log(pars$sdy)
#         #object@parameters <- pars
#         #object
#         pars
#       }#,
      transform = function(pars,...) {
        #pars <- object@parameters
        #func <- function(pars) {
        pr <- pars$r
        pq <- pars$q
        #if(is.null(pr)) if(attr(object@distance,"name") == "euclidian") pr <- 2 else pr <- 1
        #if(is.null(pq)) if(attr(object@distance,"name") == "gaussian") pq <- 2 else pq <- 1
        pars$w <- exp(pars$w)
        pars$w[pars$w == Inf] <- 1e+90
        pars$lambda <- sum(pars$w)
        pars$w <- pars$w/pars$lambda
        pars$lambda <- pars$lambda^(pr/pq)
        #object@parameters <- pars
        #object
        pars
      }
    )
    lFixList$lambda <- NULL
    lParStruct@fix <- as.logical(unlist(lFixList))
  } else {
    stop("LinearConstraints will be implemented soon!")
  }

  rpars$gamma <- log(rpars$gamma) 
  rParStruct <- TransParStruct(parameters=rpars,
    replicate=replicate,
    fixed = if(missing(rfixed)) NULL else rfixed,
    ntimes = if(missing(ntimes)) NULL else ntimes,
#     transform = function(object,...) {
#       pars <- object@parameters
#       pars$gamma <- log(pars$gamma)
#       object@parameters <- pars
#       object
#     },
    transform = function(pars,...) {
      #pars <- object@parameters
      if(!is.null(pars$gamma)) pars$gamma <- exp(pars$gamma)
      #object@parameters <- pars
      #object
      pars
    }
  )
  lmod <- new("GCMlearning",
    x=x,
    y=y,
    parStruct=lParStruct,
    nTimes=nTimes
  )
  #lmod <- runm(lmod)
  
  rmod <- new("GCMresponse",
    x = matrix(),
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

