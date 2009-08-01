# TODO:
# 1. allow matrix of means and sds
# 2. fix constructor function

setClass("GaussianMixtureResponse",
  contains="ResponseModel",
  representation(
    weights="list",
    means="list",
    sds="list"
  )
)

setMethod("getPars",signature(object="GaussianMixtureResponse"),
  function(object,which="all",unconstrained=FALSE,...) {
    if(unconstrained) {
      pars <- object@parStruct@parameters
      if(!is.null(pars$sd)) pars$sd <- log(pars$sd)
      pars <- unlist(pars)
      if(which=="free") {
        if(length(object@parStruct@fix)>0) {
          fixl <- relist(object@parStruct@fix,skeleton=object@parStruct@parameters)
          fixl$lambda <- NULL
          fix <- unlist(fixl)
          pars <- pars[!fix]
        }
      }
    } else {
      pars <- callNextMethod(object=object,which=which,unconstrained=unconstrained,...)
    }
    pars
  }
)

setMethod("setPars",signature(object="GaussianMixtureResponse"),
  function(object,pars,unconstrained=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    #object <- callNextMethod()
    parl <- callNextMethod(object=object,pars=pars,rval="parameters",...)
    if(unconstrained) {
      parl$sd <- exp(parl$sd)
    }
    switch(rval,
      object = {
        object@parStruct@parameters <- parl
        object
      },
      parameters = parl)
    #return(object)
  }
)
# sd matrix is a regression matrix for sd.
# If object@parStruct@parameters$sd usually contains a single value sd

setMethod("logLik",signature(object="GaussianMixtureResponse"),
  function(object,ntimes=NULL,discrete=FALSE,...) {
#    d <- object@x
#    for(i in 1:nrow(d)) d[i,] <- dnorm(object@r[i,],mean=object@x[i,],sd=object@sd[i,])
    #if(is.null(ntimes)) ntimes <- unlist(lapply(object@weights,nrow))
    #lt <- length(ntimes)
    #et <- cumsum(ntimes)
    #bt <- c(1,et[-lt]+1)
    LL <- vector("double")
    zwt <- 0
    nobs <- 0
    for(case in 1:object@nTimes@cases) {
      w <- object@weights[[case]]
      #if(is.null(discount))
      zw <- which(colSums(w)==0) #else zw <- which(colSums(w[,-discount])==0)
      #zw <- which(colSums(w)==0)
      #if(length(zw)>0)
      if(discrete) {
        d <- apply(as.matrix(object@y[object@nTimes@bt[case]:object@nTimes@et[case],]),1,function(x) pnorm(x+.5,mean=object@means[[case]],sd=object@sds[[case]]) - pnorm(x-.5,mean=object@means[[case]],sd=object@sds[[case]]))
      } else {
        d <- apply(as.matrix(object@y[object@nTimes@bt[case]:object@nTimes@et[case],]),1,function(x) dnorm(x,mean=object@means[[case]],sd=object@sds[[case]]))
      }
      d <- w*d
      #rmv <- unique(c(discount,zw))
      #d <- d[,-rmv]
      if(length(zw)>0) d <- d[,-zw]
      d <- colSums(d)
      miss <- is.na(d)
      LL[case] <- sum(log(d[!miss])) #+ length(zw)*default
      zwt <- zwt+length(zw)
      nobs <- nobs + length(d) - sum(miss)
    }
    out <- sum(LL)
    attr(out,"nobs") <- nobs
    #attr(out,"datapoints set to default (",default,")") <- zwt
#    if(!is.null(discount)) attr(out,"nobs") <- nobs
    out
  }
)

setMethod("predict",signature(object="GaussianMixtureResponse"),
  function(object,...) {
    res <- object@y
    for(case in 1:object@nTimes@cases) {
      if(all(dim(object@means[[case]])==dim(object@weights[[case]]))) res[object@nTimes@bt[case]:object@nTimes@et[case],] <- rowSums(object@weights[[case]]*object@means[[case]])
      if(nrow(object@means[[case]])==1) res[object@nTimes@bt[case]:object@nTimes@et[case],] <- as.vector(object@means[[case]]%*%object@weights[[case]])
    }
    res
  }
)
setMethod("runm",signature(object="GaussianMixtureResponse"),
  function(object,...) {
    for(case in 1:object@nTimes@cases) {
      if(object@parStruct@replicate) pars <- object@parStruct@parameters else pars <- object@parStruct@parameters[[case]]
      m <- pars$mean
      sds <- pars$sd
      nx <- ncol(object@weights[[case]])
      nt <- nrow(object@weights[[case]])
      if(!is.null(m)) {
        stop("means parameter not implemented yet")
      }
      if(!is.null(sds)) {
        if(length(sds)==1) object@sds[[case]] <- sds else {
          stop("multiple sds parameters not implemented yet")
        }
      }
    }
    
#    func <- function(pars,x,sd) {
#      m <- pars$mean
#      sigma <- pars$sd
#      if(!is.null(m)) {
#        if(!is.matrix(m)) m <- matrix(m,nrow=1)
#        if(ncol(m)!=ncol(x)) stop("mean should have as many columns as x")
#        if(nrow(m)==1) {
#          x <- t(matrix(m,nrow=ncol(x),ncol=nrow(x)))
#        } else {
#          if(nrow(m)!=nrow(x)) stop("mean should have as many rows as x")
#          x <- m
#        }
#      }
#      if(is.matrix(sigma)) {
      
      
#      } else {
#        if(length(sigma) > 1) {
#
#        } else {
#          if(length(object@sdmat)>0) {
#            sigma <- sigma*sdmat
#        }
      
#      }
#      if(!is.null(sig)) {
#        if(length(sig) == 1) {
#          sd <- matrix(sig,ncol=ncol(sd),nrow=nrow(sd))
#        } else {
#          if(!is.matrix(sig)) sig <- matrix(sig,nrow=1)
#          if(ncol(sig)==1) {
#            sd <- matrix(sig,ncol=ncol(sd),nrow=nrow(sd))
#          } else {
#            if(ncol(sig)!=ncol(sd)) stop("sd in parameters should have as many columns as sd in object")
#            if(nrow(sig)==1) {
#              sd <- t(matrix(sig,nrow=ncol(sd),ncol=nrow(sd)))
#            } else {
#              if(nrow(sig)!=nrow(sd)) stop("mean should have as many rows as x")
#              sd <- sig
#            }
#          }
#        }
#      }
#      return(list(x=x,sd=sd))
#    }
#    if(!is.null(ntimes)) {
#      lt <- length(ntimes)
#      et <- cumsum(ntimes)
#      bt <- c(1,et[-lt]+1)
#      for(case in 1:lt) {
#        tmp <- func(object@parStruct@parameters[[case]],object@x[bt[case]:et[case],],object@sd[bt[case]:et[case],])
#        object@x[bt[case]:et[case],] <- tmp$x
#        object@sd[bt[case]:et[case],] <- tmp$sd
#      }
#    } else {
#        tmp <- func(object@parStruct@parameters,object@x,object@sd)
#        object@x <- tmp$x
#        object@sd <- tmp$sd
#   }
    object
  }
)
setMethod("fit",signature(object="GaussianMixtureResponse"),
  function(object,method="Nelder-Mead",unconstrained=FALSE,...) {
    #if(!is.null(object@parStruct@parameters$sd)) object@parStruct@parameters$sd <- log(object@parStruct@parameters$sd)
    pstart <- getPars(object,which="free",unconstrained=unconstrained,...)
    optfun <- function(par,object,unconstrained,...) {
      #object <- setPars(object,par,unconstrained=unconstrained,...)
      object@parStruct@parameters <- setPars(object,par,rval="parameters",unconstrained=unconstrained,...)
      object <- runm(object,...)
      -logLik(object)
    }
    if(length(pstart)==1 & names(pstart)=="sd") {
      opt <- list()
      # estimate sd
      if(unconstrained) v <- v.old <- exp(pstart) else v <- v.old <- pstart
      converge <- FALSE
      while(!converge) {
        # compute responsibilities
        d <- object@weights*apply(object@y,1,function(y) dnorm(y,mean=object@x,sd=v))
        d[,colSums(object@weights)!=0] <- t(apply(d,1,"/",colSums(d)))[,colSums(object@weights)!=0]
#        d[,colSums(object@weights)==0] <- 0
        #d[object@weights==0] <- 0
        v <- sqrt(sum(d*apply(object@y,1,"-",object@x)^2)/sum(d))
        if(abs(v-v.old) < 1e-5) converge <- TRUE
        v.old <- v
      }
      if(unconstrained) opt$par <- log(v) else opt$par <- v
    } else {
      if(length(pstart)==1 & method=="Nelder-Mead") {
        if(is(object@parStruct@constraints,"BoxConstraintsList")) {
          mn <- object@parStruct@constraints@min
          mx <- object@parStruct@constraints@max
        } else {
          mn <- .Machine$double.eps
          mx <- sd(object@y)
        }
        if(unconstrained) {
          mn <- log(mn)
          mx <- log(mx)
        }
        try(opt <- optimise(f=optfun,interval=c(mn,mx),object=object,...))
        if(!is.null(opt)) opt$par <- opt$minimum
      } else {
          try(opt <- optim(pstart,fn=optfun,method=method,object=object,...))
      }
    }
    if(is.null(opt)) {
      warning("optimization failed")
    } else {
      #object <- setPars(object,opt$par,unconstrained=unconstrained,...)
      object@parStruct@parameters <- setPars(object,opt$par,rval="parameters",unconstrained=unconstrained,...)
    }
    #if(!is.null(object@parStruct@parameters$sd)) object@parStruct@parameters$sd <- exp(object@parStruct@parameters$sd)
    object <- runm(object,...)
    return(object)
  }
)

GaussianMixtureResponse <- function(formula,ncomponent=2,fixed,parStruct,data,subset,weights,ntimes=NULL,replicate=TRUE) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,remove.intercept=TRUE) else dat <- mcpl.prepare(formula,data,remove.intercept=TRUE)
  x <- dat$x
  y <- dat$y
  if(ncol(x)!=0) ncomponent <- ncol(x)
  
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$sd)) {
      parameters$sd <- 1
    } else {
      if(length(parameters$sd) > 1) {
        if(length(parameters$sd) != ncomponent) warning("sd should have length 1 or length(means). Elements will be recycled.")
        #parameters$sd <- rep(parameters$sd,length=ncomponent)
      }
    }
  } else {
    # setup a parlist
    if(length(parameters)==0) {
      nrep <- length(ntimes)
      parameters$sd <- 1
      parameters <- rep(list(parameters),nrep)
    } else warning("there is no validity check for the given parameters when combined with ntimes and replicate=FALSE \n Please make sure the supplied list is valid")
  }
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
  if(missing(parStruct)) {
    parStruct <- ParStruct(parameters=parameters,replicate=replicate,
      fixed = if(missing(fixed)) NULL else fixed,
      ntimes = if(missing(ntimes)) NULL else ntimes)
  }
  
  mod <- new("GaussianMixtureResponse",
    x = x,
    y = y,
    parStruct=parStruct,
    nTimes=nTimes
  )
}
