# learningModel
### x = cues
### y = criterion
### parameters
### contraints = list

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

setClass("ParStruct",
  representation(
    id="integer",
    fix="logical",
    constraints="ConstraintsList",
    replicate="logical" # switch to indicate identical parameters for each nrep
  )
)

makeParStruct <- function(parameters,replicate=TRUE,fixed=NULL,ntimes=NULL,
                constraints=NULL) {
    parStruct <- new("ParStruct")
    if(replicate) parStruct@replicate <- TRUE else parStruct@replicate <- FALSE
    if(is.null(ntimes) | replicate) {
      fix <- rep(FALSE,length(unlist(parameters)))
      fix <- relist(fix,skeleton=parameters)
    } else {
      fix <- rep(FALSE,length(unlist(parameters[[1]])))
      fix <- relist(fix,skeleton=parameters[[1]])
    }
    if(!is.null(fixed)) {
      for(i in 1:length(fix)) {
        if(!is.null(fixed[[names(fix)[i]]]) && fixed[[names(fix)[i]]]) fix[[i]] <- rep(TRUE,length(fix[[i]]))
      }
    }
    if(is.null(ntimes) | replicate) {
      parStruct@fix <- unlist(fix)
    } else {
      fix <- rep(list(fix,length(ntimes)))
      parStruct@fix <- unlist(fix)
    }
    parStruct
}


setClass("NTimes",
  representation(
    n="integer",
    cases="integer",
    bt="integer",
    et="integer"
  )
)

nTimes <- function(ntimes) {
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)
  out <- new("NTimes",
    n=as.integer(ntimes),
    cases=as.integer(lt),
    bt=as.integer(bt),
    et=as.integer(et)
  )
  out
}

setClass("McplBaseModel",
  representation(
    x="matrix",
    y="matrix",
    parameters="list",
    parStruct="ParStruct",
    nTimes="NTimes"
  )
)
setMethod("getPars",signature(object="McplBaseModel"),
  # Note:
  # if one element in parameters corresponding to id[i] is free, so
  #   is this parameter!
  function(object,which="all",internal=FALSE,...) {
    pars <- unlist(object@parameters)
    switch(which,
      free = {
        parid <- object@parStruct@id
        fix <- object@parStruct@fix
        if(length(parid)>0) {
          # get unique parameters
          if(length(fix)>0) parid.v <- unique(parid[!fix]) else parid.v <- unique(parid)
          parid.n <- length(parid.v)
          newpars <- vector()
          for(i in 1:parid.n) newpars[i] <- pars[which(parid==parid.v[i])[1]]
          pars <- newpars
        } else {
          if(length(fix)>0) {
            if(length(fix) == length(pars)) pars <- pars[!fix] else stop("length of fix in parStruct does not match number of parameters")
          }
        }
        if(length(pars)==0) pars <- NULL
        pars
      },
      all = pars,
      pars
    )
    return(pars)
  }
)
setMethod("setPars",signature(object="McplBaseModel"),
  function(object,pars,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    oldpars <- unlist(object@parameters)
    if(length(pars) > 0) {
      if(length(pars)!=length(oldpars)) {
        parid <- object@parStruct@id
        fix <- object@parStruct@fix
        if(length(parid)>0) {
          if(length(fix)>0) parid.v <- unique(parid[!fix]) else parid.v <- unique(parid)
          parid.n <- length(parid.v)
          if(length(pars)!=parid.n) stop("length of parid does not match length of pars")
          for(i in 1:parid.n) oldpars[which(parid==parid.v[i])] <- pars[i]
        } else {
          if(length(fix)>0) {
            if(sum(!fix)!=length(pars)) stop("parid not given and length of pars does not equal number of nonfixed parameters")
            oldpars[!fix] <- pars
          } else {
            stop("cannot work with par of this length in setPar")
          }
        }
      } else {
        oldpars <- pars
      }
      object@parameters <- relist(oldpars,skeleton=object@parameters)
    }
    switch(rval,
      object = object,
      parameters = object@parameters)
    #object <- fit(object)
    #return(object)
    #newpars <- relist(oldpars,skeleton=object@parameters)
    #newpars
  }
)
setMethod("AIC",signature(object="McplBaseModel"),
  function(object,npar,...,k=2) {
    logL <- logLik(object)
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    nobs <- attr(logL,"nobs")
    if(is.null(nobs)) nobs <- nrow(object@y)
    -2*logL + k*npar
  }
)
setMethod("AICc",signature(object="McplBaseModel"),
  function(object,npar,...) {
    logL <- logLik(object)
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    nobs <- attr(logL,"nobs")
    if(is.null(nobs)) nobs <- nrow(object@y)
    -2*logL + (2*nobs*npar)/(nobs-npar-1)
  }
)
setMethod("BIC",signature(object="McplBaseModel"),
  function(object,npar,...) {
    logL <- logLik(object)
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    nobs <- attr(logL,"nobs")
    if(is.null(nobs)) nobs <- nrow(object@y)
    -2*logL + npar*log(nobs)
  }
)
setMethod("RSquare",signature(object="McplBaseModel"),
  function(object,...) {
    p <- predict(object,type="response",...)
    SSt <- sum((object@y - mean(object@y))^2)
    SSe <- sum((object@y - p)^2)
    1-SSe/SSt
  }
)
setMethod("BIC",signature(object="missing"),
  function(object,logL,npar,nobs,...) {
    -2*logL + npar*log(nobs)
  }
)
setMethod("AIC",signature(object="missing"),
  function(object,logL,npar,nobs,...,k=2) {
    -2*logL + k*npar
  }
)
setMethod("AICc",signature(object="missing"),
  function(object,logL,npar,nobs,...) {
    -2*logL + (2*nobs*npar)/(nobs-npar-1)
  }
)

setClass("LearningModel",
  contains="McplBaseModel"
)

setClass("ResponseModel",
  contains="McplBaseModel"
)

setMethod("fit",signature(object="ResponseModel"),
  function(object,...) {
    object
  }
)

setMethod("getReplication",signature(object="McplBaseModel"),
  function(object,case,ntimes=NULL,...) {
    if(is.null(ntimes)) {
      lt <- object@nTimes@cases
      bt <- object@nTimes@bt
      et <- object@nTimes@et
    } else {
      lt <- length(ntimes)
      et <- cumsum(ntimes)
      bt <- c(1,et[-lt]+1)
    }
    x <- object@x[bt[case]:et[case],]
    y <- object@y[bt[case]:et[case],]
    if(!object@parStruct@replicate) parameters <- object@parameters[[case]] else parameters <- object@parameters
    return(list(x=x,y=y,parameters=parameters))
  }
)
setMethod("show",signature(object="McplBaseModel"),
  function(object) {
  
  }
)

setMethod("summary",signature(object="McplBaseModel"),
  function(object,fits=TRUE,...) {
    x <- object
    #cat("ResponseModel, class:",is(x),"\n\n")
    pars <- getPars(x,which="all",...)
    attr(pars,"skeleton") <- NULL
    fx <- x@parStruct@fix
    if(length(fx)>0) {
      frpars <- pars[!fx]
      fxpars <- pars[fx]
    } else {
      frpars <- pars
      fxpars <- NULL
    }
    if(length(frpars)>0) {
      cat("Estimated parameters:\n")
      print(frpars)
      cat("\n")
    }
    if(length(fxpars>0)) {
      cat("Fixed parameters:\n")
      fpars <- getPars(x,which="all",...)
      print(fpars[x@parStruct@fix])
      cat("\n")
    }
    if(fits) {
      mf <- list()
      try({
        mf$logLik=logLik(x,...)
        mf$BIC=BIC(x,...)
        mf$AICc=AICc(x,...)
      },silent=TRUE)
      if(length(mf>0)) {
        cat("Model fit:\n")
        print(unlist(mf))
        cat("\n")
      }
    }
  }
)

setMethod("estimate",signature(object="McplBaseModel"),
  function(object,method="Nelder-Mead",...) {
    MLoptfun <- function(pars,object,...) {
      object@parameters <- setPars(object,pars,...,rval="parameters",internal=TRUE)
      object <- fit(object,...)
      -logLik(object)
    }
    LSoptfun <- function(pars,object,...) {
      object@parameters <- setPars(object,pars,...,rval="parameters",internal=TRUE)
      object <- fit(object,...)
      sum((predict(object)-object@y)^2)
    }
    pars <- getPars(object,which="free",...,internal=TRUE)
    if(hasMethod("logLik",is(object))) optfun <- MLoptfun else optfun <- LSoptfun
    if(!is.null(object@parStruct@constraints)) {
      switch(is(object@parStruct@constraints),
        "LinConstraintsList" = {
          A <- object@parStruct@constraints@Amat
          b <- object@parStruct@constraints@bvec
          opt <- constrOptim(theta=pars,f=optfun,grad=NULL,ui=A,ci=b,object=object,...)
          object@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
        },
        "BoxConstraintsList" = {
          opt <- optim(pars,fn=optfun,method="L-BFGS-B",object=object,min=object@parStruct@constraints@min,max=object@parStruct@constraints@min,...)
          object@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
        },
        {
          warning("This extension of ConstraintsList is not implemented; using default optimisation.")
          opt <- optim(pars,fn=optfun,method=method,object=object,...)
          object@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
        }
       )
    } else {
      opt <- optim(pars,fn=optfun,method=method,object=object,...)
      object@parameters <- setPars(object,opt$par,...,rval="parameters",internal=TRUE)
    }
    object <- fit(object,...)
    object
  }
)
       
setMethod("summary",signature(object="LearningModel"),
  function(object,...) {
    cat("Learning Model, class:",is(object)[1],"\n\n")
    callNextMethod(object=object,...)
  }
)
setMethod("summary",signature(object="ResponseModel"),
  function(object,...) {
    cat("Response Model, class:",is(object)[1],"\n\n")
    callNextMethod(object=object,...)
  }
)
setMethod("show",signature(object="LearningModel"),
  function(object) {
    cat("Learning Model, class:",is(object)[1],"\n\n")
    callNextMethod()
  }
)
setMethod("show",signature(object="ResponseModel"),
  function(object) {
    cat("Response Model, class:",is(object)[1],"\n\n")
    callNextMethod()
  }
)
setClass("McplModel",
  representation(
    learningModel="LearningModel",
    responseModel="ResponseModel"
  )
)

setMethod("getPars",signature(object="McplModel"),
  function(object,which="all",internal=FALSE,...) {
    pars <- list()
    pars[[1]] <- getPars(object@learningModel,which=which,internal=internal,...)
    pars[[2]] <- getPars(object@responseModel,which=which,internal=internal,...)
    pars <- as.relistable(pars)
    unlist(pars)
  }
)
setMethod("setPars",signature(object="McplModel"),
  function(object,pars,parid=NULL,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    if(is.null(attr(pars,"skeleton"))) {
      parv <- getPars(object,which="free",...)
      if(length(pars)==length(parv)) {
        parl <- relist(pars,skeleton=relist(getPars(object,which="free",...))) #FIX ME!!
      } else {
        parv <- getPars(object,which="all",...)
        if(length(pars)==length(parv)) {
          parl <- relist(pars,skeleton=relist(getPars(object,which="all",...))) #FIX ME!!
        } else {
          stop("cannot relist this parameter vector; please use relistable parameter vectors.")
        }
      }
    } else {
      parl <- relist(pars)
    }
    if(length(parl)==1) length(parl) <- 2 # WATCH ME!
    switch(rval,
      object = {
        if(!is.null(parl[[1]])) object@learningModel <- setPars(object@learningModel,parl[[1]],internal=internal,rval=rval,...)
        if(!is.null(parl[[2]])) object@responseModel <- setPars(object@responseModel,parl[[2]],internal=internal,rval=rval,...)
        object
      },
      parameters = {
        outpar <- vector(mode="list",length=2)
        if(!is.null(parl[[1]])) outpar[[1]] <- setPars(object@learningModel,parl[[1]],internal=internal,rval=rval,...)
        if(!is.null(parl[[1]])) outpar[[2]] <- setPars(object@responseModel,parl[[2]],internal=internal,rval=rval,...)
        outpar
      })
    #object@learningModel <- setPars(object@learningModel,parl[[1]],rval="object",...)
    #object@responseModel <- setPars(object@responseModel,parl[[2]],...)
    #object
    #outpar <- list()
    #outpar[[1]] <- setPars(object@learningModel,parl[[1]],...)
    #outpar[[2]] <- setPars(object@responseModel,parl[[2]],...)
    #outpar
    #if(!is.relistable(object@learningModel@parameters)) stop("parameters in learning model not relistable")
    #if(!is.relistable(object@responseModel@parameters)) stop("parameters in response model not relistable")
#    oldpars <- list()
#    oldpars[[1]] <- unlist(object@learningModel@parameters)
#    oldpars[[2]] <- unlist(object@responseModel@parameters)
#    newpars <- unlist(as.relistable(oldpars))
#    if(length(pars)!=length(oldpars)) {
#      if(!is.null(parid)) {
#        if(length(parid)!=length(pars)) stop("length of parid does not match length of pars")
#        newpars[parid] <- pars
#      } else {
#        fix <- c(object@learningModel@parStruct@fix,object@responseModel@parStruct@fix)
#        if(sum(!fix)!=length(pars)) stop("parid not given and length of pars does not equal number of nonfixed parameters")
#        newpars[!fix] <- pars
#      }
#    } else {
#      skel <- attr(oldpars,"skeleton")
#      newpars <- pars
#      attr(newpars,"skeleton") <- skel
#    }
#    newpars <- relist(newpars)
#    attr(newpars[[1]],"skeleton") <- attr(oldpars[[1]],"skeleton")
#    attr(newpars[[2]],"skeleton") <- attr(oldpars[[2]],"skeleton")
#    object@learningModel@parameters <- relist(newpars[[1]])
#    object@responseModel@parameters <- relist(newpars[[2]])
  }
)
setMethod("lFr",signature(x="McplModel",y="missing"),
  function(x,y,...) {
    lFr(x@learningModel,x@responseModel,...)
  }
)
setMethod("rFl",signature(x="McplModel",y="missing"),
  function(x,y,...) {
    rFl(x@responseModel,x@learningModel,...)
  }
)
setMethod("lFr",signature(x="LearningModel",y="ResponseModel"),
  function(x,y,...) {
    y@x <- predict(x,type="link",...)
    y
  }
)
setMethod("rFl",signature(x="ResponseModel",y="LearningModel"),
  function(x,y,...) {
    y
  }
)
setMethod("fit",signature(object="McplModel"),
  function(object,lfr=TRUE,rfl=FALSE,...) {
      object@learningModel <- fit(object@learningModel,...)
      if(lfr) object@responseModel <- lFr(object,...)
      object@responseModel <- fit(object@responseModel,...)
      if(rfl) object@learningModel <- rFl(object,...)
      return(object)
  }
)

setMethod("logLik",signature(object="McplModel"),
  function(object,from=c("responseModel","learningModel"),...) {
    from <- match.arg(from)
    switch(from,
      responseModel = logLik(object@responseModel,...),
      learningModel = logLik(object@learningModel,...))
  }
)
setMethod("AIC",signature(object="McplModel"),
  function(object,npar,...,k=2) {
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    AIC(object=object@responseModel,k=2,npar=npar)
  }
)
setMethod("AICc",signature(object="McplModel"),
  function(object,npar,...) {
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    AICc(object=object@responseModel,npar=npar,...)
  }
)
setMethod("BIC",signature(object="McplModel"),
  function(object,npar,...) {
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    BIC(object=object@responseModel,npar=npar,...)
  }
)
setMethod("RSquare",signature(object="McplModel"),
  function(object,...) {
    warning("RSquare based on responseModel")
    RSquare(object=object@responseModel,...)
  }
)

setMethod("estimate",signature(object="McplModel"),
  function(object,method="Nelder-Mead",CML=FALSE,CML.method=method,...) {
    optfun <- function(pars,object,CML,CML.method,...) {
      if(CML) {
        #object@learningModel <- setPars(object@learningModel,pars,...)
        object@learningModel@parameters <- setPars(object@learningModel,pars,internal=TRUE,...,rval="parameters")
      #} else object <- setPars(object,pars,...)
      } else {
        parl <- setPars(object,pars,internal=TRUE,rval="parameters",...)
        if(!is.null(parl[[1]])) object@learningModel@parameters <- parl[[1]]
        if(!is.null(parl[[2]])) object@responseModel@parameters <- parl[[2]]
      }
      object@learningModel <- fit(object@learningModel,...)
      object@responseModel <- lFr(object,...)
      if(CML) object@responseModel <- estimate(object@responseModel,method=CML.method,...) else object@responseModel <- fit(object@responseModel,...)
      object@learningModel <- rFl(object,...)
      out <- -logLik(object,...)
      if(CML) attr(out,"rPars") <- getPars(object@responseModel,internal=TRUE,...)
      out
    }
    if(CML) {
      pars <- getPars(object@learningModel,which="free",internal=TRUE,...)
    } else pars <- getPars(object,which="free",internal=TRUE,...)
    # TODO: get and use contraints!
    opt <- optim(par=pars,fn=optfun,object=object,CML=CML,CML.method=CML.method,...)
    if(CML) {
      #object@learningModel <- setPars(object@learningModel,opt$par,...)
      object@learningModel@parameters <- setPars(object@learningModel,opt$par,internal=TRUE,...,rval="parameters")
      tmp <- optfun(pars=opt$par,object=object,CML=CML,CML.method=CML.method,...)
      #object@responseModel <- setPars(object@responseModel,attr(tmp,"rPars"),...)
      object@responseModel@parameters <- setPars(object@responseModel,attr(tmp,"rPars"),internal=TRUE,...,rval="parameters")
    #} else object <- setPars(object,opt$par,...)
    } else {
      parl <- setPars(object,opt$par,internal=TRUE,rval="parameters",...)
      if(!is.null(parl[[1]])) object@learningModel@parameters <- parl[[1]]
      if(!is.null(parl[[2]])) object@responseModel@parameters <- parl[[2]]
    }
    object <- fit(object,...)
    object
  }
)

setMethod("predict",signature(object="McplModel"),
  function(object,type="link",from=c("response","learning"),...) {
    from <- match.arg(from)
    switch(from,
      response = predict(object@responseModel,type=type,...),
      learning = predict(object@learningModel,type=type,...))
  }
)

setMethod("summary",signature(object="McplModel"),
  function(object,...) {
    cat("Mcpl Model, class:",is(object),"\n\n")
    mf <- list()
    nobs <- sum(object@responseModel@nTimes@n)
    npar <- length(unlist(getPars(object,which="free",...)))
    mf$logLik=logLik(object,...)
    mf$BIC=BIC(logL=mf$logLik,npar=npar,nobs=nobs,...)
    mf$AICc=AICc(logL=mf$logLik,npar=npar,nobs=nobs,...)
    cat("Model fit:\n")
    print(unlist(mf))
    cat("\n Submodels:\n")
    summary(object@learningModel,fits=FALSE,...)
    summary(object@responseModel,fits=FALSE,...)
  }
)


setMethod("show",signature(object="McplModel"),
  function(object) {
    cat("Mcpl Model, class:",is(object),"\n\n")
    show(object@learningModel)
    show(object@responseModel)
  }
)





