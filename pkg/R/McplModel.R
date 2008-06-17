setClass("McplModel",
  representation(
    learningModel="LearningModel",
    responseModel="ResponseModel"
  )
)

McplModel <- function(learningModel,responseModel) {
  mod <- new("McplModel",
    learningModel = learningModel,
    responseModel = responseModel)
  fit(mod)
}

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
    LL <- switch(from,
      responseModel = logLik(object@responseModel,...),
      learningModel = logLik(object@learningModel,...))
    attr(LL,"df") <- switch(from,
      responseModel = length(getPars(object,which="free",...)),
      learningModel = length(getPars(object@learningModel,which="free",...)))
    return(LL)
  }
)
setMethod("AIC",signature(object="McplModel"),
  function(object,npar,...,k=2) {
    if(missing(npar)) {
      AIC(logLik(object,...,k=k))
    } else {
      AIC(object=object@responseModel,k=k,npar=npar)
    }
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
