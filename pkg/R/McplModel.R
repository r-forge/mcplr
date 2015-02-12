setClass("McplModel",
  representation(
    learningModel="LearningModel",
    responseModel="ResponseModel"#,
    #constraints="Constraints"
  )
)

McplModel <- function(learningModel,responseModel) {
  mod <- new("McplModel",
    learningModel = learningModel,
    responseModel = responseModel)
  runm(mod)
}

setMethod("getPars",signature(object="McplModel"),
  function(object,...) {
    pars <- list()
    pars[["learning"]] <- getPars(object@learningModel,...)
    pars[["response"]] <- getPars(object@responseModel,...)
    as.relistable(pars)
  }
)

setMethod("getFreePars",signature(object="McplModel"),
  function(object,...) {
    pars <- list()
    pars[["learning"]] <- getFreePars(object@learningModel,...)
    pars[["response"]] <- getFreePars(object@responseModel,...)
    pars <- as.relistable(pars)
    unlist(pars)
  }
)

setMethod("setPars",signature(object="McplModel"),
  function(object,pars,parid=NULL,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    if(is.null(attr(pars,"skeleton"))) {
      parv <- getPars(object,which="free",internal=internal,...)
      if(length(pars)==length(parv)) {
        parl <- relist(pars,skeleton=relist(getPars(object,which="free",internal=internal,...))) #FIX ME!!
      } else {
        parv <- getPars(object,which="all",...)
        if(length(pars)==length(parv)) {
          parl <- relist(pars,skeleton=relist(getPars(object,which="all",internal=internal,...))) #FIX ME!!
        } else {
          stop("cannot relist this parameter vector; please use relistable parameter vectors.")
        }
      }
    } else {
      parl <- relist(pars)
    }
    switch(rval,
      object = {
        if(!is.null(parl[["learning"]])) object@learningModel <- setPars(object@learningModel,parl[["learning"]],internal=internal,rval=rval,...)
        if(!is.null(parl[["response"]])) object@responseModel <- setPars(object@responseModel,parl[["response"]],internal=internal,rval=rval,...)
        object
      },
      parameters = {
        outpar <- vector(mode="list")#,length=2)
        if(!is.null(parl[["learning"]])) outpar[["learning"]] <- setPars(object@learningModel,parl[["learning"]],internal=internal,rval=rval,...)
        if(!is.null(parl[["response"]])) outpar[["response"]] <- setPars(object@responseModel,parl[["response"]],internal=internal,rval=rval,...)
        outpar
      }
    )
  }
)

setMethod("setFreePars",signature(object="McplModel",pars="numeric"),
  function(object,pars,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    if(is.null(attr(pars,"skeleton"))) {
      parl <- relist(getFreePars(object,...))
    } else {
      parl <- relist(pars)
    }
    nparl <- length(parl[["learning"]])
    nparr <- length(parl[["response"]])
    if(length(pars)!=sum(c(nparl,nparr))) stop("wrong length of pars")
    if(nparl > 0) {
      object@learningModel <- setFreePars(object@learningModel,pars[1:nparl],...)
    }
    if(nparr > 0) {
      object@responseModel <- setFreePars(object@responseModel,pars[(1+nparl):(nparl+nparr)],...) 
    }
    switch(rval,
       object = {
         object
       },
       parameters = {
         object@parameters
       }
    )
  }
)

setMethod("getConstraints",signature(object="McplModel"),
  # define this for a particular McplModel to allow general constraints
  function(object,...) {
    lc <- getConstraints(object@learningModel,...)
    rc <- getConstraints(object@responseModel,...)
    ac <- combineConstraints(lc,rc)
    return(ac)
  }
)

setMethod("has.lFr",signature(object="McplModel"),
  function(object,...) {
    has.lFr(object@learningModel)
  }
)

setMethod("has.rFl",signature(object="McplModel"),
  function(object,...) {
    has.lFr(object@responseModel)
  }
)

setMethod("runm",signature(object="McplModel"),
  function(object,...) {
      if(has.runm(object@learningModel)) object@learningModel <- runm(object@learningModel,...)
      if(has.lFr(object)) object@responseModel <- lFr(object,...)
      if(has.runm(object@responseModel)) object@responseModel <- runm(object@responseModel,...)
      if(has.rFl(object)) object@learningModel <- rFl(object,...)
      return(object)
  }
)

# setMethod("logLik",signature(object="McplModel"),
#   function(object,from=c("responseModel","learningModel"),...) {
#     from <- match.arg(from)
#     LL <- switch(from,
#       responseModel = logLik(object@responseModel,...),
#       learningModel = logLik(object@learningModel,...))
#     attr(LL,"df") <- switch(from,
#       responseModel = length(getFreePars(object,...)),
#       learningModel = length(getFreePars(object@learningModel,...)))
#     return(LL)
#   }
# )

setMethod("logLik",signature=(object="McplModel"),
  function(object,discount=0,...) {
    out <- logLik(object@responseModel,discount=discount,...) 
    attr(out,"df") <- length(getFreePars(object,...))
    #class(out) <- "logLik"
    out   
  }
)


setMethod("AIC",signature(object="McplModel"),
  function(object,npar,...,k=2) {
    if(missing(npar)) {
      AIC(logLik(object,...,k=k))
    } else {
      AIC(object=object,k=k,npar=npar)
    }
  }
)
setMethod("AICc",signature(object="McplModel"),
  function(object,npar,...) {
    if(missing(npar)) npar <- length(getFreePars(object,...))
    AICc(object=object@responseModel,npar=npar,...)
  }
)
setMethod("BIC",signature(object="McplModel"),
  function(object,npar,...) {
    if(missing(npar)) npar <- length(getPars(object,which="free",...))
    BIC(object=object@responseModel,npar=npar,...)
  }
)
setMethod("Rsq",signature(object="McplModel"),
  function(object,...) {
    warning("R-Squared based on responseModel")
    Rsq(object=object@responseModel,...)
  }
)

setMethod("fit",signature(object="McplModel"),
  # 2013/05/20: delete CML methods; if needed, can be supplied in user-defined fit method...
  function(object,method="Nelder-Mead",...) {
    optfun <- function(pars,object,...) {
      # 2013/05/20: using setPars on whole object rather than below...
      #parl <- setPars(object,pars,internal=TRUE,calledBy="fit",...,rval="parameters")
      #if(!is.null(parl[[1]])) object@learningModel@parStruct@parameters <- parl[[1]]
      #if(!is.null(parl[[2]])) object@responseModel@parStruct@parameters <- parl[[2]]
      object <- setFreePars(object,pars,rval="object",...)
      object <- runm(object,...)
      
      #object@learningModel <- runm(object@learningModel,...)
      #object@responseModel <- lFr(object,...)
      #if(CML) {
      #  object@responseModel <- fit(object@responseModel,method=CML.method,...)
      #} else {
      #  object@responseModel <- runm(object@responseModel,...)
      #}
      #if(has.rFl(object)) object@learningModel <- rFl(object,...)
      out <- -logLik(object,...)
      #if(CML) attr(out,"rPars") <- getPars(object@responseModel,internal=TRUE,...)
      out
    }
    pars <- getFreePars(object,...)
    constraints <- getConstraints(object,...)
    if(is(constraints,"Unconstrained")) {
      opt <- optim(par=pars,fn=optfun,object=object,method=method,...)  
    } else if(is(constraints,"BoxContraints")) {
      opt <- optim(par=pars,fn=optfun,object=object,method="L-BFGS-B",min=constraints@min,max=constraints@max,...)  
    } else if(is(constraints,"LinearConstraints")) {
      opt <- constrOptim(theta=pars,f=optfun,grad=NULL,ui=constraints@Amat,ci=constraints@bvec,...)
    } else {
      stop("Cannot determine constraints of this MCPL model")
    }
    object <- setFreePars(object,opt$par,...,rval="object")
    #if(!is.null(parl[[1]])) object@learningModel@parStruct@parameters <- parl[[1]]
    #if(!is.null(parl[[2]])) object@responseModel@parStruct@parameters <- parl[[2]]
    #}
    object <- runm(object,...)
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
    mf$AIC=AIC(logL=mf$logLik,npar=npar,nobs=nobs,...)
    mf$AICc=AICc(logL=mf$logLik,npar=npar,nobs=nobs,...)
    cat("Model fit:\n")
    print(unlist(mf))
    cat("\n")
    #cat("\n Submodels:\n")
    #cat("\n ************** \n")
    summary(object@learningModel,fits=FALSE,...)
    #cat("\n ************** \n")
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

setMethod("simulate",signature(object="McplModel"),
  function(object,nsim=1,seed=NULL,times,...) {
    if(!missing(times)) {
      object@responseModel <- simulate(object@responseModel,nsim=nsim,seed=NULL,times,...)
      warning("simulation with a ``times'' argument may result in wrong parStructs!")
      object@learningModel@x <- object@learningModel@x[rep(times,nsim),]
      object@learningModel@y <- object@learningModel@y[rep(times,nsim),]
      object@learningModel@parStruct <- rep(object@learningModel@parStruct,times=nsim,...)
    } else {
      object@responseModel <- simulate(object@responseModel,nsim=nsim,seed=NULL,...)
      object@learningModel@x <- object@learningModel@x[rep(1:nrow(object@learningModel@x),nsim),]
      object@learningModel@y <- object@learningModel@y[rep(1:nrow(object@learningModel@y),nsim),]
      object@learningModel@parStruct <- rep(object@learningModel@parStruct,times=nsim,...)
    }
    object@learningModel@nTimes <- object@responseModel@nTimes
    object <- runm(object,...)
    object
  }
)

