# GLM Response Model
### A response model which can be fitted by glm.fit
setClass("GlmResponse",
  contains="ResponseModel",
  representation(
    family="ANY",
    formula="formula",
    data = "list",
    covariate = "logical",
    LearnModPredName = "character"
    #sigma="matrix"
  )
)

setMethod("fit",signature(object="GlmResponse"),
	function(object,...) {
    if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
      cat("herro!\n")
      if(length(object@parStruct@id)>0) stop("Constraints in GlmResponse models are currently not implemented")
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        runm <- glm.fit(x=x,y=y,family=object@family,start=pars$coefficients)
        pars$coefficients <- runm$coefficients
        if(object@family$family=="gaussian") pars$sd <- sqrt(sum((object@y - predict(object))^2)/(length(object@y)-1))
        object@parStruct@parameters[[case]] <- pars
      }
    } else {
      # y = response
      pars <- object@parStruct@parameters
      runm <- glm.fit(x=object@x,y=object@y,family=object@family,start=pars$coefficients)
      pars$coefficients <- runm$coefficients
      if(object@family$family=="gaussian") pars$sd <- sqrt(sum((object@y - predict(object))^2)/(length(object@y)-1))
      object@parStruct@parameters <- setPars(object,unlist(pars),rval="parameters",...)
    }
    #object <- setpars(object,unlist(pars))
    object <- runm(object,...)
    return(object)
	}
)
setMethod("predict","GlmResponse",
	function(object,type="link") {
    if(object@nTimes@cases > 1 & !object@parStruct@replicate) {
      mu <- vector()
      for(case in 1:object@nTimes@cases) {
        if(NCOL(object@x) > 1) mu <- rbind(mu,object@x[object@nTimes@bt[case]:object@nTimes@et[case],]%*%as.matrix(object@parStruct@parameters[[case]]$coefficients)) else mu <- c(mu,object@parStruct@parameters[[case]]$coefficients*object@x[object@nTimes@bt[case]:object@nTimes@et[case],])
      }
    } else {
    # y = response
      if(NCOL(object@x) > 1) mu <- object@x%*%as.matrix(object@parStruct@parameters$coefficients) else mu <- object@parStruct@parameters$coefficients*object@x
    }
    if(type=="link") return(mu) else {
      if(type=="response") {
        return(object@family$linkinv(mu))
      }
    }
	}
)
setMethod("logLik","GlmResponse",
  function(object) {
    nt <- NROW(object@y)
    LL <- switch(object@family$family,
      gaussian = {
        mu <- predict(object,type="response")
        sum(dnorm(x=object@y,mean=mu,sd=object@parStruct@parameters$sd,log=TRUE))
      },
      binomial = {
        p <- predict(object,type="response")
        sum(dbinom(x=object@y,size=1,prob=p,log=TRUE))
      },
      poisson = {
        lambda <- predict(object,type="response")
        sum(dpois(x=object@y,lambda=lambda,log=TRUE))
      },
      Gamma = {
        shape <- predict(object,type="response")
        sum(dgamma(x=object@y,shape=shape,log=TRUE))
      },
      stop("family",object@family$family,"not implemented (yet)")
    )
    attr(LL,"nobs") <- nt
    attr(LL,"df") <- length(getPars(object,which="free"))
    class(LL) <- "logLik"
    LL
  }
)

setMethod("simulate",signature(object="GlmResponse"),
	function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) {
      old.seed <- .Random.seed 
      set.seed(seed)
    }
    if(missing(times)) {
      pr <- predict(object,type=response)
    } else {
      pr <- predict(object,type=response)[times,]
    }
    nt <- nrow(pr)
    response <- switch(object@family$family,
      gaussian = {
        rnorm(nt*nsim,mean=pr,sd=object@parStruct@parameters$sd)
      },
      binomial = {
        if(NCOL(object@y) == 2) {
    			rbinom(nt*nsim,size=object@y[,2],prob=pr)
    		} else {
    			rbinom(nt*nsim,size=1,prob=pr)
    		}
      },
      poisson = {
        rpois(nt*nsim,lambda=pr)
      },
      Gamma = {
        rgamma(nt*nsim,shape=pr)
      },
      stop("family",object@family$family,"not implemented (yet)")
    )
		#if(nsim > 1) response <- matrix(response,ncol=nsim)
		#response <- as.matrix(response)
		
		object@y <- as.matrix(response)
		if(!missing(times)) {
      object@x <- object@x[rep(times,nsim),]
      ntim <- rep(0,length=object@nTimes@cases)
  		for(i in 1:length(ntim)) {
  		  ntim[i] <- sum(seq(object@nTimes@bt[i],object@nTimes@et[i]) %in% times)
      }
      warning("simulation with a times argument may result in wrong parStruct argument; please check parameters.")
      object@parStruct <- rep(object@parStruct,times=nsim)
    } else {
      object@x <- object@x[rep(1:nrow(object@x),nsim),]
      ntim <- object@nTimes@n
      object@parStruct <- rep(object@parStruct,times=nsim)
    }
    ntim <- rep(ntim,nsim)
    object@nTimes <- nTimes(ntim)
    if(!is.null(seed)) {
      set.seed(old.seed)
    }
		return(object)
	}
)

setMethod("lFr",signature(x="LearningModel",y="GlmResponse"),
  function(x,y,...) {
    if(y@covariate) {
      assign(as.character(y@LearnModPredName),predict(x,type="link",...))
      form <- y@formula
      form[[2]] <- NULL
      y@x <- model.matrix(object=form,data=y@data,...)
    } else {
      y@x <- predict(x,type="link",...)
    }
    #y@x <- eval(y@x.call,envir=list(parent.frame(),y@data)) #predict(x,type="link",...)
    y
  }
)

GlmResponse <- function(formula,data,family=gaussian(),learnpredname=NULL,parameters=list(),ntimes=NULL,replicate=TRUE,fixed,base=NULL,
                        parStruct,subset) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  #mf[[2]] <- eval(substitute(mf$formula, list(as.name(predname) = as.name("learningModelPrediction"))))
  mf <- eval(mf, parent.frame())
  
  if(!missing(subset)) dat <- mcpl.prepare(formula=formula,data=data,subset=subset,base=base) else
    dat <- mcpl.prepare(formula=formula,data=data,base=base)
  y <- dat$y
  x <- dat$x
  
  data <- as.data.frame(mf[,colnames(mf)!=learnpredname])
  
  parfill <- function(parameters) {
    if(family$family=="gaussian") if(is.null(parameters$sd)) parameters$sd <- 1 else parameters$sd <- parameters$sd
    if(is.null(parameters$coefficients)) parameters$coefficients <- rep(0,NCOL(x))
    parameters
  }
  
  if(is.null(ntimes) | replicate) {
    # intialize
    parameters <- parfill(parameters)
  } else {
    # setup a parlist
    nrep <- length(ntimes)
    # check structure of supplied list
    if(all(lapply(parameters,is.list)) && length(parameters)==nrep) {
      for(i in 1:nrep) {
        parameters[[i]] <- parfill(parameters[[i]])
      }
    } else {
      parameters <- parfill(parameters)
      parameters <- rep(list(parameters),nrep)
    }
  }

  if(missing(parStruct)) {
    tfix <- NULL
    if(!missing(fixed)) tfix <- fixed
    parStruct <- ParStruct(parameters,replicate=replicate,
                    fixed=tfix,ntimes=ntimes)
  }

  if(is.null(ntimes)) ntimes <- nrow(y)
  nTimes <- nTimes(ntimes)

  if(is.null(learnpredname)) {
    covariate <- FALSE
    learnpredname <- ""
  } else {
    covariate <- TRUE
  }
  mod <- new("GlmResponse",
    x = x,
    y = y,
    parStruct=parStruct,
    nTimes=nTimes,
    family=family,
    formula=formula,
    LearnModPredName=learnpredname,
    covariate=covariate,
    data=data
  )
  mod <- runm(mod)
  mod
}
