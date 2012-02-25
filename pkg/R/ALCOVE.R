setClass("ALCOVE",
  contains="McplModel")

setClass("ALCOVElearning",
  contains="LearningModel",
  representation(
    weights="list",
    output="matrix",
    humble="logical",
    exemplar.locations="list"
  )
)

setMethod("setPars",signature(object="ALCOVElearning"),
  function(object,pars,internal=FALSE,...,rval=c("object","parameters")) {
    rval <- match.arg(rval)
    sP.u <- function(pars,fixl) {
      if(!is.null(pars$eta_a) && !fixl$eta_a) pars$eta_a <- exp(pars$eta_a)
      if(!is.null(pars$eta_w) && !fixl$eta_w) pars$eta_w <- exp(pars$eta_w)
      if(!is.null(pars$r) && !fixl$r) pars$r <- exp(pars$r)
      if(!is.null(pars$q) && !fixl$q) pars$q <- exp(pars$q)
      if(!is.null(pars$spf) && !fixl$spf) pars$spf <- exp(pars$spf)
      pars
    }
    #object <- callNextMethod()
    parl <- callNextMethod(object=object,pars=pars,internal=internal,...,rval="parameters")
    if(length(object@parStruct@fix)>0) fixl <- relist(object@parStruct@fix,skeleton=object@parStruct@parameters) else fixl <- relist(rep(FALSE,length(unlist(object@parStruct@parameters))),skeleton=object@parStruct@parameters)
    if(internal && is.null(object@parStruct@constraints)) {
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) parl[[case]] <- sP.u(parl[[case]],fixl[[case]])
      } else {
        parl <- sP.u(parl,fixl)
      }
    }
    switch(rval,
      object = {
        object@parStruct@parameters <- parl
        object
      },
      parameters = parl)
    #object@parStruct@parameters
    #return(object)
  }
)

setMethod("getPars",signature(object="ALCOVElearning"),
  function(object,which="all",internal=FALSE,...) {
    gP.u <- function(pars) {
      pars$eta_a <- log(pars$eta_a)
      pars$eta_w <- log(pars$eta_w)
      pars$r <- log(pars$r)
      pars$q <- log(pars$q)
      pars$spf <- log(pars$spf)
      pars
    }
    if(internal && is.null(object@parStruct@constraints)) {
      pars <- object@parStruct@parameters
      if(object@nTimes@cases > 1 && !object@parStruct@replicate) {
        for(case in 1:object@nTimes@cases) pars[[case]] <- gP.u(pars[[case]])
      } else {
        pars <- gP.u(pars)
      }
      pars <- unlist(as.relistable(pars))
      if(which=="free") {
        if(length(object@parStruct@fix)>0) {
          pars <- pars[!object@parStruct@fix]
        }
      }
    } else {
      pars <- callNextMethod(object=object,which=which,internal=internal,...)
      #pars <- callNextMethod()
    }
    return(pars)
  }
)

setMethod("runm",signature(object="ALCOVElearning"),
  function(object,...) {
    if(object@nTimes@cases>1) {
      outp <- vector()
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        runm <- ALCOVE.runm(x=x,y=y,e=object@exemplar.locations[[case]],
          pars=pars,humble=object@humble,...)
        object@weights[[case]] <- runm$weights
        outp <- rbind(outp,runm$output)
      }
      object@output <- outp
    } else {
      runm <- ALCOVE.runm(x=object@x,y=object@y,e=object@exemplar.locations[[1]],
        pars=object@parStruct@parameters,humble=object@humble,...)
      object@weights[[1]] <- runm$weights
      object@output <- runm$output
    }
    return(object)    
  }
)

ALCOVE.runm <- function(y,x,e,pars,w.start,a.start,humble=TRUE,train.max=1,
train.min=-1,...) {
  # y: T*k (dummy) matrix for categories
  # x: T*M matrix of cues
  # e: N*M matrix with exemplar locations
  
  # M: stimulus dimension (number of attention nodes)
  # k: category dimension (number of output nodes)
  # T: number of trials
  # N: number of stimulus locations
  
  # eta_a: learning rate for attention weights (note: paper uses lambda_a) 
  # eta_w: learning rate for output weights (note: paper uses lambda_w)
  # r: 
  # q:
  # spf: specificity parameter
  
  call <- match.call()
  
  eta_w <- pars$eta_w
  eta_a <- pars$eta_a
  r <- pars$r
  q <- pars$q
  spf <- pars$spf
  
  nt <- nrow(y)
  y.n <- ncol(y)
  x.n <- ncol(x)
  e.n <- nrow(e)
  
  if(missing(w.start)) w.start <- matrix(0,nrow=e.n,ncol=y.n)
  if(missing(a.start)) a.start <- rep(1/x.n,x.n)

  a <- matrix(nrow=nt,ncol=x.n)
  w <- array(dim=c(nt,e.n,y.n))
  pred <- matrix(nrow=nt,ncol=y.n)
      
  at <- a.start
  wt <- w.start

  for(i in 1:nt) {
    # hidden node activation
    hid <- exp(-spf*(at%*%abs(t(e)-x[i,])^r)^(q/r))  # hid: 1*N matrix
    # output node activation
    pred[i,] <- out <- hid%*%wt # out: 1*y.n matrix
    
    # training signal
    if(humble) {
      train <- y[i,]*apply(rbind(out,train.max),2,max) + (1-y[i,])*apply(rbind(out,train.min),2,min)
    } else {
      train <- replace(y[i,],y[i,]==0,train.min)
      train <- replace(train,y[i,]==1,train.max)
    }
    
    # update
    dw <- eta_w*t(hid)%*%(train-out)
    da <- -eta_a*abs(t(e)-x[i,])%*%t(spf*(train-out)%*%t(wt)*hid)
    
    w[i,,] <- wt <- wt + dw
    a[i,] <- at <- as.vector(at + da)
  }
  return(list(call=call,response=y,predictor=x,output=pred,a=a,weights=w,control=list(eta_w=eta_w,eta_a=eta_a,specificity=spf)))   
}

setMethod("predict",signature(object="ALCOVElearning"),
  function(object,type=c("link","response"),...) {
    type <- match.arg(type)
    pred <- object@output
    if(type=="response") pred <- exp(pred)/rowSums(exp(pred))
    pred
  }
)

ALCOVElearning <- function(formula,parameters=list(eta_w=.05,eta_a=.05,r=1,q=1,spf=1),humble=TRUE,exemplar.locations,data,subset,fixed=list(r=TRUE,q=TRUE),parStruct,random.locations=FALSE,n.locations=10,base=NULL,ntimes=NULL,replicate=TRUE) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=TRUE) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=TRUE)
  x <- dat$x
  y <- dat$y
  
  parfill <- function(parameters) {
    if(is.null(parameters$eta_w)) parameters$eta_w <- .05
    if(is.null(parameters$eta_a)) parameters$eta_a <- .05
    if(is.null(parameters$r)) parameters$r <- 1
    if(is.null(parameters$q)) parameters$q <- 1
    if(is.null(parameters$spf)) parameters$spf <- 1
    # TODO: validate parameters
    parameters
  }
  if(is.null(ntimes) | replicate) {
    parameters <- parfill(parameters)
  } else {
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
    
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
  if(missing(parStruct)) {
    parStruct <- ParStruct(parameters=parameters,replicate=replicate,
      fixed = fixed,
      ntimes = if(missing(ntimes)) NULL else ntimes)
  }

  if(missing(exemplar.locations)) {
    exemplar.locations <- list()
    if(random.locations) {
      # randomly generate exemplar locations
      if(replicate) {
        tmp <- matrix(,nrow=n.locations,ncol=ncol(x))
        for(i in 1:ncol(tmp)) {
          tmp[,i] <- runif(n.locations,min=min(x[,i]),max=max(x[,i]))
        }
        for(case in 1:nTimes@cases) exemplar.locations[[case]] <- tmp
      } else {
        for(case in 1:nTimes@cases) {
          tmp <- matrix(,nrow=n.locations,ncol=ncol(x))
          for(i in 1:ncol(tmp)) {
            tmp[,i] <- runif(n.locations,min=min(x[nTimes@bt[case]:nTimes@et[case],i]),
              max=max(x[nTimes@bt[case]:nTimes@et[case],i]))
          }
          exemplar.locations[[case]] <- tmp
        }
      }
    } else {
      for(case in 1:nTimes@cases) {
        exemplar.locations[[case]] <- unique(x[nTimes@bt[case]:nTimes@et[case],])
      }
    }
  }
  
  mod <- new("ALCOVElearning",
    x = x,
    y = y,
    weights = list(length=nTimes@cases),
    humble = humble,
    exemplar.locations = exemplar.locations,
    parStruct=parStruct,
    nTimes=nTimes)
  mod <- runm(mod)
  mod
}

ALCOVE <- function(learning,response,parameters=list(eta_w=.05,eta_a=.05,r=1,q=1,spf=1),humble=TRUE,exemplar.locations,data,subset,fixed=list(r=TRUE,q=TRUE),random.locations=FALSE,n.locations=10,base=NULL,ntimes=NULL,replicate=TRUE) {
  
  if(!missing(subset)) {
    dat <- mcpl.prepare(learning,data,subset,base=base,remove.intercept=TRUE) 
    rdat <- mcpl.prepare(response,data,subset,base=base,remove.intercept=TRUE) 
  } else {
    dat <- mcpl.prepare(learning,data,base=base,remove.intercept=TRUE)
    rdat <- mcpl.prepare(response,data,base=base,remove.intercept=TRUE)
  }
  x <- dat$x
  y <- dat$y
  resp <- rdat$y
  
  parfill <- function(parameters) {
    if(is.null(parameters$eta_w)) parameters$eta_w <- .05
    if(is.null(parameters$eta_a)) parameters$eta_a <- .05
    if(is.null(parameters$r)) parameters$r <- 1
    if(is.null(parameters$q)) parameters$q <- 1
    if(is.null(parameters$spf)) parameters$spf <- 1
    if(is.null(parameters$beta)) parameters$beta <- 1
    # TODO: validate parameters
    parameters
  }
  if(is.null(ntimes) | replicate) {
    lpars <- parfill(parameters)
    rpars <- list(beta = lpars$beta)
    lpars$beta <- NULL
  } else {
    nrep <- length(ntimes)
    lpars <- rpars <- list()
    # check structure of supplied list
    if(all(lapply(parameters,is.list)) && length(parameters)==nrep) {
      for(i in 1:nrep) {
        lpars[[i]] <- parfill(parameters[[i]])
        rpars[[i]] <- list(beta = lpars[[i]]$beta)
        lpars[[i]]$beta <- NULL
      }
    } else {
      lpars <- parfill(parameters)
      rpars <- list(beta = lpars$beta)
      lpars$beta <- NULL
      lpars <- rep(list(lpars),nrep)
      rpars <- rep(list(rpars),nrep)
    }
  }
  
  if(!is.null(fixed)) {
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
  
  if(is.null(ntimes)) nTimes <- nTimes(nrow(y)) else nTimes <- nTimes(ntimes)
  
  #if(missing(parStruct)) {
  lParStruct <- ParStruct(parameters=lpars,replicate=replicate,
      fixed = lfixed,
      ntimes = if(missing(ntimes)) NULL else ntimes)
  rParStruct <- ParStruct(parameters=rpars,replicate=replicate,
      fixed = rfixed,
      ntimes = if(missing(ntimes)) NULL else ntimes)
  #}

  if(missing(exemplar.locations)) {
    exemplar.locations <- list()
    if(random.locations) {
      # randomly generate exemplar locations
      if(replicate) {
        tmp <- matrix(,nrow=n.locations,ncol=ncol(x))
        for(i in 1:ncol(tmp)) {
          tmp[,i] <- runif(n.locations,min=min(x[,i]),max=max(x[,i]))
        }
        for(case in 1:nTimes@cases) exemplar.locations[[case]] <- tmp
      } else {
        for(case in 1:nTimes@cases) {
          tmp <- matrix(,nrow=n.locations,ncol=ncol(x))
          for(i in 1:ncol(tmp)) {
            tmp[,i] <- runif(n.locations,min=min(x[nTimes@bt[case]:nTimes@et[case],i]),
              max=max(x[nTimes@bt[case]:nTimes@et[case],i]))
          }
          exemplar.locations[[case]] <- tmp
        }
      }
    } else {
      for(case in 1:nTimes@cases) {
        exemplar.locations[[case]] <- unique(x[nTimes@bt[case]:nTimes@et[case],])
      }
    }
  }
  
  lmod <- new("ALCOVElearning",
    x = x,
    y = y,
    weights = list(length=nTimes@cases),
    humble = humble,
    exemplar.locations = exemplar.locations,
    parStruct=lParStruct,
    nTimes=nTimes)
  lmod <- runm(lmod)
  
  rmod <- new("RatioRuleResponse",
    x = predict(lmod),
    y = resp,
    transformation = RatioRuleResponse.trans.exp,
    parStruct=rParStruct,
    nTimes=nTimes
  )
  
  tmod <- new("ALCOVE",
    learningModel = lmod,
    responseModel = rmod)
  tmod
}
