setClass("ALCOVE",
  contains="LearningModel",
  representation(
    weights="list",
    humble="logical",
    exemplar.locations="list"
  )
)

setMethod("fit",signature(object="ALCOVE"),
  function(object,...) {
    if(object@nTimes@cases>1) {
      for(case in 1:object@nTimes@cases) {
        repl <- getReplication(object,case=case,...)
        x <- repl$x
        y <- repl$y
        pars <- repl$parameters
        fit <- ALCOVE.fit(x=x,y=y,parameters=pars,distance=object@distance,similarity=object@similarity)
        object@weights[[case]] <- fit$weights
      }
    } else {
      fit <- ALOVE.fit(x=object@x,y=object@y,parameters=object@parameters,distance=object@distance,similarity=object@similarity,sampling=object@sampling)
      object@weights[[1]] <- fit$weights
    }
    return(object)    
  }
)

ALCOVE.fit <- function(y,x,e,eta_w,eta_a,r=1,q=1,spf=1,lambda=1,w.start,a.start,humble=TRUE,train.max=1,train.min=-1,...) {
    # y: T*k (dummy) matrix for categories
    # x: T*M matrix of cues
    # e: N*M matrix with exemplar locations
    
    # M: stimulus dimension (number of attention nodes)
    # k: category dimension (number of output nodes)
    # T: number of trials
    # N: number of stimulus locations
    
    # eta_a: learning rate for attention weights (note: paper uses lambda_a) 
    # eta_w: learning rate for output weights (note: paper uses lambda_w)
    
    call <- match.call()
    
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
       
    return(list(call=call,response=y,predictor=x,predict=exp(lambda*pred)/rowSums(exp(lambda*pred)),output=pred,a=a,w=w,control=list(eta_w=eta_w,eta_a=eta_a,specificity=spf)))
        
}

alcove <- function(formula,parameters=list(eta_w=.05,eta_a=.05,r=1,q=1,spf=1),humble=TRUE,exemplar.locations,random.locations=TRUE,n.locations=10,fixed,parStruct,data,subset,ntimes=NULL,replicate=TRUE,base=NULL) {
  if(!missing(subset)) dat <- mcpl.prepare(formula,data,subset,base=base,remove.intercept=TRUE) else dat <- mcpl.prepare(formula,data,base=base,remove.intercept=TRUE)
  x <- dat$x
  y <- dat$y
  
  if(is.null(ntimes) | replicate) {
    if(is.null(parameters$eta_w) parameters$eta_w <- .05
    if(is.null(parameters$eta_a) parameters$eta_a <- .05
    if(is.null(parameters$r)) parameters$r <- 1
    if(is.null(parameters$q)) parameters$q <- 1
    if(is.null(parameters$spf)) parameters$spf <- 1
  }
  
  

}