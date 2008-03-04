setClass("BDR",
  contains="LearningModel",
  representation(
    family="ANY",
    mu0="matrix",
    sigma0="matrix",
    Q="matrix")
)
setMethod("fit","BDR",
  function(object,ntimes=NULL,...) {
    if(!is.null(ntimes)) {
      lt <- length(ntimes)
      et <- cumsum(ntimes)
      bt <- c(1,et[-lt]+1)
      for(case in 1:lt) {
        x <- object@x[bt[case]:et[case],]
        y <- object@y[bt[case]:et[case],]
        eta <- object@parameters$eta[case,]
        alpha <- object@parameters$alpha[case,]
        beta <- object@parameters$beta[case,]
        ws <- object@parameters$ws[case,]
        fit <- bdr.fit(x=x,y=y,eta=eta,alpha=alpha,beta=beta,ws=ws,grad=object@gradient,window.size=object@window.size)
        object@weight[bt[case]:et[case],] <- fit$weight
      }
    } else {
      fit <- bdr.fit(x=object@x,y=object@y,eta=object@parameters$eta,alpha=object@parameters$alpha,beta=object@parameters$beta,ws=object@parameters$ws,grad=object@gradient,window.size=object@window.size)
      object@weight <- fit$weight
    }
    return(object)
  }
)
bdr