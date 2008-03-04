require(nlme)

gkfs.filter <- function(y,x,A,Q,ws,Sigma,nt=nrow(x),nx=ncol(x),fam,lt=1,bt=1,et=nt) {
  # initial filter
  # computes working y on the fly...
  # used by wkfs
  wy <- y # "working" y
  II <- diag(nx)
  w <- matrix(0,nrow=nt,ncol=nx)
  L <- P <- matrix(0,nrow=nt,ncol=length(Q))
  H <- rep(0,length=nt)
  for(case in 1:lt) {
    w[bt[case],] <- A%*%ws
    P[bt[case],] <- Pi <- A%*%Sigma%*%t(A) + Q
    for(i in bt[case]:(et[case]-1)) {
      B <- matrix(x[i,],ncol=nx)
      eta <- B%*%w[i,]
      H[i] <- 1/(B%*%Pi%*%t(B) + fam$variance(fam$linkinv(eta))/fam$mu.eta(eta)^2)
      K <- Pi%*%t(B)%*%H[i]
      L[i,] <- Li <- A%*%(II - K%*%B)
      wy[i,] <- (y[i,]-fam$linkinv(eta))/fam$mu.eta(eta) + eta
      w[i+1,] <- A%*%K%*%wy[i,] + Li%*%w[i,]
      P[i+1,] <- Pi <- Li%*%Pi%*%t(A) + Q
    }
  }
  return(list(w=w,P=P,H=H,L=L,wy=wy,like=0))
}

wkfs.filter <- function(y,wy,x,A,Q,R,ws,Sigma,nt=nrow(x),nx=ncol(x),fam,lt=1,bt=1,et=nt) {
  # used by wkfs
  like <- double(nt)
  if(fam$family == "binomial") logLike <- function(y,mu) dbinom(x=y, size=1, prob=mu, log = TRUE)
  II <- diag(nx)
  w <- matrix(0,nrow=nt,ncol=nx)
  L <- P <- matrix(0,nrow=nt,ncol=length(Q))
  H <- rep(0,length=nt)
  for(case in 1:lt) {
    w[bt[case],] <- A%*%ws
    P[bt[case],] <- Pi <- A%*%Sigma%*%t(A) + Q
    for(i in bt[case]:(et[case]-1)) {
      B <- matrix(x[i,],ncol=nx)
      eta <- B%*%w[i,]
      H[i] <- 1/(B%*%Pi%*%t(B) + R[i])
      K <- Pi%*%t(B)%*%H[i]
      L[i,] <- Li <- A%*%(II - K%*%B)
      like[i] <- logLike(y[i,],fam$linkinv(eta))
      w[i+1,] <- A%*%K%*%wy[i,] + Li%*%w[i,]
      P[i+1,] <- Pi <- Li%*%Pi%*%t(A) + Q
    }
    eta <- matrix(x[et[case],],ncol=nx)%*%w[et[case],]
    like[et[case]] <- logLike(y[et[case],],fam$linkinv(eta))
  }
  return(list(w=w,P=P,H=H,L=L,like=like))
}

wkfs.smoother <- function(wy,x,A,ws,Sigma,w,P,H,L,nt=nrow(x),nx=ncol(x),fam,lt=1,bt=1,et=nt) {
  # used by dlrm
  II <- diag(nx)
  P0s <- array(,dim=c(nx,nx,lt))
  w0s <- array(,dim=c(nx,lt))
  wks <- matrix(0,nrow=nt,ncol=nx)
  Pks <- PCks <- matrix(0,nrow=nt,ncol=ncol(P))
  u <- rep(0,nx)
  U <- matrix(0,ncol=nx,nrow=nx)
  xi <- matrix(0,ncol=nx)
  Li <- Pi <- matrix(0,ncol=nx,nrow=nx)
#  Pi <- matrix(P[i,],ncol=nx)
  for(case in 1:lt) {
    for(i in et[case]:bt[case]) {
      xi <- matrix(x[i,],ncol=nx)
      Li <- matrix(L[i,],ncol=nx)
      Pi <- matrix(P[i,],ncol=nx)
      u <- as.vector(t(xi)%*%H[i]%*%(wy[i,] - xi%*%w[i,]) + t(Li)%*%u)
      U  <- t(xi)%*%H[i]%*%xi + t(Li)%*%U%*%Li
      wks[i,] <- w[i,] + matrix(P[i,],ncol=nx)%*%u
      Pks[i,] <- as.vector(Pi - Pi%*%U%*%Pi)
      ifelse(i>bt[case],PCks[i,] <- as.vector((II - Pi%*%U)%*%Li%*%Pi),PCks[i,] <- as.vector((II - Pi%*%U)%*%A%*%Sigma))
    }
    w0s[,case] <- ws + Sigma%*%t(A)%*%u
    P0s[,,case] <- Sigma - Sigma%*%t(A)%*%U%*%A%*%Sigma
  }
  return(list(wks=wks,Pks=Pks,PCks=PCks,w0s=w0s,P0s=P0s))
}

iwkfs <- function(y,x,A,Q,ws,Sigma,nt=nrow(x),nx=ncol(x),fam,lt=1,bt=1,et=nt,maxit=100,tol=1e-4) {
  # Iterated Weighted Kalman Filter and Smoother
  #
  
  # initialization with GKFS
  filt <- gkfs.filter(y=y,x=x,A=A,Q=Q,ws=ws,Sigma=Sigma,nt=nt,nx=nx,fam=fam,lt=lt,bt=bt,et=et)
  smth <- wkfs.smoother(wy=filt$wy,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,fam=fam,lt=lt,bt=bt,et=et)
  w <- smth$wks
  j <- 0
  converge <- FALSE
  while(j <= maxit && !converge) {
    eta <- rowSums(x*smth$wks)
    wy <- (y-fam$linkinv(eta))/fam$mu.eta(eta) + eta
    R <- fam$variance(fam$linkinv(eta))/fam$mu.eta(eta)^2
    filt <- wkfs.filter(y=y,wy=wy,x=x,A=A,Q=Q,R=R,ws=ws,Sigma=Sigma,nt=nt,nx=nx,fam=fam,lt=lt,bt=bt,et=et)
    smth <- wkfs.smoother(wy=wy,x=x,A=A,ws=ws,Sigma=Sigma,w=filt$w,P=filt$P,H=filt$H,L=filt$L,nt=nt,nx=nx,fam=fam,lt=lt,bt=bt,et=et)
    if(sum(abs(w-smth$wks)) < tol) converge <- TRUE
    w <- smth$wks
    j <- j+1
  }
  return(list(filt=filt,smth=smth,iterations=j,convergence = converge))
}

iwkfs.em <- function(smth,A,nt=nrow(x),nx=ncol(x),lt=1,bt=1,et=nt) {
  # used by dglm
  a <- b <- c <- d <- matrix(0,ncol=nx,nrow=nx)
  for(case in 1:lt) {
    # used by dlrm
    a <- a + smth$P0s[,,case] + smth$w0s[,case]%*%t(smth$w0s[,case])
    b <- b + matrix(smth$PCks[bt[case],],ncol=nx) + smth$wks[bt[case],]%*%t(smth$w0s[,case])
    c <- c + matrix(smth$Pks[bt[case],],ncol=nx) + smth$wks[bt[case],]%*%t(smth$wks[bt[case],])
    for(i in (bt[case]+1):et[case]) {
      a <- a + matrix(smth$Pks[i-1,],ncol=nx) + smth$wks[i-1,]%*%t(smth$wks[i-1,])
      b <- b + matrix(smth$PCks[i,],ncol=nx) + smth$wks[i,]%*%t(smth$wks[i-1,])
      c <- c + matrix(smth$Pks[i,],ncol=nx) + smth$wks[i,]%*%t(smth$wks[i,])
    }
  }

  ws <- apply(smth$w0s,1,mean)
  
  for(case in 1:lt) {
    d <- d + as.matrix(smth$w0s[,case] - ws)%*%t(as.matrix(smth$w0s[,case] - ws))
  }
  d <- d/lt
  Sigma <- apply(smth$P0s,c(1,2),mean) + d # see Max Welling, The Kalman Filter
  
  #Sigma <- apply(smth$P0s,c(1,2),mean)

  Q <- (1/nt)*(c - b%*%t(A) - A%*%t(b) + A%*%a%*%t(A))
  Q <- (t(Q)+Q)/2 # ensure symmetry
  
  A <- b%*%solve(a) # A = Phi

  return(list(A=A,Q=Q,ws=ws,Sigma=Sigma))
}

iwkfs.opt <- function(y,x,A,Q,ws,Sigma,Q.diag=FALSE,Sigma.diag=FALSE,est.ws=TRUE,est.Sigma=TRUE,est.A=TRUE,est.Q=TRUE,method="BFGS",fam,lt=1,bt=1,et=ncol(x),maxit,tol,...) {
  # Dynamic Generalized Linear Model
  # subroutine: numerical optimization

  dimA <- attr(A,"dim")
  dimws <- attr(ws,"dim")
  nx <- ncol(x)
  nt <- nrow(x)

  func <- function(par) {
    names <- names(par)
    if(length(tmp <- grep("A",names)) > 0) {
      A <- par[tmp]
      attr(A,"dim") <- dimA
    }
    if(length(tmp <- grep("Q",names)) > 0) {
      chk <- FALSE
      try({
        if(Q.diag) Q <- as.matrix(nlme::pdDiag(par[tmp])) else Q <- as.matrix(nlme::pdSymm(par[tmp]))
        chk <- TRUE
      })
      if(!chk) return(NA)
    }
    if(length(tmp <- grep("ws",names)) > 0) {
      ws <- par[tmp]
      attr(ws,"dim") <- dimws
    }
    if(length(tmp <- grep("Sigma",names)) > 0) {
      chk <- FALSE
      try({
        if(Sigma.diag) Sigma <- as.matrix(nlme::pdDiag(par[tmp])) else Sigma <- as.matrix(nlme::pdSymm(par[tmp]))
        chk <- TRUE
      })
      if(!chk) return(NA)
    }
    mod <- iwkfs(y=y,x=x,A=A,Q=Q,ws=ws,Sigma=Sigma,nt=nrow(x),nx=ncol(x),fam=fam,lt=lt,bt=bt,et=et,maxit=maxit,tol=tol)
    -sum(mod$filt$like)
  }

  start <- list()
  if(est.A) {
    start$A <- as.numeric(A)
  }
  if(est.Q) {
    if(Q.diag) start$Q <- coef(nlme::pdDiag(Q)) else start$Q <- coef(nlme::pdSymm(Q))
  }
  if(est.ws) start$ws <- as.numeric(ws)
  if(est.Sigma) {
    if(Sigma.diag) start$Sigma <- coef(nlme::pdDiag(Sigma)) else start$Sigma <- coef(nlme::pdSymm(Sigma))
  }
  fit <- optim(unlist(start),func,method=method,...)

  names <- names(fit$par)
  if(length(tmp <- grep("A",names)) > 0) {
    A <- fit$par[tmp]
    attr(A,"dim") <- dimA
  }
  if(length(tmp <- grep("Q",names)) > 0) {
    Q <- fit$par[tmp]
    if(Q.diag) Q <- as.matrix(nlme::pdDiag(Q)) else Q <- as.matrix(nlme::pdSymm(Q))
  }
  if(length(tmp <- grep("wp",names)) > 0) {
    ws <- fit$par[tmp]
    attr(ws,"dim") <- dimws
  }
  if(length(tmp <- grep("Sigma",names)) > 0) {
    Sigma <- fit$par[tmp]
    if(Sigma.diag) Sigma <- as.matrix(nlme::pdDiag(Sigma)) else Sigma <- as.matrix(nlme::pdSymm(Sigma))
  }
  return(list(A=A,Q=Q,ws=ws,Sigma=Sigma))
}

dglm <- function(formula,family,data,maxit=100,iwkfs.maxit=50,ws,Sigma,A,Q,Q.c=NULL,Sigma.c=Q.c,ntimes=NULL,tol=1e-5,iwkfs.tol=tol,est.ws=TRUE,est.Sigma=TRUE,est.A=TRUE,est.Q=TRUE,filter.only=FALSE,verbose=FALSE,criterion="logLike",method="BFGS",switch.LL=.05) {
  # Dynamic Generalized Linear Model
  #    using Iteratively Weighted Kalman filter and smoother (IWKFS) and EM/numerical optimization
  # author: M. Speekenbrink
  # version: 0.5
  # date: 20 April 2007
  # adapted from: Fahrmeir, L. & Wagenpfeil, S. (1997). Penalized likelihood estimation and iterative Kalman smoothing for non-Gaussian dynamic regression models

  # TODO: add functionality for other models than binomial
  call <- match.call()
  criterion <- agrep(criterion,c("parameter","logLike"))
  if(length(criterion)!=1) stop("supplied value for `criterion' is unclear")

  dat <- nmcpl.prepare(formula=formula,data=data)
  x <- dat$x
  y <- matrix(dat$y,ncol=1)

  nx <- ncol(x)
  nt <- nrow(x)
  fam <- family
  
  if(is.null(ntimes)) ntimes <- nt
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)
  
  if(missing(Sigma)) Sigma <- diag(nx)
  if(missing(Q)) Q <- diag(nx)
  if(missing(A)) A <- diag(nx)
  if(missing(ws)) ws <- rep(0,nx)

  if(is.null(Q.c)) Q.c <- matrix(1,ncol=ncol(Q),nrow=nrow(Q))
  if(is.null(Sigma.c)) Sigma.c <- matrix(1,ncol=ncol(Sigma),nrow=nrow(Sigma))
  # ensure upper tri is identical to lower tri
  Q.c[upper.tri(Q.c)] <- t(Q.c)[upper.tri(Q.c)] 
  Sigma.c[upper.tri(Sigma.c)] <- t(Sigma.c)[upper.tri(Sigma.c)]
  # check whether a diagonal is wanted 
  if(sum(Q.c) == sum(diag(Q.c))) Q.diag <- TRUE else Q.diag <- FALSE
  if(sum(Sigma.c) == sum(diag(Sigma.c))) Sigma.diag <- TRUE else Sigma.diag <- FALSE
  
  # initialization
  mod <- iwkfs(y=y,x=x,A=A,Q=Q,ws=ws,Sigma=Sigma,nt=nt,nx=nx,fam=fam,lt=lt,bt=bt,et=et,maxit=iwkfs.maxit,tol=iwkfs.tol)
  LL <- LL.old <- sum(mod$filt$like)
  
  j <- 0
  k <- 0
  LL.dif <- tol+3

  if(verbose) cat("Kalman filter EM \n")
  converge <- 0
  opt.ok <- FALSE
  prev.opt <- FALSE
  while(j <= maxit & converge < 1) {
    #abs(LL.dif) > tol
    # Expectation Maximisation

    j <- j+1

    if(opt.ok) {
      if(verbose) cat("starting optim \n")
      em <- iwkfs.opt(y=y,x=x,A=A,Q=Q,ws=ws,Sigma=Sigma,Q.diag=Q.diag,Sigma.diag=Sigma.diag,est.ws=est.ws,est.Sigma=est.Sigma,est.A=est.A,est.Q=est.Q,method=method,fam=fam,lt=lt,bt=bt,et=et,maxit=iwkfs.maxit,tol=iwkfs.tol)
      prev.opt <- TRUE
#      converge <- TRUE # assume convergence
    } else {
      em <- iwkfs.em(smth=mod$smth,A=A,nt=nt,nx=nx,lt=lt,bt=bt,et=et)
      prev.opt <- FALSE
#      if(k > switch.wait) opt.ok <- TRUE #numerical optimization after em step ok
#      k <- k+1
    }

    mod <- iwkfs(y=y,x=x,
      A = if(est.A) em$A else A,
      Q = if(est.Q) abs(1-Q.c)*Q + Q.c*em$Q else Q,
      ws = if(est.ws) em$ws else ws,
      Sigma = if(est.Sigma) abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma else Sigma,
      nt=nt,nx=nx,fam=fam,lt=lt,bt=bt,et=et,maxit=iwkfs.maxit,tol=iwkfs.tol)

    LL <- sum(mod$filt$like)
    LL.dif <- LL.old - LL

    if(LL.dif <= 0) {
      # likelihood increased
      if(criterion[1]==1) {
        if(all(abs(c(ws-em$ws,Sigma-(abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma),A-em$A,Q-(abs(1-Q.c)*Q + Q.c*em$Q),R-em$R)) < tol )) converge <- TRUE
      } else {
        if(abs(LL.dif) < tol) converge <- 1
      }
    
      if(est.ws) ws <- em$ws
      if(est.Sigma) Sigma <- abs(1-Sigma.c)*Sigma + Sigma.c*em$Sigma
      if(est.A) A <- em$A # A = Phi
      if(est.Q) Q <- abs(1-Q.c)*Q + Q.c*em$Q
      LL.old <- LL
      if(abs(LL.dif) < switch.LL) opt.ok <- TRUE else opt.ok <- FALSE
    } else {
      # likelihood decreased
      LL <- LL.old
      if(!prev.opt) opt.ok <- TRUE else converge <- 2 # avoid optim loop
      #converge <- TRUE
    }

    if(verbose) cat("iteration",j,": LL =",LL,":LLdif =",LL.dif,"\n")

  }
  
  # estimate model with final parameters
  mod <- iwkfs(y=y,x=x,A=A,Q=Q,ws=ws,Sigma=Sigma,nt=nt,nx=nx,fam=fam,lt=lt,bt=bt,et=et,maxit=iwkfs.maxit,tol=iwkfs.tol)

  # Number of free parameters (EM)
  npar <- 0
  if(est.ws) npar <- npar + length(ws)
  if(est.Sigma) npar <- npar + sum(diag(Sigma.c)) + sum(Sigma.c[upper.tri(Sigma.c)])
  if(est.A) npar <- npar + length(A)
  if(est.Q) npar <- npar + sum(diag(Q.c)) + sum(Q.c[upper.tri(Q.c)])

  colnames(mod$filt$w) <- colnames(mod$smth$wks) <- dat$x.names
  return(list(call=call,response=y,predictor=x,weight=mod$smth$wks,predicted=fam$linkinv(rowSums(mod$smth$wks*x)),weight.var=mod$smth$Pks[,diag(matrix(1:ncol(mod$smth$Pks),ncol=nx))],weight.filter=mod$filt$w,A=A,Q=Q,ws=ws,Sigma=Sigma,LL=LL,npar=npar,niter=j,convergence=list(LL.dif = LL.dif,code=converge)))
}
