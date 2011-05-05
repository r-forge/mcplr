setClass("BLF",
  contains="McplModel",
  validity = function(object) {
    if(is(object@learningModel,"BLFlearn") && is(object@responseModel,"GaussianResponse")) ret <- TRUE else ret <- FALSE
  ret 
  }
)



BLF_filter <- function(y,x,A,Q,R,ws,Sigma) {
  # used by dlrm
  nt <- nrow(y)
  ny <- ncol(y)
  II <- diag(nx)
  w <- matrix(0.0,nrow=nt,ncol=nx)
  L <- P <- array(0.0,dim=c(ncol(Q),ncol(Q),nt))
  Var <- H <- onestep <- rep(0.0,length=nt)
  B <- matrix(0.0,ncol=nx)
  K <- matrix(0.0,ncol=ny,nrow=nx)
  like <- 0.0
    w[1,] <- A%*%ws
    P[,,1] <- A%*%Sigma%*%t(A) + Q
    Var[1] <- t(x[1,])%*%P[,,1]%*%x[1,] + R
    for(i in 1:(nt-1)) {
      B <- matrix(x[i,],ncol=nx)
      #t_P <- matrix(P[i,],ncol=nx)
      H[i] <- 1/Var[i]
      K <- P[,,i]%*%t(B)%*%H[i]
      L[,,i] <- A%*%(II - K%*%B)
      #L[i,] <- as.vector(t_L)

      onestep[i] <- B%*%as.matrix(w[i,])
      #if(!any(is.na(y[i,]))) like <- like + dnorm(y[i,],onestep[i],sqrt(Var[i]),log=TRUE)

      if(!any(is.na(y[i,]))) w[i+1,] <- A%*%K%*%y[i,] + L[,,i]%*%w[i,] else w[i+1,] <- A%*%w[i,]
      if(!any(is.na(y[i,]))) P[,,i+1] <- L[,,i]%*%P[,,i]%*%t(A) + Q else P[i+1,] <- A%*%P[,,i]%*%t(A) + Q
      Var[i+1] <- t(x[i+1,])%*%P[,,i+1]%*%x[i+1,] + R
    }
    B <- matrix(x[nt,],ncol=nx)
    H[nt] <- 1/Var[nt]
    K <- P[,,nt]%*%t(B)%*%H[nt]
    L[,,nt] <- A%*%(II - K%*%B)
    
    onestep[nt] <- t(x[nt,])%*%w[nt,]
  }
  return(list(w=w,P=P,H=H,L=L,onestep=onestep,like=like,var=Var))
}

