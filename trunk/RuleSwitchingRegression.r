rule.fb <- function(A,B,prior) {

    # Forward-Backward algorithm (used in Baum-Welch)
    # Returns alpha, beta, and full data likelihood
    # A = K*K matrix with transition probabilities, from row to column
    # B = T*K matrix with elements ab_{ij} = P(y_i|s_j)
    # pi = K vector with prior probabilities
    
    # NOTE: to prevent underflow, alpha and beta are scaled, using c
    
    nt <- nrow(B)
    ns <- ncol(A)
    alpha <- matrix(ncol=ns,nrow=nt)
    beta <- matrix(ncol=ns,nrow=nt)
    c <- vector(length=nt)
    
    alpha[1,] <- prior*B[1,] # initialize
    c[1] <- 1/sum(alpha[1,])
    alpha[1,] <- alpha[1,]*c[1]
    for(i in 1:(nt-1)) {
        alpha[i+1,] <- (t(A)%*%alpha[i,])*B[i+1,]
        c[i+1] <- 1/sum(alpha[i+1,])
        alpha[i+1,] <- c[i+1]*alpha[i+1,] 
    }
    
    beta[nt,] <- 1*c[nt] # initialize
    for(i in (nt-1):1) {
        beta[i,] <- (A%*%(B[i+1,]*beta[i+1,]))*c[i]
    }
    
    gamma <- alpha*beta/c
    
    xi <- array(dim=c(nt-1,ns,ns))
    for(i in 1:(nt-1)) {
            #xi[i,,] <- (gamma[i,]%*%t(B[i+1,]))*A*((1/fbo$beta[i,])%*%t(fbo$beta[i+1,]))
        xi[i,,] <- (alpha[i,]%*%t(B[i+1,]*beta[i+1,]))*A
    }
    
    #like <- as.numeric((beta[1,]*prior)%*%B[1,])
    like <- -sum(log(c))
    return(list(alpha=alpha,beta=beta,gamma=gamma,xi=xi,c=c,logLike=like))
}

rule.viterbi <- function(A,B,prior) {
    # returns the most likely state sequence
    nt <- nrow(B)
    ns <- ncol(A)
    delta <- psi <- matrix(nrow=nt,ncol=ns)
    state <- vector(length=nt)
    # initialization
    delta[1,] <- - (log(prior) + log(B[1,]))
    psi[1,] <- 0
    # recursion
    for(i in 2:nt) {
        for(j in 1:ns) {
            delta[i,j] <- min(delta[i-1,] - log(A[,j])) - log(B[i,j])
            k <- which.min(delta[i-1,] - log(A[,j]))
            if(length(k) == 0) k <- 0
            psi[i,j] <- k
        }
    }
    #trace maximum likely state
    state[nt] <- which.min(delta[nt,])
    for(i in (nt-1):1) {
        state[i] <- psi[i+1,state[i+1]]
    }
    return(state)
}

rule.BW <- function(y,X,Z,prior,A,b,tol=1e-4,maxiter=200,A.est,prior.est,b.est=seq(1,length(y)),verbose=FALSE) {
    #
    # Baum-Welch (EM) algorithm for ML estimates of the rule learning model
    # (c) 2006, M. Speekenbrink
    #
    # y = list of lenght N with response variable
    # X = list of length N, each element a T_i*k design matrix (predictors)
    # Z = S*k proportional weight matrix, 
    # prior = S vector with prior probabilities
    #
    # A.est: matrix of dim(A.est)=dim(A) with integers 0,1,... where 0 indicates fixed and 1,2,... indicate unique free parameters
    # prior.est: similar indicator vector as for A.est
    # b.est: N vector with integers 0,1,... where 0 indicates fixed and 1,2,... indicate unique free parameters
    #        defaults to (1,2,...,N)
    ni <- length(y)
    ns <- nrow(Z)
    if(length(b) < ni) b <- c(b,rep(b[length(b)],ni-length(b)))
    
    if(missing(A.est)) {
        estA <- FALSE
        A.id <- 1:length(A)
    } else {
        estA <- TRUE
        A.id <- unique(as.numeric(A.est[A.est!=0]))
        A.fix <- A[A.est==0]
    }
    if(missing(prior.est)) {
        estprior <- FALSE
        prior.id <- 1:length(prior)
    } else {
        estprior <- TRUE
        prior.id <- unique(as.numeric(prior.est[prior.est!=0]))
        prior.fix <- prior[prior.est==0]
    }
#    if(missing(b.est)) {
#        estb <- FALSE
#    } else {
        if(length(b.est)!=length(b)) stop("b.est must have length",length(b))
#        if(b.equal) warning("b.est supplied but ignored due to b.equal") 
#        estb <- TRUE
        b.id <- order(unique(as.numeric(b.est[b.est!=0])))
        b.fix <- b[b.est==0]
#    }
#    if(b.equal) {
#        estb <- TRUE
#        b.est <- rep(1,length(b))
#    }
    dat <- fbo <- B <- list()
    LL <- nt <- vector(length=ni)
    
    A_num <- array(dim=c(ni,nrow(A),ncol(A)))
    A_denom <- matrix(nrow=ni,ncol=ncol(A))
    
    for(i in 1:ni) {
        nt[i] <- nrow(X[[i]])
        dat[[i]] <- data.frame(cbind(y=y[[i]],x=as.numeric(X[[i]]%*%t(Z))))   
        py <- matrix(1/(1+exp(-b[i]*dat[[i]]$x)),ncol=ns)
        B[[i]] <- (py^y[[i]])*(1-py)^(1-y[[i]])    
        fbo[[i]] <- rule.fb(A=A,B=B[[i]],prior=prior)
        LL[i] <- fbo[[i]]$logLike
    }
    
    LL.dif <- tol+1
    j <- 0
    while(j <= maxiter && abs(LL.dif) > tol) {
        LL.old <- LL
        j <- j+1
        
        # updates
        priors <- matrix(nrow=ni,ncol=length(prior))
        for(i in 1:ni) priors[i,] <- fbo[[i]]$gamma[1,]
        prior <- colMeans(priors)
        if(estprior) {
            for(k in prior.id) prior <- replace(prior,prior.est == k,mean(prior[prior.est == k]))
            if(length(prior.fix > 0)) prior <- replace(prior,prior.est == 0,prior.fix)
            prior <- prior/sum(prior)
        }
        
        for(i in 1:ni) {
            A_num[i,,] <- apply(fbo[[i]]$xi,c(2,3),sum)
            A_denom[i,] <- colSums(fbo[[i]]$gamma[1:(nt[i]-1),])         
        }
        A <- apply(A_num,c(2,3),sum)/colSums(A_denom)
        if(estA) {
            for(k in A.id) A <- replace(A,A.est == k,mean(A[A.est == k]))
            if(length(A.fix > 0)) A <- replace(A,A.est == 0,A.fix)
            A <- A/rowSums(A)
        }
        
#        if(estb) {
            for(k in b.id) {
                tdat <- data.frame()
                gamma <- numeric(0)
                for(i in (1:ni)[b.est==k]) {
                    tdat <- rbind(tdat,dat[[i]])
                    gamma <- c(gamma,as.numeric(fbo[[i]]$gamma))
                }
                tdat$w <- gamma
                b <- replace(b,b.est==k,glm(y~x-1,data=tdat,family=binomial(),weights=w)$coefficients)
            }
            if(length(b.fix)>0) b <- replace(b,b.est==0,b.fix)
#        } else {
#            for(i in 1:ni) {
#                dat[[i]]$w <- fbo[[i]]$gamma
#                b[i] <- glm(y~x-1,data=dat[[i]],family=binomial(),weights=w)$coefficients
#            }
#        }
            
#            if(b.equal) {
#                tdat <- data.frame()
#                gamma <- numeric(0)
#                for(i in 1:ni) {
#                    tdat <- rbind(tdat,dat[[i]])
#                    gamma <- c(gamma,as.numeric(fbo[[i]]$gamma))
#                }
#                tdat$w <- gamma
#                b <- glm(y~x-1,data=tdat,family=binomial(),weights=w)$coefficients
#                b <- rep(b,ni)
#            } else {
#                for(i in 1:ni) {
#                    dat[[i]]$w <- fbo[[i]]$gamma
#                    b[i] <- glm(y~x-1,data=dat[[i]],family=binomial(),weights=w)$coefficients
#                }
#            }
#        }
        
        for(i in 1:ni) {
            py <- matrix(1/(1+exp(-b[i]*dat[[i]]$x)),ncol=ns)
            B[[i]] <- (py^y[[i]])*(1-py)^(1-y[[i]])    
            fbo[[i]] <- rule.fb(A=A,B=B[[i]],prior=prior)
            LL[i] <- fbo[[i]]$logLike
        }       
        LL.dif <- sum(LL)-sum(LL.old)
        if(verbose) cat(paste("iteration",j,": difference =",LL.dif,"\n"))
    }
    npar <- sum(c(length(A.id),length(prior.id),length(b.id)))
    return(list(A=A,B=B,prior=prior,b=b,LL=LL,df=npar,iterations=j-1,convergence=LL.dif))
}

rule.BW2 <- function(y,X,Z,prior,A,b,tol=1e-4,maxiter=200,A.est,prior.est,b.est=seq(1,length(y)),A.group=rep(1,length(y)),verbose=FALSE) {
  #
  # Baum-Welch (EM) algorithm for ML estimates of the rule learning model
  # (c) 2006, M. Speekenbrink
  #
  # y = list of lenght N with response variable
  # X = list of length N, each element a T_i*k design matrix (predictors)
  # Z = S*k proportional weight matrix,
  # prior = S vector with prior probabilities
  #
  # A.est: matrix of dim(A.est)=dim(A) with integers 0,1,... where 0 indicates fixed and 1,2,... indicate unique free parameters
  # prior.est: similar indicator vector as for A.est
  # b.est: N vector with integers 0,1,... where 0 indicates fixed and 1,2,... indicate unique free parameters
  #        defaults to (1,2,...,N)
  ni <- length(y)
  ns <- nrow(Z)
  ng <- length(unique(A.group))
  #A.group <- sort(A.group,ties="first")
  A.group <- as.numeric(factor(A.group,labels=1:ng))
  if(length(A.group)!=ni) stop("A.group must have length",ni)
  
  tmp <- list()
  if(!is.list(A)) {
    for(i in 1:ng) { tmp[[i]] <- A }
    A <- tmp
  }
  if(length(A)!=ng) stop("A must have length",ng)
  
  if(length(b) < ni) b <- c(b,rep(b[length(b)],ni-length(b)))

  if(missing(A.est)) {
    estA <- FALSE
    A.id <- 1:length(A[[1]])
  } else {
    estA <- TRUE
    A.id <- unique(as.numeric(A.est[A.est!=0]))
    A.fix <- A[A.est==0]
  }
  if(missing(prior.est)) {
    estprior <- FALSE
    prior.id <- 1:length(prior)
  } else {
    estprior <- TRUE
    prior.id <- unique(as.numeric(prior.est[prior.est!=0]))
    prior.fix <- prior[prior.est==0]
  }
  
  if(length(b.est)!=length(b)) stop("b.est must have length",length(b))
  b.id <- order(unique(as.numeric(b.est[b.est!=0])))
  b.fix <- b[b.est==0]

  dat <- fbo <- B <- list()
  LL <- nt <- vector(length=ni)

  A_num <- array(dim=c(ni,nrow(A[[1]]),ncol(A[[1]])))
  A_denom <- matrix(nrow=ni,ncol=ncol(A[[1]]))
  
  for(i in 1:ni) {
    nt[i] <- nrow(X[[i]])
    dat[[i]] <- data.frame(cbind(y=y[[i]],x=as.numeric(X[[i]]%*%t(Z))))
    py <- matrix(1/(1+exp(-b[i]*dat[[i]]$x)),ncol=ns)
    B[[i]] <- (py^y[[i]])*(1-py)^(1-y[[i]])
    fbo[[i]] <- rule.fb(A=A[[A.group[i]]],B=B[[i]],prior=prior)
    LL[i] <- fbo[[i]]$logLike
  }

  
  LL.dif <- tol+1
  j <- 0
  while(j <= maxiter && abs(LL.dif) > tol) {
    LL.old <- LL
    j <- j+1
    
    # updates
    priors <- matrix(nrow=ni,ncol=length(prior))
    for(i in 1:ni) priors[i,] <- fbo[[i]]$gamma[1,]
    prior <- colMeans(priors)
    if(estprior) {
      for(k in prior.id) prior <- replace(prior,prior.est == k,mean(prior[prior.est == k]))
      if(length(prior.fix > 0)) prior <- replace(prior,prior.est == 0,prior.fix)
      prior <- prior/sum(prior)
    }
    for(i in 1:ni) {
      A_num[i,,] <- apply(fbo[[i]]$xi,c(2,3),sum)
      A_denom[i,] <- colSums(fbo[[i]]$gamma[1:(nt[i]-1),])
    }
    for(i in 1:ng) A[[i]] <- apply(A_num[A.group==i,,],c(2,3),sum)/colSums(A_denom[A.group==i,])
    if(estA) {
      for(i in 1:ng) {
        for(k in A.id) A[[i]] <- replace(A[[i]],A.est == k,mean(A[[i]][A.est == k]))
        if(length(A.fix > 0)) A[[i]] <- replace(A[[i]],A.est == 0,A.fix)
        A[[i]] <- A[[i]]/rowSums(A[[i]])
      }
    }
    for(k in b.id) {
      tdat <- data.frame()
      gamma <- numeric(0)
      for(i in (1:ni)[b.est==k]) {
        tdat <- rbind(tdat,dat[[i]])
        gamma <- c(gamma,as.numeric(fbo[[i]]$gamma))
      }
      tdat$w <- gamma
      b <- replace(b,b.est==k,glm(y~x-1,data=tdat,family=binomial(),weights=w)$coefficients)
    }
    if(length(b.fix)>0) b <- replace(b,b.est==0,b.fix)

    for(i in 1:ni) {
      py <- matrix(1/(1+exp(-b[i]*dat[[i]]$x)),ncol=ns)
      B[[i]] <- (py^y[[i]])*(1-py)^(1-y[[i]])
      fbo[[i]] <- rule.fb(A=A[[A.group[i]]],B=B[[i]],prior=prior)
      LL[i] <- fbo[[i]]$logLike
    }

    LL.dif <- sum(LL)-sum(LL.old)
    if(verbose) cat(paste("iteration",j,": difference =",LL.dif,"\n"))
  }
  npar <- sum(c(length(A.id),length(prior.id),length(b.id)))
  return(list(A=A,B=B,prior=prior,b=b,LL=LL,df=npar,iterations=j-1,convergence=LL.dif))
}


# X,Z,prior,A,b

rule.BW.old <- function(y,X,Z,prior,A,b,tol=1e-4,maxiter=200,est.prior=TRUE) {
    #
    # Baum-Welch (EM) algorithm for ML estimates of the rule learning model
    # (c) 2006, M. Speekenbrink
    #
    # y = response variable
    # X = T*k design matrix (predictors)
    # Z = S*k proportional weight matrix, 
    # prior = S vector with prior probabilities
    #
    
    nt <- nrow(X)
    ns <- nrow(Z)
    
    dat <- data.frame(cbind(y=y,x=as.numeric(X%*%t(Z))))   
    py <- matrix(1/(1+exp(-b*dat$x)),ncol=ns)
    
    B <- (py^y)*(1-py)^(1-y)
    fbo <- rule.fb(A=A,B=B,prior=prior)
    LL <- fbo$logLike
    LL.dif <- tol+1
    j <- 0
    while(j <= maxiter && abs(LL.dif) > tol) {
        LL.old <- LL
        j <- j+1
        # updates
        if(est.prior) prior <- fbo$gamma[1,]
        A <- apply(fbo$xi,c(2,3),sum)/colSums(fbo$gamma[1:(nt-1),])
        #A <- A/rowSums(A)
        b <- glm(y~x-1,data=dat,family=binomial(),weights=as.numeric(fbo$gamma))$coefficients
        
        py <- matrix(1/(1+exp(-b*dat$x)),ncol=ns)
        B <- (py^y)*(1-py)^(1-y)
        fbo <- rule.fb(A=A,B=B,prior=prior)
        
        LL <- fbo$logLike
        LL.dif <- LL-LL.old
    
    }
    return(list(A=A,B=B,prior=prior,b=b,LL=LL,iterations=j-1,convergence=LL.dif))
}