#library(mvtnorm)
library(MASS)

# Bottleneck is the apply function...

make.positive.definite <- function(x, tol=1e-6) {
    eig <- eigen(x, symmetric=TRUE)
    rtol <- tol * eig$values[1]
    if(min(eig$values) < rtol) {
        vals <- eig$values
        vals[vals < rtol] <- rtol
        srev <- eig$vectors %*% (vals * t(eig$vectors))
        dimnames(srev) <- dimnames(x)
        attr(srev,"corrected") <- TRUE
        return(srev)
    } else {
        attr(x,"corrected") <- FALSE
        return(x)
    }
}

nmcpl.prepare <- function(formula,force.bias=FALSE,data,subset) {
    # Returns a matrix of dummy input variables and vector of dummy Y variables
    if (missing(data)) {
        data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    if(!is.factor(eval(attr(mt,"variables"),envir=mf)[[1]])) {
        stop("response variable must be a factor")
    }
    Y <- model.response(mf)
    Yn <- model.matrix(~Y)
    ifelse(dim(Yn)[2]==2,Yn <- as.numeric(Yn[,2]),stop("Multinomial logistic not implemented yet"))
    X <- if(!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else matrix(, NROW(Y), 0)
    if(colnames(X)[1] != "(Intercept)" && force.bias==FALSE) {
        Xr <- X[,2:ncol(X)] # Delete first column of contrast matrix
        attr(Xr,"assign") <- attr(X,"assign")[2:ncol(X)]
        X <- Xr
        rm(Xr)
    }
    Xn <-  matrix(as.numeric(X),nrow=dim(X)[1],ncol=dim(X)[2])
    x.names <- dimnames(Xn)[[2]]
    y.names <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
    rm(Y,X)
    return(list(x=Xn,y=Yn,x.names=x.names,y.names=y.names))
}

particleDLR.fit <- function(y,x,mean,cov,beta,npart,return.full) {
    x <- as.matrix(x)
    nt <- length(y)
    nx <- ncol(x)
    xnames <- dimnames(x)[[2]]
    weight <- matrix(0,ncol=nx,nrow=nt)
    pred <- vector()
    if (return.full==TRUE) {
        state <- array(dim=c(npart,nx,nt))
    } else {
        state <- NULL
    }
    start.state <- z <- mvrnorm(n = npart, mu=mean, Sigma=cov)
    for(i in 1:length(y)) {
        # z <- t(apply(z,1,mvrnorm,n=1,Sigma=cov))
        z <- mvrnorm(n=npart,mu=rep(0,ncol(z)),Sigma=cov) + z
        w <- logis(z%*%x[i,]) # likelihood of y=1
        w <- y[i]*w + (1-y[i])*(1-w) # importance weights
        w <- w/sum(w) # normalize importance weights
        weight[i,] <- w%*%z # compute mean weights before resampling
        rs <- sample(1:npart,size=npart,replace=TRUE,prob=w) # selection step
        if(return.full==TRUE) {
            state[,,i] <- z <- z[rs,]
        } else {
            z <- z[rs,]
        }
        pred[i] <- mean(logis(z[rs,]%*%x[i,]))
    }
    dimnames(weight)[[2]] <- xnames
    return(list(weight=weight,state=state,state.cov=cov,pred=cbind(1-pred,pred)))
}
    
particleDLR <- function(formula,npart=1000,refact=2,data,subset,mean,cov,cov.cov,sd=1,sd.cov=1,beta=1,delta=.99,force.bias=FALSE,return.full=FALSE,est.var=FALSE,est.cov=FALSE,...) {
    call <- match.call()
    tmp <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(tmp), 0)
    tmp <- tmp[c(1, m)]
    tmp[[1]] <- as.name("nmcpl.prepare")
    tmp <- eval(tmp, parent.frame())
    Xn <- tmp$x
    Yn <- tmp$y
    if(est.cov) {
        if(missing(cov.cov)) {
            cov.cov <- sqrt(sd.cov)*diag(sum(1:ncol(Xn)))
        }
    } else {
        if(est.var) {
            if(missing(cov.cov)) {
                cov.cov <- sqrt(sd.cov)*diag(ncol(Xn))
            }
        }
    }
    if(missing(mean)) mean <- rep(0,ncol(Xn))
    if(length(mean) != ncol(Xn)) stop("length of mean is ",length(mean)," but should be ",ncol(Xn))
    if(missing(cov)) {
        if(length(sd) != 1 && length(sd) != ncol(Xn)) stop("sd should have length 1 or ",ncol(Xn))
        cov <- sd^2*diag(ncol(Xn))
    }
    if(est.cov) {
        fit <- particleDLR.fit.cov(y=Yn,x=Xn,mean=mean,cov=cov,cov.cov=cov.cov,beta=beta,delta=delta,npart=npart,return.full=return.full)
    } else {
        if(est.var) {
            fit <- particleDLR.fit.var(y=Yn,x=Xn,mean=mean,cov=cov,cov.cov=cov.cov,beta=beta,delta=delta,npart=npart,refact=refact,return.full=return.full)
        } else {
            fit <- particleDLR.fit(y=Yn,x=Xn,mean=mean,cov=cov,beta=beta,npart=npart,return.full=return.full)
        }
    }
    weight <- fit$weight
    #pred <- logis(diag(Xn%*%t(weight)))
    #pred <- cbind(1-pred,pred)
    pred <- fit$pred
    colnames(pred) <- tmp$y.names
    control <- list(npart=npart,refact=refact)
    return(list(call=call,response=Yn,predictor=Xn,weight=weight,state=fit$state,state.cov=fit$state.cov,predict=pred,control=control))
}

vec2cov <- function(x,n) {
    cov <- matrix(nrow=n,ncol=n)
    diag(cov) <- x[1:n]
    cov[upper.tri(cov)] <- cov[lower.tri(cov)] <- x[(n+1):length(x)]
    return(cov)
}

part.rmvnorm <- function(x,n) {
    m <- x[1:n]
    sig <- vec2cov(x[(n+1):length(x)],n=n)
    sig <- make.positive.definite(sig)
    mvrnorm(n=1,mu=m,Sigma=sig)
}

part.rmvnorm.var <- function(x,n) {
    m <- x[1:n]
    sig <- x[(n+1):length(x)]*diag(n)
    mvrnorm(n=1,mu=m,Sigma=sig)
}

particleDLR.fit.cov <- function(y,x,mean,cov,cov.cov,beta,delta,npart,refact,return.full) {
    x <- as.matrix(x)
    nt <- length(y)
    nx <- ncol(x)
    xnames <- dimnames(x)[[2]]
    weight <- matrix(0,ncol=nx,nrow=nt)
    if (return.full==TRUE) {
        state <- array(dim=c(npart,nx,nt))
    } else {
        state <- NULL
    }

    start.state <- z <- mvrnorm(n=npart,mu=mean,Sigma=cov) # initialize
    r.start <- r <- mvrnorm(n=npart,mu=c(log(diag(cov)),cov[lower.tri(cov)]),Sigma=cov.cov)
    r.mean <- colMeans(r)
    r.cov <- cov(r)
    
    hsq <- 1-((3*delta - 1)/(2*delta))^2
    a <- sqrt(1-hsq)
    
    for(i in 1:length(y)) {
        r <- t(t(a*r) + (1-a)*r.mean)
        w <- logis(z%*%x[i,]) # likelihood of y=1
        w <- y[i]*w + (1-y[i])*(1-w) # importance weights
        #w <- w/sum(w) # normalize importance weights
        rs <- sample(1:npart,size=refact*npart,replace=TRUE,prob=w/sum(w)) # selection step
        # r <- t(apply(r[rs,],1,mvrnorm,n=1,Sigma=hsq*r.cov))
        r <- mvrnorm(n=refact*npart,mu=rep(0,ncol(r)),Sigma = hsq*r.cov) + r[rs,]
        z <- t(apply(cbind(z[rs,],cbind(exp(r[,1:nx]),r[,(nx+1):ncol(r)])),1,part.rmvnorm,n=nx))
        wt <- logis(z%*%x[i,])
        wt <- y[i]*wt + (1-y[i])*(1-wt)
        wt <- wt/w[rs] # used to be wt <- wt/w
        w <- wt/sum(wt)
        w[is.na(w)] <- 0
        weight[i,] <- w%*%z # compute mean weights before resampling
        rs <- sample(1:(refact*npart),size=npart,replace=TRUE,prob=w) # selection step
        if(return.full==TRUE) {
            state[,,i] <- z <- z[rs,]
            r <- r[rs,]
        } else {
            z <- z[rs,]
            r <- r[rs,]
        }
        r.mean <- colMeans(r)
        r.cov <- cov(r)
    }
    dimnames(weight)[[2]] <- xnames
    r.cov <- vec2cov(c(exp(r.mean[1:nx]),r.mean[(nx+1):ncol(r)]),n=nx)
    return(list(weight=weight,state=state,state.cov=r.cov))
}

particleDLR.fit.var <- function(y,x,mean,cov,cov.cov,beta,delta,npart,refact,return.full) {
    x <- as.matrix(x)
    nt <- length(y)
    nx <- ncol(x)
    xnames <- dimnames(x)[[2]]
    weight <- matrix(0,ncol=nx,nrow=nt)
    pred <- vector()
    if (return.full==TRUE) {
        state <- array(dim=c(npart,nx,nt))
    } else {
        state <- NULL
    }
    r <- z <- matrix(ncol=nx,nrow=refact*npart)
    w <- vector(length=npart,mode="double")
    wt <- vector(length=refact*npart,mode="double")
    rss <- vector(length=refact*npart,mode="integer")
    rs <- as.vector(1:npart,mode="integer")

    z[rs,] <- mvrnorm(n=npart,mu=mean,Sigma=cov) # initialize
    r[rs,] <- mvrnorm(n=npart,mu=c(log(diag(cov))),Sigma=cov.cov)
    
    r.mean <- colMeans(r[rs,])
    r.cov <- cov(r[rs,])

    hsq <- 1-((3*delta - 1)/(2*delta))^2
    a <- sqrt(1-hsq)

    for(i in 1:length(y)) {
        r[1:npart,] <- t(t(a*r[rs,]) + (1-a)*r.mean)
        w <- logis(z[rs,]%*%x[i,]) # likelihood of y=1
        if(y[i]==0) w <- 1-w
        w <- w/sum(w) # normalize importance weights
        w[is.na(w)] <- 0
        rss <- sample(1:npart,size=refact*npart,replace=TRUE,prob=w) # selection step
        #r <- t(apply(r[rss,],1,mvrnorm,n=1,Sigma=hsq*r.cov))
        r <- mvrnorm(n=refact*npart,mu=rep(0,ncol(r)),Sigma=hsq*r.cov) + r[rss,]
        
        #z <- t(apply(cbind(z[rss,],exp(r)),1,part.rmvnorm.var,n=nx))
        z <- sqrt(exp(r))*mvrnorm(n=refact*npart,mu=rep(0,ncol(r)),Sigma=diag(ncol(r))) + z[rss,]
        
        wt <- logis(z%*%x[i,])
        if(y[i]==0) wt <- 1-wt
        wt <- wt/w[rss]
        wt[is.na(wt)] <- 0
        wt <- wt/sum(wt)
        weight[i,] <- wt%*%z # compute mean weights before resampling
        rs <- sample(1:(refact*npart),size=npart,replace=TRUE,prob=wt) # selection step
        r.mean <- colMeans(r[rs,])
        r.cov <- cov(r[rs,])

        z[1:npart,] <- z[rs,]
        pred[i] <- mean(logis(z[rs,]%*%x[i,]))

        if(return.full==TRUE) {
            state[,,i] <- z[rs,]
        } 
        if(i%%5==0) cat("estimation trial",i,"complete \n")
    }
    dimnames(weight)[[2]] <- xnames
    r.cov <- exp(r.mean)*diag(nx)
    return(list(weight=weight,state=state,state.cov=r.cov,predict=cbind(1-pred,pred)))
}

particleDLR.fit.var.weight <- function(y,x,mean,cov,cov.cov,beta,delta,npart,refact,return.full) {
    x <- as.matrix(x)
    nt <- length(y)
    nx <- ncol(x)
    xnames <- dimnames(x)[[2]]
    weight <- matrix(0,ncol=nx,nrow=nt)
    if (return.full==TRUE) {
        state <- array(dim=c(npart,nx,nt))
    } else {
        state <- NULL
    }

    start.state <- z <- rmvnorm(npart,mean=mean,sigma=cov) # initialize

    r.start <- r <- rmvnorm(npart,mean=c(log(diag(cov))),sigma=cov.cov)
    r.mean <- colMeans(r)
    r.cov <- cov(r)

    hsq <- 1-((3*delta - 1)/(2*delta))^2
    a <- sqrt(1-hsq)

    for(i in 1:length(y)) {
        r <- t(t(a*r) + (1-a)*r.mean)
        l <- logis(z%*%x[i,]) # likelihood of y=1
        l <- y[i]*l + (1-y[i])*(1-l) # importance weights
        g <- w*l # normalize importance weights
        g[is.na(g)] <- 0
        g <- g/sum(g)
        rs <- sample(1:npart,size=npart,replace=TRUE,prob=g) # selection step
        r <- t(apply(r[rs,],1,rmvnorm,n=1,sigma=hsq*r.cov))
        r.mean <- colMeans(r)
        r.cov <- cov(r)
        z <- t(apply(cbind(z[rs,],cbind(exp(r))),1,part.rmvnorm.var,n=nx))
        w <- logis(z%*%x[i,])
        w <- y[i]*w + (1-y[i])*(1-w)
        w[is.na(w)] <- 0
        w <- w/l
        w <- w/sum(w)
        weight[i,] <- w%*%z # compute mean weights before resampling
        if(return.full==TRUE) {
            state[,,i] <- z
        } 
    }
    dimnames(weight)[[2]] <- xnames
    r.cov <- exp(r.mean)*diag(nx)
    return(list(weight=weight,state=state,state.cov=r.cov))
}