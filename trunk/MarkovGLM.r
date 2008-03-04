fb <- function(A,B,prior) {

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
        xi[i,,] <- (alpha[i,]%*%t(B[i+1,]*beta[i+1,]))*A
    }
    like <- -sum(log(c))
    return(list(alpha=alpha,beta=beta,gamma=gamma,xi=xi,c=c,logLike=like))
}

viterbi <- function(A,B,prior) {
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

MarkovGLM.fit <- function(A,B,prior,formula,data,family=gaussian(),A.est,maxiter=100,tol=1e-5,criterion="logLike") {
    ns <- length(prior)
    nt <- nrow(B)

    fbo <- fb(A=A,B=B,prior=prior)
    LL <- fbo$logLike

    if(missing(A.est)) {
        estA <- FALSE
    } else {
        estA <- TRUE
        A.id <- unique(as.numeric(A.est[A.est!=0]))
        A.fix <- A[A.est==0]
    }
    
    fit <- list()
    coeff <- vector()
    sd <- vector()
    for(i in 1:ns) {
        data$w <- fbo$gamma[,i] #important, since glm can't seem to find w otherwise...
        fit[[i]] <- glm(formula=formula,data=data,family=family,weights=w)
        coeff <- c(coeff,fit[[i]]$coefficients)
        sd[i] <- sqrt(sum(fbo$gamma[,i]*fit[[i]]$residuals^2/sum(fbo$gamma[,i]))) # estimated dispersion from glm seems inaccurate
        #B[,i] <- dnorm(fit[[i]]$residuals,sd=sd[i])
    }

    par.new <- c(A,prior,coeff,sd)
    j <- 0
    converge <- FALSE
    while(j <= maxiter && !converge) {
        LL.old <- LL
        par.old <- par.new
        j <- j+1
        
        # updates
        prior <- fbo$gamma[1,]

        A <-  apply(fbo$xi,c(2,3),sum)/colSums(fbo$gamma[1:(nt-1),])
        if(estA) {
            for(k in A.id) A <- replace(A,A.est==k,mean(A[A.est==k]))
            if(length(A.fix > 0)) replace(A,A.est==0,A.fix)
            A <- A/rowSums(A)
        }
        coeff <- vector()
        for(i in 1:ns) {
            data$w <- fbo$gamma[,i] # checked, should NOT be sqrt(fbo$gamma[,i])
            fit[[i]] <- glm(formula=formula,data=data,family=family,weights=w)
            coeff <- c(coeff,fit[[i]]$coefficients)
            sd[i] <- sqrt(sum(fbo$gamma[,i]*fit[[i]]$residuals^2/sum(fbo$gamma[,i])))
            B[,i] <- dnorm(fit[[i]]$residuals,sd=sd[i])
        }
        par.new <- c(A,prior,coeff,sd)
        
        fbo <- fb(A=A,B=B,prior=prior)
        LL <- fbo$logLike

        # check for convergence
        if(criterion == "param") {
            if(all(abs(par.old-par.new) < tol)) converge <- TRUE
        } else {
            if(LL-LL.old < tol) converge <- TRUE
        }
    }
    pred <- matrix(nrow=nt,ncol=ns)
    coeff <- list()
    for(i in 1:ncol(pred)) {
        pred[,i] <- fit[[i]]$fitted.values
        coeff[[i]] <- fit[[i]]$coefficients
    }
    return(list(response=fit[[1]]$y,predictor=model.matrix(attr(fit[[1]]$model,"terms"),fit[[1]]$model),predicted=pred,A=A,B=B,prior=prior,coefficients=coeff,sd=sd,LL=LL,iterations=j-1))
}
    
# Wat functies om initiele waarden voor A, B en prior te krijgen

init <- function(A,B,prior,maxiter=100,tol=1e-5) {
    j <- 0
    nt <- nrow(B)
    fbo <- fb(A=A,B=B,prior=prior)
    LL <- fbo$logLike
    LL.dif <- tol+1
    while(j <= maxiter && abs(LL.dif) > tol) {
        LL.old <- LL
        j <- j+1

        # updates
        prior <- fbo$gamma[1,]
        A <- apply(fbo$xi,c(2,3),sum)/colSums(fbo$gamma[1:(nt-1),])
        fbo <- fb(A=A,B=B,prior=prior)
        LL <- fbo$logLike
        LL.dif <- sum(LL)-sum(LL.old)
    }
    return(list(A=A,prior=prior))
}

MarkovGLM.init <- function(formula,data,family=gaussian(),maxComponents=10,nState,repFM=5,maxtry=5) {
    # use flexmix to get starting model
    if(missing(nState)) {
        mix <- stepFlexmix(formula=formula,data=data,K=1:maxComponents,nrep=repFM)
        # use BIC to determine number of components
        mix <- mix[[which.min(sapply(mix,BIC))]]
        s.id <- as.numeric(names(table(cluster(mix))))
    } else {
        ch <- FALSE
        j <- 1
        while(!ch) {
            mix <- stepFlexmix(formula=formula,data=data,K=nState,nrep=repFM)
            if(length(mix@components)==nState | j >= maxtry) ch <- TRUE
            s.id <- as.numeric(names(table(cluster(mix))))
            j <- j+1
        }
    }

    ns <- length(s.id)
    nt <- length(mix@cluster)

    # determine starting values for A and prior
    start.values <- init(A=matrix(1/ns,ncol=ns,nrow=ns),B=mix@posterior$scaled[,s.id],prior=mix@posterior$scaled[1,s.id])

    B <- matrix(nrow=nt,ncol=ns)
    for(i in 1:ns) {
        data$w <- mix@posterior$scaled[,s.id[i]]
        fit <- glm(formula=formula,data=data,family=family,weights=w)
        sd <- sqrt(sum(mix@posterior$scaled[,s.id[i]]*fit$residuals^2/sum(mix@posterior$scaled[,s.id[i]])))
        B[,i] <- dnorm(fit$residuals,sd=sd)
    }
    return(list(A=start.values$A,B=B,prior=start.values$prior))
}


BIC.markovGLM <- function(mod) {
    npar <- length(mod$A) - ncol(mod$A) + length(mod$prior) - 1 + length(unlist(mod$coefficients)) + length(mod$sd)
    nobs <- length(mod$response)
    return(-2*mod$LL + npar*log(nobs))
}

AIC.markovGLM <- function(mod) {
    npar <- length(mod$A) - ncol(mod$A) + length(mod$prior) - 1 + length(unlist(mod$coefficients)) + length(mod$sd)
    return(-2*mod$LL + 2*npar)
}

AICc.markovGLM <- function(mod) {
    npar <- length(mod$A) - ncol(mod$A) + length(mod$prior) - 1 + length(unlist(mod$coefficients)) + length(mod$sd)
    nobs <- length(mod$response)
    return(-2*mod$LL + 2*npar + (2*npar*(npar+1))/(nobs - npar - 1))
}

return.viterbi <- function(mod) {
    mod$state <- viterbi(A=mod$A,B=mod$B,prior=mod$prior)
    mod$predicted <- mod$predicted[cbind(1:nrow(mod$B),mod$state)]
    mod$weight <- t(matrix(unlist(mod$coefficients[mod$state]),nrow=length(mod$coefficients[[1]])))
    colnames(mod$weight) <- names(mod$coefficients[[1]])
    mod$sd <- mod$sd[mod$state]
    mod$LL <- sum(dnorm(mod$response,mean=mod$predicted,sd=mod$sd,log=TRUE))
    return(mod)
}

