logis <- function(x) {
    as.numeric(1/(1+exp(-x)))
}

grad.logis <- function(y,x,w) {
    if(!is.matrix(x)) {
        return(as.numeric(x*(logis(x%*%w) - y)))
        #return(as.numeric(-x*(y - logis(x%*%w))))
    #} else return(as.numeric(-t(x)%*%(y - logis(x%*%w))))
    } else return(as.numeric(t(x)%*%(logis(x%*%w) - y)))
}

onlineLR.fit <- function(y,x,ws,eta,alpha,beta,delta,method,window.size) {
    if(alpha < 0) stop("negative alpha not allowed")
    if(beta < 0) stop("negative beta not allowed")
    if(eta < 0) stop("negative eta not allowed")
    if(window.size < 0) stop("negative window.size not allowed")
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    weight <- matrix(0,ncol=ncol(x),nrow=nrow(x))
    w <- ws
    if(method=="LMS") {
        if(beta < 0 || beta > 1) stop("beta should be between 0 and 1")
        for(i in 1:length(y)) {
            j <- ifelse(i-window.size > 0, i-window.size,1)
            #w <- w - (eta/(i+1)^alpha)*(1/abs(i-j+1))*grad.logis(y=y[j:i],x=x[j:i,],w=w) + ifelse(i>2,beta*(w-as.numeric(weight[i-2,])),ifelse(i==2,beta*(w-ws),0))
            w <- w - (eta/(i)^alpha)*(1/abs(i-j+1))*grad.logis(y=y[j:i],x=x[j:i,],w=w) + ifelse(i>2,beta*(w-as.numeric(weight[i-2,])),ifelse(i==2,beta*(w-ws),0))
            weight[i,] <- w
        }
    }
    if(method=="MALR") {
        if(beta < 0 || beta > 1) stop("beta should be between 0 and 1")
        # Murata's adaptive learning rate algorithm.
        # see: Lecun, Bottou, Orr & Muller (1998) par. 4.7
        eta.o <- vector()
        for(i in 1:length(y)) {
            j <- ifelse(i-window.size > 0, i-window.size,1)
            r <- ifelse(i==1,grad.logis(y=y[j:i],x=x[j:i,],w=w),(1-delta)*r + delta*grad.logis(y=y[j:i],x=x[j:i,],w=w))
            eta.o[i] <- eta <- eta+alpha*eta*(beta*sqrt(sum(r^2))-eta)
            w <- w - eta*(1/abs(i-j+1))*grad.logis(y=y[j:i],x=x[j:i,],w=w)
            weight[i,] <- w
        }
    }
    if(method=="1SPSA") {
        # 1st-order Stochastic Perburtation Stochastic Approximation (SPSA)
        # Spall, J. C. (2000). Adaptive stochastic approximation by the simultaneous perturbation method. IEEE Transactions on automatic control, 45, 1839--1853
        nw <- length(ws)
        for(i in 1:length(y)) {
            j <- ifelse(i-window.size > 0, i-window.size,1)
            d <- 2*rbinom(nw,size=1,p=.5) - 1
            ck <- beta/(i+1)^delta
            #g <- (grad.logis(y=y[j:i],x=x[j:i,],w=(w + ck*d)) - grad.logis(y=y[j:i],x=x[j:i,],w=(w-ck*d)))/(2*ck*d)
            g <- sum((y[j:i]-logis(x[j:i,]%*%(w+ck*d)))^2 - (y[j:i]-logis(x[j:i,]%*%(w-ck*d)))^2)/(2*ck*d)
            w <- w - (eta/(i+1)^alpha)*(1/abs(i-j+1))*as.numeric(g)
            weight[i,] <- w
        }
    }
    if(method=="2SPSA") {
        # 2nd-order Stochastic Perburtation Stochastic Approximation (SPSA)
        # Spall, J. C. (2000). Adaptive stochastic approximation by the simultaneous perturbation method. IEEE Transactions on automatic control, 45, 1839--1853
        nw <- length(ws)
        for(i in 1:length(y)) {
            j <- ifelse(i-window.size > 0, i-window.size,1)
            d <- 2*rbinom(nw,size=1,p=.5) - 1
            ck <- beta/(i+1)^delta
            g <- (grad.logis(y=y[j:i],x=x[j:i,],w=(w + ck*d)) - grad.logis(y=y[j:i],x=x[j:i,],w=(w-ck*d)))
            denom <- 1/t(matrix(2*ck*d,ncol=nw))
            hess <- .5*(g*denom + t(g*denom))
            hessi <- 1/(i+1)*hess + ifelse(i>1,i/(1+i)*hessi,0)
            #g <- sum((y[j:i]-logis(x[j:i,]%*%(w+ck*d)))^2 - (y[j:i]-logis(x[j:i,]%*%(w-ck*d)))^2)/(2*ck*d)
            checker <- FALSE
            i <- 0
            ihessi <- NULL
            while(checker==FALSE) {
                ihessi <- try(solve(hessi + i*diag(nw)),silent=TRUE)
                if(ihessi) checker <- TRUE
                i <- i+.2
            }
            w <- w - (eta/(i+1)^alpha)*ihessi*grad.logis(y=y[j:i],x=x[j:i,],w=w)
            weight[i,] <- w
        }
    }
    dimnames(weight)[[2]] <- xnames
    return(list(weight=weight,if(method=="MALR") eta=eta.o))
}

onlineLR <- function(formula,eta=.1,data,subset,alpha=0,beta=0,delta=.5,ws,ws.range=.1,method="LMS",window.size=0,force.bias=FALSE,...) {
    call <- match.call()
    if(window.size < 0) stop("negative window.size not allowed")
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
        stop("response variable is not a factor")
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
    if(missing(ws)) ws <- runif(ncol(X),min=-ws.range,max=ws.range)
    if(length(ws) != ncol(X)) stop("length of ws is ",length(ws)," but should be ",ncol(X))
    rm(Y,X)
    fit <- onlineLR.fit(y=Yn,x=Xn,ws=ws,eta=eta,alpha=alpha,beta=beta,delta=delta,method=method,window.size=window.size)
    weight <- rbind(ws,fit$weight)[1:nrow(Xn),]
    wT <- fit$weight[nrow(Xn),]
    pred <- logis(diag(Xn%*%t(weight)))
    pred <- cbind(1-pred,pred) # one step ahead predictions
    colnames(pred) <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
    control <- list(eta=eta,alpha=alpha,beta=beta,ws=ws,method=method,window.size=window.size)
    if(!is.null(fit$eta.o)) control <- c(control,veta=eta.o)
    return(list(call=call,response=Yn,predictor=Xn,weight=weight,wT=wT,predict=pred,control=control))
}

onlineLR.nmcpl <- function(formula,eta=.1,data,subset,alpha=0,beta=0,delta=.5,ws,ws.range=.1,method="LMS",window.size=0,force.bias=FALSE,...) {
    # todo: subset passing to onlineLR
    call <- match.call()
    if(!missing(subset)) {
        if(!is.null(attr(data,"id"))) id <- attr(data,"id")
        data <- subset(data,subset)
        if(!is.null(id)) attr(data,"id") <- id
    } else {
        subset <- NULL
    }
    if(!is.null(attr(data,"id"))) {
        id <- data[,attr(data,"id")]
        id <- factor(id) # drop unused levels in id
        ids <- levels(id)
        ntrial <- vector(length=length(ids))
        if(length(ids) > 1) {
            if(length(eta)>1) {
                if(length(eta) < length(ids)) {
                    warning("length of eta is smaller than length of id. will fill in with last value")
                    eta <- c(eta,rep(eta[length(eta)],times=length(ids)-length(eta)))
                }
                if(length(eta) > length(ids)) {
                    warning("length of eta is greater than length of id. will cutoff last values")
                }
            } else eta <- rep(eta,length(ids))
            tdat <- data[id==ids[1],]
            ntrial[1] <- nrow(tdat)
            out <- onlineLR(formula=formula,eta=eta[1],data=tdat,alpha=alpha,beta=beta,delta=delta,ws=ws,ws.range=ws.range,method=method,window.size=window.size,force.bias=force.bias)
            for (i in 2:length(ids)) {
                tdat <- data[id==ids[i],]
                ntrial[i] <- nrow(tdat)
                fit <- onlineLR(formula=formula,eta=eta[i],data=tdat,alpha=alpha,beta=beta,delta=delta,ws=ws,ws.range=ws.range,method=method,window.size=window.size,force.bias=force.bias)
                out$weight <- rbind(out$weight,fit$weight)
                out$wT <- rbind(out$wT,fit$wT)
                out$response <- c(out$response,fit$response)
                out$predictor <- rbind(out$predictor,fit$predictor)
                out$predict <- rbind(out$predict,fit$predict)
                out$control$ws <- rbind(out$control$ws,fit$control$ws)
                if(!is.null(fit$control$veta)) out$control$veta <- c(out$control$veta,fit$control$veta)
            }
            attr(out,"id") <- id
            out$call <- call
            if(!is.null(attr(data,"design"))) attr(out,"design") <- data[,attr(data,"design")]
            attr(out,"ntrial") <- ntrial
            return(out)
        } else {
            ntrial <- nrow(data)
            out <- onlineLR(formula=formula,eta=eta,data=data,alpha=alpha,beta=beta,delta=delta,ws=ws,ws.range=ws.range,method=method,window.size=window.size,force.bias=force.bias)
            if(!is.null(attr(data,"design"))) attr(out,"design") <- data[,attr(data,"design")]
            out$call <- call
            attr(out,"ntrial") <- ntrial
            return(out)
        }
    } else {
        ntrial <- nrow(data)
        out <- onlineLR(formula=formula,eta=eta,data=data,alpha=alpha,beta=beta,delta=delta,ws=ws,ws.range=ws.range,method=method,window.size=window.size,force.bias=force.bias)
        if(!is.null(attr(data,"design"))) attr(out,"design") <- data[,attr(data,"design")]
        out$call <- call
        attr(out,"ntrial") <- ntrial
        return(out)
    }
}

onlineLR.parfit <- function(pars,...) {
    call <- match.call(expand.dots = TRUE)
    parnames <- c("alpha", "beta", "delta", "eta")
    fpar <- match(parnames, names(call), 0)
    fpar <- call[c(1, fpar)]
    fpar[[1]] <- as.name("list")
    fpar <- eval(fpar)
    return(fpar)
}

onlineLR.optim.eta <- function(pars,formula,data,response,lambda,alpha,ws) {
    eta <- pars[1]
    if(!is.null(attr(data,"id"))) {
        id <- attr(data,"id")
        if(length(levels(factor(data[,id]))) == 1) {
            fit <- onlineLR(formula=formula,data=data,eta=eta,alpha=alpha,ws=ws)
            fit$response <- model.matrix(~response-1)[,2]
            fit$weight <- lambda*fit$weight
            pred <- logis(diag(fit$predictor%*%t(fit$weight)))
            fit$predict <- cbind(1-pred,pred)
            ll <- summary.nmcplModelFit(fit)$ll
        } else {
            fit <- onlineLR.nmcpl(formula=formula,data=data,eta=eta,alpha=alpha,ws=ws)
            fit$response <- model.matrix(~response-1)[,2]
            fit$weight <- lambda*fit$weight
            pred <- logis(diag(fit$predictor%*%t(fit$weight)))
            fit$predict <- cbind(1-pred,pred)
            ll <- summary.nmcplModelFit2(fit)$ll
        }
    } else {
        fit <- onlineLR(formula=formula,data=data,eta=eta,alpha=alpha,ws=ws)
        fit$response <- model.matrix(~response-1)[,2]
        fit$weight <- lambda*fit$weight
        pred <- logis(diag(fit$predictor%*%t(fit$weight)))
        fit$predict <- cbind(1-pred,pred)
        ll <- summary.nmcplModelFit(fit)$ll
    }
    res <- sum(ll)
    if(abs(res)==Inf) res <- 1e+50
    attr(res,"ll") <- ll
    return(res)
}

onlineLR.optim.lambda <- function(pars,formula,data,response,eta,alpha,ws) {
    lambda <- pars[1]
    if(!is.null(attr(data,"id"))) {
        id <- attr(data,"id")
        if(length(levels(factor(data[,id]))) == 1) {
            fit <- onlineLR(formula=formula,data=data,eta=eta,alpha=alpha,ws=ws)
            fit$response <- model.matrix(~response-1)[,2]
            fit$weight <- lambda*fit$weight
            pred <- logis(diag(fit$predictor%*%t(fit$weight)))
            fit$predict <- cbind(1-pred,pred)
            ll <- summary.nmcplModelFit(fit)$ll
        } else {
            fit <- onlineLR.nmcpl(formula=formula,data=data,eta=eta,alpha=alpha,ws=ws)
            fit$response <- model.matrix(~response-1)[,2]
            fit$weight <- lambda*fit$weight
            pred <- logis(diag(fit$predictor%*%t(fit$weight)))
            fit$predict <- cbind(1-pred,pred)
            ll <- summary.nmcplModelFit2(fit)$ll
        }
    } else {
        fit <- onlineLR(formula=formula,data=data,eta=eta,alpha=alpha,ws=ws)
        fit$response <- model.matrix(~response-1)[,2]
        fit$weight <- lambda*fit$weight
        pred <- logis(diag(fit$predictor%*%t(fit$weight)))
        fit$predict <- cbind(1-pred,pred)
        ll <- summary.nmcplModelFit(fit)$ll
    }
    res <- sum(ll)
    if(abs(res)==Inf) res <- 1e+50
    attr(res,"ll") <- ll
    return(res)
}

onlineLR.optim <- function(par,z_formula,z_data,z_response,z_eta,z_eta.values,z_lambda,z_lambda.values,z_ws,z_alpha,z_force.bias) {

    par <- exp(par)
    if(!all(abs(par)!=Inf)) {
        res <- Inf } else {

        eta <- z_eta.values
        lambda <- z_lambda.values
        for(i in seq(length=length(par))) {
            eta <- replace(eta,z_eta==i,par[i])
            lambda <- replace(lambda,z_lambda==i,par[i])
        }
        fit <- onlineLR.nmcpl(formula=z_formula,eta=eta,data=z_data,alpha=z_alpha,ws=z_ws,force.bias=z_force.bias)
        fit$response <- model.matrix(~z_response-1)[,2]
        ntrial <- attr(fit,"ntrial")
        tlambda <- vector()
        for(i in 1:length(ntrial)) tlambda <- c(tlambda,rep(lambda[i],ntrial[i]))
        fit$predict <- fit$predict^tlambda #THIS GOES WRONG!
        fit$predict <- fit$predict/rowSums(fit$predict)
        res <- -sum(log(rowSums(model.matrix(~factor(fit$response)-1)*fit$predict)))
    }
    attr(res,"ll") <- res
    return(res)
}


onlineLR.fit.mle <- function(formula,data,response,eta,eta.start,lambda,lambda.start,alpha,ws,force.bias=FALSE,method="Nelder-Mead") {

    attr(data,"design") <- NULL
    if(!is.factor(response)) stop("response must be a factor")

    id.n <- 1
    ids <- NULL
    if(!is.null(attr(data,"id"))) {
        ids <- levels(factor(data[,attr(data,"id")]))
        id.n <- length(ids)
    }

    if(!is.numeric(eta)) stop("eta must be a numerical vector")
    if(length(eta)!=id.n) stop("eta must have length",id.n)

    if(!is.numeric(lambda)) stop("lambda must be a numerical vector")
    if(length(lambda)!=id.n) stop("lambda must have length",id.n)

    par.n <- length(unique(c(eta[eta!=0],lambda[lambda!=0])))

    eta.final <- eta.start
    lambda.final <- lambda.start

    subsets <- parsubsets(cbind(eta,lambda))
    LL <- vector()
    opt.parameters <- list()

    if(length(subsets)>1) {
        for(i in 1:length(subsets)) {
            #tmean <- mean
            #if(!is.null(mean)) {
            #    tmean <- mean[subsets[[i]],]
            #    if(!is.matrix(tmean)) tmean <- matrix(tmean,ncol=length(tmean))
            #}
            tdat <- data[data[,attr(data,"id")] %in% ids[subsets[[i]]], ]
            attr(tdat,"id") <- attr(data,"id")
            attr(tdat,"trial") <- attr(data,"trial")
            tresponse <- response[data[,attr(data,"id")] %in% ids[subsets[[i]]]]
            teta <- eta[subsets[[i]]]
            teta.start <- eta.start[subsets[[i]]]
            tlambda <- lambda[subsets[[i]]]
            tlambda.start <- lambda.start[subsets[[i]]]
            tpar.id <- unique(c(teta[teta!=0],tlambda[tlambda!=0]))
            par <- vector()
            teta[teta!=0] <- as.numeric(factor(teta[teta!=0],levels=tpar.id,labels=1:length(tpar.id)))
            tlambda[tlambda!=0] <- as.numeric(factor(tlambda[tlambda!=0],levels=tpar.id,labels=1:length(tpar.id)))
            for(j in seq(length=length(tpar.id))) par <- c(par,log(c(eta.start[eta==tpar.id[j]],lambda.start[lambda==tpar.id[j]])[1]))
            if(length(par)==1 & method=="Nelder-Mead") {
                opt <- optimize(f=onlineLR.optim,interval=c(log(1e-50),log(1e+50)),z_formula=formula,z_data=tdat,z_response=tresponse,z_eta=teta,z_eta.values=teta.start,z_lambda=tlambda,z_lambda.values=tlambda.start,z_alpha=alpha,z_ws=ws,z_force.bias=force.bias)
                opt$value <- opt$objective
                opt$par <- opt$minimum
            } else {
                opt <- optim(par,fn=onlineLR.optim,method=method,z_formula=formula,z_data=tdat,z_response=tresponse,z_eta=teta,z_eta.values=teta.start,z_lambda=tlambda,z_lambda.values=tlambda.start,z_ws=ws,z_alpha=alpha,z_force.bias=force.bias)
            }
            LL[i] <- opt$value
            opt.parameters[[i]] <- opt$par
            eta.final[subsets[[i]]] <- eta.start
            lambda.final[subsets[[i]]] <- tlambda.start
            for(j in seq(length=length(par))) {
                eta.final <- replace(eta.final,eta==tpar.id[j],exp(opt$par[j]))
                lambda.final <- replace(lambda.final,lambda==tpar.id[j],exp(opt$par[j]))
            }
            cat("estimation subset",i,"complete.\n")
        }
        LL <- sum(LL)
    } else {

        par.id <- unique(c(eta[eta!=0],lambda[lambda!=0]))

        eta[eta!=0] <- as.numeric(factor(eta[eta!=0],levels=par.id,labels=1:length(par.id)))
        lambda[lambda!=0] <- as.numeric(factor(lambda[lambda!=0],levels=par.id,labels=1:length(par.id)))
        par <- vector()
        
        for(i in seq(length=length(par.id))) par[i] <- log(c(eta.start[eta==i],lambda.start[lambda==i])[1])
        if(length(par)==1 & method=="Nelder-Mead") {
                opt <- optimize(f=onlineLR.optim,interval=c(log(1e-50),log(1e+50)),z_formula=formula,z_data=data,z_response=response,z_eta=eta,z_eta.values=eta.start,z_lambda=lambda,z_lambda.values=lambda.start,z_ws=ws,z_alpha=alpha,z_force.bias=force.bias)
                opt$value <- opt$objective
                opt$par <- opt$minimum
            } else {
                opt <- optim(par,fn=onlineLR.optim,method=method,z_formula=formula,z_data=data,z_response=response,z_eta=eta,z_eta.values=eta.start,z_lambda=lambda,z_lambda.values=lambda.start,z_ws=ws,z_alpha=alpha,z_force.bias=force.bias)
        }

        for(i in seq(length=length(par.id))) {
                eta.final <- replace(eta.final,eta==par.id[i],exp(opt$par[i]))
                lambda.final <- replace(lambda.final,lambda==par.id[i],exp(opt$par[i]))
        }
        LL <- opt$value
        opt.parameters <- opt$par
    }
    return(list(eta=eta.final,lambda=lambda.final,LL=LL,npar=par.n,BIC=2*LL + (par.n)*log(length(response)),opt.parameters=opt.parameters))
}

parsubsets <- function(x) {
    # returns number of independent rows in x
    # x: Lisrel style indicator matrix
    #    i.e. a matrix with integers denoting parameters in which (0 is fixed, and 1,2,... indicate different parameters)
    npar <- length(unique(x[x!=0]))
    set <- lapply(lapply(lapply(lapply(1:npar,"==",x),"rowSums"),">",0),"which")
    subsets <- list()
    i <- 1
    while(length(set) > 0) {
        subsets[[i]] <- set[[which.max(unlist(lapply(set,"length")))]]
        tmp <- vector(length=1)
        while(length(tmp) > 0) {
            tmp <- which(!(unlist(lapply(lapply(lapply(set,"%in%",subsets[[i]]),"!"),"all"))))
            subsets[[i]] <- unique(c(subsets[[i]],unlist(set[tmp])))
            subsets[[i]] <- subsets[[i]][order(subsets[[i]])]
            set[tmp] <- NULL
            if(length(set) == 0) break
        }
        i <- i+1
    }
    return(subsets)
}