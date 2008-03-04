mcpl.prepare <- function(formula,data,subset) {
    # Returns a matrix of numeric input variables and output variable
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
    if(is.factor(eval(attr(mt,"variables"),envir=mf)[[1]])) {
        stop("response variable is a factor")
    }
    Y <- model.response(mf)
    X <- if(!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else matrix(, NROW(Y), 0)
    Xn <- as.matrix(X)
    x.names <- dimnames(Xn)[[2]]
    #y.names <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
    return(list(x=Xn,y=Y,x.names=x.names))
}

grad.lin <- function(y,x,w) {
    if(!is.matrix(x)) {
        return(as.numeric(x*(x%*%w - y)))
    } else return(as.numeric(t(x)%*%(x%*%w - y)))
}

onlineR.fit <- function(y,x,ws,eta,alpha,beta,delta,method,window.size) {
    if(!all(alpha >= 0)) stop("negative alpha not allowed")
    if(beta < 0) stop("negative beta not allowed")
    if(!all(eta >= 0)) stop("negative eta not allowed")
    if(window.size < 0) stop("negative window.size not allowed")
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    weight <- matrix(0,ncol=ncol(x),nrow=nrow(x))
    w <- ws
    if(method=="LMS") {
        if(beta < 0 || beta > 1) stop("beta should be between 0 and 1")
        for(i in 1:length(y)) {
            j <- ifelse(i-window.size > 0, i-window.size,1)
            w <- w - (eta/(i)^alpha)*(1/abs(i-j+1))*grad.lin(y=y[j:i],x=x[j:i,],w=w) + ifelse(i>2,beta*(w-as.numeric(weight[i-2,])),ifelse(i==2,beta*(w-ws),0))
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
            r <- ifelse(i==1,grad.lin(y=y[j:i],x=x[j:i,],w=w),(1-delta)*r + delta*grad.lin(y=y[j:i],x=x[j:i,],w=w))
            eta.o[i] <- eta <- eta+alpha*eta*(beta*sqrt(sum(r^2))-eta)
            w <- w - eta*(1/abs(i-j+1))*grad.lin(y=y[j:i],x=x[j:i,],w=w)
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
            g <- sum((y[j:i]-x[j:i,]%*%(w+ck*d))^2 - (y[j:i]-x[j:i,]%*%(w-ck*d))^2)/(2*ck*d)
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
            g <- (grad.logis(y=y[j:i],x=x[j:i,],w=(w + ck*d)) - grad.lin(y=y[j:i],x=x[j:i,],w=(w-ck*d)))
            denom <- 1/t(matrix(2*ck*d,ncol=nw))
            hess <- .5*(g*denom + t(g*denom))
            hessi <- 1/(i+1)*hess + ifelse(i>1,i/(1+i)*hessi,0)
            checker <- FALSE
            i <- 0
            ihessi <- NULL
            while(checker==FALSE) {
                ihessi <- try(solve(hessi + i*diag(nw)),silent=TRUE)
                if(ihessi) checker <- TRUE
                i <- i+.2
            }
            w <- w - (eta/(i+1)^alpha)*ihessi*grad.lin(y=y[j:i],x=x[j:i,],w=w)
            weight[i,] <- w
        }
    }
    dimnames(weight)[[2]] <- xnames
    return(list(weight=weight,if(method=="MALR") eta=eta.o))
}

onlineR <- function(formula,eta=.001,data,subset,alpha=0,beta=0,delta=.5,ws,ws.range=.1,method="LMS",window.size=0,...) {
    call <- match.call()
    if(window.size < 0) stop("negative window.size not allowed")
    if (missing(data)) {
        data <- environment(formula)
    }
    
    if(missing(subset)) {
        tmp <- mcpl.prepare(formula=formula,data=data)
    } else {
        tmp <- mcpl.prepare(formula=formula,data=data,subset=subset)
    }
    
    if(length(eta)!=ncol(tmp$x)) {
        eta <- c(eta,rep(eta[length(eta)],ncol(tmp$x)-length(eta)))
    }
    if(length(alpha)!=ncol(tmp$x)) {
        alpha <- c(alpha,rep(alpha[length(alpha)],ncol(tmp$x)-length(alpha)))
    }
    if(missing(ws)) ws <- runif(ncol(tmp$x),min=-ws.range,max=ws.range)
    if(length(ws) != ncol(tmp$x)) stop("length of ws is ",length(ws)," but should be ",ncol(tmp$x))
    fit <- onlineR.fit(y=tmp$y,x=tmp$x,ws=ws,eta=eta,alpha=alpha,beta=beta,delta=delta,method=method,window.size=window.size)
    #weight <- fit$weight
    weight <- rbind(ws,fit$weight)[1:nrow(tmp$x),]
    wT <- fit$weight[nrow(tmp$x),]
    pred <- diag(tmp$x%*%t(weight))
    #pred <- cbind(1-pred,pred)
    #colnames(pred) <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
    control <- list(eta=eta,alpha=alpha,beta=beta,ws=ws,method=method,window.size=window.size)
    if(!is.null(fit$eta.o)) control <- c(control,veta=eta.o)
    return(list(call=call,response=tmp$y,predictor=tmp$x,weight=weight,wT=wT,predict=pred,control=control))
}

onlineR.mcpl <- function(formula,eta=.1,data,subset,alpha=0,beta=0,delta=.5,ws,ws.range=.1,method="LMS",window.size=0,force.bias=FALSE,...) {
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
        
        if(length(eta)==1) eta <- matrix(eta,ncol=1)
        if(length(alpha)==1) alpha <- matrix(alpha,ncol=1)
        
        if(!is.matrix(eta)) stop("eta should be a matrix")
        if(!is.matrix(alpha)) stop("alpha must be a matrix")
        
        if(nrow(eta)>1) {
            if(nrow(eta) < length(ids)) {
                warning("nrow of eta is smaller than length of id. will fill in with last given row")
                eta <- rbind(eta,matrix(rep(eta[nrow(eta),],times=length(ids)-nrow(eta)),ncol=ncol(eta),byrow=TRUE))
            }
            if(nrow(eta) > length(ids)) {
                warning("length of eta is greater than length of id. will cutoff last values")
            }
        } else eta <- matrix(rep(eta,length(ids)),ncol=ncol(eta),byrow=TRUE)
        if(nrow(alpha)>1) {
            if(nrow(alpha) < length(ids)) {
                warning("nrow of alpha is smaller than length of id. will fill in with last value")
                alpha <- rbind(alpha,matrix(rep(alpha[nrow(alpha),],times=length(ids)-nrow(alpha)),ncol=ncol(alpha),byrow=TRUE))
            }
            if(nrow(alpha) > length(ids)) {
                warning("length of alpha is greater than length of id. will cutoff last values")
            }
        } else alpha <- matrix(rep(alpha,length(ids)),ncol=ncol(alpha),byrow=TRUE)
        tdat <- data[id==ids[1],]
        ntrial[1] <- nrow(tdat)
        out <- onlineR(formula=formula,eta=eta[1,],data=tdat,alpha=alpha[1,],beta=beta,delta=delta,ws=ws,ws.range=ws.range,method=method,window.size=window.size,force.bias=force.bias)
        if(length(ids)>1) {
            for (i in 2:length(ids)) {
                tdat <- data[id==ids[i],]
                ntrial[i] <- nrow(tdat)
                fit <- onlineR(formula=formula,eta=eta[i,],data=tdat,alpha=alpha[i,],beta=beta,delta=delta,ws=ws,ws.range=ws.range,method=method,window.size=window.size,force.bias=force.bias)
                out$weight <- rbind(out$weight,fit$weight)
                out$wT <- rbind(out$wT,fit$wT)
                out$response <- c(out$response,fit$response)
                out$predictor <- rbind(out$predictor,fit$predictor)
                #out$predict <- rbind(out$predict,fit$predict)
                out$predict <- c(out$predict,fit$predict)
                out$control$ws <- rbind(out$control$ws,fit$control$ws)
                if(!is.null(fit$control$veta)) out$control$veta <- c(out$control$veta,fit$control$veta)
            }
        }
        attr(out,"id") <- id
        out$call <- call
        if(!is.null(attr(data,"design"))) attr(out,"design") <- data[,attr(data,"design")]
        return(out)
    } else {
        ntrial <- nrow(data)
        out <- onlineR(formula=formula,eta=eta,data=data,alpha=alpha,beta=beta,delta=delta,ws=ws,ws.range=ws.range,method=method,window.size=window.size,force.bias=force.bias)
        if(!is.null(attr(data,"design"))) attr(out,"design") <- data[,attr(data,"design")]
        out$call <- call
        attr(out,"ntrial") <- ntrial
        return(out)
    }
}

onlineR.optim <- function(par,z_formula,z_data,z_response,z_eta,z_eta.values,z_alpha,z_alpha.values,z_ws,z_sd = NULL) {
    par <- exp(par)
    if(!all(abs(par)!=Inf)) {
        res <- Inf } else {
        eta <- z_eta.values
        alpha <- z_alpha.values
        for(i in seq(length=length(par))) {
            eta <- replace(eta,z_eta==i,par[i])
            alpha <- replace(alpha,z_alpha==i,par[i])
        }
        if(length(eta) > 1 && !is.matrix(eta)) eta <- matrix(eta,nrow=1)
        if(length(alpha) > 1 && !is.matrix(alpha)) alpha <- matrix(alpha,nrow=1)
        fit <- onlineR.mcpl(formula=z_formula,data=z_data,eta=eta,alpha=alpha,ws=z_ws)
        #fit$response <- z_response
        #if(is.null(z_sd)) z_sd <- sqrt(sum((fit$predict - z_response)^2)/(length(z_response)-1))
        #res <- -sum(dnorm(z_response,mean=fit$predict,sd=z_sd,log=TRUE))
        res <- mean((z_response-fit$predict)^2)
        #ntrial <- attr(fit,"ntrial")
        #res <- -sum(log(rowSums(model.matrix(~factor(fit$response)-1)*fit$predict)))
    }
    #attr(res,"ll") <- res
    return(res)
}


onlineR.fit.mle <- function(formula,data,response,eta,eta.start,alpha,alpha.start,ws,lower=NULL,upper=NULL,sd=NULL,method="BFGS") {
    attr(data,"design") <- NULL
    id.n <- 1
    ids <- NULL
    if(!is.null(attr(data,"id"))) {
        ids <- levels(factor(data[,attr(data,"id")]))
        id.n <- length(ids)
    }
    
    lower <- if(!is.null(lower)) log(lower) else -Inf
    upper <- if(!is.null(upper)) log(upper) else Inf
    
    if(!is.matrix(eta)) stop("eta must be a matrix")
    if(nrow(eta)!=id.n) stop("eta must have nrow",id.n)
    if(!is.matrix(alpha)) stop("alpha must be a matrix")
    if(nrow(alpha)!=id.n) stop("alpha must have nrow",id.n)

    par.n <- length(unique(c(eta[eta!=0],alpha[alpha!=0])))
    eta.final <- eta.start
    alpha.final <- alpha.start

    subsets <- parsubsets(cbind(eta,alpha))
    MSE <- vector()
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
            teta <- eta[subsets[[i]],]
            teta.start <- eta.start[subsets[[i]],]
            talpha <- alpha[subsets[[i]],]
            talpha.start <- alpha.start[subsets[[i]],]
            tpar.id <- unique(c(teta[teta!=0],talpha[talpha!=0]))
            par <- vector()
            teta[teta!=0] <- as.numeric(factor(teta[teta!=0],levels=tpar.id,labels=1:length(tpar.id)))
            talpha[talpha!=0] <- as.numeric(factor(talpha[talpha!=0],levels=tpar.id,labels=1:length(tpar.id)))
            for(j in seq(length=length(tpar.id))) par <- c(par,log(c(eta.start[eta==tpar.id[j]],alpha.start[alpha==tpar.id[j]])[1]))
            #if(length(par)==1) {
                #if(!all(abs(c(lower,upper)) != Inf))
            #    opt <- optimize(f=onlineR.optim,interval=c(log(1e-50),log(1e+50)),z_formula=formula,z_data=tdat,z_response=tresponse,z_eta=teta,z_eta.values=teta.start,z_alpha=talpha,z_alpha.values=talpha.start,z_ws=ws)
            #    opt$value <- opt$objective
            #    opt$par <- opt$minimum
            #} else {
            if(!all(abs(c(lower,upper)) == Inf)) {
                opt <- optim(par,fn=onlineR.optim,lower=lower,upper=upper,z_formula=formula,z_data=tdat,z_response=tresponse,z_eta=teta,z_eta.values=teta.start,z_alpha=talpha,z_alpha.values=talpha.start,z_ws=ws,z_sd=sd,method="L-BFGS-B")
            } else {
                opt <- optim(par,fn=onlineR.optim,z_formula=formula,z_data=tdat,z_response=tresponse,z_eta=teta,z_eta.values=teta.start,z_alpha=talpha,z_alpha.values=talpha.start,z_ws=ws,z_sd=sd,method=method)
            }
            #}
            MSE[i] <- opt$value
            opt.parameters[[i]] <- opt$par
            eta.final[subsets[[i]],] <- teta.start
            alpha.final[subsets[[i]],] <- talpha.start
            for(j in seq(length=length(par))) {
                eta.final <- replace(eta.final,eta==tpar.id[j],exp(opt$par[j]))
                lambda.alpha <- replace(alpha.final,alpha==tpar.id[j],exp(opt$par[j]))
            }
            cat("estimation subset",i,"complete.\n")
        }
        MSE <- mean(MSE)
    } else {

        par.id <- unique(c(eta[eta!=0],alpha[alpha!=0]))

        eta[eta!=0] <- as.numeric(factor(eta[eta!=0],levels=par.id,labels=1:length(par.id)))
        alpha[alpha!=0] <- as.numeric(factor(alpha[alpha!=0],levels=par.id,labels=1:length(par.id)))
        par <- vector()
        
        for(i in seq(length=length(par.id))) par[i] <- log(c(eta.start[eta==i],alpha.start[alpha==i])[1])
        if(!all(abs(c(lower,upper)) == Inf)) {
            opt <- optim(par,fn=onlineR.optim,lower=lower,upper=upper,z_formula=formula,z_data=data,z_response=response,z_eta=eta,z_eta.values=eta.start,z_alpha=alpha,z_alpha.values=alpha.start,z_ws=ws,z_sd=sd,method="L-BFGS-B")
        } else {
            opt <- optim(par,fn=onlineR.optim,z_formula=formula,z_data=data,z_response=response,z_eta=eta,z_eta.values=eta.start,z_alpha=alpha,z_alpha.values=alpha.start,z_ws=ws,z_sd=sd,method=method)
        }
        for(i in seq(length=length(par.id))) {
            eta.final <- replace(eta.final,eta==par.id[i],exp(opt$par[i]))
            alpha.final <- replace(alpha.final,alpha==par.id[i],exp(opt$par[i]))
        }
        MSE <- opt$value
        opt.parameters <- opt$par
    }
    return(list(eta=eta.final,alpha=alpha.final,MSE=MSE,npar=par.n,opt.parameters=opt.parameters))
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