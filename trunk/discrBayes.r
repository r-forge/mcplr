discrBayes.barpoints <- function(x) {
    # x = vector with grid points
    n <- length(x)
    midpoint <- x[2:n]+(x[1:(n-1)]-x[2:n])/2
    return(c(x[1]-abs(x[1]-midpoint[1]),midpoint,x[n]+abs(x[n]-midpoint[n-1])))
}

discrBayes.dif <- function(x) {
    # x = vector with grid points
    n <- length(x)
    x <- x[2:n]-x[1:(n-1)]
    x[x < 1e-10] <- 0 # round very small difference to avoid floating point errors
    return(x)
}

discrBayes.nprior <- function(grid,mean,sd) {
    if(length(mean) < length(grid)) {
        warning("length mean not equal to length of grid, will fill in with last value")
        mean <- c(mean,rep(mean[length(mean)],length(grid)-length(mean)))
    }
    if(length(sd) < length(grid)) {
        warning("length sd not equal to length of grid, will fill in with last value")
        sd <- c(sd,rep(sd[length(sd)],length(grid)-length(sd)))
    }
    points <- lapply(grid,discrBayes.barpoints)
    p <- list()
    for(i in 1:length(grid)) {
        p[[i]] <- pnorm(points[[i]],mean=mean[i],sd=sd[i])
        p[[i]] <- discrBayes.dif(p[[i]])
        if(sum(p[[i]]) != 0) p[[i]] <- p[[i]]/sum(p[[i]]) else p[[i]] <- rep(1/length(p[[i]]),length(p[[i]]))
    }
    #p <- lapply(p,discrBayes.dif)
    #p <- lapply(p,function(x) x/sum(x)) # normalize
    
    #check <- lapply(lapply(p,complete.cases),all)
    #for(i in 1:length(grid)) if(!check[[i]]) p[[i]] <- rep(1/length(p[[i]]),length(p[[i]]))
    
    # check
    return(p)
}

#apply(expand.grid(discrBayes.nprior(grid=list(-3:3,-2:2),mean=0,sd=2)),1,prod)

discrBayes.fit.nprior <- function(pars,y,x,grid,beta,dummy=1,mean,sd) {
    if(missing(mean) && missing(sd)) {
        if(length(pars) != 2*length(grid)) stop("there must be",2*length(grid),"parameters")
        mean <- pars[1:length(grid)]
        sd <- pars[(length(grid)+1):length(pars)]
        p <- discrBayes.nprior(grid=grid,mean=mean,sd=sd)
        p <- apply(expand.grid(p),1,prod)
        p <- p/sum(p)
        fit <- discrBayes.fit(y=y,x=x,p=p,grid=grid,beta=beta)
        pred <- logis(diag(x%*%t(fit$weight)))
        ll <- -sum(log(y*pred + (1-y)*(1-pred)))
        if(abs(ll) == Inf) ll <- 1e+50
        return(ll)
    } else {
        if(missing(mean)) {
            if(length(pars) != length(grid)) stop("there must be",length(grid),"parameters")
            mean <- pars
            p <- discrBayes.nprior(grid=grid,mean=mean,sd=sd)
            p <- apply(expand.grid(p),1,prod)
            p <- p/sum(p)
            fit <- discrBayes.fit(y=y,x=x,p=p,grid=grid,beta=beta)
            pred <- logis(diag(x%*%t(fit$weight)))
            return(-sum(log(y*pred + (1-y)*(1-pred))))
        }
        if(missing(sd)) {
            if(length(pars) != length(grid)) stop("there must be",length(grid),"parameters")
            sd <- pars
            p <- discrBayes.nprior(grid=grid,mean=mean,sd=sd)
            p <- apply(expand.grid(p),1,prod)
            p <- p/sum(p)
            fit <- discrBayes.fit(y=y,x=x,p=p,grid=grid,beta=beta)
            pred <- logis(diag(x%*%t(fit$weight)))
            ll <- -sum(log(y*pred + (1-y)*(1-pred)))
            if(abs(ll) == Inf) ll <- 1e+50
            return(ll)
        }
    }
}

discrBayes.fit <- function(y,x,p,grid,beta) {
    if(!is.loaded("discrBayes")) dyn.load("C:\\Projects\\mcpl-R\\discrbayes-new")
    sw <- rowSums(as.matrix(expand.grid(grid)))
    mod <- .C("discrBayes",as.integer(y),as.integer(t(x)),as.integer(nrow(x)),as.integer(ncol(x)),as.double(unlist(grid)),as.integer(unlist(lapply(grid,length))),as.double(sw),as.integer(length(sw)),as.double(p),as.double(rep(0,length(sw))),as.double(beta),weight=as.double(x),pred=as.double(y))
    weight <- mod$weight
    pred <- mod$pred
    weight <- t(matrix(weight,nrow=ncol(x)))
    return(list(weight=weight,pred=pred))
}

discrBayes <- function(formula,grid=list(-3:3),beta=1,p,data,subset,mean=NULL,sd=NULL,force.bias=FALSE,...) {
    # p: an array with p for each grid point
    call <- match.call()
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
    nc <- ncol(Xn)
    if(length(grid)!=nc) {
        # fill grid by repeating last element in grid
        for(i in (length(grid)+1):nc) {
            grid[[i]] <- grid[[length(grid)]]
        }
    }
    if(missing(p)) p <- NULL
    if(!is.null(sd)) {
        if(length(sd)!=length(grid)) sd <- c(sd,rep(sd[length(sd)],length(grid)-length(sd)))
        if(is.null(mean)) {
            #mean <- rep(0,length(sd))
            mean <- unlist(lapply(grid,mean))
        } else {
            if(length(mean)!=length(grid)) mean <- c(mean,rep(mean[length(mean)],length(grid)-length(mean)))
        }
        p <- discrBayes.nprior(grid=grid,mean=mean,sd=sd)
        p <- apply(expand.grid(p),1,prod)
        p <- p/sum(p)
    }
    if(is.null(p)) {
        p <- array(1,dim=as.numeric(lapply(grid,length))) # initialize belief array
        p <- p/sum(p) # normalize belief array
        mean <- unlist(lapply(grid,mean))
    }
    # make sure p is ok
    p <- replace(p,is.na(p),0)
    p <- replace(p,p==Inf,1)
    p <- replace(p,p==-Inf,0)
    fit <- discrBayes.fit(y=Yn,x=Xn,p=p,grid=grid,beta=beta)
    #weight <- fit$weight
    #pred <- c(.5,fit$pred[1:(nrow(Xn)-1)]) # FIX ME!!!!
    weight <- rbind(mean,fit$weight)[1:nrow(Xn),]
    wT <- fit$weight[nrow(Xn),]
    pred <- fit$pred
    pred <- Yn*pred + (1-Yn)*(1-pred)
    pred <- cbind(1-pred,pred)
    colnames(pred) <- levels(Y)
    control <- list(grid=grid,beta=beta,mean=mean,sd=sd)
    return(list(call=call,response=Yn,predictor=Xn,weight=weight,wT=wT,predict=pred,control=control))
}

discrBayes.nmcpl <- function(formula,grid=list(-3:3),beta=1,p,data,subset,mean=NULL,sd=NULL,force.bias=FALSE,...) {
    # todo: subset passing to onlineLR
    if(!missing(subset)) {
        data <- subset(data,subset)
    } else {
        subset <- NULL
    }
    if(!is.null(attr(data,"id"))) {
        id <- factor(data[,attr(data,"id")])
        #if(!is.factor(id)) id <- as.factor(id)
        ids <- levels(id)
        ntrial <- vector(length=length(ids))
        
        if(length(ids) > 1) {
            if(!is.null(sd)) {
                if(!is.matrix(sd)) stop("sd must be a matrix")
                if(nrow(sd)!=length(ids)) for(i in (nrow(sd)+1):length(ids)) sd <- rbind(sd,sd[nrow(sd),])
            }
            if(!is.null(mean)) {
                if(!is.matrix(mean)) stop("mean must be a matrix")
                if(nrow(mean)!=length(ids)) for(i in (nrow(mean)+1):length(ids)) mean <- rbind(mean,mean[nrow(mean),])
            }

            tdat <- data[id==ids[1],]
            ntrial[1] <- nrow(tdat)

            if(!missing(sd)) {
                tsd <- matrix(sd[1,],nrow=1)
                tm <- matrix(mean[1,],nrow=1)
            } else {
                tsd <- tm <- NULL
            }
            if(!missing(p)) out <- discrBayes(formula=formula,grid=grid,data=tdat,beta=beta,p=p,mean=tm,sd=tsd,force.bias=force.bias) else out <- discrBayes(formula=formula,grid=grid,data=tdat,beta=beta,mean=tm,sd=tsd,force.bias=force.bias)
            for (i in 2:length(ids)) {
                tdat <- data[id==ids[i],]
                ntrial[i] <- nrow(tdat)
                if(!missing(sd)) {
                    tsd <- matrix(sd[1,],nrow=1)
                    tm <- matrix(mean[1,],nrow=1)
                } else {
                    tsd <- tm <- NULL
                }
                if(!missing(p)) fit <- discrBayes(formula=formula,grid=grid,data=tdat,beta=beta,p=p,mean=tm,sd=tsd,force.bias=force.bias) else fit <- discrBayes(formula=formula,grid=grid,data=tdat,beta=beta,mean=tm,sd=tsd,force.bias=force.bias)
                out$weight <- rbind(out$weight,fit$weight)
                out$wT <- rbind(out$wT,fit$wT)
                out$response <- c(out$response,fit$response)
                out$predictor <- rbind(out$predictor,fit$predictor)
                out$predict <- rbind(out$predict,fit$predict)
            }
            attr(out,"id") <- id
        } else {
        ntrial <- nrow(data)
        if(!missing(p)) out <- discrBayes(formula=formula,grid=grid,data=data,beta=beta,p=p,mean=mean,sd=sd,force.bias=force.bias) else out <- discrBayes(formula=formula,grid=grid,data=data,beta=beta,mean=mean,sd=sd,force.bias=force.bias)
        }
    } else {
        ntrial <- nrow(data)
        if(!missing(p)) out <- discrBayes(formula=formula,grid=grid,data=data,beta=beta,p=p,mean=mean,sd=sd,force.bias=force.bias) else out <- discrBayes(formula=formula,grid=grid,data=data,beta=beta,mean=mean,sd=sd,force.bias=force.bias)
    }
    if(!is.null(attr(data,"design"))) attr(out,"design") <- data[,attr(data,"design")]
    attr(out,"ntrial") <- ntrial
    return(out)
}

#discrBayes.optim.sd <- function(par,z_formula,z_data,z_response,z_grid,z_beta,z_mean,z_force.bias,z_sd,z_sd.values,z_sd.n,z_lambda,z_lambda.values,z_lambda.n) {
discrBayes.optim.sd <- function(par,z_formula,z_data,z_response,z_grid,z_beta,z_mean,z_force.bias,z_sd,z_sd.values,z_lambda,z_lambda.values) {

    par <- exp(par)
    if(!all(abs(par)!=Inf)) {
        res <- Inf } else {

        sd <- z_sd.values
        lambda <- z_lambda.values
        #for(i in seq(length=z_sd.n)) sd <- replace(sd,z_sd==i,par[i])
        #for(i in seq(length=z_lambda.n)) lambda <- replace(lambda,z_lambda==(i+z_sd.n),par[i+z_sd.n])
        for(i in seq(length=length(par))) {
            sd <- replace(sd,z_sd==i,par[i])
            lambda <- replace(lambda,z_lambda==i,par[i])
        }
        fit <- discrBayes.nmcpl(formula=z_formula,grid=z_grid,beta=z_beta,data=z_data,p=NULL,mean=z_mean,sd=sd,force.bias=z_force.bias)
        #fit$response <- model.matrix(~z_response-1)[,2]
        ntrial <- attr(fit,"ntrial")
        tlambda <- vector()
        for(i in 1:length(ntrial)) tlambda <- c(tlambda,rep(lambda[i],ntrial[i]))
        fit$predict <- fit$predict^tlambda #THIS GOES WRONG!
        fit$predict <- fit$predict/rowSums(fit$predict)
        res <- -sum(log(rowSums(model.matrix(~z_response-1)*fit$predict)))
    }
    attr(res,"ll") <- res
    return(res)
}

discrBayes.optim.sd.glm <- function(par,z_formula,z_data,z_response,z_grid,z_beta,z_mean,z_force.bias,z_sd,z_sd.values,z_lambda,z_lambda.values) {

    par <- exp(par)
    if(!all(abs(par)!=Inf)) {
        res <- Inf } else {

        sd <- z_sd.values
        #lambda <- z_lambda.values
        for(i in seq(length=length(par))) {
            sd <- replace(sd,z_sd==i,par[i])
            #lambda <- replace(lambda,z_lambda==i,par[i])
        }
        fit <- discrBayes.nmcpl(formula=z_formula,grid=z_grid,beta=z_beta,data=z_data,p=NULL,mean=z_mean,sd=sd,force.bias=z_force.bias)
        ntrial <- attr(fit,"ntrial")
        z_data$eta <- rowSums(fit$weight*fit$predictor)
        mod <- glm(formula=z_lformula,data=z_data,family=binomial())
        tlambda <- vector()
        for(i in 1:length(ntrial)) tlambda <- c(tlambda,rep(lambda[i],ntrial[i]))
        fit$predict <- fit$predict^tlambda #THIS GOES WRONG!
        fit$predict <- fit$predict/rowSums(fit$predict)
        res <- -sum(log(rowSums(model.matrix(~z_response-1)*fit$predict)))
    }
    attr(res,"ll") <- res
    return(res)
}

myLogLike <- function(fit,response,lambda) {
        ntrial <- attr(fit,"ntrial")
        tlambda <- vector()
        for(i in 1:length(ntrial)) tlambda <- c(tlambda,rep(lambda[i],ntrial[i]))
        fit$predict <- fit$predict^tlambda #THIS GOES WRONG!
        fit$predict <- fit$predict/rowSums(fit$predict)
        res <- -sum(log(rowSums(model.matrix(~response-1)*fit$predict)))
        return(res)
}

#fit <- optim(c(1,1,1,1),fn=discrBayes.optim.sd,lower=1e-50,method="L-BFGS-B",formula=e~c1+c2+c3+c4-1,data=data,response=data$r,grid=rep(list(-3:3),4),beta=1,mean=rep(0,4),force.bias=FALSE)

#fit <- optim(c(1,1,1,1),fn=discrBayes.optim.sd,lower=1e-50,method="L-BFGS-B",formula=e~c1+c2+c3+c4-1,data=dat,response=dat$r,grid=rep(list(-3:3),4),beta=1,mean=rep(0,4),force.bias=FALSE)

discrBayes.fit.mle <- function(formula,data,response,grid,beta,mean,sd,sd.start,lambda,lambda.start,force.bias=FALSE,lower=NULL,upper=NULL,est.start=FALSE,sd.values=c(.1,.5,1,2,5,10,100,7000),lambda.values=c(.2,.5,1,2,3,5),allow.optimize=FALSE,method="Nelder-Mead") {

    # est.sd: matrix with est.sd[i,j] = 0 (fixed) or j (indicator for parameter)
    # est.lambda: vector with est.sd[i,j] = 0 (fixed) or j (indicator for parameter)
    # response: a factor
    
    attr(data,"design") <- NULL
    if(!is.factor(response)) stop("response must be a factor")

    if(!is.null(lower)) lower <- log(lower) else lower <- -Inf
    if(!is.null(upper)) upper <- log(upper) else upper <- Inf
    
    id.n <- 1
    ids <- NULL
    if(!is.null(attr(data,"id"))) {
        ids <- levels(factor(data[,attr(data,"id")]))
        id.n <- length(ids)
    }

    if(!is.matrix(sd)) stop("sd must be a matrix")
    if(nrow(sd)!=id.n) stop("sd must have ",id.n,"rows")
    
    if(!is.numeric(lambda)) stop("lambda must be a numerical vector")
    if(length(lambda)!=id.n) stop("lambda must have length",id.n)
    
    par.n <- length(unique(c(sd[sd!=0],lambda[lambda!=0])))
    
    sd.final <- sd.start
    lambda.final <- lambda.start
        
    subsets <- parsubsets(cbind(sd,lambda))
    LL <- vector()
    opt.parameters <- list()
    if(length(subsets)>1) {
        for(i in 1:length(subsets)) {
            tmean <- mean
            if(!is.null(mean)) {
                tmean <- mean[subsets[[i]],]
                if(!is.matrix(tmean)) tmean <- matrix(tmean,ncol=length(tmean))
            }
            tdat <- data[data[,attr(data,"id")] %in% ids[subsets[[i]]], ]
            attr(tdat,"id") <- attr(data,"id")
            attr(tdat,"trial") <- attr(data,"trial")
            tresponse <- response[data[,attr(data,"id")] %in% ids[subsets[[i]]]]
            tsd <- sd[subsets[[i]],]
            if(!is.matrix(tsd)) tsd <- matrix(tsd,ncol=length(tsd))
            tsd.start <- sd.start[subsets[[i]],]
            if(!is.matrix(tsd.start)) tsd.start <- matrix(tsd.start,ncol=length(tsd.start))
            tlambda <- lambda[subsets[[i]]]
            tlambda.start <- lambda.start[subsets[[i]]]
            
            if(est.start) {
                tmp <- discrBayes.optim.start(formula=formula,data=tdat,response=tresponse,grid=grid,beta=beta,mean=tmean,sd=tsd,sd.start=tsd.start,lambda=tlambda,lambda.start=tlambda.start,sd.values=sd.values,lambda.values=lambda.values,force.bias=force.bias)
                tsd.start <- tmp$sd.start
                tlambda.start <- tmp$lambda.start
            }
            #parvals <- unique(c(tsd[tsd!=0],tlambda[tlambda!=0]))
            tpar.id <- unique(c(tsd[tsd!=0],tlambda[tlambda!=0]))
            par <- vector()
            tsd[tsd!=0] <- as.numeric(factor(tsd[tsd!=0],levels=tpar.id,labels=1:length(tpar.id)))
            tlambda[tlambda!=0] <- as.numeric(factor(tlambda[tlambda!=0],levels=tpar.id,labels=1:length(tpar.id)))
            #tpar.id <- unique(c(tsd[tsd!=0],tlambda[tlambda!=0]))
            for(j in seq(length=length(tpar.id))) par <- c(par,log(c(sd.start[sd==tpar.id[j]],lambda.start[lambda==tpar.id[j]])[1]))
            if(length(par)==1 & allow.optimize & method=="Nelder-Mead") {
                #cat("using optimize")
                if(upper > log(1e+50)) upper <- log(1e+50)
                if(lower < log(1e-50)) lower <- log(1e-50)
                opt <- optimize(f=discrBayes.optim.sd,interval=c(lower,upper),z_formula=formula,z_data=tdat,z_response=tresponse,z_sd=tsd,z_sd.values=tsd.start,z_lambda=tlambda,z_lambda.values=tlambda.start,z_grid=grid,z_mean=tmean,z_beta=beta,z_force.bias=force.bias)
                opt$value <- opt$objective
                opt$par <- opt$minimum
            } else {
                if(all(abs(c(lower,upper)) == Inf)) {
                    opt <- optim(par,fn=discrBayes.optim.sd,method=method,lower=lower,upper=upper,z_formula=formula,z_data=tdat,z_response=tresponse,z_sd=tsd,z_sd.values=tsd.start,z_lambda=tlambda,z_lambda.values=tlambda.start,z_grid=grid,z_mean=tmean,z_beta=beta,z_force.bias=force.bias)
                } else {
                    opt <- optim(par,fn=discrBayes.optim.sd,method="L-BFGS-B",lower=lower,upper=upper,z_formula=formula,z_data=tdat,z_response=tresponse,z_sd=tsd,z_sd.values=tsd.start,z_lambda=tlambda,z_lambda.values=tlambda.start,z_grid=grid,z_mean=tmean,z_beta=beta,z_force.bias=force.bias)
                }
            }
            LL[i] <- opt$value
            opt.parameters[[i]] <- opt$par
            sd.final[subsets[[i]],] <- sd.start[subsets[[i]],] <- tsd.start
            lambda.final[subsets[[i]]] <- lambda.start[subsets[[i]]] <- tlambda.start
            for(j in seq(length=length(par))) {
                sd.final <- replace(sd.final,sd==tpar.id[j],exp(opt$par[j]))
                lambda.final <- replace(lambda.final,lambda==tpar.id[j],exp(opt$par[j]))
            }
            cat("estimation subset",i,"complete.\n")
        }
        LL <- sum(LL)
    } else {
    
            #for(i in seq(length=sd.n)) sd.final <- replace(sd.final,sd==i,exp(opt$par[i]))
            #for(i in seq(length=lambda.n)) lambda.final <- replace(lambda.final,lambda==i,exp(opt$par[i+sd.n]))
            
    #par.sd <- seq(length=sd.n)
    #par.lambda <- seq(length=lambda.n) #+ length(par.sd)
        #sd.fixed <- sd == 0
        #lambda.fixed <- lambda == 0
        #sd.n <- length(unique(sd[!sd.fixed]))
        #lambda.n <- length(unique(lambda[!lambda.fixed]))

        par.id <- unique(c(sd[sd!=0],lambda[lambda!=0]))

    #sd[!sd.fixed] <- as.numeric(factor(sd[!sd.fixed])
        sd[sd!=0] <- as.numeric(factor(sd[sd!=0],levels=par.id,labels=1:length(par.id)))
    #lambda[!lambda.fixed] <- as.numeric(factor(lambda[!lambda.fixed]))
        lambda[lambda!=0] <- as.numeric(factor(lambda[lambda!=0],levels=par.id,labels=1:length(par.id)))
        if(est.start) {
            tmp <- discrBayes.optim.start(formula=formula,data=data,response=response,grid=grid,beta=beta,mean=mean,sd=sd,sd.start=sd.start,lambda=lambda,lambda.start=lambda.start,sd.values=sd.values,lambda.values=lambda.values,force.bias=force.bias)
            sd.start <- tmp$sd.start
            lambda.start <- tmp$lambda.start
        }
        
        par <- vector()
        #for(i in seq(length=sd.n)) par[i] <- log(sd.start[sd==i][1])
        #for(i in seq(length=lambda.n)) par[i+sd.n] <- log(lambda.start[lambda==i][1])
        for(i in seq(length=length(par.id))) par[i] <- log(c(sd.start[sd==i],lambda.start[lambda==i])[1])
        if(length(par)==1 & allow.optimize & method=="Nelder-Mead") {
                #cat("using optimize")
                if(upper > log(1e+50)) upper <- log(1e+50)
                if(lower < log(1e-50)) lower <- log(1e-50)
                opt <- optimize(f=discrBayes.optim.sd,interval=c(lower,upper),z_formula=formula,z_data=data,z_response=response,z_sd=sd,z_sd.values=sd.start,z_lambda=lambda,z_lambda.values=lambda.start,z_grid=grid,z_mean=mean,z_beta=beta,z_force.bias=force.bias)
                opt$value <- opt$objective
                opt$par <- opt$minimum
        } else {
            if(all(abs(c(lower,upper)) == Inf)) {
                opt <- optim(par,fn=discrBayes.optim.sd,method=method,lower=lower,upper=upper,z_formula=formula,z_data=data,z_response=response,z_sd=sd,z_sd.values=sd.start,z_lambda=lambda,z_lambda.values=lambda.start,z_grid=grid,z_mean=mean,z_beta=beta,z_force.bias=force.bias)
            } else {
                opt <- optim(par,fn=discrBayes.optim.sd,lower=lower,upper=upper,z_formula=formula,z_data=data,z_response=response,z_sd=sd,z_sd.values=sd.start,z_lambda=lambda,z_lambda.values=lambda.start,z_grid=grid,z_mean=mean,z_beta=beta,z_force.bias=force.bias,method="L-BFGS-B")
            }
        }

        #sd.final <- sd.start
        #lambda.final <- lambda.start
        #for(i in seq(length=sd.n)) sd.final <- replace(sd.final,sd==i,exp(opt$par[i]))
        #for(i in seq(length=lambda.n)) lambda.final <- replace(lambda.final,lambda==i,exp(opt$par[i+sd.n]))
        for(i in seq(length=length(par.id))) {
                sd.final <- replace(sd.final,sd==par.id[i],exp(opt$par[i]))
                lambda.final <- replace(lambda.final,lambda==par.id[i],exp(opt$par[i]))
        }
        LL <- opt$value
        opt.parameters <- opt$par
    }
    return(list(sd=sd.final,lambda=lambda.final,LL=LL,npar=par.n,BIC=2*LL + (par.n)*log(length(response)),sd.start=sd.start,lambda.start=lambda.start,opt.parameters=opt.parameters))
}

parsubsets <- function(x) {
    # returns number of independent rows in x
    # x: Lisrel style indicator matrix
    #    i.e. a matrix with integers denoting parameters in which (0 is fixed, and 1,2,... indicate different parameters)
    npar <- length(unique(x[x!=0]))
    if(npar < 1) return(list(1:nrow(x)))
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

getArrayIndex <- function(x,dim) {
    #return(which(array(1:prod(dim),dim=dim)==x,arr.ind=TRUE))
    rank <- length(dim)
    wh1 <- x - 1
    wh <- 1 + wh1%%dim[1]
    wh <- rep(1 + wh1%%dim[1], length = rank)
    if (rank >= 2) {
        denom <- 1
        for (i in 2:rank) {
          denom <- denom * dim[i - 1]
          nextd1 <- wh1%/%denom
          wh[i] <- 1 + nextd1%%dim[i]
        }
    }
    storage.mode(wh) <- "integer"
    return(wh)
}


discrBayes.optim.start <- function(formula,data,response,grid,beta,mean,sd,sd.start,lambda,lambda.start,sd.values=c(.1,.5,1,2,5,10,100,7000),lambda.values=c(.2,.5,1,2,3,5),force.bias=FALSE) {
    sd.id <- unique(sd[sd!=0])
    lambda.id <- unique(lambda[lambda!=0])
    par.id <- c(sd.id,lambda.id)
    sd.dim <- rep(length(sd.values),length(sd.id))
    lambda.dim <- rep(length(lambda.values),length(lambda.id))
    par.dim <- c(sd.dim,lambda.dim)
    ll <- array(dim=par.dim)
    if(prod(sd.dim) > 5000) stop("Number of models to fit VERY large! You kidding me!")
    if(length(par.id)==0) return(list(sd.start=sd.start,lambda.start=lambda.start))
    #ll <- vector()
    if(length(sd.id) > 0) {
        #pars <- list() #length=length(sd.id))
        #for(i in 1:length(sd.id)) pars[[i]] <- sd.values
        #par1 <- expand.grid(pars)
        #pars <- prod(sd.id)
        #pars <- array(1:((length(sd.values))^(length(sd.id))),dim=rep(length(sd.values),length(sd.id)))
        for(i in 1:prod(sd.dim)) {
            tsd <- sd.start
            tmp <- getArrayIndex(i,sd.dim)
            for(j in 1:length(tmp)) tsd <- replace(tsd,sd==sd.id[j],sd.values[tmp][j])
            mod <- discrBayes.nmcpl(formula=formula,data=data,grid=grid,beta=beta,mean=mean,sd=tsd)
            if(length(lambda.id) > 0) {
                for(j in 1:prod(lambda.dim)) {
                #pars <- list() #length=length(lambda.id))
                #for(i in 1:length(lambda.id)) pars[[i]] <- lambda.values
                #pars <- expand.grid(pars)
                #for(i in 1:length(mod)) {
                #for(j in 1:nrow(pars)) {
                    tlambda <- lambda.start
                    tmp2 <- getArrayIndex(j,lambda.dim)
                    for(k in 1:length(tmp2)) tlambda <- replace(tlambda,lambda==lambda.id[k],lambda.values[tmp2][k])
                    #ll <- c(ll,myLogLike(mod,response,tlambda))
                    ll[matrix(c(tmp,tmp2),nrow=1)] <- myLogLike(mod,response,tlambda)
                }
            } else {
                #ll <- c(ll,myLogLike(mod,response,lambda.start))
                ll[matrix(tmp,nrow=1)] <- myLogLike(mod,response,lambda.start)
            }
        }
     } else {
        mod <- discrBayes.nmcpl(formula=formula,data=data,grid=grid,beta=beta,mean=mean,sd=sd.start)
        if(length(lambda.id) > 0) {
            for(j in 1:prod(lambda.dim)) {
            #pars <- list() #length=length(lambda.id))
            #for(i in 1:length(lambda.id)) pars[[i]] <- lambda.values
            #pars <- expand.grid(pars)
            #for(i in 1:length(mod)) {
            #for(j in 1:nrow(pars)) {
                tlambda <- lambda.start
                tmp2 <- getArrayIndex(j,lambda.dim)
                for(k in 1:length(tmp2)) tlambda <- replace(tlambda,lambda==lambda.id[k],lambda.values[tmp2][k])
                #ll <- c(ll,myLogLike(mod,response,tlambda))
                ll[matrix(tmp2,nrow=1)] <- myLogLike(mod,response,tlambda)
            }
        } else {
            ll <- myLogLike(mod,response,lambda.start)
        }
     }
     #ll <- aperm(array(ll,dim=c(lambda.dim,sd.dim)))
     #tmp <- getArrayIndex(which.min(ll)[1],dim=c(lambda.dim,sd.dim))
     #tmp <- tmp[length(tmp):1]
     tmp <- which(ll==min(ll),arr.ind=TRUE)
     #tmp <- which(ll==min(ll),arr.ind=TRUE)
     for(i in 1:length(par.id)) {
        sd.start <- replace(sd.start,sd==par.id[i],sd.values[tmp[i]])
        lambda.start <- replace(lambda.start,lambda==par.id[i],lambda.values[tmp[i]])
     }
     #res <- vector()
     #for(i in sd.id) {
        #res[i] <- sd.values[tmp[i]]
     #   sd.start <- replace(sd.start,sd==i,sd.values[tmp[1,i]])
     #}
     #for(i in lambda.id) {
        #res[i] <- lambda.values[tmp[i]]
     #   lambda.start <- replace(lambda.start,lambda==i,lambda.values[tmp[1,i]])
     #}
     #if(!is.null(par1)) {
     #   if(!is.null(pars)) {
     #       par1 <- expand.grid(c(par1,pars))
     #   }
     #} else {
     #   par1 <- pars
     #}
     #par <- expand.grid(c(par,pars))
     #par1$ll <- ll
     #res <- par1[which.min(par1$ll),]
     
     #for(i in 1:length(sd.id)) sd.start <- replace(sd.start,sd==sd.id[i],as.numeric(res[i]))
     #for(i in sd.id) sd.start <- replace(sd.start,sd==i,as.numeric(res[i]))
     #for(i in 1:length(lambda.id)) lambda.start <- replace(lambda.start,lambda==lambda.id[i],as.numeric(res[i+length(sd.id)]))
     #for(i in lambda.id) lambda.start <- replace(lambda.start,lambda==i,as.numeric(res[i]))

     return(list(sd.start=sd.start,lambda.start=lambda.start))
}

myLogLike <- function(fit,response,lambda) {
        ntrial <- attr(fit,"ntrial")
        tlambda <- vector()
        for(i in 1:length(ntrial)) tlambda <- c(tlambda,rep(lambda[i],ntrial[i]))
        fit$predict <- fit$predict^tlambda #THIS GOES WRONG!
        fit$predict <- fit$predict/rowSums(fit$predict)
        res <- -sum(log(rowSums(model.matrix(~response-1)*fit$predict)))
        return(res)
}