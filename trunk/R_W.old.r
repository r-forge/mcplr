nmcpl.prepare <- function(formula,force.bias=FALSE,data,subset) {
    # Returns a matrix of dummy input variables and matrix of dummy Y variables
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
    Yn <- model.matrix(~Y-1)
    #ifelse(dim(Yn)[2]==2,Yn <- as.numeric(Yn[,2]),stop("Multinomial logistic not implemented yet"))
    X <- if(!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else matrix(, NROW(Y), 0)
    if(colnames(X)[1] != "(Intercept)" && force.bias==FALSE) {
        Xr <- X[,2:ncol(X)] # Delete first column of contrast matrix
        attr(Xr,"assign") <- attr(X,"assign")[2:ncol(X)]
        X <- Xr
        rm(Xr)
    }
    #Xn <-  matrix(as.numeric(X),nrow=dim(X)[1],ncol=dim(X)[2])
    x.names <- dimnames(X)[[2]]
    y.names <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
    rm(Y,X)
    return(list(x=Xn,y=Yn,x.names=x.names,y.names=y.names))
}


R_W.fit <- function(y,x,alpha,beta,lambda,ws) {
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    weight <- matrix(0,nrow=nrow(x),ncol=ncol(x))
    w <- ws
    for(i in 1:nrow(x)) {
        w <- w + y[i]*alpha*beta[1]*(lambda[1]-t(w)%*%x[i,])*x[i,] + (1-y[i])*alpha*beta[2]*(lambda[2]-t(w)%*%x[i,])*x[i,]
        weight[i,] <- w
    }
    dimnames(weight)[[2]] <- xnames
    return(list(weight=weight))
}

R_W <- function(formula,lambda=c(1,0),data,subset,force.bias=FALSE,alpha,beta,ws,ws.range=.05) {
    tmp <- nmcpl.prepare(formula=formula,force.bias=force.bias,data=data) # todo: include subset
    x <- tmp$x
    y <- tmp$y
    if(missing(ws)) ws <- runif(ncol(x),min=-ws.range,max=ws.range)
    if(length(ws) != ncol(x)) stop("length of ws is ",length(ws)," but should be ",ncol(x))
    if(missing(alpha)) alpha <- rep(.5,times=ncol(x))
    if(length(alpha) != ncol(x)) {
        warning("length of alpha is less than number of cues, will fill in with last value")
        alpha <- c(alpha,rep(alpha[length(alpha)],times = ncol(x) - length(alpha)))
    }
    if(missing(beta)) beta <- rep(.5,2)
    if(length(beta) < 2) {
        warning("length of beta is less than 2, will fill in with given value")
        beta <- c(beta,beta)
    }
    if(length(beta) > 2) {
        warning("length of beta is greater than 2, will cutoff last values")
        beta <- beta[1:2]
    }
    if(length(lambda) < 2) {
        warning("length of lambda is less than 2, will use 0 for last value")
        lambda <- c(lambda,0)
    }
    if(length(lambda) > 2) {
        warning("length of beta is greater than 2, will cutoff last values")
        lambda <- lambda[1:2]
    }
    fit <- R_W.fit(y=y,x=x,alpha=alpha,beta=beta,lambda=lambda,ws=ws)
    weight <- fit$weight
    pred <- diag(Xn%*%t(weight))
    pred <- cbind(1-pred,pred)
    colnames(pred) <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
    control <- list(alpha=alpha,beta=beta,lambda=lambda,ws=ws)
    return(list(call=call,response=Yn,predictor=Xn,weight=weight,predict=pred,control=control))
}