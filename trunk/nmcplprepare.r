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