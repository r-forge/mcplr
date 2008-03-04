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
    discr <- FALSE
    if(is.factor(eval(attr(mt,"variables"),envir=mf)[[1]])) {
        discr <- TRUE
    }
    Y <- model.response(mf)
    if(discr) Y <- model.matrix(~Y-1)
    X <- if(!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else matrix(, NROW(Y), 0)
    Xn <- as.matrix(X)
    x.names <- dimnames(Xn)[[2]]
    if(discr) y.names <- levels(eval(attr(mt,"variables"),envir=mf)[[1]])
    return(list(x=Xn,y=Y))
}
