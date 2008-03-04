mcpl <- function(formula,data,design=NULL) {
    if(is.null(design)) {
    # convert a data frame to an mcpl object
        mf <- model.frame(formula,data)
        rn <- colnames(model.response(mf))
        if(length(rn)!=2) stop("formula should be 'cbind(x,y) ~ ...'")
        attr(mf,"criterion") <- rn[1]
        attr(mf,"response") <- rn[2]
        class(mf) <- c("mcpl",class(mf)) # add class "mcpl"
    } else {
        # factor(interaction(cbind(dmcpl$id,factor(dmcpl$con))))

    }
    mf
}

rmodel <- function(x) {
    # extract response model
    form <- formula(attr(x,"terms"))
    stopifnot(inherits(form, "formula"),
            deparse(form[[1]]) == "~")
    form[[2]] <-  as.name(attr(x,"response"))
    model.frame(form,cbind(x,model.response(x)))
}

cmodel <- function(x) {
    # extract response model
    form <- formula(attr(x,"terms"))
    stopifnot(inherits(form, "formula"),
            deparse(form[[1]]) == "~")
    form[[2]] <-  as.name(attr(x,"criterion"))
    model.frame(form,cbind(x,model.response(x)))
}
