lensmodel <- function(x,family=gaussian()) {
    rmod <- rmodel(x)
    rglm <- glm(formula(attr(rmod,"terms")),family=family,data=rmod)
    cmod <- cmodel(x)
    cglm <- glm(formula(attr(cmod,"terms")),family=family,data=cmod)
    out <- list()
    if(family$family=="gaussian") {
        out$Ra <- cor(model.response(rmod),model.response(cmod))
        out$G <- cor(predict(rglm),predict(cglm))
        out$C <- cor(residuals(rglm),residuals(cglm))
        out$Re <- cor(model.response(cmod),predict(cglm))
        out$Rs <- cor(model.response(rmod),predict(rglm))
        beta.e <- cglm$coefficients*apply(model.matrix(formula(attr(cmod,"terms")),cmod),2,"sd")/sd(model.response(cmod))
        beta.s <- rglm$coefficients*apply(model.matrix(formula(attr(rmod,"terms")),rmod),2,"sd")/sd(model.response(rmod))
        out$b.e <- cglm$coefficients
        out$b.s <- rglm$coefficients
        out$validity <- beta.e[names(beta.e)!="(Intercept)"]
        out$utilization <- beta.s[names(beta.s)!="(Intercept)"]
    }
    class(out) <- c("lensmodel",class(out))
    out
}

print.lensmodel <- function(x,digits=max(3, getOption("digits") - 3)) {
    print(data.frame(r=c(x$Ra,x$G,x$C),row.names=c("Achievement:","Linear knowledge:","Unmodelled knowledge:")),digits=digits)
    print(data.frame(rbind(x$validity,x$utilization),row.names=c("validity (beta):","utilization (beta):")),digits=digits)
}
    