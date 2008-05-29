### test
library("stats4")
load("data/WPT.rda")
source("R/allGenerics.R")
source("R/classes.R")
source("R/mcplprepare.R")
source("R/relist.R")
source("R/RatioRuleResponse.R")
source("R/GaussianMixtureResponse.R")
source("R/GaussianResponse.R")
source("R/GCM.R")
source("R/SLFN.R")
source("R/RescorlaWagner.R")
source("R/RBFN.R")

# test Rescorla Wagner
dat <- subset(WPT,id=="C")
lmod <- rescorlaWagner(y~x1+x2+x3+x4,data=dat,parameters=list(alpha=.1,beta=c(1,1),lambda=c(1,-1)),intercept=FALSE,base=1,fix=list(beta=FALSE,ws=TRUE,lambda=TRUE))
rmod <- ratioRuleResponse(r~1,data=dat,base=1,ntimes=200)
rmod <- estimate(rmod)
tMod <- new("McplModel",
  learningModel = lmod,
  responseModel = rmod)
tmp <- estimate(tMod)

# test SLFN
dat <- subset(WPT,id=="C")
lmod <- slfn(y~x1+x2+x3+x4,parameters=list(eta=.3,alpha=.001),type="logistic",data=dat,intercept=FALSE,base=1)
#rbfnmod <- rbfn(y~x1+x2+x3+x4,eta=.3,location=lapply(apply(expand.grid(rep(list(c(0,1)),4)),1,list),unlist),scale=rep(list(diag(4)),16),type="logistic",data=dat,intercept=FALSE,base=1)
lmod <- fit(lmod)
lmod@parStruct@fix <- c(FALSE,FALSE,TRUE,TRUE)
rmod <- ratioRuleResponse(r~predict(lmod)-1,data=dat,base=1,ntimes=200)
rmod <- estimate(rmod)
logLik(rmod)
tMod <- new("McplModel",
  learningModel = lmod,
  responseModel = rmod)
tmp <- estimate(tMod)

# now with a redundant parameterisation
dat <- subset(WPT,id=="C")
mod <- slfn(y~x1+x2+x3+x4,parameters=list(eta=.3),type="logistic",data=dat,intercept=FALSE)
mod <- fit(mod)
lcm <- new("LuceChoiceModel")
lcm@y <- model.matrix(~dat$r-1)
lcm@x <- predict(mod)
lcm@parameters <- list(beta=1)
lcm@family <- binomial()
lcm <- estimate(lcm)
logLik(lcm)
mod@parameters <- list(eta=.3,alpha=.01,beta=.01,ws=0)
mod@parStruct@fix <- c(FALSE,FALSE,FALSE,TRUE)
tMod = new("McplModel",
  learningModel = mod,
  responseModel = lcm)
tmp <- estimate(tMod)

# test GCM
dat <- subset(WPT,id=="C")
mod <- gcm(y~x1+x2+x3+x4,data=dat)
logLik(mod)
mod <- estimate(mod,unconstrained=TRUE,method="Nelder-Mead")
logLik(mod)
logLik(setPars(mod,c(0.31243147,0.32985646,0.06847703,0.28923503,6.044744)))
# use constrained estimation
nw <- 3
A <- bdiag(list(rbind(diag(nw),rep(-1,nw)),1))
b <- c(rep(0,nw),-1-1e-10,0)
mod@parStruct@constraints <- new("LinConstraintsList",
  Amat = A, bvec=b)
mod@parameters <- list(w=rep(.25,3),lambda=1)
mod <- estimate(mod)
logLik(mod)
# with continuous data
load(paste(projectFolder,"StockMarket/R/study1/dmcpl.RData",sep=""))
dat <- subset(dmcpl,id==2)
mod <- gcm(dy~c1+c2-1,data=dat,level="interval",distance=gcm.distance("euclidian"),similarity=gcm.similarity("gaussian"),sampling=gcm.sampling("exponential"))
logLik(fit(setPars(mod,c(.5,.5,.01,1,3))))
mod <- estimate(fit(setPars(mod,c(.5,.5,.01,1,3))),unconstrained=T)
logLik(mod)
mod <- gcm(dy~c1+c2-1,data=dat,level="interval",distance=gcm.distance("euclidian"),similarity=gcm.similarity("gaussian"),sampling=gcm.sampling("exponential"))
mod@parStruct@fix[5] <- FALSE
mod <- estimate(fit(setPars(mod,c(.5,.5,.01,1,3))),unconstrained=T)

mod <- gcm(dy~c1+c2-1,data=dat,level="interval",distance=gcm.distance("euclidian"),similarity=gcm.similarity("gaussian"),sampling=gcm.sampling("exponential"))
mod <- fit(setPars(mod,c(.5,.5,.01,1,3)))
logLik(mod)
mod <- estimate(mod,unconstrained=T)
rMod <- new("GaussianMixtureModel")
rMod@means[[1]] <- t(matrix(mod@y,nrow=nrow(mod@y),ncol=ncol(mod@y)))
rMod@y <- matrix(dat$r,ncol=1)
rMod@parameters <- list(sd=10)
rMod@weights <- mod@weights
rMod@parStruct@replicate <- TRUE
tMod <- new("McplModel",
  learningModel=mod,
  responseModel=rMod)
tmp <- estimate(tMod,unconstrained=TRUE,CML=TRUE) # slow due to additional em steps
tmp2 <- estimate(tMod,unconstrained=TRUE,CML=FALSE) # faster: but, this is a relatively simple problem!

dat <- subset(dmcpl,id %in% c(2,3,4,5))
mod <- gcm(dy~c1+c2-1,data=dat,level="interval",distance=gcm.distance("euclidian"),similarity=gcm.similarity("gaussian"),sampling=gcm.sampling("exponential"),ntimes=c(300,300))
mod <- fit(setPars(mod,c(.75,.25,.005,.12,15)),ntimes=rep(300,4))
logLik(mod)
mod <- estimate(mod,unconstrained=T,ntimes=rep(300,4))

mod <- gcm(dy~c1+c2-1,data=dat,level="interval",distance=gcm.distance("euclidian"),similarity=gcm.similarity("gaussian"),sampling=gcm.sampling("uniform"),sd=3)
mod <- fit(setPars(mod,c(.5,.5,.01,3)))
logLik(mod)
mod <- estimate(mod,unconstrained=T)
logLik(mod)
rMod <- new("GaussianMixtureModel")
rMod@x <- t(matrix(mod@y,nrow=nrow(mod@y),ncol=ncol(mod@y)))
rMod@y <- matrix(dat$r,ncol=1)
rMod@parameters <- list(sd=10)
rMod@weights <- mod@weight
tMod <- new("McplModel",
  learningModel=mod,
  responseModel=rMod)
tmp <- estimate(tMod,unconstrained=TRUE)

#         w1          w2      lambda       gamma         sdy          sd
# 0.74686547  0.25313453  0.00572873  0.08722437  3.26981779 27.86448941

load(paste(projectFolder,"StockMarket/R/study1/dmcpl.RData",sep=""))
dat <- subset(dmcpl,id==2)
mod <- rbfn(dy~c1+c2,type="linear",data=dat,
  location=lapply(apply(expand.grid(rep(list(c(-75,-50,-25,0,25,50,75)),2)),1,list),unlist),
  scale=rep(list(50*diag(2)),length=49),eta=.005,beta=0)
tmp <- fit(mod)
plot(cbind(dat$dy,predict(tmp)))
cor(dat$dy,predict(tmp))

