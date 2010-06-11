require(mcplR)
data(WPT)
mod <- GCM(y~x1+x2+x3+x4,data=WPT,ntimes=c(200,200),remove.intercept=TRUE,parameters=list(w=c(.3,.2,.2,.3),lambda=2,gamma=3,r=2,q=2))
mod <- runm(mod)

dyn.load("~/Documents/RForge/mcplr/trunk/src/gcm.so")

y <- t(mod@y)
ny <- nrow(y)
x <- t(mod@x)
nx <- nrow(x)
nt <- ncol(x)
bt <- c(1,201)
et <- c(200,400)
lt <- 2
py <- matrix(0.0,ncol=nt,nrow=ny)
sim <- vector("double",length=ny)
dist <- vector("double",length=nt)
w <- c(.3,.2,.2,.3)
r <- 2
q <- 2
lambda <- -2
gamma <- 1
tmp <- .C("gcm_nominal",
	y=as.integer(y), 
	ny=as.integer(ny), 
	x=as.double(x), 
	nx=as.integer(nx), 
	bt=as.integer(bt),
	et=as.integer(et),
	lt=as.integer(lt), 
	w=as.double(w), 
	r=as.double(r), 
	g=as.double(q), 
	lambda=as.double(lambda), 
	gamma=as.double(gamma), 
	dist=as.double(dist), 
	sim=as.double(sim), 
	ypred=as.double(py)
)
t(matrix(tmp$ypred,nrow=ny)) - predict(mod,type="response")

dyn.unload("~/Documents/RForge/mcplr/trunk/src/gcm.so")

