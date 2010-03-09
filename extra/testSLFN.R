require(mcplR)
data(WPT)
## initialize model
mod <- SLFN(y~x1+x2+x3+x4,type="logistic",data=WPT,ntimes=c(200,200),parameters=list(eta=.1,alpha=0),fix=list(alpha=TRUE,beta=TRUE,ws=TRUE),remove.intercept=TRUE)
mod <- runm(mod)

y <- t(mod@y)
if(is.matrix(y)) ny <- NROW(y) else ny <- 1
x <- t(mod@x)
nx <- nrow(x)
nt <- ncol(x)
bt <- c(1,201)
et <- c(200,400)
lt <- 2
#eta <- matrix(runif((nx)*ny,max=.5,min=.1),ncol=ny)
eta <- matrix(.1,nrow=nx,ncol=ny)
ypred <- matrix(0.0,nrow=ny,ncol=nt)
w <- array(0.0,dim=c(nx,ny,nt))

dyn.load("~/Documents/RForge/mcplr/trunk/src/slfn.so")
tmp <- .C("slfn",
	y=as.double(y), 
	ny=as.integer(ny), 
	x=as.double(x), 
	nx=as.integer(nx), 
	bt=as.integer(bt),
	et=as.integer(et),
	lt=as.integer(lt), 
	eta=as.double(eta), 
	actfun = as.integer(2),
	w = as.double(w),
	ypred = as.double(ypred)
)
all(t(matrix(tmp$ypred,nrow=ny)) - predict(mod,type="response") < 1e-15)
all(aperm(array(tmp$w,dim=c(nx,ny,nt)),c(3,1,2)) - mod@weight < 1e-15)

## initialize model
mod <- SLFN(y~x1+x2+x3+x4,type="logistic",data=WPT,ntimes=c(200,200),parameters=list(eta=.1,alpha=.5),fix=list(alpha=TRUE,beta=TRUE,ws=TRUE),remove.intercept=TRUE)
mod <- runm(mod)

y <- t(mod@y)
if(is.matrix(y)) ny <- NROW(y) else ny <- 1
x <- t(mod@x)
nx <- nrow(x)
nt <- ncol(x)
bt <- c(1,201)
et <- c(200,400)
lt <- 2
#eta <- matrix(runif((nx)*ny,max=.5,min=.1),ncol=ny)
eta <- matrix(.1,nrow=nx,ncol=ny)
alpha <- .5
ypred <- matrix(0.0,nrow=ny,ncol=nt)
w <- array(0.0,dim=c(nx,ny,nt))

dyn.load("~/Documents/RForge/mcplr/trunk/src/slfn.so")
tmp <- .C("slfn_a",
	y=as.double(y), 
	ny=as.integer(ny), 
	x=as.double(x), 
	nx=as.integer(nx), 
	bt=as.integer(bt),
	et=as.integer(et),
	lt=as.integer(lt), 
	eta=as.double(eta),
	alpha=as.double(alpha),
	actfun = as.integer(2),
	w = as.double(w),
	ypred = as.double(ypred)
)
all(t(matrix(tmp$ypred,nrow=ny)) - predict(mod,type="response") < 1e-15)
all(aperm(array(tmp$w,dim=c(nx,ny,nt)),c(3,1,2)) - mod@weight < 1e-15)

## initialize model
mod <- SLFN(y~x1+x2+x3+x4,type="logistic",data=WPT,ntimes=c(200,200),parameters=list(eta=.1,alpha=0,beta=.05),fix=list(alpha=TRUE,beta=TRUE,ws=TRUE),remove.intercept=TRUE)
mod <- runm(mod)

y <- t(mod@y)
if(is.matrix(y)) ny <- NROW(y) else ny <- 1
x <- t(mod@x)
nx <- nrow(x)
nt <- ncol(x)
bt <- c(1,201)
et <- c(200,400)
lt <- 2
#eta <- matrix(runif((nx)*ny,max=.5,min=.1),ncol=ny)
eta <- matrix(.1,nrow=nx,ncol=ny)
beta <- matrix(.05,nrow=nx,ncol=ny)
ypred <- matrix(0.0,nrow=ny,ncol=nt)
w <- array(0.0,dim=c(nx,ny,nt))

dyn.load("~/Documents/RForge/mcplr/trunk/src/slfn.so")
tmp <- .C("slfn_b",
	y=as.double(y), 
	ny=as.integer(ny), 
	x=as.double(x), 
	nx=as.integer(nx), 
	bt=as.integer(bt),
	et=as.integer(et),
	lt=as.integer(lt), 
	eta=as.double(eta),
	beta=as.double(beta),
	actfun = as.integer(2),
	w = as.double(w),
	ypred = as.double(ypred)
)
all(t(matrix(tmp$ypred,nrow=ny)) - predict(mod,type="response") < 1e-15)
all(aperm(array(tmp$w,dim=c(nx,ny,nt)),c(3,1,2)) - mod@weight < 1e-14)

## initialize model
mod <- SLFN(y~x1+x2+x3+x4,type="logistic",data=WPT,ntimes=c(200,200),parameters=list(eta=.1,alpha=.5,beta=.05),fix=list(alpha=TRUE,beta=TRUE,ws=TRUE),remove.intercept=TRUE)
mod <- runm(mod)

y <- t(mod@y)
if(is.matrix(y)) ny <- NROW(y) else ny <- 1
x <- t(mod@x)
nx <- nrow(x)
nt <- ncol(x)
bt <- c(1,201)
et <- c(200,400)
lt <- 2
#eta <- matrix(runif((nx)*ny,max=.5,min=.1),ncol=ny)
eta <- matrix(.1,nrow=nx,ncol=ny)
beta <- matrix(.05,nrow=nx,ncol=ny)
ypred <- matrix(0.0,nrow=ny,ncol=nt)
w <- array(0.0,dim=c(nx,ny,nt))

dyn.load("~/Documents/RForge/mcplr/trunk/src/slfn.so")
tmp <- .C("slfn_ab",
	y=as.double(y), 
	ny=as.integer(ny), 
	x=as.double(x), 
	nx=as.integer(nx), 
	bt=as.integer(bt),
	et=as.integer(et),
	lt=as.integer(lt), 
	eta=as.double(eta),
	alpha=as.double(.5),
	beta=as.double(beta),
	actfun = as.integer(2),
	w = as.double(w),
	ypred = as.double(ypred)
)
all(t(matrix(tmp$ypred,nrow=ny)) - predict(mod,type="response") < 1e-15)
all(aperm(array(tmp$w,dim=c(nx,ny,nt)),c(3,1,2)) - mod@weight < 1e-14)

	
	
slfn <- function(y,x){
 
  nt <- 1000
  ny <- 50
  nx <- 100
	if(ny > 1) {
	  y <- matrix(runif(nt*ny),nrow=ny)
	  y <- apply(y,2,function(x) as.numeric(x == max(x)))
	} else {
	  y <- rbinom(nt,size=1,prob=.7)
	}
	if(!is.matrix(y)) y <- matrix(y,ncol=length(y))
	ypred <- matrix(0.0,nrow=nrow(y),ncol=ncol(y))
	x <- rbind(1,matrix(rnorm(nt*nx),nrow=nx))
	eta <- matrix(runif((nx+1)*ny,max=.5,min=.1),ncol=ny)
	
	dyn.load("~/Documents/RForge/mcplr/trunk/src/slfn.so")
	
	w <- array(0.0,dim=c(nx,ny,nt+1))
	
	output <- .C("slfn",
	y = as.double(y),
	ny = as.integer(ny),
	x = as.double(x),
	nx = as.integer(nx),
	bt = as.integer(1),
	et = as.integer(nt),
	lt = as.integer(1),
	eta = as.double(eta),
	actfun = as.integer(3),
	w = as.double(w),
	ypred = as.double(ypred))
	result <- list(w=array(output$w,dim=c(nx+1,ny,nt+1)),ypred=matrix(output$ypred,nrow=ny,ncol=nt))
	
	return(result)
	
	dyn.unload("~/Documents/RForge/mcplr/trunk/src/slfn.so")
}

