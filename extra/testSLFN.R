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
	nt = as.integer(nt),
	eta = as.double(eta),
	actfun = as.integer(3),
	w = as.double(w),
	ypred = as.double(ypred))
	result <- list(w=array(output$w,dim=c(nx+1,ny,nt+1)),ypred=matrix(output$ypred,nrow=ny,ncol=nt))
	
	return(result)
	
	dyn.unload("~/Documents/RForge/mcplr/trunk/src/slfn.so")
}

