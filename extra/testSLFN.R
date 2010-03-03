slfn <- function(y,x){

  ny <- 1
	if(ny > 1) {
	  y <- matrix(runif(100*ny),nrow=ny)
	  y <- apply(y,2,function(x) as.numeric(x == max(x)))
	} else {
	  y <- rbinom(100,size=1,prob=.7)
	}
	if(!is.matrix(y)) y <- matrix(y,ncol=length(y))
	x <- rbind(1,matrix(rnorm(100*3),nrow=3))
	eta <- matrix(runif(4*ny),ncol=ny)
	
	dyn.load("~/Documents/RForge/mcplr/trunk/src/slfn.so")
	
	w <- array(0.0,dim=c(nrow(x),nrow(y),ncol(x)+1))
	
	output <- .C("slfn",
	y = as.double(y),
	ny = as.integer(nrow(y)),
	x = as.double(x),
	nx = as.integer(nrow(x)),
	nt = as.integer(ncol(x)),
	eta = as.double(eta),
	actfun = as.integer(2),
	w = as.double(w),
	ypred = as.double(y))
	result <- list(w=array(output$w,dim=c(nrow(x),NROW(y),ncol(x))),ypred=matrix(output$ypred,ncol=ncol(x),nrow=nrow(y)))
	
	return(result)
	
	dyn.unload("~/Documents/RForge/mcplr/trunk/src/slfn.so")
}

