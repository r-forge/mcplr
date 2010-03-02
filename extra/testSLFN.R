slfn <- function(y,x){

	y <- rbinom(10,size=1,prob=.7)
	x <- rbind(1,matrix(rnorm(10*3),nrow=3))
	eta <- runif(4)
	
	dyn.load("~/Documents/RForge/mcplr/trunk/src/slfn.so")
	
	w <- matrix(0.0,ncol=ncol(x)+1,nrow=nrow(x))
	
	output <- .C("slfn_logistic_1y",
	y = as.integer(y),
	x = as.double(x),
	nx = as.integer(nrow(x)),
	nt = as.integer(ncol(x)),
	eta = as.double(eta),
	w = as.double(w))
	result <- matrix(output$w,ncol=ncol(w),nrow=nrow(w))
	
	return(result)
	
	dyn.unload("~/Documents/RForge/mcplr/trunk/src/slfn.so")
}

