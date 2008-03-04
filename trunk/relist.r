# relist.R -- an inverse operator to unlist
# written by Andrew Clausen <clausen@econ.upenn.edu> in 2007
# with helpful suggestions from
#	Martin Maechler <maechler@stat.math.ethz.ch>
#	Gabor Grothendieck <ggrothendieck@gmail.com>
#	Seth Falcon, sfalcon (near) fhcrc (comma) org
#
# Some functions need many parameters, which are most easily represented in
# complex structures.  Unfortunately, many mathematical functions in R,
# including optim, nlm, and grad can only operate on functions whose domain is
# a vector.  R has a function called "unlist" to convert complex objects into a
# vector representation.  This file provides an inverse operation called
# "relist" to convert vectors back to the convenient structural representation.
# Together, these functions allow structured functions to have simple
# mathematical interfaces.
#
# For example, a likelihood function for a multivariate normal model needs a
# variance-covariance matrix and a mean vector.  It would be most convenient to
# represent it as a list containing a vector and a matrix.  A typical parameter
# might look like
#
#	list(mean=c(0, 1), vcov=cbind(c(1, 1), c(1, 0)))
#
# However, optim can't operate on functions that take lists as input; it
# only likes vectors.  The solution is conversion:
#
#	initial.param <- list(mean=c(0, 1), vcov=cbind(c(1, 1), c(1, 0)))
#	initial.param <- as.relistable(initial.param)
#
#	ll <- function(param.vector)
#	{
#		param <- relist(initial.param)
#		-sum(dnorm(x, mean=param$mean, vcov=param$vcov, log=TRUE))
#		# note: dnorm doesn't do vcov... but I hope you get the point
#	}
#
#	optim(unlist(initial.param), ll)
#
# "relist" takes two parameters: skeleton and flesh.  Skeleton is a sample
# object that has the right "shape" but the wrong content.  "flesh" is a vector
# with the right content but the wrong shape.  Invoking
#
# 	relist(flesh, skeleton)
#
# will put the content of flesh on the skeleton.  You don't need to specify
# skeleton explicitly if the skeleton is stored as an attribute inside flesh.
# In particular, flesh was created from some object obj with
#
#	unlist(as.relistable(obj))
#
# then the skeleton attribute is automatically set.
#
# As long as "skeleton" has the right shape, it should be a precise inverse
# of unlist.  These equalities hold:
#
#	relist(unlist(x), skeleton) == x
#	unlist(relist(y, skeleton)) == y
#
#	x <- as.relistable(x)
#	relist(unlist(x)) == x


is.relistable <- function(obj) is(obj, "relistable")

as.relistable <- function(x)
{
	if (!inherits(x, "relistable"))
		class(x) <- c("relistable", class(x))
	x
}
 
unlist.relistable <- function(x, recursive=TRUE, use.names=TRUE)
{
	if (!recursive)
		warning("relist() requires recursively unlisted objects.")
	skeleton <- x
	class(x) <- setdiff(class(x), "relistable")
	result <- unlist(x, recursive, use.names)
	attr(result, "skeleton") <- skeleton
	result
}

relist <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
	if (is.null(skeleton))
	{
		stop(paste(c(
"The flesh argument does not contain a skeleton attribute.  ",
"Either ensure you unlist a relistable object, or specify the skeleton ",
"separately.")))
	}
	UseMethod("relist", skeleton)
}

relist.list <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
	index <- 1
	result <- skeleton
	for(i in 1:length(skeleton)) {
		size <- length(unlist(result[[i]]))
		if(size>0) result[[i]] <- relist(flesh[index:(index + size - 1)], result[[i]])
		index <- index + size
	}
	result
}

relist.numeric <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
	result <- flesh
	names(result) <- names(skeleton)
	result
}

relist.matrix <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
	if (is.numeric(skeleton[1,1]))
		return(matrix(flesh, nrow=nrow(skeleton),
			      dimnames=dimnames(skeleton)))
	n = nrow(skeleton)
	m = ncol(skeleton)
	result <- skeleton
	index <- 1
	for (j in 1:m) for (i in 1:n) 
	{
		size <- length(unlist(skeleton[[i, j]]))
		result[[i, j]] <- relist(flesh[index:(index + size - 1)],
					 skeleton[[i, j]])
		index <- index + size
	}
	result
}

relist.factor <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
	l <- levels(skeleton)
	as.factor(l[flesh])
}

relist.pdMat <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
	result <- flesh
	class(result) <- class(skeleton)
	attr(result,"names") <- NULL
	attr(result,"Dimnames") <- attr(skeleton,"Dimnames")
	result
}

relist.NULL <- function(flesh, skeleton=attr(flesh, "skeleton")) {
  result <- NULL
  result
}