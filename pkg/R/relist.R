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