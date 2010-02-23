\name{GaussianMixtureResponse}
\alias{GaussianMixtureResponse}

\title{Generate a Gaussian Mixture Response model}
\description{
Constructs a \code{\linkS4class{GaussianMixtureResponse}} object.
}
\usage{
GaussianMixtureResponse(formula,ncomponent=2,fixed,parStruct,data,subset,
  weights,ntimes=NULL,replicate=TRUE)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to
 that class): a symbolic description of the model to be fitted. For more details
 of model specification, see \code{\link{lm}} or \code{\link{glm}}.}
\item{ncomponent}{number of mixture components.}
\item{fixed}{a logical vector indicating whether model parameters are fixed or 
free}
\item{parStruct}{a ParStruct object. If supplied, the fixed argument above will 
be ignored}
\item{data}{(optional) data frame.}
\item{subset}{(optional) subset.}
\item{weights}{(optional) mixture weights.}
\item{ntimes}{an optional vector with, for each repetition in the data, the
total number of trials.}
\item{replicate}{logical to indicate whether model parameters are identical for 
each replication in ntimes.}
}
\details{The \code{GaussianMixtureResponse} function sets up a simple Normal
finite mixture model.}
\value{A (fitted) object of class \code{\linkS4class{GaussianMixtureResponse}}.}

\keyword{models}
