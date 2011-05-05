\name{MaxResponse}
\docType{class}
\alias{MaxResponse}

\title{Maximising responses with error}
\description{Construct a \code{\linkS4class{MaxResponse}} object.}

\usage{
MaxResponse(formula,data,parameters=list(beta=2.944439),ntimes=NULL,
  replicate=TRUE,fixed,parStruct,subset)
}

\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to
 that class): a symbolic description of the model to be fitted. For more details
 of model specification, see \code{\link{lm}} or \code{\link{glm}}.}
\item{data}{(optional) data frame.}
\item{parameters}{an (optional) named list with (starting) values of the 
parameters. The model has one parameter (beta), which sets the error rate 
(see details). If no values are supplied, defaults are used. }
\item{ntimes}{an optional vector with, for each repetition in the data, the
total number of trials.}
\item{replicate}{logical to indicate whether model parameters are identical for 
each replication in ntimes.}
\item{fixed}{a logical vector indicating whether model parameters are fixed or 
free}
\item{parStruct}{a ParStruct object. If supplied, the fixed argument above will 
be ignored}
\item{subset}{(optional) subset.}
}
\details{The \code{MaxResponse} function sets up a 
\code{\linkS4class{MaxResponse}} model. This model predicts categorical 
responses R = 1,...,k on the basis of predictors x[j], j = 1,...,k, as
P(R=j) = 1 - exp(beta)/(1+exp(beta))
if j = arg max x[j]
and
P(R=j) = exp(beta)/(1+exp(beta))
otherwise.
}

\value{A (fitted) object of class \code{\linkS4class{MaxResponse}}.}

\keyword{models}
