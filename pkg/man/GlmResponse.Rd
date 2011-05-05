\name{GlmResponse}
\alias{GlmResponse}

\title{Generate a Generalized Linear Model (GLM) Response model}
\description{
Constructs a \code{\linkS4class{GlmResponse}} object.
}
\usage{
GlmResponse(formula,family=gaussian(),data,covariate,
combine=function(x,z) { cbind(x,z) },parameters=list(),ntimes=NULL,
replicate=TRUE,fixed,base=NULL,parStruct,subset)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to
 that class): a symbolic description of the model to be fitted. For more details
 of model specification, see \code{\link{lm}} or \code{\link{glm}}.}
\item{family}{a \code{\link{family}} object. For more details, see 
 \code{\link{glm}}.}
\item{data}{(optional) data frame.}
\item{covariate}{(optional) a matrix or a formula specifying additional 
covariates. See details.}
\item{combine}{a function taking two arguments which produces a new predictor
matrix for the GlmModel. The first argument is the prediction from a 
LearningModel and the second argument is matrix specified through the covariate
argument.}
\item{parameters}{an (optional) named list with (starting) values of the 
parameters. If no values are supplied, defaults are used. }
\item{ntimes}{an optional vector with, for each repetition in the data, the
total number of trials.}
\item{replicate}{logical to indicate whether model parameters are identical for 
each replication in ntimes.}
\item{fixed}{a logical vector indicating whether model parameters are fixed or 
free}
\item{base}{assign one of the levels of the criterion variable the role of base 
category.}
\item{parStruct}{a ParStruct object. If supplied, the fixed argument above will 
be ignored}
\item{subset}{(optional) subset.}
}
\details{The \code{GlmResponse} function sets up a GLM to use as ResponseModel.
If the returned \code{\linkS4class{GlmResponse}} is to be used as part of an
\code{\linkS4class{McplModel}}, then the formula should contain an entry
for the predictions of the learningModel in that McplModel. When fitting the
McplModel, the predictions of the learningModel at the current iteration will be 
used. It is therefore necessary to give the name of the predictions of the 
learning model in the formula, so that this variable can be replaced with
the current variable in the estimation. Ideally, the variable specified by
\code{learnmmodpred} should not be contained in the data.frame specified 
by \code{data}.}

\value{A (fitted) object of class \code{\linkS4class{GlmResponse}}.}

\keyword{models}

\author{Maarten Speekenbrink}
