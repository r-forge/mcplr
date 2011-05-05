\name{ParStruct}
\alias{ParStruct}
\title{Parameter structure objects}
\description{Construct a Construct a \code{\linkS4class{ParStruct}} object.}
\usage{
ParStruct(parameters,replicate=TRUE,fixed=NULL,ntimes=NULL,constraints=NULL)
}
\arguments{
\item{parameters}{A (named) list with parameter values}
\item{replicate}{logical. Indicates whether parameters are replicated for 
repeated data structures. If FALSE, parameters should be an unnamed list of
length equal to the number of replications (participants) in the data. Each 
element of this list should be a named list with parameter values.}
\item{fixed}{A logical vector indicating whether parameter values are fixed
(TRUE) or free (FALSE).}
\item{ntimes}{A vector specifying the number of cases (observations) for each 
replication (participant) in the data. Note: only the length of the vector is used
in this constructor function.}
\item{constraints}{A \code{"\linkS4class{ConstraintsList}"} object, specifying
parameter constraints.}
}
\details{This constructor function creates \code{"\linkS4class{ParStruct}"} 
objects.}
\value{A (fitted) object of class \code{\linkS4class{ParStruct}}.}
\author{Maarten Speekenbrink}
