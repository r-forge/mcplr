\name{ParStruct-class}

\docType{class}

\alias{ParStruct-class}
\alias{getPars,ParStruct-method}
\alias{setPars,ParStruct-method}
\alias{rep,ParStruct-method}



\title{Class "ParStruct"}

\description{A parameter structure class.}

\section{Slots}{

	\describe{

	\item{\code{parameters}:}{A (named) list with parameter values.}

	\item{\code{id}:}{A numerical vector with indexes for each parameter, used
  to set equality constraints for parameters}
  
	\item{\code{fix}:}{A logical vector indicating whether each parameter is
  fixed (TRUE) or estimable (FALSE) }
  
  \item{\code{constraints}}{An object of class \code{ConstraintsList}}
  
  \item{\code{replicate}}{Logical, indicates whether parameter values are 
  replicated for each replication in NTimes (TRUE), or not (FALSE)}
  
  \item{\code{transform}}{Function to transform parameters to a more suitable 
  metric for optimisation}
  
  \item{\code{inv.transform}}{Function implementing the inverse transformation
  to put parameters back to the normal metric}

	}
}

\section{Methods}{
  \describe{
    \item{getPars}{\code{signature(object = "ParStruct")}: Get current 
    parameter values.}
    \item{rep}{\code{signature(object = "ParStruct")}: Repeat suitable 
    elements.}
    \item{setPars}{\code{signature(object = "ParStruct")}: Set parameter 
    values.}
    \item{summary}{\code{signature(object = "McplModel")}: Generate object 
    summary.}
  }
}
\author{Maarten Speekenbrink}

\keyword{classes}
