\name{MaxResponse-class}
\docType{class}
\alias{MaxResponse-class}
\alias{predict,MaxResponse-method}
\alias{logLik,MaxResponse-method}
\title{Class "MaxResponse"}
\description{A class to model maximising responses with error.}
\section{Objects from the Class}{
Objects can be created by the constructor function \code{\link{MaxResponse}} or 
directly by calls of the form \code{new("MaxResponse", ...)}. For details on 
MaxResponse models, see the documentation for \code{\link{MaxResponse}}.
}
\section{Slots}{
	 \describe{
	  \item{\code{x}:}{A \code{"matrix"}, containing values of the response 
	  options.}
    \item{\code{y}:}{A \code{"matrix"}, containing dummy coded response 
    variable.}  
    \item{\code{parStruct}:}{A \code{"\linkS4class{ParStruct}"} object, 
    specifying parameters and (optional) constraints.}
    \item{\code{nTimes}:}{A \code{"\linkS4class{NTimes}"} object, specifying the
    repeated structure of the data.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ResponseModel}"}, directly.
Class \code{"\linkS4class{McplBaseModel}"}, by class "ResponseModel", 
distance 2.
}
\section{Methods}{
  \describe{
    \item{logLik}{\code{signature(object = "MaxResponse")}: ... }
    \item{predict}{\code{signature(object = "MaxResponse")}: ... }
	 }
}
\author{Maarten Speekenbrink}
\keyword{classes}
