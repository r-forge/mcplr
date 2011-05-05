\name{GlmResponse-class}
\docType{class}
\alias{GlmResponse-class}
\alias{fit,GlmResponse-method}
\alias{predict,GlmResponse-method}
\alias{logLik,GlmResponse-method}
\title{Class "GlmResponse"}
\description{Use a Generalized Linear Model to model responses.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GlmResponse", ...)}.
}
\section{Slots}{
	 \describe{
	  \item{\code{family}:}{A \code{"family"} object. See \code{\link[stats]{glm}}}
	  \item{\code{z}:}{A \code{"matrix"} containing the values of the covariate(s).}
	  \item{\code{combine}:}{A \code{"function"} to combine predictions from a \code{"\linkS4class{LearningModel}"} and covariate(s).}
	  \item{\code{x}:}{A \code{"matrix"}, containing values of the response options.}
    \item{\code{y}:}{A \code{"matrix"}, containing dummy coded response variable.}  
    \item{\code{parStruct}:}{A \code{"ParStruct"}, specifying parameters and (optional) constraints.}
    \item{\code{nTimes}:}{Object of class \code{"NTimes"}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ResponseModel}"}, directly.
Class \code{"\linkS4class{McplBaseModel}"}, by class "ResponseModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{fit}{\code{signature(object = "GlmResponse")}: ... }
    \item{predict}{\code{signature(object = "GlmResponse")}: ... }
	 }
}
\author{Maarten Speekenbrink}
\keyword{classes}
