\name{GlmResponse-class}
\docType{class}
\alias{GlmResponse-class}
\alias{estimate,GlmResponse-method}
\alias{predict,GlmResponse-method}

\title{Class "GlmResponse"}
\description{Use a Generalized Linear Model to model responses.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GlmResponse", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{x}:}{A \code{"matrix"}, containing values of the response options.}
    \item{\code{y}:}{Object of class \code{"matrix"}, containing dummy coded response variable.}
    \item{\code{parameters}:}{A (named) \code{"list"} with parameters. }
    \item{\code{parStruct}:}{Object of class \code{"ParStruct"}, specifying parameter constraints.}
    \item{\code{nTimes}:}{Object of class \code{"NTimes"}.}
    \item{\code{family}:}{A family object. See \code{\link[stats]{glm}}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ResponseModel}"}, directly.
Class \code{"\linkS4class{McplBaseModel}"}, by class "ResponseModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(object = "GlmResponse")}: ... }
    \item{predict}{\code{signature(object = "GlmResponse")}: ... }
	 }
}

\author{Maarten Speekenbrink}

\keyword{classes}
