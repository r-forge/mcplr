\name{GaussianResponse-class}
\docType{class}
\alias{GaussianResponse-class}
\alias{estimate,GaussianResponse-method}
\alias{logLik,GaussianResponse-method}
\alias{predict,GaussianResponse-method}

\title{Class "GaussianResponse"}
\description{A simple Gaussian response model.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GaussianResponse", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{x}:}{A \code{"matrix"}, containing values of the response options.}
    \item{\code{y}:}{Object of class \code{"matrix"}, containing dummy coded response variable.}
    \item{\code{parameters}:}{A (named) \code{"list"} with parameters. }
    \item{\code{parStruct}:}{Object of class \code{"ParStruct"}, specifying parameter constraints.}
    \item{\code{nTimes}:}{Object of class \code{"NTimes"}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ResponseModel}"}, directly.
Class \code{"\linkS4class{McplBaseModel}"}, by class "ResponseModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(object = "GaussianResponse")}: ... }
    \item{logLik}{\code{signature(object = "GaussianResponse")}: ... }
    \item{predict}{\code{signature(object = "GaussianResponse")}: ... }
	 }
}
%\references{Luce, R. D. (1959) Individual choice behavior: A theoretical analysis. New York: John Wiley and Sons.}

\author{Maarten Speekenbrink}

% \note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\keyword{classes}
