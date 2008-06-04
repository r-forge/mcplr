\name{GaussianMixtureResponse}
\docType{class}
\alias{GaussianMixtureResponse-class}
\alias{estimate,GaussianMixtureResponse-method}
\alias{logLik,GaussianMixtureResponse-method}
\alias{predict,GaussianMixtureResponse-method}

\title{Class "GaussianMixtureResponse"}
\description{Model responses according to a Gaussian mixture.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GaussianMixtureResponse", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{weights}:}{A \code{"list"}, with for each repetition in nTimes, a T*K matrix of mixture weights (T=total number of trials, k = number of mixture components).}
    \item{\code{means}:}{A \code{"list"}, with for each repetition in nTimes, a matrix of mixture means. When component means vary over trials, this should be a T*K matrix (T=total number of trials, k = number of mixture components), else a 1*k matrix.}
    \item{\code{sds}:}{A \code{"list"}, with for each repetition in nTimes, the standard deviation of the mixture components. Currently, only components with identical means are supported, and the values in this slot are set by a single "sd" parameter (see "parameters" slot).}
    \item{\code{x}:}{A \code{"matrix"}, containing values of the response options.}
    \item{\code{y}:}{Object of class \code{"matrix"}, containing dummy coded response variable.}
    \item{\code{parameters}:}{A (named) \code{"list"} with parameters.}
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
    \item{estimate}{\code{signature(object = "GaussianMixtureResponse")}: ... }
    \item{logLik}{\code{signature(object = "GaussianMixtureResponse")}: ... }
    \item{predict}{\code{signature(object = "GaussianMixtureResponse")}: ... }
	 }
}

\author{Maarten Speekenbrink}

% \note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("GaussianMixtureResponse")
}
\keyword{classes}
