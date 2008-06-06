\name{RatioRuleResponse-class}
\docType{class}
\alias{RatioRuleResponse-class}
\alias{estimate,RatioRuleResponse-method}
\alias{logLik,RatioRuleResponse-method}
\alias{predict,RatioRuleResponse-method}

\title{Class "RatioRuleResponse"}
\description{Responses according to Luce's ratio response rule.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RatioRuleResponse", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{family}:}{Object of S3 class "family". See \code{\link{family}} for details.}
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
    \item{estimate}{\code{signature(object = "RatioRuleResponse")}: ... }
    \item{logLik}{\code{signature(object = "RatioRuleResponse")}: ... }
    \item{predict}{\code{signature(object = "RatioRuleResponse")}: ... }
	 }
}
\references{
Luce, R. D. (1959) \emph{Individual choice behavior: A theoretical analysis}. New York: John
Wiley and Sons.}

\author{Maarten Speekenbrink}

% \note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}

\keyword{classes}
