\name{RatioRuleResponse}
\docType{class}
\alias{RatioRuleResponse}
\alias{ExpRatioRuleResponse}
\alias{PowRatioRuleResponse}
\alias{RatioRuleResponse-class}
\alias{fit,RatioRuleResponse-method}
\alias{logLik,RatioRuleResponse-method}
\alias{predict,RatioRuleResponse-method}

\title{Class "RatioRuleResponse"}
\description{Responses according to Luce's ratio response rule.}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{\link{RatioRuleResponse}} 
or by \code{new("RatioRuleResponse", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{x}:}{A \code{"matrix"}, containing values of the response options.}
    \item{\code{y}:}{Object of class \code{"matrix"}, containing dummy coded response variable.}
    \item{\code{parStruct}:}{Object of class \code{"ParStruct"}, specifying parameters and (optional) constraints.}
    \item{\code{nTimes}:}{Object of class \code{"NTimes"}.}
  }
}
\section{Extends}{
Class \code{\linkS4class{ResponseModel}}, directly.
Class \code{\linkS4class{McplBaseModel}}, by class "ResponseModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{fit}{\code{signature(object = "RatioRuleResponse")}: ... }
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
