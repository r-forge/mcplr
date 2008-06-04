\name{estimate-methods}
\docType{methods}
\alias{estimate-methods}
\alias{estimate,ANY-method}
%\alias{estimate,GaussianMixtureResponse-method}
%\alias{estimate,GaussianResponse-method}
%\alias{estimate,GCM-method}
%\alias{estimate,McplBaseModel-method}
%\alias{estimate,McplModel-method}
%\alias{estimate,RatioRuleResponse-method}
\title{Estimate model parameters.}
\description{
The \code{estimate} method estimates an objects parameters. For objects deriving 
from McplBaseModel, default estimation is by Maximum Likelihood if a 
\code{logLik} method is defined for the object. If not, Least Squares estimation
is tried.
}
\section{Methods}{
\describe{

\item{object = "ANY"}{Generic function.}

}}
\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
