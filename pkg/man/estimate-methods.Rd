\name{fit-methods}
\docType{methods}
\alias{fit-methods}
\alias{fit,ANY-method}
%\alias{fit,GaussianMixtureResponse-method}
%\alias{fit,GaussianResponse-method}
%\alias{fit,GCM-method}
%\alias{fit,McplBaseModel-method}
%\alias{fit,McplModel-method}
%\alias{fit,RatioRuleResponse-method}
\title{Fit a model (estimate model parameters).}
\description{
The \code{fit} method fits the object to data (i.e., estimates an object's
parameters). For objects deriving from McplBaseModel, default estimation is by
Maximum Likelihood if a \code{logLik} method is defined for the object. If not,
Least Squares estimation is attempted.
}
\section{Methods}{
\describe{

\item{object = "ANY"}{Generic function.}

}}
\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
