\name{McplBaseModel-class}

\docType{class}

\alias{McplBaseModel-class}
\alias{LearningModel-class}
\alias{ResponseModel-class}

%\alias{getReplication}
\alias{AIC,McplBaseModel-method}
\alias{AICc,McplBaseModel-method}
\alias{BIC,McplBaseModel-method}
\alias{fit}
\alias{fit,McplBaseModel-method}
\alias{runm}
\alias{runm,McplBaseModel-method}
\alias{runm,ResponseModel-method}
\alias{runm,LearningModel-method}
\alias{getPars}
\alias{getPars,McplBaseModel-method}
\alias{getPars,McplBaseModel,ANY-method}
\alias{getReplication,McplBaseModel-method}
\alias{logLik,McplBaseModel-method}
\alias{ntimes,McplBaseModel-method}
\alias{Rsq,McplBaseModel-method}
\alias{show,McplBaseModel-method}
\alias{show,ResponseModel-method}
\alias{show,LearningModel-method}
\alias{setPars}
\alias{setPars,McplBaseModel-method}
\alias{summary,McplBaseModel-method}
\alias{summary,ResponseModel-method}
\alias{summary,LearningModel-method}
%\alias{ntimes}

\title{Class "McplBaseModel"}

\description{Basic class underlying \code{\link{LearningModel-class}} and \code{\link{ResponseModel-class}}.}

\section{Slots}{

	\describe{

	\item{\code{x}:}{Matrix with predictor variables.}

	\item{\code{y}}{Matrix with criterion/response variable.}

	\item{\code{parameters}:}{(named) list with parameter values.}

	\item{\code{parStruct}:}{A \code{\link{ParStruct-class}} object.}

	\item{\code{nTimes}:}{An \code{\link{NTimes-class}} object.}

	}
}

\section{Methods}{
  \describe{
    \item{AIC}{\code{signature(object = "McplBaseModel")}: Akaike Information Criterion.}
    \item{AICc}{\code{signature(object = "McplBaseModel")}: Corrected Akaike Information Criterion.}
    \item{BIC}{\code{signature(object = "McplBaseModel")}: Bayesian Information Criterion/Schwartz Information Criterion.}
    \item{fit}{\code{signature(object = "McplBaseModel")}: Estimate model parameters (by Maximum Likelihood).}
    \item{runm}{\code{signature(object = "McplBaseModel")}: Run model (calculate model states using current parameters).}
    \item{getPars}{\code{signature(object = "McplBaseModel")}: Get current parameter values.}
    \item{logLik}{\code{signature(object = "McplBaseModel")}: Log-likelihood for current parameter values.}
    \item{RSquare}{\code{signature(object = "McplBaseModel")}: Squared correlation coefficient ("proportion variance explained")}
    \item{show}{\code{signature(object = "McplBaseModel")}: Display object briefly.}
    \item{summary}{\code{signature(object = "McplBaseModel")}: Generate object summary.}
  }
}

\author{Maarten Speekenbrink}

\keyword{classes}
