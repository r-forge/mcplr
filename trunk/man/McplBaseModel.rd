\name{McplBaseModel-class}

\docType{class}

\alias{McplBaseModel-class}
\alias{LearningModel-class}
\alias{ResponseModel-class}

\alias{getReplication}
\alias{getReplication,mcplr-method}

\alias{ntimes}
\alias{ntimes,depmix-method}

\title{Class "McplBaseModel"}

\description{Basic class underlying \code{LearningModel} and \code{ResponseModel}.}

\section{Slots}{

	\describe{

	\item{\code{x}:}{Matrix with predictor variables.}

	\item{\code{y}}{Matrix with criterion/response variable.}

	\item{\code{parameters}:}{(named) list with parameter values.}

	\item{\code{parStruct}:}{A \code{\link{ParStruct}} object.}

	\item{\code{nTimes}:}{An \code{\link{NTimes}} object.}

	}
}

\section{Methods}{
  \describe{
    \item{AIC}{\code{signature(object = "McplBaseModel")}: Akaike Information Criterion.}
    \item{AICc}{\code{signature(object = "McplBaseModel")}: Corrected Akaike Information Criterion.}
    \item{BIC}{\code{signature(object = "McplBaseModel")}: Bayesian Information Criterion/Schwartz Information Criterion.}
    \item{estimate}{\code{signature(object = "McplBaseModel")}: Estimate model parameters (by Maximum Likelihood).}
    \item{fit}{\code{signature(object = "McplBaseModel")}: Estimate (trial dependent) model states given current parameters.}
    \item{getPars}{\code{signature(object = "McplBaseModel")}: Get current parameter values.}
    \item{logLik}{\code{signature(object = "McplBaseModel")}: Log-likelihood for current parameter values.}
    \item{RSquare}{\code{signature(object = "McplBaseModel")}: Squared correlation coefficient ("proportion variance explained")}
    \item{show}{\code{signature(object = "McplBaseModel")}: Display object briefly.}
    \item{summary}{\code{signature(object = "McplBaseModel")}: Generate object summary.}
  }
}

\author{Maarten Speekenbrink}

\keyword{classes}
