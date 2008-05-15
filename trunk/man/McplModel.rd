\name{McplModel-class}

\docType{class}

\alias{McplModel-class}


\alias{getReplication}
\alias{getReplication,mcplr-method}

\alias{ntimes}
\alias{ntimes,depmix-method}

\title{Class "McplModel"}

\description{An \code{\link{mcplR}} model.}

\section{Slots}{

	\describe{

	\item{\code{learningModel}:}{A \code{\link{LearningModel}} object.}

	\item{\code{learningModel}:}{A \code{\link{ResponseModel}} object.}

	}
}

\section{Methods}{
  \describe{
    \item{AIC}{\code{signature(object = "McplModel")}: Akaike Information Criterion- based on response model.}
    \item{AICc}{\code{signature(object = "McplModel")}: Corrected Akaike Information Criterion - based on response model.}
    \item{BIC}{\code{signature(object = "McplModel")}: Bayesian Information Criterion/Schwartz Information Criterion - based on response model.}
    \item{estimate}{\code{signature(object = "McplModel")}: Estimate model parameters (by Maximum Likelihood).}
    \item{fit}{\code{signature(object = "McplModel")}: Estimate (trial dependent) model states given current parameters.}
    \item{getPars}{\code{signature(object = "McplModel")}: Get current parameter values.}
    \item{logLik}{\code{signature(object = "McplModel")}: Log-likelihood for current parameter values - based on response model.}
    \item{RSquare}{\code{signature(object = "McplModel")}: Squared correlation coefficient - a.k.a "proportion variance explained" - based on response model.}
    \item{show}{\code{signature(object = "McplModel")}: Display object briefly.}
    \item{summary}{\code{signature(object = "McplModel")}: Generate object summary.}
  }
}
\author{Maatren Speekenbrink}

\keyword{classes}
