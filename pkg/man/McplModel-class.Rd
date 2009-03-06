\name{McplModel-class}

\docType{class}

\alias{McplModel-class}
\alias{getReplication,McplModel-method}
\alias{ntimes,McplModel-method}
\alias{AIC,McplModel-method}
\alias{AICc,McplModel-method}
\alias{BIC,McplModel-method}
\alias{estimate,McplModel-method}
\alias{fit,McplModel-method}
\alias{getPars,McplModel-method}
\alias{logLik,McplModel-method}
\alias{Rsq,McplModel-method}
\alias{predict,McplModel-method}
\alias{setPars,McplModel-method}
\alias{show,McplModel-method}
\alias{summary,McplModel-method}



\title{Class "McplModel"}

\description{An MCPL model, containing a learning and response model.}

\section{Slots}{

	\describe{

	\item{\code{learningModel}:}{A \code{\link{LearningModel-class}} object.}

	\item{\code{learningModel}:}{A \code{\link{ResponseModel-class}} object.}

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
    \item{setPars}{\code{signature(object = "McplModel")}: Set parameter values.}
    \item{show}{\code{signature(object = "McplModel")}: Display object briefly.}
    \item{summary}{\code{signature(object = "McplModel")}: Generate object summary.}
  }
}
\author{Maarten Speekenbrink}

\keyword{classes}
