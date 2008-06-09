\name{AICc}
\docType{genericFunction}
\alias{AICc}
\alias{AICc,ANY-method}
\alias{AICc,logLik-method}
\title{Corrected Akaike Information Criterion}
\description{
  This generic function calculates the corrected Akaike information criterion,
  also known as Schwarz's Bayesian criterion (SBC), for one or several
  fitted model objects for which a log-likelihood value can be obtained,
  according to the formula \eqn{-2 \mbox{log-likelihood} + 2 n_{par}
  (n{obs}/(n_{obs} - n_{par} -2))}{-2*log-likelihood +
  2*npar*(nobs/(nobs-npar-2))}, where \eqn{n_{par}}{npar}  represents the
  number of parameters and \eqn{n_{obs}}{nobs} the number of
  observations in the fitted model.
}
\usage{
AICc(object, \dots)
}
\arguments{
  \item{object}{An object of a suitable class for the AICc to be
    calculated - usually a \code{"logLik"} object or an object for which
    a \code{\link[stats4:logLik-methods]{logLik}} method exists.
  }
  \item{\dots}{Some methods for this generic function may take
    additional, optional arguments.  At present none do.}
}
\value{
  Returns a numeric value with the corresponding AICc.
}
\references{
  Hurvich, C.M. and Tsai, C.L. (1989) "Regression and Time Series Model Selection in Small Samples", \emph{Biometrika},
\bold{76}, 297-307.
}
\author{Maarten Speekenbrink}
%\note{}
\seealso{\code{\link{logLik-methods}}, \code{\link{AIC-methods}}}
\examples{

}
% \keyword{models}

