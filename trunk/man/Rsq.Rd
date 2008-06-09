\name{Rsq}
\docType{genericFunction}
\alias{Rsq}
\alias{Rsq,ANY-method}
\title{R-Squared}
\description{
  This generic function calculates the squared multiple correlation, or
  'proportion variance explained', according to the formula
  \deqn{R^2 = 1 - \frac{SS_{obj}}{SS_{tot}}}{R^2 = 1 - SS(obj)/SS(tot)},
  where \eqn{SS_{obj} = \sum_{i=1}^N (y_i - \hat{y}_y)^2}{SS(obj) = sum (y -
  pred(y))^2} is the squared error of predictions \eqn{\hat{y}}{pred(y)}, and
  \eqn{SS_{tot} = \sum_{i=1}^N (y_i - \bar{y}_y)^2}{SS(tot) = sum (y -
  mean(y))^2} the total sum of squares.
}
\usage{
Rsq(object, \dots)
}
\arguments{
  \item{object}{An object of a suitable class for the R-Squared to be
    calculated - an object with a predict method, and if y is missing, a slot
    \code{y}.
  }
  \item{\dots}{Some methods for this generic function accept an argument y,
    providing the observed data. At present, no other additional arguments are
    used}
}
\value{
  Returns a numeric value with the corresponding R-Squared.
}

\author{Maarten Speekenbrink}
%\note{}
%\examples{

%}
% \keyword{models}

