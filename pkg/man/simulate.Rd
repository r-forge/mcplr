\name{simulate}
\docType{genericFunction}
\alias{simulate}
\alias{simulate,ANY-method}
\title{simulate}
\description{
  This generic function simulates observations from a model.
}
\usage{
simulate(object, nsim=1, seed=NULL, \dots)
}
\arguments{
  \item{object}{An object of a suitable class for simulation.}
  \item{nsim}{Number of times to simulate the (whole) model.}
  \item{seed}{This can be used to set the seed of the random number generator}
  \item{...}{Some methods for this generic function may accept additional
  argumments}
}
\value{
  Returns a new object of the same class as the object in the argument.
}

\author{Maarten Speekenbrink}
%\note{}
%\examples{

%}
% \keyword{models}

