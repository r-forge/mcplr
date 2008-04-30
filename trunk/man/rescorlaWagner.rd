\name{gcm}
\alias{gcm}
\title{Rescorla-Wagner Model}
\description{
Constructs a Rescorla-Wagner learning model.
}
\usage{
rescorlaWagner(formula,parameters=list(alpha=.1,beta=c(1,1),lambda=c(1,0),ws=0),
data,intercept=TRUE,base=NULL,ntimes=NULL,replicate=TRUE,
fixed,parStruct,subset)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to 
 that class): a symbolic description of the model to be fitted. For more details
 of model speecification, see \code{lm} or \code{glm}.}
\item{parameters}{a list with (starting) values of the parameters. Missing 
  parameters in the list will be set to default values.}
\item{data}{(optional) data frame for evaluation of the formula}
\item{intercept}{should an intercept be included in the predictors? If set to 
FALSE, an intercept will be removed from the model.frame. This is especially
useful when the predictors are categorical variables.}
\item{base}{which level of the criterion variable is considered the base 
category? Defaults to the first level.}
\item{ntimes}{an optional vector with, for each repetition in the data, the 
total number of trials.}
\item{replicate}{are the repeated series true replications, i.e., are the model
parameters considered identical for each series?}
}
\details{The Rescorla-Wagner model (Rescorla & Wagner, 1972) is an associative
learning model. 
}
\value{A (fitted) object of class \code{LearningModel}}
\references{Rescorla, Wagner (1972). }
\examples{
## open weather prediction data
data(WP)
## initialize model
mod <- gcm(y~x1+x2+x3+x4-1,data=WP)
## estimate free parameters
mod <- estimate(mod,unconstrained=T)
summary(mod,unconstrained=T)
}
\keyword{models}