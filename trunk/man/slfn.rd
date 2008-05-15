\name{slfn}
\alias{slfn}
\title{Single Layer Feedforward Network}
\description{
Constructs a general artificial neural network with a single input and
output layer.
}
\usage{
formula,parameters=list(eta=.01,alpha=0,beta=0,ws=0),
type=c("linear","logistic"),data,window.size=0,intercept=TRUE,base=NULL,
ntimes=NULL,replicate=T,subset)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to
 that class): a symbolic description of the model to be fitted. For more details
 of model specification, see \code{lm} or \code{glm}.}
\item{parameters}{a list with starting values for the parameters. See details.}
\item{type}{the name of the activation function. Currently, only linear and
logistic activation functions are implemented. See Example on how to use other
activation functions.}
\item{data}{(optional) data frame.}
\item{window.size}{an integer >= 0 specifying the number of previous data points
used to compute the gradient (in addition to the current data point).}
\item{intercept}{logical. If set to FALSE, the intercept term will be removed
if included in the model through the model formula. If the
formula specifies to remove the intercept (via -1), setting this to TRUE will
not add an intercept.}
\item{base}{if the criterion (rhs of formula) is a factor, an overparametrized
dummy coding will be used by default. That is, for a criterion with n levels, a
dummy matrix will be used with n columns. By setting base to an integer k,
1 <= k <= n, column k will be removed from the matrix.}

\details{The \code{slfn} function sets up a simple ANN useful for deriving
online model predictions etc.}
\value{A (fitted) object of class \code{SLFN} extending \code{LearningModel}}
\references{Gluck, M & Bower, G. (1988). }
\examples{
## open weather prediction data
data(WP)
## initialize model
mod <- slfn(y~x1+x2+x3+x4-1,type="logistic",data=WP)
## estimate free parameters
mod <- estimate(mod,unconstrained=T)
summary(mod,unconstrained=T)

## TODO: add other activation function

##
}
\keyword{models}