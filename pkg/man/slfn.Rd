\name{SLFN}
\alias{SLFN}
\alias{SLFN-class}
\alias{fit,SLFN-method}
\alias{runm,SLFN-method}
\alias{getPars,SLFN-method}
\alias{predict,SLFN-method}
\alias{setPars,SLFN-method}

\title{Single Layer Feedforward Network}
\description{
Constructs a SLFN, an artificial neural network with a single input and
output layer, and no hidden layer.
}
\usage{
SLFN(formula,parameters=list(eta=.01,alpha=0,beta=0,ws=0),
type=c("linear","logistic"),data,subset,fixed,parStruct,window.size=0,
remove.intercept=FALSE,base=NULL,ntimes=NULL,replicate=T)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to
 that class): a symbolic description of the model to be fitted. For more details
 of model specification, see \code{lm} or \code{glm}.}
\item{parameters}{a named list with (starting) values for the parameters. 
See details.}
\item{type}{the name of the activation function. Currently, only linear and
logistic activation functions are implemented. See Example on how to use other
activation functions.}
\item{data}{(optional) data frame.}
\item{subset}{(optional) subset.}
\item{fixed}{(optional) logical vector indicating whether parameters are fixed 
(TRUE) or freely estimable (FALSE).}
\item{parStruct}{(optional) ParStruct object. Note that if parStruct is given, 
the `fixed' argument above will be ignored.}
\item{window.size}{an integer >= 0 specifying the number of previous data points
used to compute the gradient (in addition to the current data point).}
\item{remove.intercept}{should an intercept be included in the predictors? If 
set to TRUE, an intercept term (if present) will be removed from x. This can 
be useful when the predictors are categorical variables, where removing the 
intercept in the formula will result in an extra dummy variable.}
\item{base}{if the criterion (rhs of formula) is a factor, an overparametrized
dummy coding will be used by default. That is, for a criterion with n levels, a
dummy matrix will be used with n columns. By setting base to an integer k,
1 <= k <= n, column k will be removed from the matrix.}
\item{ntimes}{an optional vector with, for each repetition in the data, the 
total number of trials.}
\item{replicate}{are the repeated series true replications, i.e., are the model
parameters considered identical for each series?}
}
\details{The \code{SLFN} function sets up a simple ANN useful for deriving
online model predictions etc.}
\value{A (fitted) object of class \code{SLFN} extending \code{LearningModel}}
\references{
Bishop, C.M. (1995). \emph{Neural Networks for Pattern Recognition}. Oxford, UK: Oxford University Press.

Gluck, M. A., & Bower, G. H. (1988). From conditioning to category learning: 
An adaptive network model. \emph{Journal of Experimental Psychology: General}, 
\emph{117}, 227--247.
}
\examples{
## open weather prediction data
data(WPT)
## initialize model
mod <- SLFN(y~x1+x2+x3+x4-1,type="logistic",data=WPT,ntimes=c(200,200),fix=list(beta=TRUE,ws=TRUE))
## estimate free parameters
\dontrun{mod <- fit(mod)}
summary(mod)

## TODO: add other activation function

##
}
\keyword{models}