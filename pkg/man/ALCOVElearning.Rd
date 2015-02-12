\name{ALCOVElearning}
\alias{ALCOVElearning}
\alias{ALCOVElearning-class}
\alias{runm,ALCOVElearning-method}
\alias{fit,ALCOVElearning-method}
\alias{getPars,ALCOVElearning-method}
\alias{getPars,ALCOVElearning,ANY-method}
\alias{predict,ALCOVElearning-method}
\alias{setPars,ALCOVElearning-method}

\title{Attention Learning COVEring map model (ALCOVE)}
\description{
Constructs the learning submodel of an Attention Learning COVEring map model.
}
\usage{
ALCOVElearning(formula,parameters=list(eta_w=.05,eta_a=.05,r=1,q=1,spf=1),humble=TRUE,
exemplar.locations,data,subset,fixed,parStruct,
random.locations=FALSE,n.locations=10,base=NULL,ntimes,replicate=TRUE)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to
 that class): a symbolic description of the model to be fitted. For more details
 of model speecification, see \code{\link{lm}} or \code{\link{glm}}.}
\item{parameters}{an (optional) list with (starting) values of the parameters.
 If no values are supplied, defaults are used. }
\item{humble}{logical. If TRUE, humble teaching signal is used.}
\item{exemplar.locations}{(optional) list with exemplar node locations. See 
details.}
\item{data}{(optional) data frame for evaluation of the formula.}
\item{subset}{(optional) subset of the data.}
\item{fixed}{(optional) logical vector indicating whether parameters are fixed 
(TRUE) or freely estimable (FALSE).}
\item{parStruct}{(optional) ParStruct object. Note that if parStruct is given, 
the `fixed' argument above will be ignored.}
\item{random.locations}{If no exemplar.locations are given, should they be 
determined randomly? If FALSE (default), then unique values of training cues 
are used as exemplar locations.}
\item{n.locations}{Number of randomly generated exemplar locations, if 
random.locations = TRUE} 
\item{base}{which level of the criterion variable is considered the base 
category? Defaults to the first level.}
\item{ntimes}{an optional vector with, for each repetition in the data, the 
total number of trials.}
\item{replicate}{are the repeated series true replications, i.e., are the model
parameters considered identical for each series?}
}

\details{ALCOVE (Kruschke, 1992) is based on the \code{\link{GCM}} model, but
has a mechanism to learn the attention weights. It is formulated as an ANN.
}
\value{A (fitted) object of class \code{ALCOVElearning}}
\references{
Kruschke, J. K. (1992). ALCOVE: An exemplar-based connectionist model of 
category learning. \emph{Psychological Review}, \emph{99}, 22-44.
}
\author{Maarten Speekenbrink}
\examples{
## open weather prediction data
data(WPT)
controls <- subset(WPT,id \%in\% paste("C",1:16,sep=""))
## initialize model
mod <- ALCOVElearning(y~x1+x2+x3+x4,data=controls,
  fix=list(r=TRUE,q=TRUE),ntimes=rep(200,16))
## estimate free parameters (discounting first 5 trials)
\dontrun{mod <- fit(mod,discount=5)}
summary(mod)
}

