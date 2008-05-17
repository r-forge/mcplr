\name{gcm}
\alias{gcm}
\title{Generalized Context Model}
\description{
Constructs a Generalized Context Model.
}
\usage{
gcm(formula,level=c("nominal","interval"),distance=gcm.distance("cityblock"),
 similarity=gcm.similarity("exponential"),sampling=gcm.sampling("uniform"),
 parameters=list(),fixed,parStruct,data,subset,ntimes=NULL,replicate=TRUE,
 base=NULL)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to 
 that class): a symbolic description of the model to be fitted. For more details
 of model speecification, see \code{lm} or \code{glm}.}
\item{level}{the measurement level of the dependent variable, either nominal
 or interval.}
\item{distance}{a function which returns a T*T matrix with distances between the
 cues.}
\item{similarity}{a function which converts the distance matrix to a similarity
  matrix.}
\item{sampling}{a function which returns a T*T matrix with sampling weights. 
 See \code{gcm.sampling} for more details.}
\item{parameters}{an (optional) list with (starting) values of the parameters. 
 If no values are supplied, defaults are used. }
}
\details{The Generalized Context Model (Nosofsky, 1986) is an exemplar model.
It predicts the value of a criterion based on the similarity of a probe cue to
stored cues.

The GCM can be seen as a mixture model. Each encountered exemplar adds a
new component to the mixture. The mixture proportions are defined by the
similarities. See the package manual for more information.

The model implemented by \code{gcm} extends the original GCM (Nosofsky, 1986) 
by allowing (1) a continuous criterion, and (2) memory decay of exemplars.
}
\value{A (fitted) object of class \code{GcmInterval}, \code{GcmNominal}, 
  \code{GcmUnconstrainedInterval} or \code{GcmUnconstrainedNominal}}
\references{Nosofsky, R. (1986). Speekenbrink, M. \& Shanks, D.R. }
\examples{
## open weather prediction data
data(WP)
## initialize model
mod <- gcm(y~x1+x2+x3+x4-1,data=WP)
## estimate free parameters
mod <- estimate(mod)
summary(mod)
}
\keyword{models}