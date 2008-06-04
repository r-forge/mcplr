\name{GCM}
\alias{GCM}
\alias{GCM-class}
\alias{GCMinterval-class}
\alias{GCMnominal-class}

\title{Generalized Context Model}
\description{
Constructs a Generalized Context Model.
}
\usage{
GCM(formula,level=c("nominal","interval"),distance=gcm.distance("cityblock"),
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
cues. See details.}
\item{similarity}{a function which converts the T*T distance matrix to a T*T 
similarity matrix. See datails.}
\item{sampling}{a function which returns a T*T matrix with sampling weights. 
See details.}
\item{parameters}{an (optional) named list with (starting) values of the 
parameters. If no values are supplied, defaults are used. }
\item{fixed}{a logical vector indicating whether model parameters are fixed or 
free}
\item{parStruct}{a ParStruct object. If supplied, the fixed argument above will 
be ignored}
\item{data}{data frame containing the variables listed in the formula argument.}
\item{subset}{subset of data to fit model to.}
\item{ntimes}{a vector with number of observations for each replication.}
\iten{replicate}{logical to indicate whether model parameters are identical for 
each replication in ntimes.}
}
\details{The Generalized Context Model (Nosofsky, 1986) is an exemplar model.
It predicts the value of a criterion at t based on the similarity of a probe cue 
x[t] to stored cues x[t-k], k = 1,...,t-1.

The GCM can be seen as a mixture model. Each encountered exemplar adds a
new component to the mixture. The mixture proportions are defined by the
similarities. See the package manual for more information.

The model implemented by \code{GCM} extends the original GCM (Nosofsky, 1986) 
by allowing (1) a continuous criterion, and (2) memory decay of exemplars.
}
\value{A (fitted) object of class \code{GCMinterval} or \code{GCMnominal}.}

\references{
Nosofsky, R.M. (1986). Attention, Similarity, and the 
Identification-Categorization Relationship. \emph{Journal of Experimental 
Psychology: General}, \emph{115}, 39-57.

Speekenbrink, M. \& Shanks, D.R. (2008). Learning in a changing environment. 
Manuscript submitted for publication.}
\examples{
## open weather prediction data
data(WPT)
## initialize model
mod <- gcm(y~x1+x2+x3+x4-1,data=WPT,ntimes=c(200,200))
## estimate free parameters
mod <- estimate(mod)
summary(mod)
}
\keyword{models}